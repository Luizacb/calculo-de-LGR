"""
Núcleo matemático do LGR (Lugar Geométrico das Raízes).

Este módulo NÃO depende de Streamlit. Toda a lógica de cálculo aqui
reproduz, termo a termo, o algoritmo do projeto de referência
(igorservo159/calc_lgr, app_lgr.py), preservando:

  - as mesmas bibliotecas numéricas (numpy.roots / numpy.convolve / etc.)
  - as mesmas tolerâncias numéricas
  - a mesma ordem de operações

Qualquer alteração de comportamento em relação ao original está
sinalizada explicitamente com um comentário `# DESVIO DO ORIGINAL:`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import sympy

# Tolerâncias numéricas (idênticas ao original — não alterar sem validar)
TOL_COEF_ZERO = 1e-12       # coeficiente considerado nulo
TOL_RAIZ_REAL = 1e-8        # |Im(raiz)| abaixo disso -> raiz "real" (segmentos/eixo real)
TOL_RAIZ_REAL_BK = 1e-6     # idem, usado em breakaway / cruzamento jw
TOL_ANGULO_LGR = 5.0        # graus de tolerância no critério do ângulo (Passo 11)
TOL_DEDUP_CRUZAMENTO = 1e-4  # dedup de cruzamentos jw (K e omega)
TOL_INTEIRO = 1e-10         # para simbolizar coeficientes como inteiros na tabela de Routh


# ============================================================
# Entrada / parsing
# ============================================================

def parse_coefs(texto: str) -> Optional[np.ndarray]:
    """Converte string de coeficientes ('1 4 0') em array numpy, em ordem
    decrescente de potência de s. Retorna None se inválido ou vazio."""
    try:
        vals = [float(x) for x in texto.strip().split()]
        if len(vals) == 0:
            return None
        return np.array(vals)
    except Exception:
        return None


# ============================================================
# Formatação (LaTeX / texto) — não afeta valores numéricos
# ============================================================

def poly_latex(coefs, var: str = "s") -> str:
    """Converte array de coeficientes numpy para string LaTeX."""
    grau = len(coefs) - 1
    partes = []
    for k, c in enumerate(coefs):
        exp = grau - k
        if abs(c) < TOL_COEF_ZERO:
            continue
        ac = abs(c)
        cs = f"{ac:g}" if not (exp != 0 and abs(ac - 1) < TOL_COEF_ZERO) else ""
        if exp == 0:
            termo = cs
        elif exp == 1:
            termo = f"{cs}{var}" if cs else var
        else:
            termo = f"{cs}{var}^{{{exp}}}" if cs else f"{var}^{{{exp}}}"
        if not partes:
            partes.append(f"-{termo}" if c < 0 else termo)
        else:
            partes.append(f" - {termo}" if c < 0 else f" + {termo}")
    return "".join(partes) if partes else "0"


def cx_latex(c: complex, dec: int = 4) -> str:
    """Formata número complexo para LaTeX (a + bj)."""
    r = round(c.real, dec)
    i = round(c.imag, dec)
    if abs(i) < 1e-10:
        return f"{r:g}"
    if abs(r) < 1e-10:
        if abs(abs(i) - 1) < 1e-10:
            return "j" if i > 0 else "-j"
        return f"{i:g}j"
    sinal = "+" if i >= 0 else "-"
    return f"{r:g} {sinal} {abs(i):g}j"


def formatar_complexo(c: complex) -> str:
    """Formato texto simples a+bj (para exibição fora do LaTeX)."""
    r = round(c.real, 4)
    i = round(c.imag, 4)
    if abs(i) < 1e-10:
        return f"{r}"
    s = "+" if i >= 0 else "-"
    return f"{r} {s} {abs(i)}j"


def fatorado_latex(raizes, var: str = "s") -> str:
    """Forma fatorada do polinômio em LaTeX, agrupando raízes repetidas."""
    if len(raizes) == 0:
        return "1"
    fatores = []
    for r in raizes:
        if abs(r.imag) < 1e-8:
            rv = round(r.real, 4)
            if abs(rv) < 1e-8:
                fatores.append(var)
            elif rv > 0:
                fatores.append(f"({var} - {rv:g})")
            else:
                fatores.append(f"({var} + {abs(rv):g})")
        else:
            fatores.append(rf"\left({var} - ({cx_latex(r)})\right)")
    vistos, contagem = [], {}
    for f in fatores:
        if f not in contagem:
            contagem[f] = 0
            vistos.append(f)
        contagem[f] += 1
    partes = [f if contagem[f] == 1 else f"{f}^{{{contagem[f]}}}" for f in vistos]
    return "".join(partes)


# ============================================================
# Passo 1 — Equação característica
# ============================================================

def montar_malha_aberta(nG, dG, nH, dH):
    """G(s)H(s) = K * num/den. Equação característica: den + K*num = 0."""
    num = np.convolve(nG, nH)
    den = np.convolve(dG, dH)
    n = max(len(num), len(den))
    num = np.pad(num, (n - len(num), 0))
    den = np.pad(den, (n - len(den), 0))
    return num, den


# ============================================================
# Passo 4 — Segmentos do eixo real
# ============================================================

def segmentos_eixo_real(zeros, polos):
    """Segmentos do eixo real que pertencem ao LGR: à esquerda de um número
    ímpar de polos+zeros reais."""
    reais = [p.real for p in polos if abs(p.imag) < TOL_RAIZ_REAL]
    reais += [z.real for z in zeros if abs(z.imag) < TOL_RAIZ_REAL]
    if not reais:
        return []
    fronteiras = sorted(set(round(r, 8) for r in reais), reverse=True)
    segs = []
    for i in range(len(fronteiras) - 1):
        meio = (fronteiras[i] + fronteiras[i + 1]) / 2
        cont = sum(1 for r in reais if r > meio + 1e-10)
        if cont % 2 == 1:
            segs.append((fronteiras[i + 1], fronteiras[i]))
    if len(reais) % 2 == 1:
        segs.append((-np.inf, fronteiras[-1]))
    return segs


# ============================================================
# Passo 7 — Assíntotas
# ============================================================

def assintotas(zeros, polos):
    """Retorna (sigma_a, [angulos_graus]) ou (None, []) se n_p == n_z."""
    n_p, n_z = len(polos), len(zeros)
    diff = n_p - n_z
    if diff == 0:
        return None, []
    sigma = (np.sum(polos).real - np.sum(zeros).real) / diff
    angs = [(2 * q + 1) * 180.0 / diff for q in range(diff)]
    return sigma, angs


# ============================================================
# Passo 8 — Breakaway / break-in
# ============================================================

def breakaway(num, den, polos, zeros):
    """Pontos de saída/entrada via dK/ds = 0 -> D'(s)N(s) - D(s)N'(s) = 0."""
    dN = np.polyder(num)
    dD = np.polyder(den)
    eq = np.polysub(np.convolve(num, dD), np.convolve(den, dN))
    raizes = np.roots(eq)
    reais_pz = [p.real for p in polos if abs(p.imag) < TOL_RAIZ_REAL]
    reais_pz += [z.real for z in zeros if abs(z.imag) < TOL_RAIZ_REAL]
    pts = []
    for r in raizes:
        vn = np.polyval(num, r)
        Kv = -np.polyval(den, r) / vn if abs(vn) > 1e-12 else None
        if abs(r.imag) < TOL_RAIZ_REAL_BK:
            rr = r.real
            cont = sum(1 for x in reais_pz if x > rr + 1e-10)
            no_lgr = cont % 2 == 1
            Kr = Kv.real if Kv is not None else np.inf
            if no_lgr and Kr > 0:
                pts.append((rr, Kr))
        else:
            if Kv is not None and abs(Kv.imag) < TOL_RAIZ_REAL_BK and Kv.real > 0:
                pts.append((complex(r), Kv.real))
    return pts


def equacao_breakaway_latex(num, den):
    """Retorna (N, N', D, D', eq_exibicao) usados na exibição do Passo 8."""
    dN = np.polyder(num)
    dD = np.polyder(den)
    eq_display = np.polysub(np.convolve(num, dD), np.convolve(den, dN))
    return dN, dD, eq_display


# ============================================================
# Passo 9 — Cruzamento com o eixo imaginário (Routh + s=jw)
# ============================================================

def separar_re_im_jw(coefs):
    """Separa um polinômio em Re(coefs) e Im(coefs) ao substituir s=j*omega."""
    grau = len(coefs) - 1
    re, im = {}, {}
    for k, c in enumerate(coefs):
        pot = grau - k
        r = pot % 4
        if r == 0:
            re[pot] = re.get(pot, 0) + c
        elif r == 1:
            im[pot] = im.get(pot, 0) + c
        elif r == 2:
            re[pot] = re.get(pot, 0) - c
        else:
            im[pot] = im.get(pot, 0) - c

    def montar(d):
        if not d:
            return np.array([0.0])
        g = max(d.keys())
        arr = np.zeros(g + 1)
        for p, v in d.items():
            arr[g - p] = v
        return arr

    return montar(re), montar(im)


def cruzamento_jw(den, num):
    """Cruzamentos com o eixo imaginário (K>0). Retorna (lista[(K,omega)], info)."""
    Re_D, Im_D = separar_re_im_jw(den)
    Re_N, Im_N = separar_re_im_jw(num)
    cross = np.polysub(np.convolve(Re_D, Im_N), np.convolve(Im_D, Re_N))
    while len(cross) > 1 and abs(cross[0]) < TOL_COEF_ZERO:
        cross = cross[1:]
    info_jw = {"Re_D": Re_D, "Im_D": Im_D, "Re_N": Re_N, "Im_N": Im_N, "cross": cross}
    if len(cross) <= 1:
        return [], info_jw
    ws = np.roots(cross)
    resultado = []
    for w in ws:
        if abs(w.imag) > 1e-6 or w.real < 1e-8:
            continue
        omega = w.real
        ImN = np.polyval(Im_N, omega)
        ImD = np.polyval(Im_D, omega)
        ReN = np.polyval(Re_N, omega)
        ReD = np.polyval(Re_D, omega)
        if abs(ImN) > 1e-12:
            K = -ImD / ImN
        elif abs(ReN) > 1e-12:
            K = -ReD / ReN
        else:
            continue
        if K > 1e-10:
            if not any(abs(K - kk) < TOL_DEDUP_CRUZAMENTO and abs(omega - ww) < TOL_DEDUP_CRUZAMENTO
                       for kk, ww in resultado):
                resultado.append((K, omega))
    return resultado, info_jw


def tabela_routh(den, num):
    """Tabela de Routh-Hurwitz simbólica (K como parâmetro positivo)."""
    K = sympy.Symbol('K', positive=True)
    n = max(len(den), len(num))
    dp = np.pad(den, (n - len(den), 0))
    np_ = np.pad(num, (n - len(num), 0))

    def sym(c):
        r = round(c)
        return sympy.Integer(r) if abs(c - r) < TOL_INTEIRO else sympy.nsimplify(c, rational=True)

    coefs = [sym(d) + K * sym(nn) for d, nn in zip(dp, np_)]
    grau = len(coefs) - 1
    cols = (grau + 2) // 2

    tab = [[sympy.S.Zero] * cols for _ in range(grau + 1)]
    for j in range(cols):
        if 2 * j < len(coefs):
            tab[0][j] = coefs[2 * j]
        if 2 * j + 1 < len(coefs):
            tab[1][j] = coefs[2 * j + 1]

    for i in range(2, grau + 1):
        piv = tab[i - 1][0]
        if piv == 0:
            # DESVIO DO ORIGINAL (documentado, não corrigido silenciosamente):
            # o original interrompe aqui sem aplicar o método do polinômio
            # auxiliar para o caso de pivô nulo / linha de zeros. Mantido
            # idêntico ao original para equivalência numérica; a interface
            # sinaliza esse caso ao usuário (ver README / auditoria).
            break
        for j in range(cols - 1):
            n_ = tab[i - 1][0] * tab[i - 2][j + 1] - tab[i - 2][0] * tab[i - 1][j + 1]
            tab[i][j] = sympy.simplify(n_ / piv)

    conds = []
    k_crit = set()
    for i in range(grau + 1):
        e = sympy.simplify(tab[i][0])
        if e.has(K):
            try:
                c = sympy.solve(e > 0, K)
                cond_ltx = r"\forall\; K > 0" if (c is True or c == sympy.S.true) else sympy.latex(c)
            except Exception:
                cond_ltx = None
            conds.append((grau - i, e, cond_ltx))
            sols = sympy.solve(sympy.Eq(e, 0), K)
            for s_val in sols:
                if s_val.is_real and s_val > 0:
                    k_crit.add(float(s_val))

    return {
        'tab': tab, 'grau': grau, 'cols': cols,
        'conds': conds, 'k_crits': sorted(k_crit), 'K_sym': K,
    }


# ============================================================
# Traçado completo do LGR (para o gráfico)
# ============================================================

def ordenar_raizes(prev, curr):
    """Reordena `curr` pelo casamento guloso de menor distância com `prev`,
    para que os ramos do LGR fiquem contínuos entre valores de K vizinhos."""
    n = len(curr)
    if n == 0:
        return curr
    out = np.zeros(n, dtype=complex)
    usado = set()
    for i in range(n):
        melhor, dist_min = None, np.inf
        for j in range(n):
            if j not in usado:
                d = abs(prev[i] - curr[j])
                if d < dist_min:
                    dist_min, melhor = d, j
        out[i] = curr[melhor]
        usado.add(melhor)
    return out


def calcular_lgr(num, den, Kmax=None):
    """Varre K e retorna (Ks, raizes[K, ramo]) usados no gráfico completo."""
    n_p = len(den) - 1
    if Kmax is None:
        Kmax = 100.0
        for kt in [100, 500, 1000, 5000]:
            poly = np.polyadd(den, kt * num)
            rr = np.roots(poly)
            if np.max(np.abs(rr)) > 50:
                Kmax = kt
                break
        else:
            Kmax = 1000.0
    k1 = np.linspace(0, 0.1, 200)
    k2 = np.logspace(-1, np.log10(Kmax), 4800)
    Ks = np.unique(np.concatenate([k1, k2]))
    Ks.sort()
    raizes = np.zeros((len(Ks), n_p), dtype=complex)
    for i, k in enumerate(Ks):
        poly = np.polyadd(den, k * num)
        r = np.roots(poly)
        if i > 0:
            r = ordenar_raizes(raizes[i - 1], r)
        raizes[i] = r
    return Ks, raizes


# ============================================================
# Passos 11/12 — Critério do ângulo e do módulo (ponto de teste)
# ============================================================

def testar_angulo(s_teste, zeros, polos):
    """Retorna (ang, ang_normalizado, pertence_ao_lgr, K_no_ponto|None)."""
    ap = sum(np.degrees(np.angle(s_teste - p)) for p in polos)
    az = sum(np.degrees(np.angle(s_teste - z)) for z in zeros)
    ang = az - ap
    ang_norm = ((ang + 180) % 360) - 180
    pertence = abs(abs(ang_norm) - 180) < TOL_ANGULO_LGR
    K_val = None
    if pertence:
        prod_p = np.prod([abs(s_teste - p) for p in polos]) if len(polos) > 0 else 1.0
        prod_z = np.prod([abs(s_teste - z) for z in zeros]) if len(zeros) > 0 else 1.0
        K_val = prod_p / prod_z if prod_z > 1e-12 else np.inf
    return ang, ang_norm, pertence, K_val


def calcular_K_ponto(s_pt, zeros, polos):
    """K = produto das distâncias aos polos / produto das distâncias aos zeros."""
    prod_p = np.prod([abs(s_pt - p) for p in polos]) if len(polos) > 0 else 1.0
    prod_z = np.prod([abs(s_pt - z) for z in zeros]) if len(zeros) > 0 else 1.0
    return prod_p / prod_z if prod_z > 1e-12 else np.inf


# ============================================================
# Limites de gráfico (geometria, não afeta valores numéricos do LGR)
# ============================================================

def limites_grafico(polos, zeros, extra_real=None):
    todos = np.concatenate([polos, zeros]) if len(zeros) > 0 else polos
    spread = max(np.ptp(todos.real), np.ptp(todos.imag), 1.0)
    marg = max(1.0, 0.5 * spread)
    x_pts = list(todos.real)
    if extra_real:
        x_pts.extend(extra_real)
    xl = (min(x_pts) - marg, max(x_pts) + marg)
    yl = max(abs(todos.imag).max() + marg * 0.5, marg)
    return xl, yl, marg


# ============================================================
# Estrutura de resultado agregada (conveniência para a UI)
# ============================================================

@dataclass
class ResultadoLGR:
    nG: np.ndarray
    dG: np.ndarray
    nH: np.ndarray
    dH: np.ndarray
    num: np.ndarray
    den: np.ndarray
    zeros: np.ndarray
    polos: np.ndarray
    segmentos: list = field(default_factory=list)
    sigma_a: Optional[float] = None
    angulos_assintotas: list = field(default_factory=list)
    breakaway_pts: list = field(default_factory=list)
    cruzamentos_jw: list = field(default_factory=list)
    info_jw: dict = field(default_factory=dict)
    routh: dict = field(default_factory=dict)
    Ks_lgr: np.ndarray = None
    raizes_lgr: np.ndarray = None


def calcular_tudo(txt_nG: str, txt_dG: str, txt_nH: str, txt_dH: str) -> Optional[ResultadoLGR]:
    """Ponto de entrada único: executa todos os cálculos a partir das 4
    strings de coeficientes. Retorna None se alguma entrada for inválida."""
    nG, dG, nH, dH = (parse_coefs(t) for t in (txt_nG, txt_dG, txt_nH, txt_dH))
    if any(x is None for x in (nG, dG, nH, dH)):
        return None

    num, den = montar_malha_aberta(nG, dG, nH, dH)
    zeros = np.roots(num)
    polos = np.roots(den)

    segs = segmentos_eixo_real(zeros, polos)
    sigma_a, angs = assintotas(zeros, polos)
    bk_pts = breakaway(num, den, polos, zeros)
    routh = tabela_routh(den, num)
    cruzs, info_jw = cruzamento_jw(den, num)
    Ks_lgr, raizes_lgr = calcular_lgr(num, den)

    return ResultadoLGR(
        nG=nG, dG=dG, nH=nH, dH=dH, num=num, den=den, zeros=zeros, polos=polos,
        segmentos=segs, sigma_a=sigma_a, angulos_assintotas=angs,
        breakaway_pts=bk_pts, cruzamentos_jw=cruzs, info_jw=info_jw,
        routh=routh, Ks_lgr=Ks_lgr, raizes_lgr=raizes_lgr,
    )
