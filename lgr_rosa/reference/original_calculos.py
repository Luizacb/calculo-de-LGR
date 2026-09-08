import numpy as np
import matplotlib.pyplot as plt
import sympy

# ============================================================
# Formatação LaTeX
# ============================================================

def poly_latex(coefs, var="s"):
    """Converte array de coeficientes numpy para string LaTeX."""
    grau = len(coefs) - 1
    partes = []
    for k, c in enumerate(coefs):
        exp = grau - k
        if abs(c) < 1e-12:
            continue
        ac = abs(c)
        if exp == 0:
            cs = f"{ac:g}"
        elif abs(ac - 1) < 1e-12:
            cs = ""
        else:
            cs = f"{ac:g}"
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


def cx_latex(c, dec=4):
    """Formata número complexo para LaTeX."""
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


def fatorado_latex(raizes, var="s"):
    """Forma fatorada do polinômio em LaTeX."""
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
    vistos = []
    contagem = {}
    for f in fatores:
        if f not in contagem:
            contagem[f] = 0
            vistos.append(f)
        contagem[f] += 1
    partes = []
    for f in vistos:
        n = contagem[f]
        partes.append(f if n == 1 else f"{f}^{{{n}}}")
    return "".join(partes)


def formatar_complexo(c):
    """Formato texto simples (uso interno)."""
    r = round(c.real, 4)
    i = round(c.imag, 4)
    if abs(i) < 1e-10:
        return f"{r}"
    s = "+" if i >= 0 else "-"
    return f"{r} {s} {abs(i)}j"


# ============================================================
# Funções computacionais
# ============================================================

def parse_coefs(texto):
    try:
        vals = [float(x) for x in texto.strip().split()]
        if len(vals) == 0:
            return None
        return np.array(vals)
    except Exception:
        return None


def separar_jw(coefs):
    grau = len(coefs) - 1
    re = {}
    im = {}
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


def fazer_passo1(nG, dG, nH, dH):
    num = np.convolve(nG, nH)
    den = np.convolve(dG, dH)
    n = max(len(num), len(den))
    num = np.pad(num, (n - len(num), 0))
    den = np.pad(den, (n - len(den), 0))
    return num, den


def achar_segmentos_eixo_real(zeros, polos):
    reais = []
    for p in polos:
        if abs(p.imag) < 1e-8:
            reais.append(p.real)
    for z in zeros:
        if abs(z.imag) < 1e-8:
            reais.append(z.real)
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


def calcular_assintotas(zeros, polos):
    np_ = len(polos)
    nz = len(zeros)
    diff = np_ - nz
    if diff == 0:
        return None, []
    sigma = (np.sum(polos).real - np.sum(zeros).real) / diff
    angs = [(2 * q + 1) * 180.0 / diff for q in range(diff)]
    return sigma, angs


def achar_breakaway(num, den, polos, zeros):
    dN = np.polyder(num)
    dD = np.polyder(den)
    eq = np.polysub(np.convolve(num, dD), np.convolve(den, dN))
    raizes = np.roots(eq)
    reais_pz = []
    for p in polos:
        if abs(p.imag) < 1e-8:
            reais_pz.append(p.real)
    for z in zeros:
        if abs(z.imag) < 1e-8:
            reais_pz.append(z.real)
    pts = []
    for r in raizes:
        vn = np.polyval(num, r)
        Kv = -np.polyval(den, r) / vn if abs(vn) > 1e-12 else None
        if abs(r.imag) < 1e-6:
            rr = r.real
            cont = sum(1 for x in reais_pz if x > rr + 1e-10)
            no_lgr = cont % 2 == 1
            Kr = Kv.real if Kv is not None else np.inf
            if no_lgr and Kr > 0:
                pts.append((rr, Kr))
        else:
            if Kv is not None and abs(Kv.imag) < 1e-6 and Kv.real > 0:
                pts.append((complex(r), Kv.real))
    return pts


def cruzamento_jw(den, num):
    Re_D, Im_D = separar_jw(den)
    Re_N, Im_N = separar_jw(num)
    cross = np.polysub(np.convolve(Re_D, Im_N), np.convolve(Im_D, Re_N))
    while len(cross) > 1 and abs(cross[0]) < 1e-12:
        cross = cross[1:]
    info_jw = {
        "Re_D": Re_D, "Im_D": Im_D,
        "Re_N": Re_N, "Im_N": Im_N,
        "cross": cross,
    }
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
            if not any(abs(K - kk) < 1e-4 and abs(omega - ww) < 1e-4
                       for kk, ww in resultado):
                resultado.append((K, omega))
    return resultado, info_jw


def tabela_routh(den, num):
    K = sympy.Symbol('K', positive=True)
    n = max(len(den), len(num))
    dp = np.pad(den, (n - len(den), 0))
    np_ = np.pad(num, (n - len(num), 0))

    def sym(c):
        r = round(c)
        return sympy.Integer(r) if abs(c - r) < 1e-10 else sympy.nsimplify(c, rational=True)

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
                if c is True or c == sympy.S.true:
                    cond_ltx = r"\forall\; K > 0"
                else:
                    cond_ltx = sympy.latex(c)
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


def ordenar_raizes(prev, curr):
    n = len(curr)
    if n == 0:
        return curr
    out = np.zeros(n, dtype=complex)
    usado = set()
    for i in range(n):
        melhor = None
        dist_min = np.inf
        for j in range(n):
            if j not in usado:
                d = abs(prev[i] - curr[j])
                if d < dist_min:
                    dist_min = d
                    melhor = j
        out[i] = curr[melhor]
        usado.add(melhor)
    return out


def calcular_lgr(num, den, Kmax=None):
    np_ = len(den) - 1
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
    raizes = np.zeros((len(Ks), np_), dtype=complex)
    for i, k in enumerate(Ks):
        poly = np.polyadd(den, k * num)
        r = np.roots(poly)
        if i > 0:
            r = ordenar_raizes(raizes[i - 1], r)
        raizes[i] = r
    return Ks, raizes


def testar_angulo(s_teste, zeros, polos):
    ap = sum(np.degrees(np.angle(s_teste - p)) for p in polos)
    az = sum(np.degrees(np.angle(s_teste - z)) for z in zeros)
    ang = az - ap
    ang_norm = ((ang + 180) % 360) - 180
    pertence = abs(abs(ang_norm) - 180) < 5.0
    K_val = None
    if pertence:
        prod_p = np.prod([abs(s_teste - p) for p in polos]) if len(polos) > 0 else 1.0
        prod_z = np.prod([abs(s_teste - z) for z in zeros]) if len(zeros) > 0 else 1.0
        K_val = prod_p / prod_z if prod_z > 1e-12 else np.inf
    return ang, ang_norm, pertence, K_val


def calcular_K_ponto(s_pt, zeros, polos):
    prod_p = np.prod([abs(s_pt - p) for p in polos]) if len(polos) > 0 else 1.0
    prod_z = np.prod([abs(s_pt - z) for z in zeros]) if len(zeros) > 0 else 1.0
    return prod_p / prod_z if prod_z > 1e-12 else np.inf


# ============================================================
# Funções de plotagem
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


def desenhar_polos_zeros(ax, polos, zeros):
    ax.plot(polos.real, polos.imag, "rx", ms=10, mew=2, label="Polos", zorder=5)
    if len(zeros) > 0:
        ax.plot(zeros.real, zeros.imag, "go", ms=8, mew=2, fillstyle="none",
                label="Zeros", zorder=5)


def finalizar_grafico(ax, xl, yl, titulo="", limites_usr=None):
    if limites_usr is not None:
        ax.set_xlim(limites_usr[0], limites_usr[1])
        ax.set_ylim(limites_usr[2], limites_usr[3])
    else:
        ax.set_xlim(xl)
        ax.set_ylim(-yl, yl)
    ax.axhline(0, color='k', lw=0.5, alpha=0.3)
    ax.axvline(0, color='k', lw=0.5, alpha=0.3)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(r"Real ($\sigma$)")
    ax.set_ylabel(r"Imaginario ($j \omega$)")
    ax.set_title(titulo)
    ax.legend(fontsize=9)


def desenhar_segmentos(ax, segs, xl):
    for i, (a, b) in enumerate(segs):
        a_plot = max(a, xl[0] - 5) if np.isfinite(a) else xl[0] - 5
        b_plot = min(b, xl[1] + 5) if np.isfinite(b) else xl[1] + 5
        lbl = "Segmento LGR" if i == 0 else ""
        ax.plot([a_plot, b_plot], [0, 0], 'b-', linewidth=4, alpha=0.6,
                solid_capstyle='round', label=lbl)


def desenhar_assintotas(ax, sigma_a, angs, xl, yl):
    line_len = max(abs(xl[0]), abs(xl[1]), yl) * 2
    for i, ang_deg in enumerate(angs):
        ang_rad = np.radians(ang_deg)
        dx = line_len * np.cos(ang_rad)
        dy = line_len * np.sin(ang_rad)
        lbl = "Assintotas" if i == 0 else ""
        ax.plot([sigma_a, sigma_a + dx], [0, dy], '--', color='darkorange',
                linewidth=1.5, alpha=0.7, label=lbl)


def desenhar_lgr_fundo(ax, todas_raizes, xl, yl):
    for j in range(todas_raizes.shape[1]):
        ramo = todas_raizes[:, j]
        ax.plot(ramo.real, ramo.imag, '-', color='gray', linewidth=1.5, alpha=0.4)

