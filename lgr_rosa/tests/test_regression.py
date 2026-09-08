"""
Testes de regressão: comparam core.calculos (nova implementação) com
reference.original_calculos (funções extraídas verbatim do app_lgr.py
original) para os MESMOS conjuntos de entrada.

Critério de sucesso: resultados numericamente equivalentes dentro de
tolerância razoável (1e-6 a 1e-9 conforme o cálculo).
"""
import sys
import os
import math

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "core"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import calculos as novo               # noqa: E402
import original_calculos as orig      # noqa: E402


# ------------------------------------------------------------------
# Conjuntos de teste (nG, dG, nH, dH) -> cobrem: polos reais, polos
# complexos, polo na origem, ordem alta, sem zeros finitos, com zeros,
# realimentação H(s) != 1, e um caso "extremo" de alta ordem.
# ------------------------------------------------------------------
CASOS = {
    "default_original": ("1 2", "1 4 0", "1", "1 1"),
    "classico_3_polos_reais": ("1", "1 6 11 6", "1", "1"),          # p=-1,-2,-3
    "com_zero_finito": ("1 3", "1 5 6", "1", "1"),                   # z=-3, p=-2,-3
    "polos_complexos": ("1", "1 2 5", "1", "1"),                     # p = -1 ± 2j
    "polo_na_origem": ("1", "1 2 0", "1", "1"),                      # p=0, -2
    "realimentacao_nao_unitaria": ("1", "1 2 0", "1 1", "1 5"),
    "ordem_alta": ("1 1", "1 6 11 6 0", "1", "1 2"),
    "sem_zeros_np_igual_nz_menos1": ("1 4", "1 1 0", "1", "1"),
    "quarta_ordem_complexo": ("1", "1 4 8 8 4", "1", "1"),
}


def _rot_key(c: complex, dec=6):
    """Chave de ordenação estável e independente de convenção de sinal de
    zero, para comparar conjuntos de raízes sem depender da ordem retornada
    por np.roots (que já deve ser idêntica, mas isso deixa o teste robusto)."""
    return (round(c.real, dec), round(c.imag, dec))


def _raizes_equivalentes(a, b, tol=1e-6):
    if len(a) != len(b):
        return False
    a_sorted = sorted(a, key=_rot_key)
    b_sorted = sorted(b, key=_rot_key)
    return all(abs(x - y) < tol for x, y in zip(a_sorted, b_sorted))


@pytest.mark.parametrize("nome", CASOS.keys())
def test_polos_zeros_equivalentes(nome):
    txt_nG, txt_dG, txt_nH, txt_dH = CASOS[nome]
    nG = novo.parse_coefs(txt_nG)
    dG = novo.parse_coefs(txt_dG)
    nH = novo.parse_coefs(txt_nH)
    dH = novo.parse_coefs(txt_dH)

    num_n, den_n = novo.montar_malha_aberta(nG, dG, nH, dH)
    num_o, den_o = orig.fazer_passo1(nG, dG, nH, dH)
    assert np.allclose(num_n, num_o)
    assert np.allclose(den_n, den_o)

    zeros_n, polos_n = np.roots(num_n), np.roots(den_n)
    zeros_o, polos_o = np.roots(num_o), np.roots(den_o)
    assert _raizes_equivalentes(zeros_n, zeros_o)
    assert _raizes_equivalentes(polos_n, polos_o)


@pytest.mark.parametrize("nome", CASOS.keys())
def test_segmentos_eixo_real(nome):
    txt_nG, txt_dG, txt_nH, txt_dH = CASOS[nome]
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)
    zeros, polos = np.roots(num), np.roots(den)

    segs_n = novo.segmentos_eixo_real(zeros, polos)
    segs_o = orig.achar_segmentos_eixo_real(zeros, polos)
    assert len(segs_n) == len(segs_o)
    for (a1, b1), (a2, b2) in zip(segs_n, segs_o):
        assert math.isclose(a1, a2, abs_tol=1e-6) or (math.isinf(a1) and math.isinf(a2))
        assert math.isclose(b1, b2, abs_tol=1e-6) or (math.isinf(b1) and math.isinf(b2))


@pytest.mark.parametrize("nome", CASOS.keys())
def test_assintotas(nome):
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)
    zeros, polos = np.roots(num), np.roots(den)

    sigma_n, angs_n = novo.assintotas(zeros, polos)
    sigma_o, angs_o = orig.calcular_assintotas(zeros, polos)
    if sigma_o is None:
        assert sigma_n is None
        assert angs_n == []
    else:
        assert math.isclose(sigma_n, sigma_o, abs_tol=1e-9)
        assert len(angs_n) == len(angs_o)
        assert all(math.isclose(a, b, abs_tol=1e-9) for a, b in zip(sorted(angs_n), sorted(angs_o)))


@pytest.mark.parametrize("nome", CASOS.keys())
def test_breakaway(nome):
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)
    zeros, polos = np.roots(num), np.roots(den)

    pts_n = novo.breakaway(num, den, polos, zeros)
    pts_o = orig.achar_breakaway(num, den, polos, zeros)
    assert len(pts_n) == len(pts_o)
    reais_n = sorted(p for p, k in pts_n if not isinstance(p, complex))
    reais_o = sorted(p for p, k in pts_o if not isinstance(p, complex))
    assert np.allclose(reais_n, reais_o, atol=1e-6)


@pytest.mark.parametrize("nome", CASOS.keys())
def test_cruzamento_jw(nome):
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)

    cruzs_n, info_n = novo.cruzamento_jw(den, num)
    cruzs_o, info_o = orig.cruzamento_jw(den, num)
    assert np.allclose(info_n["cross"], info_o["cross"])
    assert len(cruzs_n) == len(cruzs_o)
    kn = sorted(k for k, w in cruzs_n)
    ko = sorted(k for k, w in cruzs_o)
    assert np.allclose(kn, ko, atol=1e-4)


@pytest.mark.parametrize("nome", CASOS.keys())
def test_tabela_routh_k_criticos(nome):
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)

    r_n = novo.tabela_routh(den, num)
    r_o = orig.tabela_routh(den, num)
    assert r_n["grau"] == r_o["grau"]
    assert np.allclose(sorted(r_n["k_crits"]), sorted(r_o["k_crits"]), atol=1e-6)


@pytest.mark.parametrize("nome", CASOS.keys())
def test_lgr_completo_ramos(nome):
    """Compara os ramos do LGR (mesma varredura de K, mesmo Kmax automático,
    mesmo pareamento guloso de raízes)."""
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)

    Ks_n, raizes_n = novo.calcular_lgr(num, den)
    Ks_o, raizes_o = orig.calcular_lgr(num, den)
    assert np.allclose(Ks_n, Ks_o)
    assert np.allclose(raizes_n, raizes_o)


# ------------------------------------------------------------------
# Pontos de teste para os critérios de ângulo/módulo (Passos 11/12)
# ------------------------------------------------------------------
PONTOS_TESTE = [
    complex(-1, 0),
    complex(-2, 1),
    complex(0.5, 0.5),
    complex(-3, 2),
    complex(-1, -1),
]


@pytest.mark.parametrize("nome", CASOS.keys())
@pytest.mark.parametrize("s0", PONTOS_TESTE)
def test_criterio_angulo_e_modulo(nome, s0):
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)
    zeros, polos = np.roots(num), np.roots(den)

    ang_n, angn_n, pert_n, k_n = novo.testar_angulo(s0, zeros, polos)
    ang_o, angn_o, pert_o, k_o = orig.testar_angulo(s0, zeros, polos)
    assert math.isclose(ang_n, ang_o, abs_tol=1e-9)
    assert math.isclose(angn_n, angn_o, abs_tol=1e-9)
    assert pert_n == pert_o

    kp_n = novo.calcular_K_ponto(s0, zeros, polos)
    kp_o = orig.calcular_K_ponto(s0, zeros, polos)
    if math.isinf(kp_o):
        assert math.isinf(kp_n)
    else:
        assert math.isclose(kp_n, kp_o, rel_tol=1e-9)


# ------------------------------------------------------------------
# Formatação (LaTeX) — precisa reproduzir texto idêntico
# ------------------------------------------------------------------
@pytest.mark.parametrize("nome", CASOS.keys())
def test_formatacao_latex_identica(nome):
    nG, dG, nH, dH = (novo.parse_coefs(t) for t in CASOS[nome])
    num, den = novo.montar_malha_aberta(nG, dG, nH, dH)
    zeros, polos = np.roots(num), np.roots(den)

    assert novo.poly_latex(num) == orig.poly_latex(num)
    assert novo.poly_latex(den) == orig.poly_latex(den)
    assert novo.fatorado_latex(zeros) == orig.fatorado_latex(zeros)
    assert novo.fatorado_latex(polos) == orig.fatorado_latex(polos)
    for p in polos:
        assert novo.cx_latex(p) == orig.cx_latex(p)


# ------------------------------------------------------------------
# Casos de erro / entrada inválida
# ------------------------------------------------------------------
@pytest.mark.parametrize("txt", ["", "abc", "1 2 x", "   "])
def test_parse_coefs_entrada_invalida(txt):
    assert novo.parse_coefs(txt) is None
    assert orig.parse_coefs(txt) is None


def test_calcular_tudo_entrada_invalida_retorna_none():
    assert novo.calcular_tudo("1 2", "abc", "1", "1 1") is None


def test_ponto_coincide_com_zero_da_K_infinito():
    """Reproduz o caso de erro: K não pode ser calculado se s0 == zero."""
    nG, dG, nH, dH = "1 2", "1 4 0", "1", "1 1"
    nG_a, dG_a, nH_a, dH_a = (novo.parse_coefs(t) for t in (nG, dG, nH, dH))
    num, den = novo.montar_malha_aberta(nG_a, dG_a, nH_a, dH_a)
    zeros = np.roots(num)
    s0 = zeros[0]  # coincide exatamente com um zero
    k_n = novo.calcular_K_ponto(s0, zeros, np.roots(den))
    k_o = orig.calcular_K_ponto(s0, zeros, np.roots(den))
    assert math.isinf(k_n) and math.isinf(k_o)


# ------------------------------------------------------------------
# Caso extremo: sistema já estável (sem cruzamento jw, sem K crítico)
# ------------------------------------------------------------------
def test_caso_sem_cruzamento_imaginario():
    nG, dG, nH, dH = "1", "1 3 2", "1", "1"  # p=-1,-2, sem zeros -> estável p/ todo K>0
    nG_a, dG_a, nH_a, dH_a = (novo.parse_coefs(t) for t in (nG, dG, nH, dH))
    num, den = novo.montar_malha_aberta(nG_a, dG_a, nH_a, dH_a)
    cruzs_n, _ = novo.cruzamento_jw(den, num)
    cruzs_o, _ = orig.cruzamento_jw(den, num)
    assert cruzs_n == [] and cruzs_o == []


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
