"""
Calculadora de LGR (Lugar Geométrico das Raízes) — versão rosa.

Interface própria e reorganizada, otimizada para uso rápido durante prova.
Toda a matemática vem de core/calculos.py (validada por testes de regressão
contra o projeto original em tests/test_regression.py).
"""
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

from core import calculos as calc

st.set_page_config(page_title="LGR Rosa · Calculadora de Sistemas de Controle",
                    layout="wide", page_icon="🌸")

# ============================================================
# Tema rosa
# ============================================================
st.markdown("""
<style>
:root {
    --pink-primary: #d6336c;
    --pink-primary-dark: #a61e4d;
    --pink-light: #fff0f6;
    --pink-card: #ffe3ef;
    --pink-border: #f5c2dd;
}
.stApp { background-color: var(--pink-light); }
h1, h2, h3 { color: var(--pink-primary-dark) !important; }
.katex-display { text-align: left !important; }

div.stButton > button:first-child {
    background-color: var(--pink-primary);
    color: white;
    font-weight: 700;
    border: none;
    border-radius: 10px;
    padding: 0.7rem 2rem;
    font-size: 1.05rem;
    width: 100%;
}
div.stButton > button:first-child:hover { background-color: var(--pink-primary-dark); color: white; }

[data-testid="stExpander"] {
    background-color: white;
    border-radius: 10px;
    border: 1px solid var(--pink-border);
    margin-bottom: 0.5rem;
}
[data-testid="stVerticalBlockBorderWrapper"] {
    border-radius: 12px;
}
.card {
    background: white;
    border-radius: 12px;
    padding: 1.1rem 1.3rem;
    border: 1px solid var(--pink-border);
    margin-bottom: 0.8rem;
}
.destaque {
    background: white;
    border-radius: 12px;
    padding: 1rem 1.3rem;
    border-left: 6px solid var(--pink-primary);
    margin-bottom: 0.6rem;
}
.destaque .label { font-size: 0.8rem; color: #a6336c; font-weight: 600; text-transform: uppercase; }
.destaque .valor { font-size: 1.4rem; color: #3b0a20; font-weight: 700; }
.stTabs [data-baseweb="tab-list"] { gap: 4px; }
.stTabs [data-baseweb="tab"] {
    background-color: var(--pink-card);
    border-radius: 8px 8px 0 0;
    padding: 8px 16px;
}
.stTabs [aria-selected="true"] { background-color: var(--pink-primary) !important; color: white !important; }
</style>
""", unsafe_allow_html=True)

st.title("🌸 LGR Rosa")
st.caption("Calculadora de Lugar Geométrico das Raízes — apoio rápido para provas de Sistemas de Controle")

# ============================================================
# Entrada do sistema
# ============================================================
with st.container(border=True):
    st.subheader("1 · Entrada do sistema")
    st.latex(r"G(s) = K \cdot \frac{N_G(s)}{D_G(s)} \qquad H(s) = \frac{N_H(s)}{D_H(s)}")

    cg1, cg2, ch1, ch2 = st.columns(4)
    with cg1:
        txt_nG = st.text_input("Numerador G(s)", value="1 2",
                                placeholder="ex: 1 2",
                                help="Coeficientes em ordem decrescente de s, separados por espaço")
    with cg2:
        txt_dG = st.text_input("Denominador G(s)", value="1 4 0",
                                placeholder="ex: 1 4 0",
                                help="Coeficientes em ordem decrescente de s")
    with ch1:
        txt_nH = st.text_input("Numerador H(s)", value="1",
                                placeholder="ex: 1",
                                help="Coeficientes em ordem decrescente de s")
    with ch2:
        txt_dH = st.text_input("Denominador H(s)", value="1 1",
                                placeholder="ex: 1 1",
                                help="Coeficientes em ordem decrescente de s")

with st.container(border=True):
    st.subheader("2 · Configurações")
    cc1, cc2, cc3 = st.columns([1, 1, 2])
    with cc1:
        sr = st.number_input("Ponto de teste — parte real", value=0.0, format="%.4f")
    with cc2:
        si = st.number_input("Ponto de teste — parte imaginária", value=0.0, format="%.4f")
    with cc3:
        usar_limites = st.checkbox("Definir limites dos gráficos manualmente", value=False)
        lim_usr = None
        if usar_limites:
            cl1, cl2, cl3, cl4 = st.columns(4)
            with cl1:
                x_min_usr = st.number_input("x min", value=-10.0, format="%.2f")
            with cl2:
                x_max_usr = st.number_input("x max", value=2.0, format="%.2f")
            with cl3:
                y_min_usr = st.number_input("y min", value=-10.0, format="%.2f")
            with cl4:
                y_max_usr = st.number_input("y max", value=10.0, format="%.2f")
            lim_usr = (x_min_usr, x_max_usr, y_min_usr, y_max_usr)

calcular = st.button("🌸 Calcular LGR", type="primary")
if calcular:
    st.session_state["calcular_lgr"] = True

if not st.session_state.get("calcular_lgr", False):
    st.info("Preencha os coeficientes acima e clique em **Calcular LGR** para ver os resultados.")
    st.stop()

resultado = calc.calcular_tudo(txt_nG, txt_dG, txt_nH, txt_dH)
if resultado is None:
    st.error("⚠️ Confira os coeficientes — algum dos quatro campos não é uma lista válida de números.")
    st.stop()

# ============================================================
# Cálculos derivados para exibição (mesmas funções puras do core)
# ============================================================
zeros, polos = resultado.zeros, resultado.polos
num, den = resultado.num, resultado.den
segs = resultado.segmentos
sigma_a, angs = resultado.sigma_a, resultado.angulos_assintotas
bk_pts = resultado.breakaway_pts
routh = resultado.routh
cruzs, info_jw = resultado.cruzamentos_jw, resultado.info_jw
Ks_lgr, todas_raizes = resultado.Ks_lgr, resultado.raizes_lgr

xl, yl, marg = calc.limites_grafico(polos, zeros, [sigma_a] if sigma_a is not None else None)
s_test = complex(sr, si)
ang, ang_n, pert, _ = calc.testar_angulo(s_test, zeros, polos)
K_ponto = calc.calcular_K_ponto(s_test, zeros, polos)


def desenhar_pz(ax):
    ax.plot(polos.real, polos.imag, "rx", ms=10, mew=2, label="Polos", zorder=5)
    if len(zeros) > 0:
        ax.plot(zeros.real, zeros.imag, "go", ms=8, mew=2, fillstyle="none", label="Zeros", zorder=5)


def finalizar(ax, titulo, xlim=xl, ylim=yl):
    if lim_usr is not None:
        ax.set_xlim(lim_usr[0], lim_usr[1])
        ax.set_ylim(lim_usr[2], lim_usr[3])
    else:
        ax.set_xlim(xlim)
        ax.set_ylim(-ylim, ylim)
    ax.axhline(0, color='k', lw=0.5, alpha=0.3)
    ax.axvline(0, color='k', lw=0.5, alpha=0.3)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(r"Real ($\sigma$)")
    ax.set_ylabel(r"Imaginário ($j\omega$)")
    ax.set_title(titulo)
    ax.legend(fontsize=9)


# ============================================================
# Resultado em destaque
# ============================================================
st.markdown("---")
st.subheader("3 · Resultado em destaque")

d1, d2, d3, d4 = st.columns(4)
with d1:
    st.markdown(f"""<div class="destaque"><div class="label">Polos / Zeros</div>
    <div class="valor">{len(polos)} / {len(zeros)}</div></div>""", unsafe_allow_html=True)
with d2:
    txt_sigma = f"{sigma_a:.4f}" if sigma_a is not None else "—"
    st.markdown(f"""<div class="destaque"><div class="label">Centroide σₐ</div>
    <div class="valor">{txt_sigma}</div></div>""", unsafe_allow_html=True)
with d3:
    k_crits = routh['k_crits']
    txt_kcrit = ", ".join(f"{k:.4f}" for k in k_crits) if k_crits else "nenhum"
    st.markdown(f"""<div class="destaque"><div class="label">K crítico (Routh)</div>
    <div class="valor">{txt_kcrit}</div></div>""", unsafe_allow_html=True)
with d4:
    txt_cruz = ", ".join(f"±{w:.3f}j (K={k:.3f})" for k, w in cruzs) if cruzs else "sem cruzamento"
    st.markdown(f"""<div class="destaque"><div class="label">Cruzamento eixo jω</div>
    <div class="valor" style="font-size:1.05rem">{txt_cruz}</div></div>""", unsafe_allow_html=True)

d5, d6 = st.columns(2)
with d5:
    veredito = "✅ pertence ao LGR" if pert else "❌ não pertence ao LGR"
    st.markdown(f"""<div class="destaque"><div class="label">Ponto de teste s₀ = {calc.cx_latex(s_test)}</div>
    <div class="valor" style="font-size:1.1rem">{veredito} — Δθ = {ang_n % 360:.2f}°</div></div>""",
                unsafe_allow_html=True)
with d6:
    txt_k = f"{K_ponto:.6f}" if np.isfinite(K_ponto) else "∞ (coincide com zero)"
    st.markdown(f"""<div class="destaque"><div class="label">K no ponto de teste</div>
    <div class="valor">{txt_k}</div></div>""", unsafe_allow_html=True)

st.code(
    "Polos: " + "; ".join(calc.formatar_complexo(p) for p in polos) + "\n"
    "Zeros: " + ("; ".join(calc.formatar_complexo(z) for z in zeros) if len(zeros) else "nenhum") + "\n"
    f"Equação característica: {calc.poly_latex(den)} + K({calc.poly_latex(num)}) = 0",
    language="text",
)

# ============================================================
# Detalhes por seção (abas)
# ============================================================
st.markdown("---")
st.subheader("4 · Detalhamento")

abas = st.tabs([
    "① Polos & Zeros", "② Eixo real & simetria", "③ Assíntotas",
    "④ Breakaway", "⑤ Routh & cruzamento jω", "⑥ Ângulos partida/chegada",
    "⑦ Ponto de teste", "⑧ Gráfico completo",
])

# --- Aba 1: Polos & Zeros ---
with abas[0]:
    nG_l, dG_l = calc.poly_latex(resultado.nG), calc.poly_latex(resultado.dG)
    nH_l, dH_l = calc.poly_latex(resultado.nH), calc.poly_latex(resultado.dH)
    num_l, den_l = calc.poly_latex(num), calc.poly_latex(den)
    st.latex(rf"G(s) = K \cdot \frac{{{nG_l}}}{{{dG_l}}} \qquad H(s) = \frac{{{nH_l}}}{{{dH_l}}}")
    st.latex(rf"G(s)H(s) = K \cdot \frac{{{num_l}}}{{{den_l}}} \;\;\Longrightarrow\;\; "
             rf"{den_l} + K\left({num_l}\right) = 0")
    st.latex(rf"P(s) = \frac{{{calc.fatorado_latex(zeros)}}}{{{calc.fatorado_latex(polos)}}}")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(f"**Polos** (n_p = {len(polos)})")
        for i, p in enumerate(polos):
            st.latex(rf"p_{{{i+1}}} = {calc.cx_latex(p)}")
    with c2:
        st.markdown(f"**Zeros** (n_z = {len(zeros)})")
        if len(zeros) > 0:
            for i, z in enumerate(zeros):
                st.latex(rf"z_{{{i+1}}} = {calc.cx_latex(z)}")
        else:
            st.markdown("*Nenhum zero finito*")

    fig, ax = plt.subplots(figsize=(9, 5))
    desenhar_pz(ax)
    for i, p in enumerate(polos):
        ax.annotate(rf"$p_{{{i+1}}}$", (p.real, p.imag), xytext=(8, 8),
                    textcoords="offset points", fontsize=9, color='red')
    for i, z in enumerate(zeros):
        ax.annotate(rf"$z_{{{i+1}}}$", (z.real, z.imag), xytext=(8, 8),
                    textcoords="offset points", fontsize=9, color='green')
    finalizar(ax, "Polos e Zeros no plano s")
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

# --- Aba 2: Eixo real & simetria ---
with abas[1]:
    st.markdown("**Regra:** pertencem ao LGR os segmentos do eixo real à esquerda de um número "
                "ímpar de polos e zeros reais.")
    if segs:
        for a, b in segs:
            ea = f"{a:.4f}" if np.isfinite(a) else r"-\infty"
            eb = f"{b:.4f}" if np.isfinite(b) else r"+\infty"
            st.latex(rf"\left[{ea}\;,\; {eb}\right]")
    else:
        st.info("Nenhum segmento no eixo real pertence ao LGR.")

    st.markdown("O LGR é **simétrico em relação ao eixo real**, pois raízes complexas de "
                "polinômios com coeficientes reais sempre ocorrem em pares conjugados.")
    st.latex(rf"n_p = {len(polos)}, \quad n_z = {len(zeros)}, \quad "
             rf"L_s = \max(n_p, n_z) = {max(len(polos), len(zeros))}")

    fig, ax = plt.subplots(figsize=(9, 5))
    desenhar_pz(ax)
    for i, (a, b) in enumerate(segs):
        a_plot = max(a, xl[0] - 5) if np.isfinite(a) else xl[0] - 5
        b_plot = min(b, xl[1] + 5) if np.isfinite(b) else xl[1] + 5
        ax.plot([a_plot, b_plot], [0, 0], 'b-', linewidth=4, alpha=0.6,
                solid_capstyle='round', label="Segmento LGR" if i == 0 else "")
    finalizar(ax, "Segmentos do eixo real pertencentes ao LGR")
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

# --- Aba 3: Assíntotas ---
with abas[2]:
    if sigma_a is not None:
        na = len(polos) - len(zeros)
        st.latex(rf"n_a = n_p - n_z = {len(polos)} - {len(zeros)} = {na}")
        st.latex(r"\sigma_a = \frac{\sum \operatorname{Re}(p_i) - \sum \operatorname{Re}(z_j)}{n_p - n_z}"
                  rf"= {sigma_a:.4f}")
        st.markdown("**Ângulos das assíntotas:**")
        for q, a_deg in enumerate(angs):
            st.latex(rf"q={q}: \quad \phi_a = \frac{{(2\cdot{q}+1)\cdot180^\circ}}{{{na}}} = {a_deg:.1f}^\circ")

        fig, ax = plt.subplots(figsize=(9, 5))
        desenhar_pz(ax)
        for i, (a, b) in enumerate(segs):
            a_plot = max(a, xl[0] - 5) if np.isfinite(a) else xl[0] - 5
            b_plot = min(b, xl[1] + 5) if np.isfinite(b) else xl[1] + 5
            ax.plot([a_plot, b_plot], [0, 0], 'b-', linewidth=4, alpha=0.6, solid_capstyle='round')
        line_len = max(abs(xl[0]), abs(xl[1]), yl) * 2
        for i, ang_deg in enumerate(angs):
            rad = np.radians(ang_deg)
            dx, dy = line_len * np.cos(rad), line_len * np.sin(rad)
            ax.plot([sigma_a, sigma_a + dx], [0, dy], '--', color='darkorange', linewidth=1.5,
                    alpha=0.7, label="Assíntotas" if i == 0 else "")
        ax.plot(sigma_a, 0, "k+", ms=12, mew=2, label=rf"Centroide ($\sigma_a={sigma_a:.2f}$)")
        finalizar(ax, "LGR — Assíntotas")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
    else:
        st.latex(r"n_p = n_z \;\Rightarrow\; \text{sem assíntotas}")

# --- Aba 4: Breakaway ---
with abas[3]:
    dN, dD, bk_eq_display = calc.equacao_breakaway_latex(num, den)
    st.latex(rf"K = -\frac{{D(s)}}{{N(s)}} \;\;\Rightarrow\;\; "
             rf"\frac{{dK}}{{ds}}=0 \;\;\Rightarrow\;\; D'(s)N(s) - D(s)N'(s) = 0")
    st.latex(rf"{calc.poly_latex(bk_eq_display)} = 0")

    if bk_pts:
        st.markdown("**Pontos de breakaway/break-in válidos** (pertencem ao LGR, K>0):")
        for s_bk, k_bk in bk_pts:
            if isinstance(s_bk, complex):
                st.latex(rf"s = {calc.cx_latex(s_bk)}, \quad K = {k_bk:.4f}")
            else:
                st.success(rf"s = {s_bk:.4f}, \quad K = {k_bk:.4f}")
    else:
        st.info("Nenhum ponto de breakaway/break-in válido encontrado.")

    fig, ax = plt.subplots(figsize=(9, 5))
    desenhar_pz(ax)
    for i, (a, b) in enumerate(segs):
        a_plot = max(a, xl[0] - 5) if np.isfinite(a) else xl[0] - 5
        b_plot = min(b, xl[1] + 5) if np.isfinite(b) else xl[1] + 5
        ax.plot([a_plot, b_plot], [0, 0], 'b-', linewidth=4, alpha=0.6, solid_capstyle='round')
    bk_validos = [(s_bk, k_bk) for s_bk, k_bk in bk_pts if not isinstance(s_bk, complex)]
    if bk_validos:
        bk_x = [s_bk for s_bk, _ in bk_validos]
        ax.plot(bk_x, [0] * len(bk_x), 'md', ms=10, mew=2, label="Breakaway/Break-in", zorder=6)
    finalizar(ax, "Pontos de saída/entrada (breakaway)")
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

# --- Aba 5: Routh & cruzamento jw ---
with abas[4]:
    import sympy
    tab, grau, cols, K_sym = routh['tab'], routh['grau'], routh['cols'], routh['K_sym']
    st.markdown("**Tabela de Routh-Hurwitz** (a partir de $D(s) + K\\cdot N(s) = 0$):")
    linhas_ltx = []
    for i in range(grau + 1):
        exp = grau - i
        cells = [sympy.latex(sympy.simplify(tab[i][j])) for j in range(cols)]
        linhas_ltx.append(f"s^{{{exp}}} & " + " & ".join(cells))
    tabela_ltx = (r"\begin{array}{c|" + "c" * cols + "}\n \\hline\n"
                  + (r" \\" + "\n").join(linhas_ltx) + r" \\ \hline" + "\n\\end{array}")
    st.latex(tabela_ltx)

    if routh['conds']:
        st.markdown("**Condições de estabilidade** (primeira coluna > 0):")
        for exp, expr, cond_ltx in routh['conds']:
            if cond_ltx is not None:
                st.latex(rf"s^{{{exp}}}: \quad {sympy.latex(expr)} > 0 \;\Rightarrow\; {cond_ltx}")
            else:
                st.latex(rf"s^{{{exp}}}: \quad {sympy.latex(expr)} > 0")
    if routh['k_crits']:
        st.markdown("**Valores críticos de K:**")
        for kc in routh['k_crits']:
            st.latex(rf"K_{{\text{{crit}}}} = {kc:.4f}")

    st.markdown("---")
    st.markdown("**Método alternativo — substituição $s=j\\omega$:**")
    st.latex(rf"{calc.poly_latex(info_jw['cross'], var=r'\omega')} = 0")
    if cruzs:
        for Kc, wc in cruzs:
            st.latex(rf"\omega = {wc:.4f} \;\Rightarrow\; K = {Kc:.4f} \;\Rightarrow\; "
                     rf"s = \pm\,{wc:.4f}\,j")
    else:
        st.info("O LGR não cruza o eixo imaginário para K > 0.")

    fig, ax = plt.subplots(figsize=(9, 5))
    for j in range(todas_raizes.shape[1]):
        ramo = todas_raizes[:, j]
        ax.plot(ramo.real, ramo.imag, '-', color='gray', linewidth=1.5, alpha=0.4)
    desenhar_pz(ax)
    for Kc, wc in cruzs:
        ax.plot(0, wc, 's', ms=10, color='deeppink', markeredgecolor='#a61e4d', mew=2, zorder=6,
                label=rf"$j\omega={wc:.2f}j\;(K={Kc:.2f})$")
        ax.plot(0, -wc, 's', ms=10, color='deeppink', markeredgecolor='#a61e4d', mew=2, zorder=6)
    finalizar(ax, "LGR — Cruzamento com o eixo imaginário")
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

# --- Aba 6: Ângulos de partida/chegada ---
with abas[5]:
    polos_cx = [p for p in polos if p.imag > 1e-8]
    zeros_cx = [z for z in zeros if z.imag > 1e-8]
    if not (polos_cx or zeros_cx):
        st.info("Sem polos/zeros complexos — este passo não se aplica.")
    else:
        angulos_partida, angulos_chegada = {}, {}
        if polos_cx:
            st.markdown("**Ângulos de partida (polos complexos):**")
            st.latex(r"\theta_d = 180^\circ - \sum_{j\neq k}\angle(p_k-p_j) + \sum_j \angle(p_k-z_j)")
            for pk in polos_cx:
                ang_polos = [np.degrees(np.angle(pk - pj)) for pj in polos if abs(pj - pk) > 1e-10]
                ang_zeros = [np.degrees(np.angle(pk - zj)) for zj in zeros]
                theta = 180.0 - sum(ang_polos) + sum(ang_zeros)
                theta = ((theta + 180) % 360) - 180
                st.latex(rf"p_k={calc.cx_latex(pk)}: \quad \theta_d = {theta % 360:.2f}^\circ "
                         rf"\quad(\text{{conjugado: }}{(-theta) % 360:.2f}^\circ)")
                angulos_partida[pk] = theta
        if zeros_cx:
            st.markdown("**Ângulos de chegada (zeros complexos):**")
            st.latex(r"\theta_a = 180^\circ - \sum_{j\neq k}\angle(z_k-z_j) + \sum_j \angle(z_k-p_j)")
            for zk in zeros_cx:
                ang_zeros = [np.degrees(np.angle(zk - zj)) for zj in zeros if abs(zj - zk) > 1e-10]
                ang_polos = [np.degrees(np.angle(zk - pj)) for pj in polos]
                theta = 180.0 - sum(ang_zeros) + sum(ang_polos)
                theta = ((theta + 180) % 360) - 180
                st.latex(rf"z_k={calc.cx_latex(zk)}: \quad \theta_a = {theta % 360:.2f}^\circ "
                         rf"\quad(\text{{conjugado: }}{(-theta) % 360:.2f}^\circ)")
                angulos_chegada[zk] = theta

        fig, ax = plt.subplots(figsize=(9, 6))
        for j in range(todas_raizes.shape[1]):
            ramo = todas_raizes[:, j]
            ax.plot(ramo.real, ramo.imag, '-', color='gray', linewidth=1.5, alpha=0.4)
        desenhar_pz(ax)
        todos_arr = np.concatenate([polos, zeros]) if len(zeros) > 0 else polos
        spread = max(np.ptp(todos_arr.real), np.ptp(todos_arr.imag), 1.0)
        arrow_len = spread * 0.3
        for pk, theta in angulos_partida.items():
            rad = np.radians(theta)
            dx, dy = arrow_len * np.cos(rad), arrow_len * np.sin(rad)
            ax.annotate('', xy=(pk.real + dx, pk.imag + dy), xytext=(pk.real, pk.imag),
                        arrowprops=dict(arrowstyle='->', color='darkred', lw=2))
            rad_c = np.radians(-theta)
            dxc, dyc = arrow_len * np.cos(rad_c), arrow_len * np.sin(rad_c)
            ax.annotate('', xy=(pk.real + dxc, -pk.imag + dyc), xytext=(pk.real, -pk.imag),
                        arrowprops=dict(arrowstyle='->', color='darkred', lw=2))
        for zk, theta in angulos_chegada.items():
            rad = np.radians(theta)
            dx, dy = arrow_len * np.cos(rad), arrow_len * np.sin(rad)
            ax.annotate('', xy=(zk.real + dx, zk.imag + dy), xytext=(zk.real, zk.imag),
                        arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
        finalizar(ax, "Ângulos de partida/chegada",
                  xlim=(xl[0] - marg * 0.5, xl[1] + marg * 0.5), ylim=yl + marg * 0.5)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

# --- Aba 7: Ponto de teste (critérios de ângulo e módulo) ---
with abas[6]:
    st.latex(r"\sum\angle(s_0-z_j) - \sum\angle(s_0-p_i) = \pm180^\circ(2q+1)")
    st.markdown(f"**Ponto de teste:** $s_0 = {calc.cx_latex(s_test)}$")

    theta_poles = [np.degrees(np.angle(s_test - p)) for p in polos]
    for i, (p, a) in enumerate(zip(polos, theta_poles)):
        st.latex(rf"\theta_{{{i+1}}} = \angle(s_0 - p_{{{i+1}}}) = {a:.2f}^\circ")
    soma_theta = sum(theta_poles)

    if len(zeros) > 0:
        phi_zeros = [np.degrees(np.angle(s_test - z)) for z in zeros]
        for i, (z, a) in enumerate(zip(zeros, phi_zeros)):
            st.latex(rf"\phi_{{{i+1}}} = \angle(s_0 - z_{{{i+1}}}) = {a:.2f}^\circ")
        soma_phi = sum(phi_zeros)
    else:
        soma_phi = 0.0

    delta = soma_theta - soma_phi
    st.latex(rf"\Delta\theta = \sum\theta_i - \sum\phi_j = {delta:.2f}^\circ \quad"
             rf"(\text{{normalizado: }}{ang_n % 360:.2f}^\circ)")
    if pert:
        st.success(rf"O ponto **pertence** ao LGR (Δθ ≈ 180°, tolerância {calc.TOL_ANGULO_LGR}°)")
    else:
        st.warning(rf"O ponto **não pertence** ao LGR (Δθ ≠ 180°, tolerância {calc.TOL_ANGULO_LGR}°)")

    st.markdown("---")
    st.latex(r"K = \frac{\prod_i |s_0 - p_i|}{\prod_j |s_0 - z_j|}")
    prod_p = np.prod([abs(s_test - p) for p in polos]) if len(polos) > 0 else 1.0
    prod_z = np.prod([abs(s_test - z) for z in zeros]) if len(zeros) > 0 else 1.0
    if prod_z > 1e-12:
        st.latex(rf"K = \frac{{{prod_p:.4f}}}{{{prod_z:.4f}}} = {prod_p/prod_z:.6f}")
    else:
        st.error("Não é possível calcular K: o ponto coincide com um zero.")

    fig, ax = plt.subplots(figsize=(9, 5))
    desenhar_pz(ax)
    cor = 'limegreen' if pert else 'red'
    marcador = '*' if pert else 'X'
    ax.plot(s_test.real, s_test.imag, marcador, ms=14, color=cor, markeredgecolor='black',
            label=rf"$s_0={calc.cx_latex(s_test)}$", zorder=6)
    for p in polos:
        ax.plot([p.real, s_test.real], [p.imag, s_test.imag], ':', color='red', alpha=0.4)
    for z in zeros:
        ax.plot([z.real, s_test.real], [z.imag, s_test.imag], ':', color='green', alpha=0.4)
    all_x = list(polos.real) + [s_test.real] + (list(zeros.real) if len(zeros) else [])
    all_y = list(abs(polos.imag)) + [abs(s_test.imag)] + (list(abs(zeros.imag)) if len(zeros) else [])
    finalizar(ax, "Critério de Ângulo", xlim=(min(all_x) - 2, max(all_x) + 2), ylim=max(all_y) + 2)
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

# --- Aba 8: Gráfico completo ---
with abas[7]:
    fig, ax = plt.subplots(figsize=(9, 6))
    for j in range(todas_raizes.shape[1]):
        ramo = todas_raizes[:, j]
        ax.plot(ramo.real, ramo.imag, '-', color='#d6336c', linewidth=2.2, alpha=0.85)
    desenhar_pz(ax)
    if sigma_a is not None:
        ax.plot(sigma_a, 0, "k+", ms=12, mew=2, label=rf"Centroide ($\sigma_a={sigma_a:.2f}$)")
    finalizar(ax, "Lugar Geométrico das Raízes — completo")
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)
