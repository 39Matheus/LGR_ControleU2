import warnings

import numpy as np
import sympy as sp
import streamlit as st

from calculations import AnalisadorLGR

warnings.filterwarnings("ignore")

st.set_page_config(page_title="Calculadora LGR Completa", layout="wide")


# -----------------------------------------------------------------------------
# Funções de interface
# -----------------------------------------------------------------------------
def latex(expr):
    return sp.latex(sp.expand(expr))


def latex_factor(expr):
    return sp.latex(sp.factor(expr))


def fmt_complex(z, digits=4):
    z = complex(z)
    if abs(z.imag) < 10 ** (-digits):
        return f"{z.real:.{digits}f}"
    sign = "+" if z.imag >= 0 else "-"
    return f"{z.real:.{digits}f} {sign} {abs(z.imag):.{digits}f}j"


def fmt_interval(a, b):
    left = "-∞" if np.isneginf(a) else f"{a:.4f}"
    right = "+∞" if np.isposinf(b) else f"{b:.4f}"
    return f"({left}, {right})"


def render_vector_table(vectors):
    rows = []
    for i, v in enumerate(vectors, 1):
        rows.append(
            {
                "Vetor": i,
                "Origem": fmt_complex(v.origem),
                "Destino": fmt_complex(v.destino),
                "ΔRe": f"{v.dx:.4f}",
                "ΔIm": f"{v.dy:.4f}",
                "|vetor|": f"{v.magnitude:.4f}",
                "Ângulo (°)": f"{v.angulo_deg:.4f}",
            }
        )
    st.dataframe(rows, use_container_width=True, hide_index=True)




def render_vector_details(vectors, nome):
    st.write(f"**Cálculo trigonométrico — {nome}:**")
    for i, v in enumerate(vectors, 1):
        st.latex(
            rf"v_{i}=({v.dx:.4f})+j({v.dy:.4f}),\qquad "
            rf"|v_{i}|=\sqrt{{({v.dx:.4f})^2+({v.dy:.4f})^2}}={v.magnitude:.4f}"
        )
        st.latex(
            rf"\theta_{i}=\operatorname{{atan2}}({v.dy:.4f},{v.dx:.4f})="
            rf"{v.angulo_deg:.4f}^\circ"
        )

def render_routh_table(table):
    max_cols = max(len(row) for _, row in table)
    headers = ["Linha"] + [f"Coluna {j + 1}" for j in range(max_cols)]
    markdown = "| " + " | ".join(headers) + " |\n"
    markdown += "| " + " | ".join(["---"] * len(headers)) + " |\n"
    for power, values in table:
        cells = [f"$s^{{{int(power)}}}$"]
        for j in range(max_cols):
            value = values[j] if j < len(values) else 0
            cells.append(f"${sp.latex(sp.sympify(value))}$")
        markdown += "| " + " | ".join(cells) + " |\n"
    st.markdown(markdown)


# -----------------------------------------------------------------------------
# Entrada
# -----------------------------------------------------------------------------
st.title("Lugar Geométrico das Raízes — Resolução em 12 Passos")
st.markdown(
    "Este aplicativo foi reorganizado para separar a **interface Streamlit** dos "
    "**cálculos**. A saída agora mostra as contas intermediárias para facilitar a "
    "reprodução manual da resolução em prova."
)

with st.form("lgr_form"):
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### Planta $G(s)$")
        num_g_input = st.text_input("Numerador $G(s)$", "1, 2")
        den_g_input = st.text_input("Denominador $G(s)$", "1, 4, 0")

    with col2:
        st.markdown("### Sensor $H(s)$")
        num_h_input = st.text_input("Numerador $H(s)$", "1")
        den_h_input = st.text_input("Denominador $H(s)$", "1")

    st.markdown("---")
    st.markdown("### Teste de um ponto no LGR — Passos 11 e 12")
    col_pt1, col_pt2 = st.columns(2)
    with col_pt1:
        teste_real = st.number_input("Parte Real de $s_t$", value=-2.4, format="%0.4f")
    with col_pt2:
        teste_imag = st.number_input("Parte Imaginária de $s_t$", value=0.0, format="%0.4f")

    submit_button = st.form_submit_button(label="Calcular e Gerar Relatório")


if submit_button:
    try:
        num_G = [float(x.strip()) for x in num_g_input.split(",")]
        den_G = [float(x.strip()) for x in den_g_input.split(",")]
        num_H = [float(x.strip()) for x in num_h_input.split(",")]
        den_H = [float(x.strip()) for x in den_h_input.split(",")]
        ponto_teste = complex(teste_real, teste_imag)

        analisador = AnalisadorLGR(num_G, den_G, num_H, den_H)
        p1 = analisador.dados_passo1()
        analisador.calcular_polos_zeros()
        p4 = analisador.calcular_segmentos_eixo_real()
        p7 = analisador.calcular_assintotas()
        p8 = analisador.calcular_pontos_saida_chegada()
        p9 = analisador.calcular_cruzamento_eixo_imaginario()
        p10 = analisador.analisar_angulo_partida_chegada()
        ptest = analisador.analisar_ponto_teste(ponto_teste)

        st.markdown("---")
        st.header("Relatório Passo-a-Passo")

        # ------------------------------------------------------------------
        # Passo 1
        # ------------------------------------------------------------------
        st.subheader("Passo 1 — Escrever o polinômio característico")
        st.latex(r"1 + G(s)H(s) = 0")
        st.latex(
            rf"1 + K P(s) = 1 + K\left({latex_factor(p1['gh'])}\right) = 0"
        )
        st.markdown("**Parte racional sem o ganho $K$:**")
        st.latex(
            rf"P(s)=\frac{{{latex(p1['num_expr'])}}}{{{latex(p1['den_expr'])}}}"
        )
        st.markdown("**Forma fatorada de $P(s)$:**")
        st.latex(
            rf"P(s)={latex_factor(p1['gh'])}=C\frac{{\prod_j(s-z_j)}}{{\prod_i(s-p_i)}},\qquad C={latex(analisador.ganho_constante)}"
        )
        st.markdown("**Polinômio característico da malha fechada:**")
        st.latex(rf"\Phi(s,K) = {latex(p1['char_expr'])}=0")
        st.markdown(r"**Coeficientes de $\Phi(s,K)$, do maior grau para o menor:**")
        st.write([sp.latex(c) for c in p1["char_poly_coeffs"]])
        st.markdown(r"**Derivada de $\Phi(s,K)$ em relação a $s$:**")
        st.latex(rf"\frac{{\partial\Phi}}{{\partial s}}={latex(p1['char_derivative'])}")

        # ------------------------------------------------------------------
        # Passo 2
        # ------------------------------------------------------------------
        st.subheader("Passo 2 — Fatorar $P(s)$ e encontrar pólos e zeros")
        st.latex(rf"P(s)={latex_factor(p1['gh'])}")
        st.markdown("**Pólos de malha aberta:**")
        for p in analisador.polos:
            st.latex(rf"s_p={fmt_complex(p)}")
        st.markdown("**Zeros de malha aberta:**")
        if analisador.zeros:
            for z in analisador.zeros:
                st.latex(rf"s_z={fmt_complex(z)}")
        else:
            st.write("Nenhum zero finito.")
        st.write(f"$n_P={analisador.np}$ e $n_Z={analisador.nz}$.")
        grupos_polos = analisador._group_points(analisador.polos)
        grupos_mult = [g for g in grupos_polos if g["multiplicity"] > 1]
        if grupos_mult:
            for g in grupos_mult:
                st.info(
                    f"Pólo múltiplo em {fmt_complex(g['point'])}: multiplicidade "
                    f"{g['multiplicity']}. No gráfico, o marcador desse pólo recebe menor opacidade para facilitar a visualização dos ramos."
                )

        # ------------------------------------------------------------------
        # Passo 3
        # ------------------------------------------------------------------
        st.subheader("Passo 3 — Assinalar pólos (×) e zeros (○)")
        st.info("Os pólos são marcados com × e os zeros com ○ no plano-s.")
        st.plotly_chart(analisador.fig, use_container_width=True, key="fig_p3")

        # ------------------------------------------------------------------
        # Passo 4
        # ------------------------------------------------------------------
        st.subheader("Passo 4 — Segmentos do eixo real que pertencem ao LGR")
        st.markdown(
            "Regra usada no material: um ponto pertence ao LGR quando há um "
            "**número ímpar de pólos e zeros reais à sua direita**."
        )
        for item in p4["testes"]:
            st.latex(
                rf"s_t={item['ponto_teste']:.4f}:\quad "
                rf"N_{{direita}}={item['elementos_direita']}\;(" + item["paridade"] + ")"
            )
            if item["pertence"]:
                st.success(f"Intervalo {fmt_interval(item['esquerda'], item['direita'])}: pertence ao LGR.")
            else:
                st.write(f"Intervalo {fmt_interval(item['esquerda'], item['direita'])}: não pertence ao LGR.")
        st.markdown("**Segmentos obtidos:**")
        for seg in p4["segmentos"]:
            st.latex(rf"{fmt_interval(seg[0], seg[1])}")

        # ------------------------------------------------------------------
        # Passo 5
        # ------------------------------------------------------------------
        st.subheader("Passo 5 — Número de lugares separados")
        ls = analisador.np
        st.latex(rf"LS=n_P={ls}")
        st.write(f"Logo, o LGR possui **{ls} ramos**.")

        # ------------------------------------------------------------------
        # Passo 6
        # ------------------------------------------------------------------
        st.subheader("Passo 6 — Simetria")
        st.write("Como os coeficientes são reais, o LGR é simétrico em relação ao eixo real.")
        st.latex(r"s=\sigma+j\omega \Longrightarrow \sigma-j\omega")

        # ------------------------------------------------------------------
        # Passo 7
        # ------------------------------------------------------------------
        st.subheader("Passo 7 — Assíntotas")
        if p7["numero"] > 0:
            st.latex(rf"N_A=n_P-n_Z={p7['numero']}")
            st.latex(
                rf"\sigma_A=\frac{{\sum p_i-\sum z_i}}{{n_P-n_Z}}="
                rf"\frac{{{latex(sp.expand(p7['soma_polos']))}-{latex(sp.expand(p7['soma_zeros']))}}}{{{p7['numero']}}}="
                rf"{p7['centro']:.4f}"
            )
            st.markdown("**Ângulos:**")
            for q, angle in enumerate(p7["angulos"]):
                st.latex(rf"\phi_{{A,{q}}}=\frac{{(2q+1)180^\circ}}{{n_P-n_Z}}={angle:.4f}^\circ")
        else:
            st.info(r"Não existem assíntotas para $n_P\le n_Z$.")

        # ------------------------------------------------------------------
        # Passo 8
        # ------------------------------------------------------------------
        st.subheader("Passo 8 — Ponto(s) de saída/chegada no eixo real")
        st.markdown("A resolução segue exatamente a ideia do material: primeiro isolamos $K$ como função de $s$ e depois impomos $dK/ds=0$.")
        st.latex(rf"K(s)=-\frac{{D(s)}}{{N(s)}}=-\frac{{{latex(p8['K_den'])}}}{{{latex(p8['K_num'])}}}")
        st.markdown("**Derivando:**")
        st.latex(rf"\frac{{dK}}{{ds}}={latex(p8['dK_ds'])}")
        st.markdown("**Forma equivalente para zerar o numerador:**")
        st.latex(
            rf"\frac{{dK}}{{ds}}=0\Longleftrightarrow D'(s)N(s)-D(s)N'(s)=0"
        )
        st.latex(
            rf"D'(s)={latex(p8['Dp'])},\qquad N'(s)={latex(p8['Np'])}"
        )
        st.latex(rf"D'(s)N(s)-D(s)N'(s)={latex(p8['equacao_numerador'])}=0")

        if p8["candidatos"]:
            for cand in p8["candidatos"]:
                kval = cand["K"]
                if cand["real"]:
                    ktxt = "indefinido" if kval is None else f"{kval:.6f}"
                    st.write(
                        f"Candidato: $s={cand['s'].real:.6f}$, $K={ktxt}$. "
                        + ("VÁLIDO para $K>0$ e para o segmento real do LGR." if cand["valido"] else "Ignorado (não satisfaz simultaneamente as condições do LGR).")
                    )
                else:
                    st.write(f"Candidato complexo ignorado: $s={fmt_complex(cand['s'])}$.")
        else:
            st.write("Nenhuma raiz da equação de saída/chegada foi encontrada.")

        # ------------------------------------------------------------------
        # Passo 9
        # ------------------------------------------------------------------
        st.subheader("Passo 9 — Cruzamento do eixo imaginário via Routh-Hurwitz")
        st.markdown("**Tabela de Routh-Hurwitz:**")
        render_routh_table(p9["routh"]["table"])
        st.markdown("**Primeira coluna da tabela:**")
        for power, row in p9["routh"]["table"]:
            st.latex(rf"s^{{{int(power)}}}:\quad {latex(row[0])}")

        st.markdown(r"**Condição de estabilidade:** os elementos da primeira coluna devem permanecer positivos (para coeficiente líder positivo).")
        try:
            conds = [sp.Gt(sp.sympify(expr), 0) for expr in p9["routh"]["first_column"]]
            faixa = sp.reduce_inequalities(conds, analisador.K)
            st.latex(rf"\text{{Faixa de estabilidade: }} {sp.latex(faixa)}")
        except Exception:
            st.write("A desigualdade simbólica não pôde ser reduzida automaticamente; use os elementos da primeira coluna acima.")

        if p9["routh"]["used_epsilon"]:
            st.warning(r"Foi necessário introduzir $\epsilon>0$ no cálculo da tabela de Routh porque apareceu um primeiro elemento nulo em uma linha.")

        if p9["routh"]["K_candidates"]:
            st.markdown("**Valores de $K$ candidatos obtidos anulando a primeira coluna:**")
            for kval in p9["routh"]["K_candidates"]:
                st.latex(rf"K={kval:.6f}")

        # Método pedido pelo professor: usar explicitamente a equação da linha s²
        # para as substituições e a determinação do ganho no cruzamento.
        s2 = p9.get("s2_method", {})
        if s2.get("disponivel"):
            st.markdown("**Determinação do cruzamento usando a equação da linha $s^2$:**")
            a, b = s2["linha_s2"]
            st.latex(rf"\text{{Linha }}s^2:\quad {latex(a)}s^2+{latex(b)}=0")
            st.latex(rf"s=j\omega\Longrightarrow s^2=-\omega^2\Longrightarrow {latex(s2['aux_jw'])}=0")
            if s2.get("k_expression") is not None:
                st.latex(rf"\text{{Da equação de }}s^2:\quad K={latex(s2['k_expression'])}")
            st.latex(rf"Re[\Phi(j\omega,K)]={latex(s2['char_re'])}")
            st.latex(rf"Im[\Phi(j\omega,K)]={latex(s2['char_im'])}")
            st.latex(rf"\frac{{Im[\Phi(j\omega,K)]}}{{\omega}}=0\quad(\omega\ne0)\Longrightarrow {latex(s2['omega_equation'])}=0")
            for omega in s2.get("omega_values", []):
                st.latex(rf"\omega={omega:.6f}")
            for calc in s2.get("calculos", []):
                if s2.get("k_expression") is not None:
                    k_sub = sp.N(s2["k_expression"].subs(analisador.w, calc["w"]))
                    st.latex(
                        rf"\omega={calc['w']:.6f}\Longrightarrow "
                        rf"K={latex(s2['k_expression'])}\Big|_{{\omega={calc['w']:.6f}}}="
                        rf"{float(k_sub):.6f}"
                    )
                else:
                    st.latex(
                        rf"\omega={calc['w']:.6f}\Longrightarrow "
                        rf"{latex(s2['aux_jw'])}=0\Longrightarrow K={calc['K']:.6f}"
                    )
                st.latex(
                    rf"\text{{Verificação: }}Re[\Phi]={calc['real_residual']:.3e},\quad "
                    rf"Im[\Phi]={calc['imag_residual']:.3e}"
                )
        else:
            st.info(r"A equação da linha $s^2$ não está disponível para este grau do polinômio; foi usado o procedimento geral de cruzamento.")

        st.markdown(r"**Polinômio característico em $s=j\omega$:**")
        st.latex(rf"D(j\omega)={latex(p9['D_jw'])}")
        st.latex(rf"N(j\omega)={latex(p9['N_jw'])}")
        st.latex(
            rf"Re[D(j\omega)]={latex(p9['re_D'])},\quad Im[D(j\omega)]={latex(p9['im_D'])}"
        )
        st.latex(
            rf"Re[N(j\omega)]={latex(p9['re_N'])},\quad Im[N(j\omega)]={latex(p9['im_N'])}"
        )
        st.markdown("**Eliminando $K$ entre as duas equações:**")
        st.latex(rf"{latex(p9['eliminacao'])}=0")

        if p9["cruzamentos"]:
            for c in p9["cruzamentos"]:
                st.success(rf"Cruzamento: $s=\pm j{c['w']:.6f}$ para $K={c['K']:.6f}$.")
                if p9["routh"]["crossing_candidates"]:
                    for detail in p9["routh"]["crossing_candidates"]:
                        if abs(detail["K"] - c["K"]) < 1e-5 and detail.get("auxiliary"):
                            aux = detail["auxiliary"]
                            st.markdown("**Polinômio auxiliar da linha de Routh:**")
                            st.latex(rf"A(s)={latex(aux['polynomial'])}=0")
                            roots_txt = ", ".join(fmt_complex(r) for r in aux["roots"])
                            st.latex(rf"s={roots_txt}")
        else:
            st.info("Não foi identificado cruzamento do eixo imaginário para $K>0$.")

        # ------------------------------------------------------------------
        # Passo 10
        # ------------------------------------------------------------------
        st.subheader("Passo 10 — Ângulos de partida e chegada")
        if not p10["resultados"]:
            st.write("Não há pólos/zeros complexos que exijam ângulo de partida/chegada.")
        for item in p10["resultados"]:
            if item["tipo"] == "partida":
                st.markdown(f"**Ângulo de partida em $s={fmt_complex(item['ponto'])}$**")
                if item["vetores_polos"]:
                    st.write("Ângulos provenientes dos demais pólos:")
                    render_vector_table(item["vetores_polos"])
                    render_vector_details(item["vetores_polos"], "demais pólos")
                if item["vetores_zeros"]:
                    st.write("Ângulos provenientes dos zeros:")
                    render_vector_table(item["vetores_zeros"])
                    render_vector_details(item["vetores_zeros"], "zeros")
                st.latex(rf"\sum\theta_i={item['soma_polos']:.4f}^\circ,\qquad \sum\phi_j={item['soma_zeros']:.4f}^\circ")
                st.latex(rf"\theta_{{partida}}=180^\circ-\sum\theta_i+\sum\phi_j={item['angulo']:.4f}^\circ")
            else:
                st.markdown(f"**Ângulo de chegada em $s={fmt_complex(item['ponto'])}$**")
                if item["vetores_zeros"]:
                    st.write("Ângulos provenientes dos demais zeros:")
                    render_vector_table(item["vetores_zeros"])
                    render_vector_details(item["vetores_zeros"], "demais zeros")
                if item["vetores_polos"]:
                    st.write("Ângulos provenientes dos pólos:")
                    render_vector_table(item["vetores_polos"])
                    render_vector_details(item["vetores_polos"], "pólos")
                st.latex(rf"\sum\phi_i={item['soma_zeros']:.4f}^\circ,\qquad \sum\theta_j={item['soma_polos']:.4f}^\circ")
                st.latex(rf"\theta_{{chegada}}=180^\circ-\sum\phi_i+\sum\theta_j={item['angulo']:.4f}^\circ")

        # ------------------------------------------------------------------
        # Passo 11
        # ------------------------------------------------------------------
        st.subheader(f"Passo 11 — Teste da condição de ângulo no ponto $s_t={ponto_teste}$")
        st.markdown("Para deixar a conta reproduzível na folha, o app explicita cada vetor usando trigonometria:")
        st.latex(r"\theta=\operatorname{atan2}(\Delta Im,\Delta Re)")
        st.latex(r"|v|=\sqrt{(\Delta Re)^2+(\Delta Im)^2}")

        if ptest["vetores_polos"]:
            st.write("**Vetores do ponto de teste até os pólos:**")
            render_vector_table(ptest["vetores_polos"])
            render_vector_details(ptest["vetores_polos"], "pólos do ponto de teste")
        if ptest["vetores_zeros"]:
            st.write("**Vetores do ponto de teste até os zeros:**")
            render_vector_table(ptest["vetores_zeros"])
            render_vector_details(ptest["vetores_zeros"], "zeros do ponto de teste")

        st.latex(
            rf"\arg(C)={ptest['fase_constante_deg']:.4f}^\circ,\qquad "
            rf"\sum\theta_i={ptest['soma_angulos_polos']:.4f}^\circ,\qquad "
            rf"\sum\phi_j={ptest['soma_angulos_zeros']:.4f}^\circ"
        )
        st.latex(
            rf"\angle P(s_t)=\arg(C)+\sum\phi_j-\sum\theta_i="
            rf"{ptest['fase_bruta']:.4f}^\circ\equiv {ptest['fase_mod']:.4f}^\circ\pmod{{360^\circ}}"
        )

        if ptest["pertence"]:
            st.success("O ponto satisfaz a condição de ângulo e pertence ao LGR para $K>0$.")
        else:
            st.error("O ponto não satisfaz a condição de ângulo para $K>0$.")

        # ------------------------------------------------------------------
        # Passo 12
        # ------------------------------------------------------------------
        st.subheader("Passo 12 — Determinar o ganho $K$ pelo critério de módulo")
        if ptest["pertence"]:
            st.latex(r"|K P(s_t)|=1\Longrightarrow K=\frac{\prod_i|s_t-p_i|}{|C|\prod_j|s_t-z_j|}")
            st.latex(
                rf"K=\frac{{{ptest['produto_mod_polos']:.6f}}}{{"
                rf"|{latex(ptest['ganho_constante'])}|\times {ptest['produto_mod_zeros']:.6f}}}="
                rf"{ptest['K']:.6f}"
            )
            st.markdown(
                rf"Aqui, $C={latex(ptest['ganho_constante'])}$ é o ganho constante da forma fatorada de $P(s)$."
            )
            st.info(f"Logo, para $s_t={ponto_teste}$, o ganho é **K = {ptest['K']:.6f}**.")
        else:
            st.write("Como a condição de ângulo não foi satisfeita, o cálculo de $K$ é descartado para $K>0$.")

        # ------------------------------------------------------------------
        # Gráfico final
        # ------------------------------------------------------------------
        st.subheader("Gráfico Final — LGR com os elementos calculados")
        analisador.adicionar_elementos_geometricos()
        analisador.adicionar_marcadores_especiais()
        analisador.adicionar_ponto_teste(ponto_teste)
        analisador.calcular_lgr_exato()
        st.plotly_chart(analisador.fig, use_container_width=True, key="fig_final")

        # ------------------------------------------------------------------
        # Exportação da resolução
        # ------------------------------------------------------------------
        st.markdown("---")
        st.subheader("Exportar resolução")
        st.markdown(
            "Use o botão abaixo para gerar um PDF contendo os dados informados e "
            "toda a resolução calculada nos 12 passos, incluindo tabelas, contas "
            "intermediárias e o diagrama final do LGR."
        )

        pdf_resolucao = analisador.gerar_pdf_resolucao(
            p1=p1,
            p4=p4,
            p7=p7,
            p8=p8,
            p9=p9,
            p10=p10,
            ptest=ptest,
            ponto_teste=ponto_teste,
        )

        col_pdf1, col_pdf2 = st.columns(2)

        with col_pdf1:
            st.download_button(
                label="Baixar resolução completa em PDF",
                data=pdf_resolucao,
                file_name="resolucao_LGR_completa.pdf",
                mime="application/pdf",
                on_click="ignore",
                use_container_width=True,
            )

        pdf_prova = analisador.gerar_pdf_modo_prova(
            p1=p1,
            p4=p4,
            p7=p7,
            p8=p8,
            p9=p9,
            p10=p10,
            ptest=ptest,
            ponto_teste=ponto_teste,
        )

        with col_pdf2:
            st.download_button(
                label="Baixar PDF — Modo Prova",
                data=pdf_prova,
                file_name="resolucao_LGR_modo_prova.pdf",
                mime="application/pdf",
                on_click="ignore",
                use_container_width=True,
            )

    except Exception as exc:
        st.error("Erro na execução. Verifique os coeficientes informados.")
        st.exception(exc)
