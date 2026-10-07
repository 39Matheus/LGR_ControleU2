import numpy as np
import sympy as sp
import streamlit as st

from core.lgr import AnalisadorLGR
from reports.pdf import build_lgr_pdf
from .renderers import fmt_complex


def _interval(a,b):
    left="-∞" if np.isneginf(a) else f"{a:.5f}"
    right="+∞" if np.isposinf(b) else f"{b:.5f}"
    return f"({left}, {right})"


def render_lgr_results(data,point):
    a=AnalisadorLGR(
        data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )
    p1=a.dados_passo1()
    pole_zero=a.calcular_polos_zeros()
    p4=a.calcular_segmentos_eixo_real()
    p7=a.calcular_assintotas()
    p8=a.calcular_pontos_saida_chegada()
    p9=a.calcular_cruzamento_eixo_imaginario()
    p10=a.analisar_angulo_partida_chegada()
    pt=a.analisar_ponto_teste(point)

    st.subheader("1 — Polinômio característico")
    st.latex(rf"\Phi(s,K)={sp.latex(p1['char_expr'])}=0")

    st.subheader("2–3 — Polos e zeros da malha aberta")
    st.caption(
        "São as raízes do denominador e do numerador de G(s)H(s) após "
        "cancelamentos exatos de fatores comuns."
    )
    st.write("Polos:",", ".join(fmt_complex(p) for p in a.polos) or "nenhum")
    st.write("Zeros:",", ".join(fmt_complex(z) for z in a.zeros) or "nenhum")
    st.plotly_chart(a.fig,use_container_width=True,key="lgr_open_loop_points")

    st.subheader("4 — Eixo real")
    for x in p4["testes"]:
        st.write(
            f"{_interval(x['esquerda'],x['direita'])}: "
            f"N_direita={x['elementos_direita']} → "
            f"{'LGR' if x['pertence'] else 'não pertence'}"
        )

    st.subheader("5–6 — Ramos e simetria")
    st.latex(rf"LS=n_P={a.np}")
    st.write("O LGR é simétrico em relação ao eixo real.")

    st.subheader("7 — Assíntotas")
    st.write("Centro:",p7.get("centro"),"Ângulos:",p7.get("angulos",[]))

    st.subheader("8 — Saída/chegada")
    st.latex(r"K(s)=-D(s)/N(s),\quad dK/ds=0")
    for c in p8["candidatos"]:
        if c.get("real"):
            st.write(
                f"s={c['s'].real:.6f}, K={c.get('K')} — "
                f"{'válido' if c.get('valido') else 'descartado'}"
            )

    st.subheader("9 — Cruzamento imaginário")
    for c in p9["cruzamentos"]:
        st.latex(rf"s=\pm j{c['w']:.6f},\quad K={c['K']:.6f}")
    if not p9["cruzamentos"]:
        st.write("Nenhum para K>0.")

    st.subheader("10 — Partida/chegada")
    for item in p10["resultados"]:
        st.write(
            f"{item['tipo']}: {fmt_complex(item['ponto'])} → "
            f"{item['angulo']:.5f}°"
        )

    st.subheader("11 — Condição de ângulo")
    st.latex(rf"\angle P(s_t)={pt['fase_mod']:.6f}^\circ")
    if pt["pertence"]:
        st.success("O ponto pertence ao LGR.")
    else:
        st.error("O ponto não pertence ao LGR para K>0.")

    st.subheader("12 — Critério de módulo")
    if pt["pertence"]:
        st.latex(rf"K={pt['K']:.8g}")
    else:
        st.write("K descartado: condição de ângulo não satisfeita.")

    a.adicionar_elementos_geometricos()
    a.adicionar_marcadores_especiais()
    a.adicionar_ponto_teste(point)
    a.calcular_lgr_exato()
    st.subheader("LGR final")
    st.plotly_chart(a.fig,use_container_width=True,key="lgr_final")

    pdf=build_lgr_pdf(a,p1,p4,p7,p8,p9,p10,pt,point)
    st.download_button(
        "Baixar LGR — Modo Prova",
        pdf,
        "resolucao_LGR_modo_prova.pdf",
        "application/pdf",
    )
