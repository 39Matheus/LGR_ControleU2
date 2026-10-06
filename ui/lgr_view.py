import numpy as np
import sympy as sp
import streamlit as st
from core.lgr import AnalisadorLGR
from reports.pdf import build_lgr_pdf
from .renderers import fmt_complex

def _coeffs(text):
    return [float(x.strip().replace(",",".")) for x in text.split(";") if x.strip()]

def _interval(a,b):
    left="-∞" if np.isneginf(a) else f"{a:.5f}"
    right="+∞" if np.isposinf(b) else f"{b:.5f}"
    return f"({left}, {right})"

def render_lgr_view():
    st.header("Análise do LGR — 12 passos")
    with st.form("lgr_form"):
        c1,c2=st.columns(2)
        ng=c1.text_input("Numerador G(s)","1;2"); dg=c1.text_input("Denominador G(s)","1;4;0")
        nh=c2.text_input("Numerador H(s)","1"); dh=c2.text_input("Denominador H(s)","1")
        a,b=st.columns(2); xr=a.number_input("Re(s_t)",value=-2.0); xi=b.number_input("Im(s_t)",value=2.0)
        run=st.form_submit_button("Calcular LGR",type="primary")
    if not run:
        return
    try:
        point=complex(xr,xi)
        a=AnalisadorLGR(_coeffs(ng),_coeffs(dg),_coeffs(nh),_coeffs(dh))
        p1=a.dados_passo1()
        a.calcular_polos_zeros()
        p4=a.calcular_segmentos_eixo_real()
        p7=a.calcular_assintotas()
        p8=a.calcular_pontos_saida_chegada()
        p9=a.calcular_cruzamento_eixo_imaginario()
        p10=a.analisar_angulo_partida_chegada()
        pt=a.analisar_ponto_teste(point)

        st.subheader("1 — Polinômio característico")
        st.latex(rf"\Phi(s,K)={sp.latex(p1['char_expr'])}=0")

        st.subheader("2–3 — Polos e zeros")
        st.write("Polos:",", ".join(fmt_complex(p) for p in a.polos))
        st.write("Zeros:",", ".join(fmt_complex(z) for z in a.zeros) or "nenhum")
        st.plotly_chart(a.fig,use_container_width=True)

        st.subheader("4 — Eixo real")
        for x in p4["testes"]:
            st.write(f"{_interval(x['esquerda'],x['direita'])}: N_direita={x['elementos_direita']} → {'LGR' if x['pertence'] else 'não pertence'}")

        st.subheader("5–6 — Ramos e simetria")
        st.latex(rf"LS=n_P={a.np}")
        st.write("O LGR é simétrico em relação ao eixo real.")

        st.subheader("7 — Assíntotas")
        st.write("Centro:",p7.get("centro"),"Ângulos:",p7.get("angulos",[]))

        st.subheader("8 — Saída/chegada")
        st.latex(r"K(s)=-D(s)/N(s),\quad dK/ds=0")
        for c in p8["candidatos"]:
            if c.get("real"):
                st.write(f"s={c['s'].real:.6f}, K={c.get('K')} — {'válido' if c.get('valido') else 'descartado'}")

        st.subheader("9 — Cruzamento imaginário")
        for c in p9["cruzamentos"]:
            st.latex(rf"s=\pm j{c['w']:.6f},\quad K={c['K']:.6f}")
        if not p9["cruzamentos"]:
            st.write("Nenhum para K>0.")

        st.subheader("10 — Partida/chegada")
        for item in p10["resultados"]:
            st.write(f"{item['tipo']}: {fmt_complex(item['ponto'])} → {item['angulo']:.5f}°")

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
        st.plotly_chart(a.fig,use_container_width=True)

        pdf=build_lgr_pdf(a,p1,p4,p7,p8,p9,p10,pt,point)
        st.download_button("Baixar LGR — Modo Prova",pdf,"resolucao_LGR_modo_prova.pdf","application/pdf")
    except Exception as exc:
        st.error("Erro na execução.")
        st.exception(exc)
