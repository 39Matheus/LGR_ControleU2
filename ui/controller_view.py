import sympy as sp
import streamlit as st

from core.controllers import ControllerDesigner
from core.presets import EXERCISES
from core.specs import DesignSpecs
from reports.pdf import build_controller_pdf
from .renderers import render_controller_design

def _parse(text):
    vals=[float(x.strip().replace(",",".")) for x in text.split(";") if x.strip()]
    if not vals:
        raise ValueError("Coeficientes vazios.")
    if abs(vals[0])<1e-15:
        raise ValueError("O primeiro coeficiente não pode ser zero. Remova zeros à esquerda.")
    return vals

def _poly_latex(coeffs):
    s=sp.symbols("s")
    degree=len(coeffs)-1
    expr=sum(sp.Float(c)*s**(degree-i) for i,c in enumerate(coeffs))
    return sp.latex(sp.expand(expr))

def _show_interpreted_system(data):
    st.markdown("**Sistema interpretado pelo aplicativo:**")
    st.latex(
        rf"G(s)=\frac{{{_poly_latex(data['num_g'])}}}"
        rf"{{{_poly_latex(data['den_g'])}}}"
    )
    st.latex(
        rf"H(s)=\frac{{{_poly_latex(data['num_h'])}}}"
        rf"{{{_poly_latex(data['den_h'])}}}"
    )

def _solve(data):
    d=ControllerDesigner(data["num_g"],data["den_g"],data["num_h"],data["den_h"])
    return d.design(data["controller"],data["specs"]())

def _friendly_error(exc,data=None):
    message=str(exc)
    st.error(f"Falha no projeto: {message}")
    if data is not None:
        _show_interpreted_system(data)
    if "Fase incompatível" in message or "zero fora" in message:
        st.warning(
            "A função de transferência informada produz uma geometria diferente da esperada "
            "para esse controlador e polo desejado. Confira principalmente coeficientes zero "
            "de potências ausentes. Exemplo: s³+4s²+4s deve ser digitado como 1;4;4;0, "
            "e não 1;4;4."
        )

def render_controller_view():
    st.header("Projeto de controlador pelo LGR")
    st.caption("PD, PI e PID (zeros reais iguais), no roteiro da apostila: polo → ângulo → zero(s) → módulo → verificação.")
    source=st.radio("Entrada",["Exercícios da lista U2","Problema personalizado"],horizontal=True)

    if source=="Exercícios da lista U2":
        name=st.selectbox("Exercício",list(EXERCISES))
        data=EXERCISES[name]
        st.info(data["description"])
        if st.button("Resolver exercício",type="primary",use_container_width=True):
            try:
                st.session_state["last_design"]=(_solve(data),data)
            except Exception as exc:
                _friendly_error(exc,data)
    else:
        with st.form("custom_ctrl"):
            kind=st.selectbox("Controlador",["PD","PI","PID"])
            st.caption(
                "Coeficientes em ordem decrescente, separados por ponto e vírgula. "
                "Inclua 0 para toda potência ausente. Ex.: s³+4s²+4s → 1;4;4;0."
            )
            c1,c2=st.columns(2)
            ng=c1.text_input("Numerador G(s)","1")
            dg=c1.text_input("Denominador G(s)","1;3;2")
            nh=c2.text_input("Numerador H(s)","1")
            dh=c2.text_input("Denominador H(s)","1")

            sk=st.radio("Especificação",["Mp + ts","ξ + ωn","Polos desejados"],horizontal=True)
            if sk=="Mp + ts":
                a,b,c=st.columns(3)
                mp=a.number_input("Mp máx (%)",0.01,99.0,10.0)
                ts=b.number_input("ts máx (s)",0.001,value=4.0)
                band=c.selectbox("Critério ts",[5,2])
            elif sk=="ξ + ωn":
                a,b=st.columns(2)
                zeta=a.number_input("ξ",0.001,0.999,0.7)
                wn=b.number_input("ωn",0.001,value=0.5)
            else:
                a,b=st.columns(2)
                pr=a.number_input("Re(s_d)",value=-4.0)
                pi=b.number_input("Im(s_d)",0.001,value=4.0)

            go=st.form_submit_button("Projetar",type="primary")

            if go:
                data=None
                try:
                    if sk=="Mp + ts":
                        specs=DesignSpecs.from_mp_ts(mp,ts,band/100)
                    elif sk=="ξ + ωn":
                        specs=DesignSpecs.from_zeta_wn(zeta,wn)
                    else:
                        specs=DesignSpecs.from_pole(pr,pi)

                    data={
                        "controller":kind,
                        "num_g":_parse(ng),
                        "den_g":_parse(dg),
                        "num_h":_parse(nh),
                        "den_h":_parse(dh),
                        "specs":lambda s=specs:s,
                    }
                    st.session_state["last_interpreted_system"]=data
                    st.session_state["last_design"]=(_solve(data),data)
                    st.session_state.pop("last_error",None)
                except Exception as exc:
                    st.session_state.pop("last_design",None)
                    st.session_state["last_error"]=(str(exc),data)

        if "last_interpreted_system" in st.session_state:
            with st.expander("Conferir função de transferência interpretada",expanded=False):
                _show_interpreted_system(st.session_state["last_interpreted_system"])
        if "last_error" in st.session_state:
            msg,data=st.session_state.pop("last_error")
            _friendly_error(ValueError(msg),data)

    if "last_design" in st.session_state:
        design,data=st.session_state["last_design"]
        st.markdown("---")
        render_controller_design(design,data)
        pdf=build_controller_pdf(design,data["num_g"],data["den_g"],data["num_h"],data["den_h"])
        st.download_button(
            "Baixar resolução — Modo Prova",
            pdf,
            f"projeto_{design.controller_type}_modo_prova.pdf",
            "application/pdf",
            use_container_width=True,
        )
