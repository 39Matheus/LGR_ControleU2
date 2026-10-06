import streamlit as st
from core.controllers import ControllerDesigner
from core.presets import EXERCISES
from core.specs import DesignSpecs
from reports.pdf import build_controller_pdf
from .renderers import render_controller_design

def _parse(text):
    vals=[float(x.strip().replace(",",".")) for x in text.split(";") if x.strip()]
    if not vals: raise ValueError("Coeficientes vazios.")
    return vals

def _solve(data):
    d=ControllerDesigner(data["num_g"],data["den_g"],data["num_h"],data["den_h"])
    return d.design(data["controller"],data["specs"]())

def render_controller_view():
    st.header("Projeto de controlador pelo LGR")
    st.caption("PD, PI e PID (zeros reais iguais), no roteiro da apostila: polo → ângulo → zero(s) → módulo → verificação.")
    source=st.radio("Entrada",["Exercícios da lista U2","Problema personalizado"],horizontal=True)
    if source=="Exercícios da lista U2":
        name=st.selectbox("Exercício",list(EXERCISES))
        data=EXERCISES[name]; st.info(data["description"])
        if st.button("Resolver exercício",type="primary",use_container_width=True):
            try: st.session_state["last_design"]=(_solve(data),data)
            except Exception as exc: st.error(f"Falha no projeto: {exc}")
    else:
        with st.form("custom_ctrl"):
            kind=st.selectbox("Controlador",["PD","PI","PID"])
            st.caption("Coeficientes em ordem decrescente, separados por ponto e vírgula.")
            c1,c2=st.columns(2)
            ng=c1.text_input("Numerador G(s)","1"); dg=c1.text_input("Denominador G(s)","1;3;2")
            nh=c2.text_input("Numerador H(s)","1"); dh=c2.text_input("Denominador H(s)","1")
            sk=st.radio("Especificação",["Mp + ts","ξ + ωn","Polos desejados"],horizontal=True)
            if sk=="Mp + ts":
                a,b,c=st.columns(3); mp=a.number_input("Mp máx (%)",0.01,99.0,10.0); ts=b.number_input("ts máx (s)",0.001,value=4.0); band=c.selectbox("Critério ts",[5,2])
            elif sk=="ξ + ωn":
                a,b=st.columns(2); zeta=a.number_input("ξ",0.001,0.999,0.7); wn=b.number_input("ωn",0.001,value=0.5)
            else:
                a,b=st.columns(2); pr=a.number_input("Re(s_d)",value=-4.0); pi=b.number_input("Im(s_d)",0.001,value=4.0)
            go=st.form_submit_button("Projetar",type="primary")
            if go:
                try:
                    if sk=="Mp + ts": specs=DesignSpecs.from_mp_ts(mp,ts,band/100)
                    elif sk=="ξ + ωn": specs=DesignSpecs.from_zeta_wn(zeta,wn)
                    else: specs=DesignSpecs.from_pole(pr,pi)
                    data={"controller":kind,"num_g":_parse(ng),"den_g":_parse(dg),"num_h":_parse(nh),"den_h":_parse(dh),"specs":lambda s=specs:s}
                    st.session_state["last_design"]=(_solve(data),data)
                except Exception as exc: st.error(f"Falha no projeto: {exc}")
    if "last_design" in st.session_state:
        design,data=st.session_state["last_design"]; st.markdown("---"); render_controller_design(design)
        pdf=build_controller_pdf(design,data["num_g"],data["den_g"],data["num_h"],data["den_h"])
        st.download_button("Baixar resolução — Modo Prova",pdf,f"projeto_{design.controller_type}_modo_prova.pdf","application/pdf",use_container_width=True)
