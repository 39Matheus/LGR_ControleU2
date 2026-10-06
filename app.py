import warnings
import streamlit as st
from ui.controller_view import render_controller_view
from ui.lgr_view import render_lgr_view

warnings.filterwarnings("ignore")
st.set_page_config(page_title="Calculadora LGR — Unidade 2",layout="wide")
st.title("Calculadora LGR — Projeto de Sistemas de Controle")
st.caption("Escopo: conteúdo até Projeto de Controladores pelo LGR. O capítulo 5 (aproximação discreta) não está incluído.")

mode=st.radio("Modo",["Projeto de controlador","Análise do LGR (12 passos)"],horizontal=True)
if mode=="Projeto de controlador":
    render_controller_view()
else:
    render_lgr_view()
