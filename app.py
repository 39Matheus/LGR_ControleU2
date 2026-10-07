import warnings

import streamlit as st

from ui.main_view import render_main_view

warnings.filterwarnings("ignore")

st.set_page_config(page_title="Calculadora LGR — Unidade 2",layout="wide")
st.title("Calculadora LGR — Projeto de Sistemas de Controle")
st.caption(
    "Escopo: conteúdo até Projeto de Controladores pelo LGR. "
    "O capítulo 5 (aproximação discreta) não está incluído."
)

render_main_view()
