import warnings

import streamlit as st

from ui.main_view import render_main_view

warnings.filterwarnings("ignore")

st.set_page_config(page_title="Calculadora LGR — Unidade 2",layout="wide")
st.title("Calculadora LGR — Projeto de Sistemas de Controle")
st.caption(
    "Escopo: análise e projeto pelo LGR + aproximação discreta de funções "
    "de transferência contínuas (Forward/Euler, Backward, Tustin e invariância ao degrau)."
)

render_main_view()
