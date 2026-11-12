import warnings

import streamlit as st

from ui.main_view import render_main_view
from ui.state_space_view import render_state_space_view
from ui.state_space_u3l2_view import render_state_space_u3l2_view

warnings.filterwarnings("ignore")

st.set_page_config(page_title="Calculadora de Controle — U2 e U3", layout="wide")
st.title("Calculadora — Projeto de Sistemas de Controle")
st.caption("Unidade 2: LGR e discretização · Unidade 3: Lista 1, espaço de estados · Lista 2, seguidores de referência.")

module = st.radio(
    "Módulo de estudo",
    ["Unidade 3 — Lista 2 (seguidores de referência)", "Unidade 3 — Lista 1 (espaço de estados)", "Unidade 2 — LGR e discretização"],
    horizontal=True,
    key="study_module",
)
if "Lista 2" in module:
    render_state_space_u3l2_view()
elif "Lista 1" in module:
    render_state_space_view()
else:
    render_main_view()