import streamlit as st

from core.controllers import ControllerDesigner
from reports.pdf import build_controller_pdf
from .renderers import render_controller_design


def solve_controller(data):
    designer=ControllerDesigner(
        data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )
    return designer.design(data["controller"],data["specs"])


def render_controller_result(design,data):
    render_controller_design(design,data)
    pdf=build_controller_pdf(
        design,data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )
    st.download_button(
        "Baixar resolução — Modo Prova",
        pdf,
        f"projeto_{design.controller_type}_modo_prova.pdf",
        "application/pdf",
        use_container_width=True,
    )
