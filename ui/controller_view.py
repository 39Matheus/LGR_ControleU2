import re

import streamlit as st

from core.controllers import ControllerDesigner
from reports.html import build_controller_html
from reports.pdf import build_controller_pdf
from .renderers import render_controller_design


def solve_controller(data):
    designer=ControllerDesigner(
        data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )
    return designer.design(
        data["controller"],
        data["specs"],
        auto_refine=bool(data.get("auto_refine",False)),
    )


def _safe_filename(name,default):
    name=(name or "").strip()
    if not name:
        name=default
    name=re.sub(r'[<>:"/\\|?*]+',"_",name)
    name=name.rstrip(". ")
    return name or default


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

    st.markdown("---")
    st.subheader("Baixar resolução completa")
    st.caption(
        "O relatório completo mantém o desenvolvimento detalhado e inclui os "
        "gráficos Plotly. Ele é salvo em HTML para preservar a apresentação e "
        "a interatividade dos gráficos."
    )
    default_name=f"resolucao_completa_{design.controller_type}"
    filename=st.text_input(
        "Nome do arquivo",
        value=default_name,
        key="full_report_filename",
        help="A extensão .html é adicionada automaticamente.",
    )
    clean_name=_safe_filename(filename,default_name)
    if clean_name.lower().endswith(".html"):
        clean_name=clean_name[:-5]
    full_report=build_controller_html(design,data)
    st.download_button(
        "Baixar resolução completa",
        data=full_report,
        file_name=f"{clean_name}.html",
        mime="text/html",
        use_container_width=True,
        type="primary",
    )
