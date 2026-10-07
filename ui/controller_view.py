import re

import streamlit as st

from core.controllers import ControllerDesigner
from reports.full_pdf import build_controller_full_pdf
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
        "Escolha HTML para manter os gráficos interativos ou PDF para uma "
        "versão completa e portátil com os mesmos passos e gráficos estáticos."
    )
    c1,c2=st.columns([1,2])
    with c1:
        report_format=st.selectbox(
            "Formato",
            ["HTML (interativo)","PDF (completo)"],
            key="full_report_format",
        )
    with c2:
        default_name=f"resolucao_completa_{design.controller_type}"
        filename=st.text_input(
            "Nome do arquivo",
            value=default_name,
            key="full_report_filename",
            help="A extensão correta é adicionada automaticamente.",
        )

    clean_name=_safe_filename(filename,default_name)
    for ext in (".html",".pdf"):
        if clean_name.lower().endswith(ext):
            clean_name=clean_name[:-len(ext)]

    if report_format=="HTML (interativo)":
        st.caption(
            "HTML: tema claro fixo, equações formatadas e gráficos Plotly interativos."
        )
        full_report=build_controller_html(design,data)
        extension="html"
        mime="text/html"
    else:
        st.caption(
            "PDF: mantém toda a sequência da resolução e inclui os gráficos "
            "como imagens estáticas para impressão/arquivo."
        )
        full_report=build_controller_full_pdf(design,data)
        extension="pdf"
        mime="application/pdf"

    st.download_button(
        f"Baixar resolução completa em {extension.upper()}",
        data=full_report,
        file_name=f"{clean_name}.{extension}",
        mime=mime,
        use_container_width=True,
        type="primary",
    )
