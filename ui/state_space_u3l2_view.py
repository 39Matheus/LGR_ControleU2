"""Tela U3L2: 4 questões completas, modo individual editável e PDF/HTML/ZIP."""
import re

import streamlit as st
import sympy as sp

from core.state_space import S, controllable_canonical
from core.state_space_steps import canonical_steps
from core.state_space_u3l2 import (U3L2_PRESETS, observable_canonical,
                                   design_step_servo)
from core.state_space_u3l2_questions import (
    QUESTION_TITLES, solve_full_question, individual_document, _input_steps,
    observable_canonical_steps, servo_steps,
)
from reports.state_space_export import (build_state_space_pdf, build_state_space_html,
                                        build_all_state_space_zip)


def _render_steps(steps):
    extra=[]
    for s in steps:
        if s.kind=="heading":
            st.markdown("#### "+s.value)
        elif s.kind=="math":
            st.latex(s.value)
        elif s.kind=="success":
            st.success(s.value)
        elif s.kind=="warning":
            st.warning(s.value)
        elif s.kind=="detail_math":
            extra.append(s.value)
        else:
            st.caption(s.value)
    if extra:
        with st.expander("Contas complementares"):
            for value in extra:
                st.latex(value)


@st.cache_data(show_spinner=False)
def _complete_question(number):
    return solve_full_question(number)


@st.cache_data(show_spinner=False)
def _pdf(document):
    return build_state_space_pdf(document)


@st.cache_data(show_spinner=False)
def _html(document):
    return build_state_space_html(document)


@st.cache_data(show_spinner=False)
def _all_zip():
    return build_all_state_space_zip(solve_full_question)


def _complete_mode():
    st.caption("Uma questão inteira por preset: itens a) e b) aparecem juntos; na questão 2 há apenas o projeto pedido.")
    number=st.selectbox("Questão completa",list(QUESTION_TITLES),
                        format_func=lambda n:QUESTION_TITLES[n],key="u3l2_full_question")
    try:
        doc=_complete_question(number)
        pdf=_pdf(doc)
        html=_html(doc)
    except (ValueError, ArithmeticError, sp.SympifyError) as exc:
        st.error(f"Falha ao preparar a questão completa: {exc}")
        return
    a,b=st.columns(2)
    a.download_button("Baixar PDF da questão",data=pdf,mime="application/pdf",
                      file_name=f"U3L2_Questao_{number}.pdf",use_container_width=True,
                      key="u3l2_full_pdf")
    b.download_button("Baixar HTML da questão",data=html,mime="text/html",
                      file_name=f"U3L2_Questao_{number}.html",use_container_width=True,
                      key="u3l2_full_html")
    for sec in doc.sections:
        st.divider()
        st.markdown("### "+sec.title)
        _render_steps(sec.steps)
    st.divider()
    with st.expander("Pacote completo para entrega — 4 PDFs + 4 HTMLs"):
        st.caption("Os oito arquivos são separados por questão e reunidos em um único ZIP.")
        if st.button("Preparar pacote U3L2",key="u3l2_prepare_zip"):
            with st.spinner("Montando as quatro questões..."):
                st.session_state["u3l2_zip"]= _all_zip()
        if "u3l2_zip" in st.session_state:
            st.download_button("Baixar ZIP completo (8 arquivos)",
                               data=st.session_state["u3l2_zip"],
                               file_name="U3L2_Questoes_1_a_4_PDF_HTML.zip",
                               mime="application/zip",key="u3l2_zip_download",
                               use_container_width=True)


def _mat_txt(rows):
    return "\n".join("; ".join(map(str,row)) for row in rows)


def _load_individual_preset():
    number=st.session_state["u3l2_individual_question"]
    p=U3L2_PRESETS[number]
    st.session_state["u3l2_input_type"]=p["representation"]
    if p["representation"]=="G(s)":
        st.session_state["u3l2_num"]="; ".join(map(str,p["num"]))
        st.session_state["u3l2_den"]="; ".join(map(str,p["den"]))
        st.session_state["u3l2_canonical"] = "Observável" if p["canonical"] == "observavel" else "Controlável"
    else:
        for label in ("A", "B", "C"):
            st.session_state["u3l2_" + label] = _mat_txt(p[label])
    st.session_state["u3l2_poles"] = "; ".join(str(x).replace("I", "i") for x in p["poles"])
    st.session_state.pop("u3l2_individual_result", None)


def _scalar(value):
    return sp.Rational(value.strip().replace(",", ".").replace("−", "-"))


def _coefficients(text):
    vals = [_scalar(x) for x in text.split(";") if x.strip()]
    if not vals:
        raise ValueError("Forneça os coeficientes separados por ';'.")
    return vals


def _matrix(text):
    lines = [x for x in text.splitlines() if x.strip()]
    vals = [[_scalar(x) for x in line.split(";")] for line in lines]
    if not vals or any(len(row) != len(vals[0]) for row in vals):
        raise ValueError("As linhas das matrizes devem conter a mesma quantidade de colunas.")
    return sp.Matrix(vals)


def _poles(text):
    parts = [x.strip().replace("−", "-") for x in text.split(";") if x.strip()]
    if not parts:
        raise ValueError("Informe os polos desejados.")
    result = []
    for part in parts:
        if not re.fullmatch(r"[\d\s\.,\+\-\*/()ijIJ]+", part):
            raise ValueError(f"Polo inválido: {part!r}")
        expr = part.replace(",", ".").replace("i", "I").replace("j", "I")
        result.append(sp.nsimplify(sp.sympify(expr, locals={"I": sp.I}), rational=True))
    return result


def _signature():
    return tuple(str(st.session_state.get(x, "")) for x in (
        "u3l2_individual_question", "u3l2_item", "u3l2_input_type", "u3l2_canonical",
        "u3l2_num", "u3l2_den", "u3l2_A", "u3l2_B", "u3l2_C", "u3l2_poles",
    ))


def _individual_mode():
    if "u3l2_individual_question" not in st.session_state:
        st.session_state["u3l2_individual_question"] = 1
        _load_individual_preset()
    number = st.selectbox("Questão original", list(QUESTION_TITLES),
                          format_func=lambda n: QUESTION_TITLES[n],
                          on_change=_load_individual_preset, key="u3l2_individual_question")
    if number == 2:
        item = "Projeto do seguidor de degrau"
    else:
        item = st.radio("Item a resolver", ["a) Forma canônica", "b) Seguidor de degrau"],
                        horizontal=True, key="u3l2_item")
    st.radio("Representação de entrada", ["G(s)", "Matrizes A, B, C"],
             horizontal=True, key="u3l2_input_type")
    if st.session_state["u3l2_input_type"] == "G(s)":
        a, b = st.columns(2)
        a.text_input("Numerador (coeficientes decrescentes, separados por ;)", key="u3l2_num")
        b.text_input("Denominador (incluir zeros ausentes)", key="u3l2_den")
        st.radio("Realização da FT", ["Controlável", "Observável"], horizontal=True, key="u3l2_canonical")
        try:
            num = _coefficients(st.session_state["u3l2_num"])
            den = _coefficients(st.session_state["u3l2_den"])
            numerator = sp.Poly.from_list(num, gens=S).as_expr()
            denominator = sp.Poly.from_list(den, gens=S).as_expr()
            st.latex("G(s)=" + sp.latex(numerator/denominator))
        except (ValueError, TypeError, sp.SympifyError) as exc:
            st.warning(f"Entrada parcial: {exc}")
    else:
        st.caption("Uma linha por linha, coeficientes separados por ponto e vírgula.")
        st.text_area("Matriz A", key="u3l2_A", height=105)
        a, b = st.columns(2)
        a.text_area("Vetor B", key="u3l2_B", height=95)
        b.text_area("Linha C", key="u3l2_C", height=95)
        try:
            A, B, C = [_matrix(st.session_state.get("u3l2_"+v, "")) for v in ("A", "B", "C")]
            x, y, z = st.columns([2, 1, 2])
            x.latex("A=" + sp.latex(A))
            y.latex("B=" + sp.latex(B))
            z.latex("C=" + sp.latex(C))
        except (ValueError, TypeError, sp.SympifyError):
            pass
    if number == 2 or item.startswith("b)"):
        st.text_input("Polos desejados, incluindo multiplicidades (;)", key="u3l2_poles",
                      help="Ex.: -1+i/10; -1-i/10; -5; -12")
    if st.button("Resolver item U3L2", type="primary", use_container_width=True):
        try:
            canonical = None
            if st.session_state["u3l2_input_type"] == "G(s)":
                function = (observable_canonical if st.session_state["u3l2_canonical"] == "Observável"
                            else controllable_canonical)
                canonical = function(_coefficients(st.session_state["u3l2_num"]),
                                     _coefficients(st.session_state["u3l2_den"]))
                A, B, C = canonical.A, canonical.B, canonical.C
            else:
                A, B, C = [_matrix(st.session_state["u3l2_"+v]) for v in ("A", "B", "C")]
            if number != 2 and item.startswith("a)"):
                if canonical is None:
                    raise ValueError("Para obter a forma canônica, informe uma função G(s).")
                steps = _input_steps(number, canonical, None)
                if st.session_state["u3l2_canonical"] == "Observável":
                    steps += observable_canonical_steps(canonical)
                else:
                    steps += canonical_steps(canonical)
            else:
                servo = design_step_servo(A, B, C, _poles(st.session_state["u3l2_poles"]))
                steps = [*_input_steps(number, canonical, servo), *servo_steps(servo)]
            document = individual_document(f"U3L2 — Questão {number} — {item}", steps)
            st.session_state["u3l2_individual_result"] = (document, _signature())
        except (ValueError, TypeError, ArithmeticError, sp.SympifyError) as exc:
            st.session_state.pop("u3l2_individual_result", None)
            st.error(f"Falha na resolução: {exc}")
    saved = st.session_state.get("u3l2_individual_result")
    if not saved:
        return
    document, signature = saved
    if signature != _signature():
        st.info("Os dados foram alterados. Execute novamente para atualizar os cálculos.")
        return
    st.divider()
    st.subheader("Resolução individual")
    for section in document.sections:
        _render_steps(section.steps)
    a, b = st.columns(2)
    a.download_button("Baixar PDF do item", data=_pdf(document),
                      file_name="U3L2_Item_Individual.pdf", mime="application/pdf",use_container_width=True)
    b.download_button("Baixar HTML do item", data=_html(document),
                      file_name="U3L2_Item_Individual.html", mime="text/html",use_container_width=True)


def render_state_space_u3l2_view():
    st.header("Unidade 3 · Lista 2 — Seguidores de referência")
    st.caption("Apostila §7.10: e=y-r, sistema aumentado, Ackermann, lei integral, polos e erro nulo em degrau.")
    mode = st.radio("Modo de resolução", ["Questão completa (recomendado)", "Item individual (editável)"],
                    key="u3l2_mode", horizontal=True)
    if mode.startswith("Questão completa"):
        _complete_mode()
    else:
        _individual_mode()