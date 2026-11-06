"""Interface U3L1: contas essenciais reproduzíveis, sem repetir as matrizes-base."""

import re

import streamlit as st
import sympy as sp

from core.state_space import (
    U3L1_SYSTEMS, analyze_system, controllable_canonical,
    place_observer, place_state_feedback,
)
from core.state_space_steps import build_resolution
from core.state_space_questions import QUESTION_TITLES, solve_full_question, individual_document
from reports.state_space_export import (
    build_state_space_pdf, build_state_space_html, build_all_state_space_zip,
)

ANALYZE = "Estabilidade, controlabilidade e observabilidade"
FEEDBACK = "Realimentação de estados — K"
OBSERVER = "Observador de estados — L"
CANONICAL = "Forma canônica controlável"

QUESTIONS = {
    "Q1 (Sistema 1) — itens a, b, c": ("Sistema 1 (Q1 e Q2)", ANALYZE),
    "Q1 (Sistema 2) — itens a, b, c": ("Sistema 2 (Q1 e Q3)", ANALYZE),
    "Q2a — Realimentação do Sistema 1": ("Sistema 1 (Q1 e Q2)", FEEDBACK),
    "Q2b — Observador do Sistema 1": ("Sistema 1 (Q1 e Q2)", OBSERVER),
    "Q3c — Realimentação do Sistema 2": ("Sistema 2 (Q1 e Q3)", FEEDBACK),
    "Q3d — Observador do Sistema 2": ("Sistema 2 (Q1 e Q3)", OBSERVER),
    "Q4a — Forma canônica controlável": ("Sistema da Q4 (G(s))", CANONICAL),
    "Q4b — Observador (−1, −1, −2)": ("Sistema da Q4 (G(s))", OBSERVER),
    "Q4c — Realimentação (−4, −4, −5)": ("Sistema da Q4 (G(s))", FEEDBACK),
}


def _matrix_text(value):
    return "\n".join("; ".join(str(v) for v in row) for row in value)


def _poles_text(values):
    return "; ".join(str(v).replace("I", "i") for v in values)


def _select_preset():
    name, operation = QUESTIONS[st.session_state["u3_question"]]
    p = U3L1_SYSTEMS[name]
    st.session_state["u3_operation"] = operation
    st.session_state["u3_representation"] = "G(s)" if "num" in p else "Matrizes A, B, C"
    if "num" in p:
        st.session_state["u3_num"] = "; ".join(map(str, p["num"]))
        st.session_state["u3_den"] = "; ".join(map(str, p["den"]))
    else:
        for key in ("A", "B", "C"):
            st.session_state["u3_" + key] = _matrix_text(p[key])
    desired = p["observer_poles"] if operation == OBSERVER else p["feedback_poles"]
    st.session_state["u3_poles"] = _poles_text(desired)
    st.session_state.pop("u3_result", None)


def _scalar(text):
    clean = text.strip().replace(",", ".").replace("−", "-")
    try:
        return sp.Rational(clean)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"Número inválido: {text!r}. Use inteiros, decimais ou frações (1/3).") from exc


def _coefficients(text):
    parts = [v.strip() for v in text.split(";") if v.strip()]
    if not parts:
        raise ValueError("Informe os coeficientes separados por ponto e vírgula.")
    return [_scalar(v) for v in parts]


def _matrix(text, label):
    rows = [row.strip() for row in text.splitlines() if row.strip()]
    if not rows:
        raise ValueError(f"Matriz {label}: campo vazio.")
    data = [[_scalar(x) for x in row.split(";")] for row in rows]
    if any(len(row) != len(data[0]) for row in data):
        raise ValueError(f"Matriz {label}: todas as linhas devem ter o mesmo número de colunas.")
    return sp.Matrix(data)


def _poles(text):
    parts = [x.strip() for x in text.replace("−", "-").split(";") if x.strip()]
    if not parts:
        raise ValueError("Informe os polos desejados separados por ';'.")
    result = []
    for part in parts:
        try:
            if not re.fullmatch(r"[\d\s\.\+\-\*/()ijIJ]+", part):
                raise ValueError("Caracteres não permitidos no polo.")
            expr = part.replace("j", "I").replace("i", "I")
            value = sp.nsimplify(sp.sympify(expr, locals={"I": sp.I}), rational=True)
            if value.free_symbols or not value.is_number:
                raise ValueError("Polo não numérico.")
            result.append(value)
        except (ValueError, TypeError, SyntaxError, sp.SympifyError) as exc:
            raise ValueError(f"Polo inválido: {part!r}. Exemplo: -2+i; -2-i; -10.") from exc
    return result


def _signature():
    """Evita mostrar uma solução anterior após alterações dos parâmetros."""
    names = ("u3_question", "u3_operation", "u3_representation", "u3_A", "u3_B",
             "u3_C", "u3_num", "u3_den", "u3_poles")
    return tuple(str(st.session_state.get(name, "")) for name in names)


def _show_input_preview():
    """A, B, C (ou G) aparecem uma única vez, antes da resolução."""
    with st.expander("Conferir dados interpretados", expanded=True):
        if st.session_state["u3_representation"] == "G(s)":
            num = _coefficients(st.session_state["u3_num"])
            den = _coefficients(st.session_state["u3_den"])
            numerator = sp.Poly.from_list(num, gens=sp.Symbol("s")).as_expr()
            denominator = sp.Poly.from_list(den, gens=sp.Symbol("s")).as_expr()
            if denominator == 0:
                raise ValueError("Denominador de G(s) não pode ser zero.")
            st.latex(r"G(s)=\frac{" + sp.latex(numerator) + "}{" + sp.latex(denominator) + "}")
        else:
            A, B, C = (_matrix(st.session_state["u3_"+key], key) for key in ("A", "B", "C"))
            # Confirmação sem repetições nas etapas de Ackermann ou no determinante.
            c1, c2, c3 = st.columns([2, 1, 2])
            c1.latex("A=" + sp.latex(A))
            c2.latex("B=" + sp.latex(B))
            c3.latex("C=" + sp.latex(C))
            if A.rows != A.cols or B.shape != (A.rows, 1) or C.shape != (1, A.rows):
                st.warning("Confira as dimensões: A deve ser n×n, B deve ser n×1 e C deve ser 1×n.")


def _render_steps(steps):
    details = []
    for step in steps:
        if step.kind == "heading":
            st.markdown("#### " + step.value)
        elif step.kind == "math":
            st.latex(step.value)
        elif step.kind == "success":
            st.success(step.value)
        elif step.kind == "warning":
            st.warning(step.value)
        elif step.kind == "detail_math":
            details.append(step.value)
        else:
            st.caption(step.value)
    if details:
        with st.expander("Contas complementares (opcionais)"):
            for expression in details:
                st.latex(expression)


def _operation_key(operation):
    return {ANALYZE: "analysis", FEEDBACK: "K", OBSERVER: "L", CANONICAL: "canonical"}[operation]



@st.cache_data(show_spinner=False)
def _complete_question(number):
    return solve_full_question(number)


@st.cache_data(show_spinner=False)
def _pdf_bytes(document):
    return build_state_space_pdf(document)


@st.cache_data(show_spinner=False)
def _html_bytes(document):
    return build_state_space_html(document)


def _render_complete_question_mode():
    st.caption(
        "Um preset por questão: Q1 reúne os dois sistemas; Q2 e Q3 reúnem "
        "realimentação e observador; Q4 reúne forma canônica, observador e realimentação. "
        "Os dados originais aparecem uma vez por sistema."
    )
    number = st.selectbox(
        "Selecione a questão completa",
        list(QUESTION_TITLES),
        format_func=lambda n: QUESTION_TITLES[n],
        key="u3_full_question",
    )
    try:
        document = _complete_question(number)
    except (ValueError, ArithmeticError, sp.SympifyError) as exc:
        st.error(f"Não foi possível resolver a questão completa: {exc}")
        return

    st.subheader("Resolução completa — fórmula, operações e conferência")
    st.caption("Os mesmos passos são preservados nos arquivos PDF e HTML.")
    col_pdf, col_html = st.columns(2)
    col_pdf.download_button(
        "Baixar PDF da questão",
        data=_pdf_bytes(document),
        file_name=f"U3L1_Questao_{number}.pdf",
        mime="application/pdf",
        key=f"u3_pdf_full_{number}",
        use_container_width=True,
    )
    col_html.download_button(
        "Baixar HTML da questão",
        data=_html_bytes(document),
        file_name=f"U3L1_Questao_{number}.html",
        mime="text/html",
        key=f"u3_html_full_{number}",
        use_container_width=True,
    )
    for section in document.sections:
        st.markdown("---")
        st.markdown("### " + section.title)
        _render_steps(section.steps)
    st.divider()
    with st.expander("Pacote para entrega — 4 PDFs + 4 HTMLs"):
        st.caption(
            "Gera oito arquivos com os mesmos cálculos apresentados, "
            "organizados por questão e reunidos em um único ZIP."
        )
        if st.button("Preparar pacote completo", key="u3_prepare_zip"):
            with st.spinner("Gerando os relatórios das quatro questões..."):
                st.session_state["u3_complete_zip"] = build_all_state_space_zip(_complete_question)
        if "u3_complete_zip" in st.session_state:
            st.download_button(
                "Baixar ZIP (4 PDFs + 4 HTMLs)",
                data=st.session_state["u3_complete_zip"],
                file_name="U3L1_Questoes_1_a_4_PDF_HTML.zip",
                mime="application/zip",
                use_container_width=True,
                key="u3_zip_download",
            )


def render_state_space_view():
    st.header("Unidade 3 · Lista 1 — Espaço de estados")
    st.caption("Resolução para reproduzir na prova: fórmula → substituição principal → resultado → conferência. Apostila: capítulo 7.")
    scope=st.radio(
        "Modo de resolução",
        ["Questão completa (recomendado)", "Item individual (editável)"],
        horizontal=True,
        key="u3_scope",
    )
    if scope.startswith("Questão completa"):
        _render_complete_question_mode()
        return
    if "u3_question" not in st.session_state:
        st.session_state["u3_question"] = next(iter(QUESTIONS))
        _select_preset()

    st.selectbox("Questão da lista", list(QUESTIONS), key="u3_question", on_change=_select_preset)
    st.radio("Operação", [ANALYZE, FEEDBACK, OBSERVER, CANONICAL], key="u3_operation")
    st.radio("Representação da entrada", ["Matrizes A, B, C", "G(s)"], key="u3_representation", horizontal=True)

    if st.session_state["u3_representation"] == "G(s)":
        st.caption("Coeficientes em ordem decrescente; incluir zeros em potências ausentes. Ex.: 3s²+4s−2 → 3;4;-2.")
        a, b = st.columns(2)
        a.text_input("Numerador de G(s)", key="u3_num")
        b.text_input("Denominador de G(s)", key="u3_den")
    else:
        st.caption("Uma linha por linha de texto, com valores separados por ';'. Aceita inteiros, decimais e frações.")
        st.text_area("A", key="u3_A", height=110)
        a, b = st.columns(2)
        a.text_area("B (coluna)", key="u3_B", height=90)
        b.text_area("C (linha)", key="u3_C", height=90)

    if st.session_state["u3_operation"] in (FEEDBACK, OBSERVER):
        st.text_input("Polos desejados (separados por ';')", key="u3_poles",
                      help="Ex.: -2+i; -2-i; -10. Repita os polos com multiplicidade.")
    try:
        _show_input_preview()
    except (ValueError, TypeError, sp.SympifyError) as exc:
        st.warning(f"Ainda não foi possível interpretar os dados: {exc}")

    if st.button("Resolver item da U3L1", type="primary", use_container_width=True):
        try:
            canonical = None
            if st.session_state["u3_representation"] == "G(s)":
                canonical = controllable_canonical(_coefficients(st.session_state["u3_num"]),
                                                   _coefficients(st.session_state["u3_den"]))
                A, B, C = canonical.A, canonical.B, canonical.C
            else:
                A, B, C = (_matrix(st.session_state["u3_"+key], key) for key in ("A", "B", "C"))
            analysis = analyze_system(A, B, C)
            operation = st.session_state["u3_operation"]
            placement = None
            if operation == FEEDBACK:
                placement = place_state_feedback(A, B, C, _poles(st.session_state["u3_poles"]))
            elif operation == OBSERVER:
                placement = place_observer(A, B, C, _poles(st.session_state["u3_poles"]))
            if operation == CANONICAL and canonical is None:
                raise ValueError("Para obter a forma canônica, selecione G(s) e informe seu numerador e denominador.")
            st.session_state["u3_result"] = (analysis, placement, canonical, operation, _signature())
        except (ValueError, ArithmeticError, sp.SympifyError) as exc:
            st.session_state.pop("u3_result", None)
            st.error(f"Não foi possível resolver: {exc}")

    result = st.session_state.get("u3_result")
    if result is None:
        return
    analysis, placement, canonical, operation, saved_signature = result
    if saved_signature != _signature():
        st.info("Os parâmetros foram alterados. Clique em 'Resolver item da U3L1' para atualizar as contas.")
        return

    mode = _operation_key(operation)
    # A entrada está confirmada no painel anterior; a resolução não a repete.
    steps = build_resolution(analysis, placement=placement, canonical=canonical,
                             operation=mode, include_input=False)
    st.divider()
    st.subheader("Resolução — contas principais")
    _render_steps(steps)

    if operation == CANONICAL:
        with st.expander("Opcional: conferir estabilidade, controlabilidade e observabilidade da realização"):
            extra = build_resolution(analysis, operation="analysis", include_input=False)
            _render_steps(extra)

    # O download inclui também os dados originais, uma única vez, e as contas da tela.
    export_steps = build_resolution(analysis, placement=placement, canonical=canonical,
                                    operation=mode, include_input=True)
    document=individual_document(st.session_state["u3_question"], export_steps)
    c_pdf,c_html=st.columns(2)
    c_pdf.download_button(
        "Baixar PDF do item", data=_pdf_bytes(document),
        file_name="U3L1_Item_Individual.pdf", mime="application/pdf",
        use_container_width=True, key="u3_individual_pdf",
    )
    c_html.download_button(
        "Baixar HTML do item", data=_html_bytes(document),
        file_name="U3L1_Item_Individual.html", mime="text/html",
        use_container_width=True, key="u3_individual_html",
    )