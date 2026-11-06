"""Interface Streamlit dedicada à Lista U3.1, sem alterar a calculadora da U2."""

import re

import streamlit as st
import sympy as sp

from core.state_space import (
    U3L1_SYSTEMS, analyze_system, controllable_canonical,
    place_observer, place_state_feedback,
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
    name, op = QUESTIONS[st.session_state["u3_question"]]
    preset = U3L1_SYSTEMS[name]
    st.session_state["u3_operation"] = op
    st.session_state["u3_representation"] = "G(s)" if "num" in preset else "Matrizes A, B, C"
    if "num" in preset:
        st.session_state["u3_num"] = "; ".join(map(str, preset["num"]))
        st.session_state["u3_den"] = "; ".join(map(str, preset["den"]))
    else:
        for key in ("A", "B", "C"):
            st.session_state["u3_" + key] = _matrix_text(preset[key])
    poles = preset["observer_poles"] if op == OBSERVER else preset["feedback_poles"]
    st.session_state["u3_poles"] = _poles_text(poles)


def _scalar(text):
    clean = text.strip().replace(",", ".").replace("−", "-")
    try:
        return sp.Rational(clean)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"Número inválido: {text!r}. Use inteiros, decimais ou frações (1/3).") from exc


def _coefficients(text):
    values = [v.strip() for v in text.split(";") if v.strip()]
    if not values:
        raise ValueError("Informe os coeficientes separados por ponto e vírgula.")
    return [_scalar(v) for v in values]


def _matrix(text, label):
    rows = [row.strip() for row in text.splitlines() if row.strip()]
    if not rows:
        raise ValueError(f"Matriz {label}: campo vazio.")
    values = [[_scalar(x) for x in row.split(";")] for row in rows]
    if any(len(row) != len(values[0]) for row in values):
        raise ValueError(f"Matriz {label}: todas as linhas devem ter o mesmo número de colunas.")
    return sp.Matrix(values)


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


def _eigen_text(eigenvalues):
    values = []
    for z in eigenvalues:
        if abs(z.imag) < 1e-8:
            values.append(f"{z.real:.6g}")
        else:
            values.append(f"{z.real:.6g}{z.imag:+.6g}i")
    return "; ".join(values)


def _line(name, expr):
    st.latex(rf"{name}={sp.latex(expr)}")


def _show_analysis(result):
    st.markdown("#### 1. Polinômio característico e estabilidade")
    st.latex(r"\Delta(s)=\det(sI-A)=" + sp.latex(result.characteristic))
    st.write("Autovalores de A:", _eigen_text(result.eigenvalues))
    if result.stable:
        st.success("Estável: todos os autovalores possuem parte real negativa.")
    else:
        st.warning("Não é assintoticamente estável: existe autovalor com parte real não negativa.")

    st.markdown("#### 2. Matriz de controlabilidade")
    st.latex(r"U=[B\quad AB\quad\ldots\quad A^{n-1}B]=" + sp.latex(result.U))
    st.write(f"posto(U) = {result.rank_U} de n = {result.A.rows}; det(U) = {result.U.det()}.")
    if result.controllable:
        st.success("Sistema controlável.")
    else:
        st.warning("Sistema não controlável: não é possível posicionar arbitrariamente todos os polos por realimentação.")

    st.markdown("#### 3. Matriz de observabilidade")
    st.latex(r"V=\begin{bmatrix}C\\CA\\\vdots\\CA^{n-1}\end{bmatrix}=" + sp.latex(result.V))
    st.write(f"posto(V) = {result.rank_V} de n = {result.A.rows}; det(V) = {result.V.det()}.")
    if result.observable:
        st.success("Sistema observável.")
    else:
        st.warning("Sistema não observável: não é possível posicionar arbitrariamente todos os polos do observador.")


def _show_placement(result, analysis):
    st.markdown("#### 1. Verificação da condição de existência")
    is_feedback = result.mode == "K"
    if is_feedback:
        st.write(f"posto(U) = {analysis.rank_U} = n: sistema controlável.")
    else:
        st.write(f"posto(V) = {analysis.rank_V} = n: sistema observável.")
    st.markdown("#### 2. Polinômio desejado e q(A)")
    st.latex(r"\Delta_d(s)=\prod_{j=1}^{n}(s-\lambda_j)=" + sp.latex(sp.factor(result.desired)))
    _line("q(A)", result.qA)
    st.markdown("#### 3. Fórmula de Ackermann (convenção da apostila)")
    if is_feedback:
        st.latex(r"K=-\begin{bmatrix}0&\cdots&0&1\end{bmatrix}U^{-1}q(A)")
        _line(r"U^{-1}", result.inverse_matrix)
        _line("K", result.gain)
        st.markdown("#### 4. Conferência de A + BK")
        _line("A+BK", result.result_A)
    else:
        st.latex(r"L=q(A)V^{-1}\begin{bmatrix}0&\cdots&0&1\end{bmatrix}^{T}")
        _line(r"V^{-1}", result.inverse_matrix)
        _line("L", result.gain)
        st.markdown("#### 4. Conferência de A − LC")
        _line("A-LC", result.result_A)
    st.latex(r"\det(sI-A_{\mathrm{final}})=" + sp.latex(result.result_characteristic))
    st.write("Polos efetivamente obtidos:", _eigen_text(result.result_eigenvalues))
    st.success("Polinômio e polos conferidos exatamente com os valores desejados.")


def _report(analysis, placement=None, canonical=None):
    lines = ["# Resolução — Lista 1 da Unidade 3", "",
             "Metodologia: apostila, seções 7.3 a 7.8. Convenção: u=Kx, A_mf=A+BK; erro do observador: A-LC.", ""]
    if canonical is not None:
        lines.extend(["## Forma canônica controlável",
                      f"G(s) = {canonical.original_transfer}",
                      f"A = {canonical.A}; B = {canonical.B}; C = {canonical.C}; D = {canonical.D}",
                      f"G_reconstruída(s) = {canonical.reconstructed_transfer}", ""])
    lines.extend(["## Verificação do sistema", f"A = {analysis.A}; B = {analysis.B}; C = {analysis.C}",
                  f"Δ(s) = {analysis.characteristic}",
                  f"Autovalores: {_eigen_text(analysis.eigenvalues)}",
                  f"Estabilidade assintótica: {'sim' if analysis.stable else 'não'}",
                  f"U = {analysis.U}; posto(U) = {analysis.rank_U}; controlável: {analysis.controllable}",
                  f"V = {analysis.V}; posto(V) = {analysis.rank_V}; observável: {analysis.observable}", ""])
    if placement is not None:
        lines.extend(["## Projeto por Ackermann", f"Δ_d(s) = {placement.desired}",
                      f"q(A) = {placement.qA}", f"Matriz inversa usada = {placement.inverse_matrix}",
                      f"{placement.mode} = {placement.gain}",
                      f"Matriz dinâmica final = {placement.result_A}",
                      f"Δ_final(s) = {placement.result_characteristic}",
                      f"Polos finais: {_eigen_text(placement.result_eigenvalues)}"])
    return "\n\n".join(lines)


def render_state_space_view():
    st.header("Unidade 3 · Lista 1 — Espaço de estados")
    st.caption("Questões 1 a 4; conforme as fórmulas de estabilidade, U, V, Ackermann e forma canônica do capítulo 7 da apostila.")
    if "u3_question" not in st.session_state:
        st.session_state["u3_question"] = next(iter(QUESTIONS))
        _select_preset()

    st.selectbox("Selecione o item da lista", list(QUESTIONS), key="u3_question", on_change=_select_preset)
    st.radio("Operação", [ANALYZE, FEEDBACK, OBSERVER, CANONICAL], key="u3_operation")
    st.radio("Representação de entrada", ["Matrizes A, B, C", "G(s)"], horizontal=True, key="u3_representation")

    if st.session_state["u3_representation"] == "G(s)":
        st.caption("Coeficientes em ordem decrescente, incluindo zeros ausentes. Ex.: 3s²+4s−2 → 3;4;-2.")
        a, b = st.columns(2)
        a.text_input("Numerador de G(s)", key="u3_num")
        b.text_input("Denominador de G(s)", key="u3_den")
    else:
        st.caption("Uma linha da matriz por linha de texto; elementos separados por ';'. Aceita inteiros, decimais e frações.")
        st.text_area("A", key="u3_A", height=110)
        a, b = st.columns(2)
        a.text_area("B (vetor coluna)", key="u3_B", height=90)
        b.text_area("C (vetor linha)", key="u3_C", height=90)

    if st.session_state["u3_operation"] in (FEEDBACK, OBSERVER):
        st.text_input("Polos desejados (separados por ';')", key="u3_poles",
                      help="Ex.: -2+i; -2-i; -10. Escreva os polos repetidos duas vezes.")
    if st.button("Resolver item da U3L1", type="primary", use_container_width=True):
        try:
            canonical = None
            if st.session_state["u3_representation"] == "G(s)":
                canonical = controllable_canonical(_coefficients(st.session_state["u3_num"]),
                                                   _coefficients(st.session_state["u3_den"]))
                A, B, C = canonical.A, canonical.B, canonical.C
            else:
                A, B, C = (_matrix(st.session_state["u3_" + key], key) for key in ("A", "B", "C"))
            analysis = analyze_system(A, B, C)
            operation = st.session_state["u3_operation"]
            placement = None
            if operation == FEEDBACK:
                placement = place_state_feedback(A, B, C, _poles(st.session_state["u3_poles"]))
            elif operation == OBSERVER:
                placement = place_observer(A, B, C, _poles(st.session_state["u3_poles"]))
            if operation == CANONICAL and canonical is None:
                raise ValueError("Para obter a forma canônica, escolha a representação G(s) e informe seu numerador/denominador.")
            st.session_state["u3_result"] = (analysis, placement, canonical, operation)
        except (ValueError, ArithmeticError, sp.SympifyError) as exc:
            st.session_state.pop("u3_result", None)
            st.error(f"Não foi possível resolver: {exc}")

    if "u3_result" not in st.session_state:
        return
    analysis, placement, canonical, operation = st.session_state["u3_result"]
    st.divider()
    st.subheader("Resolução passo a passo")
    if canonical is not None:
        st.markdown("#### Forma canônica controlável — equação (7.12)")
        _line("G(s)", canonical.original_transfer)
        _line("A_c", canonical.A)
        _line("B_c", canonical.B)
        _line("C_c", canonical.C)
        _line("D", canonical.D)
        st.latex(r"G(s)=C_c(sI-A_c)^{-1}B_c+D="+sp.latex(canonical.reconstructed_transfer))
        st.success("Realização conferida pela função de transferência reconstruída.")
    if operation == ANALYZE:
        _show_analysis(analysis)
    elif operation in (FEEDBACK, OBSERVER):
        _show_placement(placement, analysis)
    elif operation == CANONICAL:
        with st.expander("Verificar também controlabilidade e observabilidade da realização"):
            _show_analysis(analysis)
    st.download_button("Baixar resolução (.md)", data=_report(analysis, placement, canonical),
                       file_name="resolucao_U3L1.md", mime="text/markdown")