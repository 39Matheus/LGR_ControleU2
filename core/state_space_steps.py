"""Resoluções da U3L1: fórmulas + operações intermediárias + conclusão.

Esta camada é pura (sem Streamlit) e serve tanto à tela quanto à exportação.
As matrizes de entrada são exibidas apenas uma vez na resolução.
"""
from dataclasses import dataclass

import sympy as sp

from .state_space import S


@dataclass(frozen=True)
class ResolutionStep:
    kind: str  # heading | math | note | success | warning | detail_math
    value: str


def _heading(value):
    return ResolutionStep("heading", value)


def _math(value):
    return ResolutionStep("math", value)


def _note(value):
    return ResolutionStep("note", value)


def _good(value):
    return ResolutionStep("success", value)


def _warn(value):
    return ResolutionStep("warning", value)


def _detail(value):
    return ResolutionStep("detail_math", value)


def _tex(value):
    return sp.latex(value)


def _eq(name, value):
    return _math(rf"{name}={_tex(value)}")


def _poly_at_matrix(poly, n):
    coeffs = sp.Poly(poly, S).all_coeffs()
    terms = []
    for i, coefficient in enumerate(coeffs):
        power = n - i
        if coefficient == 0:
            continue
        base = "I" if power == 0 else "A" if power == 1 else rf"A^{{{power}}}"
        if coefficient == 1:
            terms.append(base)
        elif coefficient == -1:
            terms.append("-" + base)
        else:
            terms.append(_tex(coefficient) + base)
    expression = "+".join(terms).replace("+-", "-")
    return rf"q(A)={expression}"


def _eigen_tex(values):
    results = []
    for z in values:
        r, im = z.real, z.imag
        if abs(im) <= 1e-9:
            results.append(f"{r:.5g}")
        else:
            sign = "+" if im >= 0 else "-"
            results.append(rf"{r:.5g}{sign}{abs(im):.5g}i")
    return r"\left\{ " + r",\; ".join(results) + r" \right\}"


def _cofactor_expansion(M):
    """Primeira linha: usa menores 2×2 já avaliados, sem duplicar A."""
    if M.rows != 3:
        return None
    pieces = []
    for j in range(3):
        coeff = (-1)**j * M[0, j]
        if coeff == 0:
            continue
        minor = sp.factor(M.minor_submatrix(0, j).det())
        pieces.append(sp.factor(coeff*minor))
    if not pieces:
        return "0"
    return sp.latex(sp.Add(*pieces, evaluate=False))


def _rank_steps(symbol, matrix, rank, n, property_name):
    det = sp.factor(matrix.det())  # SISO quadrado (n×n)
    steps = [
        _math(rf"\det({symbol})={_tex(det)}"),
        _math(rf"\operatorname{{posto}}({symbol})={rank}{'=' if rank == n else '<'}n={n}"),
    ]
    if rank == n:
        steps.append(_good(f"O sistema é {property_name}: posto completo ({rank})."))
    else:
        steps.append(_warn(f"O sistema não é {property_name}: posto {rank} inferior à ordem {n}."))
        steps.append(_detail(rf"\operatorname{{rref}}({symbol})={_tex(matrix.rref()[0])}"))
    return steps


def _input_steps(analysis, canonical):
    if canonical is not None:
        return [_heading("Dados da função de transferência"), _eq("G(s)", canonical.original_transfer)]
    return [
        _heading("Conferência dos dados de entrada"),
        _math(rf"A={_tex(analysis.A)},\quad B={_tex(analysis.B)},\quad C={_tex(analysis.C)}"),
    ]


def canonical_steps(realization):
    """Mostra a leitura de α e β e a realização (7.12) uma única vez."""
    n = realization.A.rows
    a = [-realization.A[n-1, n-i] for i in range(1,n+1)]
    b = [realization.C[0,n-i] for i in range(1,n+1)]
    den_poly = S**n + sum(a[i]*S**(n-1-i) for i in range(n))
    num_poly = sum(b[i]*S**(n-1-i) for i in range(n))
    alpha_tex = ",\\;".join(rf"\alpha_{{{i+1}}}={_tex(v)}" for i,v in enumerate(a))
    beta_tex = ",\\;".join(rf"\beta_{{{i+1}}}={_tex(v)}" for i,v in enumerate(b))
    result = [
        _heading("Forma canônica controlável — equação (7.12)"),
        _note("Ler os coeficientes do denominador mônico e do numerador em ordem decrescente."),
        _math(rf"D(s)={_tex(sp.expand(den_poly))}\quad\Longrightarrow\quad {alpha_tex}"),
        _math(rf"N(s)={_tex(sp.expand(num_poly))}\quad\Longrightarrow\quad {beta_tex}"),
        _note("Preencher a última linha de A com os α negativos; B termina em 1; C recebe os β em ordem inversa."),
        _math(rf"A_c={_tex(realization.A)},\quad B_c={_tex(realization.B)},\quad C_c={_tex(realization.C)}"),
        _eq("d", realization.D),
        _math(rf"C_c(sI-A_c)^{{-1}}B_c+d={_tex(realization.reconstructed_transfer)}=G(s)"),
        _good("A realização reproduz a função de transferência dada."),
    ]
    return result


def analysis_steps(analysis):
    n=analysis.A.rows
    M=S*sp.eye(n)-analysis.A
    steps = [
        _heading("1. Polinômio característico e estabilidade — seção 7.3"),
        _math(rf"\Delta(s)=\det(sI-A)=\det\left({_tex(M)}\right)"),
    ]
    cofactor = _cofactor_expansion(M)
    if cofactor is not None:
        steps.append(_math(rf"\Delta(s)={cofactor}"))
    steps += [
        _eq(r"\Delta(s)", sp.expand(analysis.characteristic)),
    ]
    factored=sp.factor(analysis.characteristic)
    if sp.expand(factored-analysis.characteristic)==0 and factored!=sp.expand(analysis.characteristic):
        steps.append(_math(rf"\Delta(s)=0\ \Longleftrightarrow\ {_tex(factored)}=0"))
    steps += [
        _math(rf"\Delta(s)=0\ \Rightarrow\ \lambda(A)\approx{_eigen_tex(analysis.eigenvalues)}"),
        _good("Estável: todas as partes reais dos polos são negativas.") if analysis.stable
        else _warn("Não é assintoticamente estável: há polo com parte real não negativa."),
        _heading("2. Controlabilidade — seção 7.4"),
    ]
    if n>1:
        product_tex = r",\quad ".join(rf"A^{{{k}}}B={_tex(analysis.A**k*analysis.B)}" if k>1 else rf"AB={_tex(analysis.A*analysis.B)}" for k in range(1,n))
        steps.append(_math(product_tex))
    steps.append(_math(rf"U=[B\quad AB\quad\cdots\quad A^{{{n-1}}}B]={_tex(analysis.U)}"))
    steps.extend(_rank_steps("U",analysis.U,analysis.rank_U,n,"controlável"))
    steps.append(_heading("3. Observabilidade — seção 7.5"))
    if n>1:
        product_tex = r",\quad ".join(rf"CA^{{{k}}}={_tex(analysis.C*analysis.A**k)}" if k>1 else rf"CA={_tex(analysis.C*analysis.A)}" for k in range(1,n))
        steps.append(_math(product_tex))
    steps.append(_math(rf"V=\begin{{bmatrix}}C\\CA\\\vdots\\CA^{{{n-1}}}\end{{bmatrix}}={_tex(analysis.V)}"))
    steps.extend(_rank_steps("V",analysis.V,analysis.rank_V,n,"observável"))
    return steps


def placement_steps(analysis, placement):
    if placement is None:
        return []
    is_feedback=placement.mode=="K"
    n=analysis.A.rows
    target_matrix = analysis.U if is_feedback else analysis.V
    symbol = "U" if is_feedback else "V"
    rank = analysis.rank_U if is_feedback else analysis.rank_V
    goal = "realimentação de estados" if is_feedback else "observador de estados"
    method_section = "7.7 — equação (7.15)" if is_feedback else "7.8 — equação (7.18)"
    factor=sp.factor(placement.desired)
    terms=_poly_at_matrix(placement.desired, n)
    output=[
        _heading(f"1. Verificar a possibilidade de projetar {goal}"),
        _math(rf"{symbol}={_tex(target_matrix)}"),
        _math(rf"\det({symbol})={_tex(target_matrix.det())}\neq0\quad\Longrightarrow\quad\operatorname{{posto}}({symbol})={rank}=n={n}"),
        _note("Posto completo: a alocação de todos os polos é possível."),
        _heading("2. Formar o polinômio desejado"),
        _math(rf"\Delta_d(s)={_tex(factor)}={_tex(sp.expand(placement.desired))}"),
        _heading("3. Calcular o polinômio matricial q(A)"),
        _math(terms),
        _eq("q(A)",placement.qA),
        _heading(f"4. Aplicar Ackermann — apostila, seção {method_section}"),
    ]
    e=sp.zeros(n,1);e[n-1]=1
    if is_feedback:
        row=e.T
        reduced=sp.simplify(row*placement.inverse_matrix)
        output.extend([
            _math(r"K=-e_n^TU^{-1}q(A),\qquad e_n^T=" +_tex(row)),
            _math(rf"U^{{-1}}={_tex(placement.inverse_matrix)}"),
            _math(rf"e_n^TU^{{-1}}={_tex(reduced)}"),
            _math(rf"K=-\left(e_n^TU^{{-1}}\right)q(A)={_tex(placement.gain)}"),
            _note("Convenção da apostila: u(t)=Kx(t)+r(t), então A_mf=A+BK."),
            _heading("5. Conferir o sistema realimentado"),
            _math(rf"A_{{mf}}=A+BK={_tex(placement.result_A)}"),
            _math(rf"\det(sI-A_{{mf}})={_tex(sp.expand(placement.result_characteristic))}=\Delta_d(s)"),
        ])
    else:
        reduced=sp.simplify(placement.inverse_matrix*e)
        output.extend([
            _math(r"L=q(A)V^{-1}e_n,\qquad e_n=" +_tex(e)),
            _math(rf"V^{{-1}}={_tex(placement.inverse_matrix)}"),
            _math(rf"V^{{-1}}e_n={_tex(reduced)}"),
            _math(rf"L=q(A)\left(V^{{-1}}e_n\right)={_tex(placement.gain)}"),
            _note("A dinâmica do erro do observador é A-LC."),
            _heading("5. Conferir os polos do observador"),
            _math(rf"A_o=A-LC={_tex(placement.result_A)}"),
            _math(rf"\det(sI-A_o)={_tex(sp.expand(placement.result_characteristic))}=\Delta_d(s)"),
        ])
    output.append(_math(rf"\lambda_{{obtidos}}={_eigen_tex(placement.result_eigenvalues)}"))
    output.append(_good("Conferido: os polos obtidos são os solicitados."))
    return output


def build_resolution(analysis, placement=None, canonical=None, operation="analysis", include_input=True):
    """Sequência didática para exibição e download, sem repetir matrizes base."""
    steps=_input_steps(analysis,canonical) if include_input else []
    if canonical is not None:
        steps.extend(canonical_steps(canonical))
    if operation=="analysis":
        steps.extend(analysis_steps(analysis))
    elif operation in ("K","L"):
        steps.extend(placement_steps(analysis,placement))
    elif operation=="canonical":
        steps.append(_note("A realização pode ser analisada também quanto à controlabilidade e observabilidade."))
    else:
        raise ValueError(f"Operação desconhecida: {operation}")
    return steps


def resolution_markdown(steps):
    """Exportação com os mesmos passos visíveis na tela."""
    output=["# Resolução — Unidade 3, Lista 1", "Metodologia: apostila, capítulo 7. Convenção: u=Kx; A_mf=A+BK e A_o=A-LC."]
    for step in steps:
        if step.kind=="heading":
            output.append("\n## "+step.value)
        elif step.kind in ("math","detail_math"):
            output.append("\n$$\n"+step.value+"\n$$")
        else:
            output.append("\n"+step.value)
    return "\n\n".join(output)+"\n"