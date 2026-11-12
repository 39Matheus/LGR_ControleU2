"""Quatro respostas integrais da U3L2, com passos principais iguais na tela/PDF/HTML."""
import sympy as sp

from .state_space import S
from .state_space_steps import ResolutionStep, canonical_steps
from .state_space_questions import QuestionDocument, QuestionSection
from .state_space_u3l2 import (U3L2_PRESETS, StepServo, observable_canonical,
                              solve_preset_servo)

QUESTION_TITLES = {
    1: "Questão 1 — Forma canônica controlável e seguidor de degrau",
    2: "Questão 2 — Seguidor de degrau no sistema de estados",
    3: "Questão 3 — Forma canônica controlável e polos em -1",
    4: "Questão 4 — Forma canônica observável e seguidor de degrau",
}


def _step(kind, value):
    return ResolutionStep(kind, str(value))


def _m(tex):
    return _step("math", tex)


def _header(title):
    return _step("heading", title)


def _remark(text):
    return _step("note", text)


def _ok(text):
    return _step("success", text)


def _section(title, steps):
    return QuestionSection(title, tuple(steps))


def _polynomial_matrix_tex(polynomial, n):
    coeff = sp.Poly(polynomial, S).all_coeffs()
    tokens = []
    for power in range(n, -1, -1):
        factor = coeff[n-power]
        if factor == 0:
            continue
        name = "I" if power == 0 else "A_a" if power == 1 else rf"A_a^{{{power}}}"
        if factor == 1:
            tokens.append(name)
        elif factor == -1:
            tokens.append("-" + name)
        else:
            tokens.append(sp.latex(factor)+name)
    return "+".join(tokens).replace("+-", "-")


def _poles_tex(poles):
    return r"\left\{" + ",\\;".join(sp.latex(p) for p in poles) + r"\right\}"


def _input_steps(number, realization, servo):
    if realization is not None:
        return [
            _header("Conferência da função de transferência do enunciado"),
            _m(r"G(s)="+sp.latex(realization.original_transfer)),
        ]
    return [
        _header("Conferência das matrizes A, B e C do enunciado"),
        _m(r"A="+sp.latex(servo.A)+r",\quad B="+sp.latex(servo.B)+r",\quad C="+sp.latex(servo.C)),
    ]


def observable_canonical_steps(realization):
    n = realization.A.rows
    coefficients_a = [-realization.A[n-1-i, n-1] for i in range(n)]
    coefficients_b = [realization.B[n-1-i, 0] for i in range(n)]
    denominator = sp.Poly.from_list([1, *coefficients_a], gens=S).as_expr()
    numerator = sum(coefficients_b[j] * S**(n-j-1) for j in range(n))
    alpha_str = r",\quad ".join(rf"\alpha_{{{i+1}}}={sp.latex(v)}" for i,v in enumerate(coefficients_a))
    beta_str = r",\quad ".join(rf"\beta_{{{i+1}}}={sp.latex(v)}" for i,v in enumerate(coefficients_b))
    return [
        _header("Forma canônica observável — apostila, equação (7.11)"),
        _remark("Identificar os coeficientes do denominador mônico e do numerador."),
        _m("D(s)="+sp.latex(denominator)+r"\quad\Longrightarrow\quad "+alpha_str),
        _m("N(s)="+sp.latex(numerator)+r"\quad\Longrightarrow\quad "+beta_str),
        _remark("Na forma observável: última coluna de A com os -α, B com os β em ordem inversa e C=[0 ... 1]."),
        _m(r"A_o="+sp.latex(realization.A)+r",\quad B_o="+sp.latex(realization.B)+r",\quad C_o="+sp.latex(realization.C)),
        _m(r"d="+sp.latex(realization.D)),
        _m(r"C_o(sI-A_o)^{-1}B_o+d="+sp.latex(realization.reconstructed_transfer)+"=G(s)"),
        _ok("Conferido: a realização observável reproduz exatamente a FT do enunciado."),
    ]


def servo_steps(design: StepServo):
    n = design.A.rows
    N = n+1
    desired_factor = sp.prod(S - p for p in design.desired_poles)
    output = [
        _header("1. Modelo interno do degrau — apostila, seção 7.10"),
        _m(r"e(t)=y(t)-r(t),\quad \dot e(t)=C\dot x(t),\quad z=\dot x,\quad w=\dot u"),
        _m(r"\frac{d}{dt}\begin{bmatrix}e\\z\end{bmatrix}="
           r"\underbrace{\begin{bmatrix}0&C\\0&A\end{bmatrix}}_{A_a}"
           r"\begin{bmatrix}e\\z\end{bmatrix}+"
           r"\underbrace{\begin{bmatrix}0\\B\end{bmatrix}}_{B_a}w"),
        _m("A_a="+sp.latex(design.augmented_A)+r",\quad B_a="+sp.latex(design.augmented_B)),
        _remark(f"Como a planta tem {n} estados, o seguidor exige {N} polos (1 integrador extra)."),
        _header("2. Matriz de controlabilidade do sistema aumentado"),
    ]
    # Produtos intermediários, mas jamais substitui A e B originais repetidamente.
    for power in range(1, N):
        label = "A_aB_a" if power == 1 else rf"A_a^{{{power}}}B_a"
        output.append(_m(label+"="+sp.latex(design.augmented_A**power * design.augmented_B)))
    output += [
        _m(r"U_a=[B_a\quad A_a B_a\quad\cdots\quad A_a^{"+str(n)+r"}B_a]="+sp.latex(design.augmented_U)),
        _m(r"\det(U_a)="+sp.latex(sp.factor(design.augmented_U.det()))),
        _m(r"\operatorname{posto}(U_a)="+str(design.augmented_rank)+"=N="+str(N)),
        _ok("O sistema aumentado é controlável; é possível alocar todos os polos solicitados."),
        _header("3. Polinômio desejado e polinômio matricial"),
        _m(r"\lambda_d="+_poles_tex(design.desired_poles)),
        _m(r"\Delta_d(s)="+sp.latex(desired_factor)+"="+sp.latex(sp.expand(design.desired_polynomial))),
        _m("q(A_a)="+_polynomial_matrix_tex(design.desired_polynomial, N)),
        _m("q(A_a)="+sp.latex(design.q_augmented)),
        _header("4. Fórmula de Ackermann — convenção u=Kx"),
        _m(r"K_a=-\begin{bmatrix}0&\cdots&0&1\end{bmatrix}U_a^{-1}q(A_a)"),
        _m(r"U_a^{-1}="+sp.latex(design.augmented_U_inverse)),
        _m(r"e_N^TU_a^{-1}="+sp.latex(design.ackermann_row)),
        _m(r"K_a=-\left(e_N^TU_a^{-1}\right)q(A_a)="+sp.latex(design.gain)),
        _m(r"K_a=[k_1\;k_2],\quad k_1="+sp.latex(design.k1)+r",\quad k_2="+sp.latex(design.k2)),
        _header("5. Lei de controle e conferência"),
        _m(r"w=k_1e+k_2z,\qquad u(t)=k_1\int_0^t (y(\tau)-r(\tau))d\tau+k_2x(t)"),
        _m(r"u(t)="+sp.latex(design.k1)+r"\int_0^t (y(\tau)-r(\tau))d\tau+"
           +sp.latex(design.k2)+"x(t)"),
        _m(r"A_{mf}=\begin{bmatrix}A+Bk_2&Bk_1\\C&0\end{bmatrix}="+sp.latex(design.closed_loop_A)),
        _m(r"B_r="+sp.latex(design.reference_B)+r",\quad C_r="+sp.latex(design.output_C)),
        _m(r"\det(sI-A_{mf})="+sp.latex(design.closed_loop_polynomial)+r"=\Delta_d(s)"),
        _m(r"G_{yr}(s)=C_r(sI-A_{mf})^{-1}B_r="+sp.latex(design.transfer_reference)),
        _m(r"\lim_{s\to0}G_{yr}(s)="+sp.latex(design.steady_state_gain)),
        _ok("Polos conferidos exatamente e erro em regime permanente nulo para referência degrau."),
    ]
    return output


def solve_full_question(number):
    if number not in QUESTION_TITLES:
        raise ValueError("Selecione uma questão de 1 a 4 da U3L2.")
    realization, servo = solve_preset_servo(number)
    data = _input_steps(number, realization, servo)
    sections = []
    if number in (1, 3, 4):
        if number == 4:
            canonical = observable_canonical_steps(realization)
            title = "a) Forma canônica observável"
        else:
            canonical = canonical_steps(realization)
            title = "a) Forma canônica controlável"
        sections.append(_section(title, [*data, *canonical]))
        sections.append(_section("b) Projeto do seguidor de referência — degrau", servo_steps(servo)))
    else:
        sections.append(_section("Projeto do seguidor de referência — degrau", [*data, *servo_steps(servo)]))
    return QuestionDocument(QUESTION_TITLES[number], tuple(sections), number=number, list_label="U3L2")


def individual_document(title, steps):
    return QuestionDocument(title, (_section("Resolução individual", steps),), number=0, list_label="U3L2")