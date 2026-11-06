"""Unidade 3, Lista 1: espaço de estados pela convenção da apostila.

A realimentação é u(t)=Kx(t)+r(t), logo A_mf=A+B*K.
O observador tem erro de estimação regido por A-L*C.
Os cálculos usam SymPy para preservar frações e verificações exatas.
"""
from dataclasses import dataclass

import sympy as sp

S = sp.Symbol("s")


@dataclass(frozen=True)
class StateAnalysis:
    A: sp.Matrix
    B: sp.Matrix
    C: sp.Matrix
    characteristic: sp.Expr
    eigenvalues: tuple
    U: sp.Matrix
    V: sp.Matrix
    rank_U: int
    rank_V: int
    stable: bool
    controllable: bool
    observable: bool


@dataclass(frozen=True)
class PolePlacement:
    mode: str
    desired: sp.Expr
    qA: sp.Matrix
    gain: sp.Matrix
    result_A: sp.Matrix
    result_characteristic: sp.Expr
    result_eigenvalues: tuple
    inverse_matrix: sp.Matrix


@dataclass(frozen=True)
class CanonicalRealization:
    A: sp.Matrix
    B: sp.Matrix
    C: sp.Matrix
    D: sp.Expr
    original_transfer: sp.Expr
    reconstructed_transfer: sp.Expr


def as_state_matrices(A, B, C):
    A, B, C = sp.Matrix(A), sp.Matrix(B), sp.Matrix(C)
    if A.rows != A.cols or A.rows == 0:
        raise ValueError("A deve ser uma matriz quadrada não vazia.")
    n = A.rows
    if B.shape != (n, 1):
        raise ValueError(f"B deve ser coluna {n}×1 (um sinal de entrada).")
    if C.shape != (1, n):
        raise ValueError(f"C deve ser linha 1×{n} (uma saída).")
    if any(x.free_symbols for M in (A, B, C) for x in M):
        raise ValueError("Informe somente entradas numéricas nas matrizes A, B e C.")
    return A, B, C


def _numerical_eigenvalues(matrix):
    result = []
    for eig, multiplicity in matrix.eigenvals().items():
        result.extend([complex(sp.N(eig, 14))] * int(multiplicity))
    return tuple(sorted(result, key=lambda z: (round(z.real, 10), round(z.imag, 10))))


def analyze_system(A, B, C):
    A, B, C = as_state_matrices(A, B, C)
    n = A.rows
    U = sp.Matrix.hstack(*(A**j * B for j in range(n)))
    V = sp.Matrix.vstack(*(C * A**j for j in range(n)))
    eigen = _numerical_eigenvalues(A)
    return StateAnalysis(
        A=A, B=B, C=C,
        characteristic=sp.factor(A.charpoly(S).as_expr()),
        eigenvalues=eigen, U=U, V=V,
        rank_U=U.rank(), rank_V=V.rank(),
        stable=all(z.real < -1e-9 for z in eigen),
        controllable=U.rank() == n, observable=V.rank() == n,
    )


def desired_polynomial(poles, n):
    if len(poles) != n:
        raise ValueError(f"Informe exatamente {n} polos, incluindo multiplicidades.")
    values = [sp.sympify(pole) for pole in poles]
    if any(v.free_symbols for v in values):
        raise ValueError("Os polos devem ser numéricos.")
    poly = sp.Poly(sp.expand(sp.prod(S - v for v in values)), S)
    if any(sp.simplify(sp.im(coef)) != 0 for coef in poly.all_coeffs()):
        raise ValueError("Para matrizes reais, informe polos reais ou pares conjugados.")
    return sp.expand(poly.as_expr())


def _evaluate_polynomial_at_matrix(poly, A):
    coeff = sp.Poly(poly, S).all_coeffs()
    return sum((coef * A**(len(coeff) - i - 1) for i, coef in enumerate(coeff)), sp.zeros(A.rows))


def place_state_feedback(A, B, C, poles):
    analysis = analyze_system(A, B, C)
    A, B, C = analysis.A, analysis.B, analysis.C
    if not analysis.controllable:
        raise ValueError(f"Não é possível alocar arbitrariamente os polos: posto(U)={analysis.rank_U} < n={A.rows}.")
    desired = desired_polynomial(poles, A.rows)
    qA = _evaluate_polynomial_at_matrix(desired, A)
    row = sp.zeros(1, A.rows)
    row[0, A.rows - 1] = 1
    inverse = analysis.U.inv()
    gain = sp.simplify(-row * inverse * qA)  # u=Kx => A+BK
    result_A = A + B * gain
    actual = sp.expand(result_A.charpoly(S).as_expr())
    if sp.simplify(actual - desired) != 0:
        raise ArithmeticError("Falha na verificação do polinômio de malha fechada.")
    return PolePlacement("K", desired, qA, gain, result_A, actual, _numerical_eigenvalues(result_A), inverse)


def place_observer(A, B, C, poles):
    analysis = analyze_system(A, B, C)
    A, B, C = analysis.A, analysis.B, analysis.C
    if not analysis.observable:
        raise ValueError(f"Não é possível alocar arbitrariamente os polos do observador: posto(V)={analysis.rank_V} < n={A.rows}.")
    desired = desired_polynomial(poles, A.rows)
    qA = _evaluate_polynomial_at_matrix(desired, A)
    col = sp.zeros(A.rows, 1)
    col[A.rows - 1, 0] = 1
    inverse = analysis.V.inv()
    gain = sp.simplify(qA * inverse * col)  # A-LC
    result_A = A - gain * C
    actual = sp.expand(result_A.charpoly(S).as_expr())
    if sp.simplify(actual - desired) != 0:
        raise ArithmeticError("Falha na verificação do polinômio do observador.")
    return PolePlacement("L", desired, qA, gain, result_A, actual, _numerical_eigenvalues(result_A), inverse)


def controllable_canonical(num, den):
    """A forma companheira da eq. (7.12), sem cancelar polos/zeros da FT."""
    num_coeff = [sp.sympify(v) for v in num]
    den_coeff = [sp.sympify(v) for v in den]
    if not num_coeff or not den_coeff or den_coeff[0] == 0:
        raise ValueError("Numerador e denominador precisam de coeficientes válidos.")
    N = sp.Poly.from_list(num_coeff, gens=S)
    D = sp.Poly.from_list(den_coeff, gens=S)
    n = D.degree()
    if n < 1 or N.degree() > n:
        raise ValueError("A função deve ser própria e ter denominador de grau pelo menos 1.")
    if any(v.free_symbols for v in num_coeff + den_coeff):
        raise ValueError("Informe apenas coeficientes numéricos.")
    N = sp.Poly(N.as_expr() / D.LC(), S)
    D = sp.Poly(D.as_expr() / D.LC(), S)
    direct, remainder = sp.div(N, D)
    d = sp.simplify(direct.as_expr())
    alphas = [D.nth(k) for k in range(n)]
    betas = [remainder.nth(k) for k in range(n)]
    A = sp.zeros(n)
    for i in range(n - 1):
        A[i, i + 1] = 1
    for j in range(n):
        A[n - 1, j] = -alphas[j]
    B = sp.zeros(n, 1)
    B[n - 1, 0] = 1
    C = sp.Matrix([betas])
    original = sp.factor(N.as_expr() / D.as_expr())
    rebuilt = sp.factor((C * (S * sp.eye(n) - A).inv() * B)[0] + d)
    if sp.simplify(original - rebuilt) != 0:
        raise ArithmeticError("A forma canônica não reproduziu G(s).")
    return CanonicalRealization(A, B, C, d, original, rebuilt)


U3L1_SYSTEMS = {
    "Sistema 1 (Q1 e Q2)": {
        "A": [[-4, 3, 0], [-1, -4, -1], [0, 0, -2]],
        "B": [[0], [0], [3]], "C": [[1, 0, 0]],
        "feedback_poles": [-2+sp.I, -2-sp.I, -10],
        "observer_poles": [-1+sp.I, -1-sp.I, -10],
    },
    "Sistema 2 (Q1 e Q3)": {
        "A": [[-5, 1, 0], [2, 7, -1], [0, -3, 6]],
        "B": [[0], [0], [1]], "C": [[1, -2, 3]],
        "feedback_poles": [-2+sp.I, -2-sp.I, -10],
        "observer_poles": [-1+sp.I, -1-sp.I, -10],
    },
    "Sistema da Q4 (G(s))": {
        "num": [3, 4, -2], "den": [1, 3, 7, 5],
        "feedback_poles": [-4, -4, -5],
        "observer_poles": [-1, -1, -2],
    },
}