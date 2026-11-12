"""Unidade 3, Lista 2: seguidor de referência do tipo degrau (apostila, §7.10).

Convenção do professor: e=y-r, w=du/dt=k1*e+k2*z, z=dx/dt,
logo u=k1*integral(e)+k2*x. Usa A_a=[0 C;0 A], B_a=[0;B] e
K_a=-e_last^T U_a^(-1) q(A_a) (malha fechada A_a+B_a*K_a).
Os coeficientes são exatos (SymPy), sem arredondamento dos polos.
"""
from dataclasses import dataclass
from typing import Tuple

import sympy as sp

from .state_space import (
    S, CanonicalRealization, as_state_matrices, controllable_canonical,
    desired_polynomial,
)


def observable_canonical(num, den):
    """Forma canônica observável da equação (7.11), dual da controlável (7.12)."""
    controllable = controllable_canonical(num, den)
    realization = CanonicalRealization(
        A=controllable.A.T,
        B=controllable.C.T,
        C=controllable.B.T,
        D=controllable.D,
        original_transfer=controllable.original_transfer,
        reconstructed_transfer=controllable.original_transfer,
    )
    transfer = sp.factor(
        (realization.C * (S * sp.eye(realization.A.rows) - realization.A).inv() * realization.B)[0]
        + realization.D
    )
    if sp.simplify(transfer - realization.original_transfer) != 0:
        raise ArithmeticError("A forma canônica observável não reproduziu G(s).")
    return realization


@dataclass(frozen=True)
class StepServo:
    A: sp.Matrix
    B: sp.Matrix
    C: sp.Matrix
    augmented_A: sp.Matrix
    augmented_B: sp.Matrix
    augmented_U: sp.Matrix
    augmented_U_inverse: sp.Matrix
    augmented_rank: int
    desired_poles: Tuple[sp.Expr, ...]
    desired_polynomial: sp.Expr
    q_augmented: sp.Matrix
    ackermann_row: sp.Matrix
    gain: sp.Matrix
    k1: sp.Expr
    k2: sp.Matrix
    closed_loop_A: sp.Matrix
    reference_B: sp.Matrix
    output_C: sp.Matrix
    closed_loop_polynomial: sp.Expr
    transfer_reference: sp.Expr
    steady_state_gain: sp.Expr


def design_step_servo(A, B, C, poles):
    """Projeta o seguidor de degrau e valida todas as identidades da §7.10."""
    A, B, C = as_state_matrices(A, B, C)
    n = A.rows
    if len(poles) != n + 1:
        raise ValueError(f"O seguidor deve ter {n+1} polos (a planta tem {n} estados mais o integrador).")
    desired = desired_polynomial(poles, n + 1)
    augmented_A = sp.zeros(n + 1)
    augmented_A[0, 1:] = C
    augmented_A[1:, 1:] = A
    augmented_B = sp.zeros(n + 1, 1)
    augmented_B[1:, 0] = B
    U = sp.Matrix.hstack(*(augmented_A**j * augmented_B for j in range(n + 1)))
    rank = U.rank()
    if rank != n + 1:
        raise ValueError(
            f"Seguidor não alocável: posto(U_a)={rank} < {n+1}; "
            "o sistema aumentado não é controlável."
        )
    inverted = U.inv()
    poly = sp.Poly(desired, S)
    qA = sum((coef * augmented_A**power for power, coef in
               enumerate(reversed(poly.all_coeffs()))), sp.zeros(n + 1))
    e = sp.zeros(1, n + 1)
    e[0, n] = 1
    row = sp.simplify(e * inverted)
    gain = sp.simplify(-row * qA)
    k1 = gain[0, 0]
    k2 = gain[:, 1:]
    F = sp.Matrix.vstack(
        sp.Matrix.hstack(A + B*k2, B*k1),
        sp.Matrix.hstack(C, sp.zeros(1, 1)),
    )
    Br = sp.zeros(n + 1, 1)
    Br[n, 0] = -1
    Cr = sp.Matrix.hstack(C, sp.zeros(1, 1))
    characteristic = sp.expand(F.charpoly(S).as_expr())
    if sp.simplify(characteristic - desired) != 0:
        raise ArithmeticError("A alocação não produziu o polinômio desejado.")
    Gyr = sp.factor((Cr * (S*sp.eye(n+1) - F).inv() * Br)[0])
    dc = sp.simplify(sp.limit(Gyr, S, 0))
    if sp.simplify(dc - 1) != 0:
        raise ArithmeticError("O seguidor não apresenta ganho estacionário unitário para degrau.")
    return StepServo(
        A, B, C, augmented_A, augmented_B, U, inverted, rank,
        tuple(sp.sympify(p) for p in poles), desired, qA, row, gain,
        k1, k2, F, Br, Cr, characteristic, Gyr, dc,
    )


U3L2_PRESETS = {
    1: {
        "representation": "G(s)", "num": (3, 4, -2), "den": (1, 3, 7, 5),
        "canonical": "controlavel",
        "poles": (-1+sp.I/10, -1-sp.I/10, -5, -12),
    },
    2: {
        "representation": "Matrizes A, B, C",
        "A": ((-6, 1, 0), (-11, 0, 1), (-6, 0, 0)),
        "B": ((2,), (6,), (2,)), "C": ((1, 0, 0),),
        "poles": (-1+sp.I/10, -1-sp.I/10, -5, -12),
    },
    3: {
        "representation": "G(s)", "num": (2, 1), "den": (1, 3, 0),
        "canonical": "controlavel", "poles": (-1, -1, -1),
    },
    4: {
        "representation": "G(s)", "num": (3,), "den": (1, 2, -4),
        "canonical": "observavel", "poles": (-2, -2, -3),
    },
}


def preset_realization(number):
    p = U3L2_PRESETS[number]
    if p["representation"] != "G(s)":
        return None
    return (observable_canonical if p["canonical"] == "observavel" else controllable_canonical)(
        p["num"], p["den"]
    )


def solve_preset_servo(number):
    p = U3L2_PRESETS[number]
    realization = preset_realization(number)
    if realization is not None:
        A, B, C = realization.A, realization.B, realization.C
    else:
        A, B, C = p["A"], p["B"], p["C"]
    return realization, design_step_servo(A, B, C, p["poles"])