"""Resoluções por questão completa da Lista U3.1.

Q1 possui dois sistemas; Q2 e Q3 possuem dois projetos no respectivo
sistema; Q4 possui uma realização e dois projetos. A fonte é a lista,
cap. 7 da apostila e a mesma convenção u=Kx da calculadora.
"""
from dataclasses import dataclass
from typing import Tuple

import sympy as sp

from .state_space import (U3L1_SYSTEMS, analyze_system, controllable_canonical,
                          place_observer, place_state_feedback)
from .state_space_steps import ResolutionStep, build_resolution


QUESTION_TITLES = {
    1: "Questão 1 — Estabilidade, controlabilidade e observabilidade",
    2: "Questão 2 — Realimentação e observador (Sistema 1)",
    3: "Questão 3 — Realimentação e observador (Sistema 2)",
    4: "Questão 4 — Forma canônica e alocação de polos",
}


@dataclass(frozen=True)
class QuestionSection:
    title: str
    steps: Tuple[ResolutionStep, ...]


@dataclass(frozen=True)
class QuestionDocument:
    title: str
    sections: Tuple[QuestionSection, ...]
    number: int = 0


def _section(title, steps):
    return QuestionSection(title, tuple(steps))


def _data_once(analysis):
    return (
        ResolutionStep("heading", "Conferência dos dados do sistema"),
        ResolutionStep("math", rf"A={sp.latex(analysis.A)},\quad B={sp.latex(analysis.B)},\quad C={sp.latex(analysis.C)}"),
    )


def solve_full_question(number):
    """Gera todas as partes do enunciado, sem duplicar os dados do sistema."""
    if number not in QUESTION_TITLES:
        raise ValueError("Questão inválida: escolha 1, 2, 3 ou 4.")

    if number == 1:
        sections = []
        for i, key in enumerate(("Sistema 1 (Q1 e Q2)", "Sistema 2 (Q1 e Q3)"), start=1):
            preset = U3L1_SYSTEMS[key]
            analysis = analyze_system(preset["A"], preset["B"], preset["C"])
            sections.append(_section(
                f"Sistema {i} — itens a), b) e c)",
                build_resolution(analysis, operation="analysis", include_input=True),
            ))
    elif number in (2, 3):
        key = "Sistema 1 (Q1 e Q2)" if number == 2 else "Sistema 2 (Q1 e Q3)"
        preset = U3L1_SYSTEMS[key]
        analysis = analyze_system(preset["A"], preset["B"], preset["C"])
        feedback = place_state_feedback(preset["A"], preset["B"], preset["C"], preset["feedback_poles"])
        observer = place_observer(preset["A"], preset["B"], preset["C"], preset["observer_poles"])
        # Dados uma única vez antes das contas das duas partes.
        sections = [
            _section(f"Sistema {number-1} — dados", _data_once(analysis)),
            _section("a) Realimentação de estados" if number == 2 else "c) Realimentação de estados",
                     build_resolution(analysis, feedback, operation="K", include_input=False)),
            _section("b) Observador de estados" if number == 2 else "d) Observador de estados",
                     build_resolution(analysis, observer, operation="L", include_input=False)),
        ]
    else:
        preset = U3L1_SYSTEMS["Sistema da Q4 (G(s))"]
        canonical = controllable_canonical(preset["num"], preset["den"])
        analysis = analyze_system(canonical.A, canonical.B, canonical.C)
        observer = place_observer(canonical.A, canonical.B, canonical.C, preset["observer_poles"])
        feedback = place_state_feedback(canonical.A, canonical.B, canonical.C, preset["feedback_poles"])
        sections = [
            _section("a) Forma canônica controlável",
                     build_resolution(analysis, canonical=canonical, operation="canonical", include_input=True)),
            _section("b) Observador de estados (−1, −1, −2)",
                     build_resolution(analysis, observer, operation="L", include_input=False)),
            _section("c) Realimentação de estados (−4, −4, −5)",
                     build_resolution(analysis, feedback, operation="K", include_input=False)),
        ]
    return QuestionDocument(QUESTION_TITLES[number], tuple(sections), number)


def individual_document(title, steps):
    return QuestionDocument(title, (_section("Resolução individual", steps),), number=0)