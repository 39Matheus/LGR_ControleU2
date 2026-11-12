"""U3L2: metodologia do seguidor de degrau, 4 documentos e regressão U3L1."""
from io import BytesIO
from zipfile import ZipFile

import pytest
import sympy as sp

from core.state_space import S
from core.state_space_u3l2 import (
    U3L2_PRESETS, design_step_servo, observable_canonical,
    solve_preset_servo,
)
from core.state_space_u3l2_questions import QUESTION_TITLES, solve_full_question
from core.state_space_questions import solve_full_question as solve_u3l1
from reports.state_space_export import (build_state_space_html, build_state_space_pdf,
                                        build_all_state_space_zip)


@pytest.mark.parametrize("number, expected_gain", [
    (1, [sp.Rational(303, 10), -sp.Rational(25337, 100), -sp.Rational(17891, 100), -16]),
    (2, [-sp.Rational(303, 10), sp.Rational(3373, 40), -sp.Rational(6617, 200), sp.Rational(843, 100)]),
    (3, [-1, -1, 0]),
    (4, [-4, -sp.Rational(5, 3), -sp.Rational(10, 3)]),
])
def test_preset_servo_exact_gain_and_integrator(number, expected_gain):
    realization, servo = solve_preset_servo(number)
    n = servo.A.rows
    assert len(servo.desired_poles) == n+1
    assert servo.augmented_U.shape == (n+1, n+1)
    assert servo.augmented_rank == n+1
    assert list(servo.gain) == expected_gain
    assert sp.simplify(servo.closed_loop_polynomial-servo.desired_polynomial) == 0
    assert sp.simplify(servo.steady_state_gain-1) == 0
    assert servo.reference_B[-1] == -1
    assert servo.closed_loop_A[:n, :n] == servo.A + servo.B*servo.k2
    assert servo.closed_loop_A[:n, n] == servo.B*servo.k1


def test_q1_and_q3_controllable_realizations():
    for number in (1,3):
        realization, servo = solve_preset_servo(number)
        assert realization.B[-1,0] == 1
        assert sp.simplify(
            realization.C*(S*sp.eye(realization.A.rows)-realization.A).inv()*realization.B
            - sp.Matrix([[realization.original_transfer]])
        ) == sp.zeros(1,1)
    _, q3 = solve_preset_servo(3)
    assert sp.factor(q3.desired_polynomial) == (S+1)**3


def test_q2_enunciado_matrices_are_correct():
    p = U3L2_PRESETS[2]
    assert sp.Matrix(p['A']) == sp.Matrix([[-6,1,0],[-11,0,1],[-6,0,0]])
    assert sp.Matrix(p['B']) == sp.Matrix([2,6,2])
    assert sp.Matrix(p['C']) == sp.Matrix([[1,0,0]])


def test_q4_observable_realization_matches_apostila_7_11():
    realization = observable_canonical([3],[1,2,-4])
    assert realization.A == sp.Matrix([[0,4],[1,-2]])
    assert realization.B == sp.Matrix([3,0])
    assert realization.C == sp.Matrix([[0,1]])
    assert sp.simplify(realization.original_transfer-3/(S**2+2*S-4)) == 0
    _, servo = solve_preset_servo(4)
    assert sp.factor(servo.desired_polynomial) == (S+2)**2*(S+3)


@pytest.mark.parametrize('question, expected_sections', [(1,2),(2,1),(3,2),(4,2)])
def test_complete_questions_sections_and_input_once(question,expected_sections):
    doc=solve_full_question(question)
    assert doc.list_label == 'U3L2'
    assert doc.title == QUESTION_TITLES[question]
    assert len(doc.sections)==expected_sections
    all_steps=[st for section in doc.sections for st in section.steps]
    assert sum('Conferência' in s.value and s.kind=='heading' for s in all_steps)==1
    text='\n'.join(s.value for s in all_steps)
    assert 'U_a' in text and 'q(A_a)' in text and 'k_1' in text
    if question == 4:
        assert 'observável' in doc.sections[0].title


def test_reject_noncontrollable_augmented_system_and_bad_poles():
    with pytest.raises(ValueError,match='posto\\(U_a\\)'):
        design_step_servo([[-1]],[[1]],[[0]],[-1,-2])
    with pytest.raises(ValueError,match='3 polos'):
        design_step_servo(sp.diag(-1,-2),sp.Matrix([1,1]),sp.Matrix([[1,1]]),[-2,-3])
    with pytest.raises(ValueError,match='pares conjugados'):
        design_step_servo([[-1]],[[1]],[[1]],[-1+sp.I,-2])


def test_u3l2_pdf_html_and_zip_while_preserving_u3l1():
    question=solve_full_question(4)
    pdf=build_state_space_pdf(question)
    html=build_state_space_html(question).decode('utf-8')
    assert pdf.startswith(b'%PDF-') and len(pdf)>1200
    assert 'MathJax' in html and 'U3L2' in html
    assert 'observável' in html and 'seguidor' in html
    with ZipFile(BytesIO(build_all_state_space_zip(solve_full_question))) as folder:
        names=sorted(folder.namelist())
        assert names==sorted(f'U3L2_Questao_{n}.{ext}' for n in range(1,5) for ext in ('pdf','html'))
        for n in range(1,5):
            assert folder.read(f'U3L2_Questao_{n}.pdf').startswith(b'%PDF-')
            assert b'MathJax' in folder.read(f'U3L2_Questao_{n}.html')
    existing=solve_u3l1(1)
    assert existing.list_label=='U3L1'
    assert 'U3L1' in build_state_space_html(existing).decode('utf-8')