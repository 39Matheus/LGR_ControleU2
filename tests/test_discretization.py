import math

from core.controllers import ControllerDesigner
from core.discretization import (
    discretize_project,discretize_transfer,difference_equation_latex,
)
from core.presets import get_exercise
from ui.discretization_view import coefficient_rows


def _design(name):
    x=get_exercise(name)
    d=ControllerDesigner(
        x["num_g"],x["den_g"],x["num_h"],x["den_h"]
    ).design(
        x["controller"],x["specs"](),
        controller_config=x.get("controller_config"),
    )
    data=dict(x)
    data["specs"]=x["specs"]()
    return d,data


def test_forward_integrator_matches_chapter_formula():
    r=discretize_transfer([1],[1,0],2.0,"Euler / Forward")
    # 1/s -> T/(z-1)
    assert len(r.discrete_num)==1
    assert math.isclose(r.discrete_num[0],2.0,rel_tol=1e-10)
    assert len(r.discrete_den)==2
    assert math.isclose(r.discrete_den[0],1.0,rel_tol=1e-10)
    assert math.isclose(r.discrete_den[1],-1.0,rel_tol=1e-10)


def test_backward_integrator_matches_chapter_formula():
    r=discretize_transfer([1],[1,0],2.0,"Backward")
    # 1/s -> 2z/(z-1)
    assert len(r.discrete_num)==2
    assert math.isclose(r.discrete_num[0],2.0,rel_tol=1e-10)
    assert math.isclose(r.discrete_num[1],0.0,abs_tol=1e-12)
    assert r.discrete_den==[1.0,-1.0]


def test_tustin_integrator_matches_chapter_formula():
    r=discretize_transfer([1],[1,0],2.0,"Tustin")
    # 1/s -> (z+1)/(z-1) para T=2
    assert len(r.discrete_num)==2
    assert math.isclose(r.discrete_num[0],1.0,rel_tol=1e-10)
    assert math.isclose(r.discrete_num[1],1.0,rel_tol=1e-10)
    assert r.discrete_den==[1.0,-1.0]


def test_second_list_q1_tustin_controller():
    d,data=_design("Lista U2.2 — Questão 1 — PID + Tustin")
    r=discretize_project(d,data)
    assert r.method=="Tustin / Bilinear"
    assert math.isclose(r.sample_time,2.0)
    expected_num=[2.931910246997871,0.8352268820672055,0.05948373975994125]
    expected_den=[1.0,0.0,-1.0]
    assert len(r.discrete_num)==3
    for got,exp in zip(r.discrete_num,expected_num):
        assert math.isclose(got,exp,rel_tol=2e-6,abs_tol=2e-8)
    for got,exp in zip(r.discrete_den,expected_den):
        assert math.isclose(got,exp,rel_tol=2e-6,abs_tol=2e-8)


def test_second_list_q2_euler_complete_loop():
    d,data=_design("Lista U2.2 — Questão 2 — Genérico a/(s+b) + Euler")
    r=discretize_project(d,data)
    assert r.method=="Euler / Forward"
    expected_num=[7.599809885931558,15.199619771863116,0.0]
    expected_den=[1.0,5.806083650190114,8.224334600760457,5.806083650190114,7.224334600760457]
    for got,exp in zip(r.discrete_num,expected_num):
        assert math.isclose(got,exp,rel_tol=2e-6,abs_tol=2e-8)
    for got,exp in zip(r.discrete_den,expected_den):
        assert math.isclose(got,exp,rel_tol=2e-6,abs_tol=2e-8)
    assert "y[k]" in r.difference_equation_latex


def test_second_list_q2_display_is_normalized_and_readable():
    d,data=_design("Lista U2.2 — Questão 2 — Genérico a/(s+b) + Euler")
    r=discretize_project(d,data)
    assert r.discrete_latex.startswith(r"\frac{7.599")
    assert "10^{" not in r.discrete_latex
    assert "10^{" not in r.continuous_latex
    assert "5.8060837" in r.discrete_latex
    assert "8.2243346" in r.discrete_latex


def test_coefficient_table_aligns_by_power_and_hides_numeric_noise():
    rows=coefficient_rows(
        [7.599809886,15.19961977,1.77635684e-15],
        [1.0,5.80608365,8.224334601,5.80608365,7.224334601],
    )
    assert [r["potência"] for r in rows]==["z^4","z^3","z^2","z","1"]
    assert [r["numerador"] for r in rows]==[
        "0","0","7.599809886","15.19961977","0"
    ]
    assert [r["denominador"] for r in rows]==[
        "1","5.80608365","8.224334601","5.80608365","7.224334601"
    ]


def test_forward_pid_improper_difference_equation_is_explicitly_noncausal():
    r=discretize_transfer(
        [0.53904178,1.4362133,0.95665522],
        [1.0,0.0],
        2.0,
        "Euler / Forward",
        target_label="Gc(s)",
    )
    assert "inf" not in r.difference_equation_latex.lower()
    assert r"u[k+1]" in r.difference_equation_latex
    assert r"y[k-1]" in r.difference_equation_latex
    assert any("imprópria" in note for note in r.notes)


def test_second_list_q1_preset_uses_tustin_ui_label():
    x=get_exercise("Lista U2.2 — Questão 1 — PID + Tustin")
    assert x["discretization_method"]=="Tustin / Bilinear"
    assert math.isclose(x["sample_time"],2.0)
    assert x["discretization_target"]=="Controlador Gc(s)"
