import math

from core.controllers import ControllerDesigner
from core.discretization import (
    discretize_project,discretize_transfer,difference_equation_latex,
)
from core.presets import get_exercise


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
