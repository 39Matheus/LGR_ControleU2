import math
from core.controllers import ControllerDesigner
from core.presets import EXERCISES

def solve(name):
    x=EXERCISES[name]
    return ControllerDesigner(x["num_g"],x["den_g"],x["num_h"],x["den_h"]).design(x["controller"],x["specs"]())

def has_pole(poles,target,tol=1e-5):
    return any(abs(p-target)<tol for p in poles)

def test_q2_pd_reference():
    d=solve("Questão 2 — PD (ξ e ωn)")
    assert abs(d.desired_pole-complex(-0.35,0.3570714214))<1e-8
    assert math.isclose(d.zero_parameter,2.038857142857,rel_tol=1e-7)
    assert math.isclose(d.kc,7000.0,rel_tol=1e-7)
    assert math.isclose(d.kp,14272.0,rel_tol=1e-7)
    assert math.isclose(d.kd,7000.0,rel_tol=1e-7)
    assert has_pole(d.closed_loop_poles,d.desired_pole)

def test_q3_pi_requested_poles():
    d=solve("Questão 3 — PI (polos desejados)")
    assert math.isclose(d.zero_parameter,24/7,rel_tol=1e-7)
    assert math.isclose(d.kc,7.0,rel_tol=1e-7)
    assert math.isclose(d.kp,7.0,rel_tol=1e-7)
    assert math.isclose(d.ki,24.0,rel_tol=1e-7)
    assert has_pole(d.closed_loop_poles,complex(-4,4))
    assert has_pole(d.closed_loop_poles,complex(-4,-4))

def test_q4_pid_equal_zeros():
    d=solve("Questão 4 — PID (zeros iguais)")
    assert len(d.zero_locations)==2
    assert math.isclose(d.zero_locations[0],d.zero_locations[1],rel_tol=1e-12)
    assert has_pole(d.closed_loop_poles,d.desired_pole,tol=2e-5)

def test_q1_full_response_after_refinement():
    d=solve("Questão 1 — PD (Mp e ts 5%)")
    assert d.metrics and d.metrics.stable
    assert d.metrics.overshoot_percent is None or d.metrics.overshoot_percent<=10.001
    assert d.metrics.settling_time is None or d.metrics.settling_time<4.0
