import math

from core.controllers import ControllerDesigner
from core.geometry import angle_breakdown,magnitude_breakdown
from core.presets import EXERCISES,get_exercise


def solve(name,auto_refine=False):
    x=get_exercise(name)
    return ControllerDesigner(
        x["num_g"],x["den_g"],x["num_h"],x["den_h"]
    ).design(
        x["controller"],x["specs"](),
        auto_refine=auto_refine,
        controller_config=x.get("controller_config"),
    )


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


def test_pid_distinct_zeros_is_supported():
    x=get_exercise("Lista U2.2 — Questão 1 — PID + Tustin")
    d=ControllerDesigner(
        x["num_g"],x["den_g"],x["num_h"],x["den_h"]
    ).design(
        "PID",x["specs"](),
        controller_config={"pid_zero_mode":"distinct","fixed_zero":1.0},
    )
    assert len(d.zero_locations)==2
    assert not math.isclose(d.zero_locations[0],d.zero_locations[1],rel_tol=1e-6)
    assert math.isclose(d.zero_locations[0],-1.0,rel_tol=1e-10)
    assert math.isclose(d.zero_locations[1],-1.70624842595,rel_tol=1e-6)
    assert has_pole(d.closed_loop_poles,d.desired_pole,tol=2e-5)


def test_q1_uses_boundary_pole_by_default():
    d=solve("Questão 1 — PD (Mp e ts 5%)")
    assert math.isclose(d.desired_pole.real,-0.75,rel_tol=1e-10)
    assert not d.specification_summary["adjusted"]
    assert d.refinement_history==[]


def test_q1_full_response_after_optional_refinement():
    d=solve("Questão 1 — PD (Mp e ts 5%)",auto_refine=True)
    assert d.metrics and d.metrics.stable
    assert d.metrics.overshoot_percent is None or d.metrics.overshoot_percent<=10.001
    assert d.metrics.settling_time is None or d.metrics.settling_time<4.0


def test_second_list_q1_pid_reference_boundary_design():
    d=solve("Lista U2.2 — Questão 1 — PID + Tustin")
    assert math.isclose(d.desired_pole.real,-1.0,rel_tol=1e-10)
    assert math.isclose(d.desired_pole.imag,1.36437635384,rel_tol=1e-8)
    assert math.isclose(d.zero_locations[0],-1.33219104446,rel_tol=1e-7)
    assert math.isclose(d.zero_locations[1],d.zero_locations[0],rel_tol=1e-12)
    assert math.isclose(d.kc,0.53904177617,rel_tol=1e-7)
    assert has_pole(d.closed_loop_poles,d.desired_pole,tol=2e-5)


def test_second_list_q2_generic_pole_controller():
    d=solve("Lista U2.2 — Questão 2 — Genérico a/(s+b) + Euler")
    assert d.controller_type=="Genérico"
    assert d.solved_singularity_kind=="pole"
    assert math.isclose(d.solved_singularity_parameter,2.80608365019,rel_tol=1e-7)
    assert math.isclose(d.kc,3.79990494297,rel_tol=1e-7)
    assert len(d.zero_locations)==0
    assert len(d.pole_locations)==1
    assert math.isclose(d.pole_locations[0],-2.80608365019,rel_tol=1e-7)
    assert has_pole(d.closed_loop_poles,complex(-2.5,2.0),tol=2e-5)
    assert has_pole(d.closed_loop_poles,complex(-2.5,-2.0),tol=2e-5)


def test_magnitude_breakdown_matches_manual_gain_separation():
    x=get_exercise("Questão 2 — PD (ξ e ωn)")
    d=solve("Questão 2 — PD (ξ e ωn)")
    m=magnitude_breakdown(d,x["num_g"],x["den_g"],x["num_h"],x["den_h"])
    assert math.isclose(m["kg"],1/10000,rel_tol=1e-12)
    assert math.isclose(m["kt"],0.7,rel_tol=1e-7)
    assert math.isclose(m["kc"],d.kc,rel_tol=1e-10)


def test_angle_breakdown_exposes_labels_for_every_preset():
    for name,x in EXERCISES.items():
        d=solve(name)
        info=angle_breakdown(d,x["num_g"],x["den_g"],x["num_h"],x["den_h"])
        assert "plant_ids" in info
        assert set(info["plant_ids"])=={"zero","pole"}
        for contrib in info["plant_zeros"]:
            assert contrib.label in info["plant_ids"]["zero"]
        for contrib in info["plant_poles"]:
            assert contrib.label in info["plant_ids"]["pole"]
