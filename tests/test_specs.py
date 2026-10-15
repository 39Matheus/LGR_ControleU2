import math
from core.specs import DesignSpecs,resolve_specs,zeta_from_overshoot

def test_zeta_10_percent():
    assert math.isclose(zeta_from_overshoot(10),0.5911550338,rel_tol=1e-8)

def test_ts_5_percent_rule():
    p,zeta,wn,sigma=resolve_specs(DesignSpecs.from_mp_ts(10,4,0.05))
    assert math.isclose(sigma,0.75,rel_tol=1e-12)
    assert math.isclose(p.real,-0.75,rel_tol=1e-12)
    assert p.imag>0


def test_zeta_wn_direct_pole_formula():
    p,zeta,wn,sigma=resolve_specs(DesignSpecs.from_zeta_wn(0.7,0.5))
    assert math.isclose(sigma,0.35,rel_tol=1e-12)
    assert math.isclose(p.real,-0.35,rel_tol=1e-12)
    assert math.isclose(p.imag,0.5*math.sqrt(1-0.7**2),rel_tol=1e-12)
