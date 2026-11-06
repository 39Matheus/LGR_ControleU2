import pytest
import sympy as sp
from core.state_space import (S, U3L1_SYSTEMS, analyze_system, controllable_canonical,
                              place_observer, place_state_feedback, desired_polynomial)


def test_u3l1_q1_system_1():
    p = U3L1_SYSTEMS['Sistema 1 (Q1 e Q2)']
    a = analyze_system(p['A'], p['B'], p['C'])
    assert sp.expand(a.characteristic) == S**3 + 10*S**2 + 35*S + 38
    assert a.U == sp.Matrix([[0,0,-9],[0,-3,18],[3,-6,12]])
    assert a.V == sp.Matrix([[1,0,0],[-4,3,0],[13,-24,-3]])
    assert a.stable and a.controllable and a.observable


def test_u3l1_q1_system_2():
    p = U3L1_SYSTEMS['Sistema 2 (Q1 e Q3)']
    a = analyze_system(p['A'], p['B'], p['C'])
    assert sp.expand(a.characteristic) == S**3-8*S**2-28*S+207
    assert a.rank_U == a.rank_V == 3
    assert not a.stable


@pytest.mark.parametrize('key, K, L', [
    ('Sistema 1 (Q1 e Q2)', [sp.Rational(8,3),-sp.Rational(22,3),-sp.Rational(4,3)],
     [2,-sp.Rational(25,3),-sp.Rational(16,3)]),
    ('Sistema 2 (Q1 e Q3)', [72,117,-22],
     [-sp.Rational(752,4827),-sp.Rational(41341,4827),sp.Rational(4870,4827)]),
])
def test_u3l1_feedback_and_observer(key, K, L):
    p = U3L1_SYSTEMS[key]
    fk = place_state_feedback(p['A'],p['B'],p['C'],p['feedback_poles'])
    ob = place_observer(p['A'],p['B'],p['C'],p['observer_poles'])
    assert list(fk.gain) == K
    assert list(ob.gain) == L
    assert fk.result_characteristic == S**3+14*S**2+45*S+50
    assert ob.result_characteristic == S**3+12*S**2+22*S+20


def test_u3l1_q4():
    p = U3L1_SYSTEMS['Sistema da Q4 (G(s))']
    f = controllable_canonical(p['num'],p['den'])
    assert f.A == sp.Matrix([[0,1,0],[0,0,1],[-5,-7,-3]])
    assert f.B == sp.Matrix([0,0,1])
    assert f.C == sp.Matrix([[-2,4,3]])
    assert f.D == 0
    assert sp.simplify(f.original_transfer-f.reconstructed_transfer)==0
    ana = analyze_system(f.A,f.B,f.C)
    assert ana.controllable and ana.observable
    fk=place_state_feedback(f.A,f.B,f.C,p['feedback_poles'])
    ob=place_observer(f.A,f.B,f.C,p['observer_poles'])
    assert fk.gain == sp.Matrix([[-75,-49,-10]])
    assert ob.gain == sp.Matrix([-sp.Rational(23,241),sp.Rational(75,241),-sp.Rational(35,241)])


def test_noncontrollable_and_nonobservable_fail_gracefully():
    A=sp.diag(-1,-2); B=sp.Matrix([1,0]); C=sp.Matrix([[1,0]])
    a=analyze_system(A,B,C)
    assert not a.controllable and not a.observable
    with pytest.raises(ValueError,match='posto\\(U\\)'):
        place_state_feedback(A,B,C,[-3,-4])
    with pytest.raises(ValueError,match='posto\\(V\\)'):
        place_observer(A,B,C,[-3,-4])


def test_unconjugated_poles_and_bad_dimensions():
    with pytest.raises(ValueError,match='pares conjugados'):
        desired_polynomial([-1+sp.I,-2,-3],3)
    with pytest.raises(ValueError,match='quadrada'):
        analyze_system([[1,2]],[[1]],[[1]])