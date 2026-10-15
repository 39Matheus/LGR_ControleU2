import math
from core.lgr import AnalisadorLGR

def test_second_order_regression():
    a=AnalisadorLGR([1],[1,4,0],[1],[1])
    a.calcular_polos_zeros()
    seg=a.calcular_segmentos_eixo_real()
    ass=a.calcular_assintotas()
    test=a.analisar_ponto_teste(complex(-2,2))
    assert len(a.polos)==2
    assert any(abs(p)<1e-8 for p in a.polos)
    assert any(abs(p+4)<1e-8 for p in a.polos)
    assert any(abs(x[0]+4)<1e-8 and abs(x[1])<1e-8 for x in seg["segmentos"])
    assert math.isclose(ass["centro"],-2.0,abs_tol=1e-8)
    assert test["pertence"]
    assert math.isclose(test["K"],8.0,rel_tol=1e-7)


def test_question_3_open_loop_cancellation_and_repeated_pole():
    a=AnalisadorLGR([5,25,20],[1,4,4],[0.2],[1,1])
    result=a.calcular_polos_zeros()

    assert len(a.zeros)==1
    assert abs(a.zeros[0]+4)<1e-8
    assert len(a.polos)==2
    assert all(abs(p+2)<1e-8 for p in a.polos)

    # G(s)H(s) = (s+4)/(s+2)^2 após cancelar o fator (s+1).
    assert str(result["P_expr"]) in {"(s + 4)/(s**2 + 4*s + 4)","(s + 4)/(s + 2)**2"}


def test_question_3_lgr_can_be_computed_with_repeated_pole():
    a=AnalisadorLGR([5,25,20],[1,4,4],[0.2],[1,1])
    a.calcular_polos_zeros()
    rlist,klist=a.calcular_lgr_exato()
    assert rlist.shape[1]==2
    assert len(klist)==rlist.shape[0]
