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
