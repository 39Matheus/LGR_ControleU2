from core.controllers import ControllerDesigner
from core.presets import EXERCISES
from reports.html import build_controller_html


def _q1():
    x=EXERCISES["Questão 1 — PD (Mp e ts 5%)"]
    d=ControllerDesigner(
        x["num_g"],x["den_g"],x["num_h"],x["den_h"]
    ).design(x["controller"],x["specs"](),auto_refine=False)
    data={
        "controller":x["controller"],
        "num_g":x["num_g"],
        "den_g":x["den_g"],
        "num_h":x["num_h"],
        "den_h":x["den_h"],
        "specs":x["specs"](),
        "auto_refine":False,
    }
    return d,data


def test_complete_html_report_contains_full_derivation_and_plots():
    design,data=_q1()
    report=build_controller_html(design,data).decode("utf-8")

    assert "1 — Obtenção do polo desejado" in report
    assert "2 — Condição de ângulo" in report
    assert "3 — Condição de módulo" in report
    assert "4 — Controlador" in report
    assert "5 — Representação visual do projeto" in report
    assert "6 — Verificação numérica" in report
    assert "plotly" in report.lower()
    assert r"\,\cdot\," in report
    assert r"\cdotA" not in report
