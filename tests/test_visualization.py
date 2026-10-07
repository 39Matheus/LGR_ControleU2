from core.visualization import _focus_ranges,focus_description


def test_focus_ranges_default_is_proportional_20_x_10_y():
    xr,yr=_focus_ranges([complex(-20,2),complex(1,15)])
    assert xr==[-24.0,24.0]
    assert yr==[-16.5,16.5]
    assert focus_description()=="automático proporcional: +20% em x e +10% em y"


def test_focus_ranges_accepts_custom_proportional_margins():
    focus={
        "mode":"proportional",
        "x_percent":5.0,
        "y_percent":30.0,
    }
    xr,yr=_focus_ranges([complex(-20,2),complex(1,15)],focus=focus)
    assert xr==[-21.0,21.0]
    assert yr==[-19.5,19.5]


def test_focus_ranges_accepts_custom_linear_margins():
    focus={
        "mode":"linear",
        "x_padding":2.0,
        "y_padding":4.0,
    }
    xr,yr=_focus_ranges([complex(-20,2),complex(1,15)],focus=focus)
    assert xr==[-22.0,22.0]
    assert yr==[-19.0,19.0]


def test_focus_ranges_keeps_minimum_ten():
    xr,yr=_focus_ranges([complex(-2,1),complex(1,3)])
    assert xr==[-10.0,10.0]
    assert yr==[-10.0,10.0]


def test_legacy_padding_argument_remains_linear_and_symmetric():
    xr,yr=_focus_ranges([complex(-20,2),complex(1,15)],padding=7.0)
    assert xr==[-27.0,27.0]
    assert yr==[-22.0,22.0]
