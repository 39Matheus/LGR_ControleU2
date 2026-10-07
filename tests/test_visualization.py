from core.visualization import _focus_ranges


def test_focus_ranges_default_padding_is_seven():
    xr,yr=_focus_ranges([complex(-8,2),complex(1,4)])
    assert xr==[-15.0,15.0]
    assert yr==[-11.0,11.0]


def test_focus_ranges_accepts_custom_padding():
    xr,yr=_focus_ranges([complex(-8,2),complex(1,4)],padding=2.0)
    assert xr==[-10.0,10.0]
    assert yr==[-10.0,10.0]
