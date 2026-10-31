import numpy as np

from digsim import state as st
from digsim.enums import SMALL_INTESTINE, Component, Segment


def test_segment_numbers_match_pascal_header():
    assert (Segment.MOUTH, Segment.CORPUS, Segment.ANTRUM) == (2, 3, 4)
    assert Segment.DUODENUM_1 == 5 and Segment.JEJUNUM_1 == 7 and Segment.ILEUM_1 == 10
    assert Segment.ILEUM_6 == 15 and Segment.COLON == 16


def test_small_intestine_slice_covers_duodenum_to_ileum():
    rows = np.arange(17)[SMALL_INTESTINE]
    assert rows[0] == Segment.DUODENUM_1 and rows[-1] == Segment.ILEUM_6 and len(rows) == 11


def test_component_numbers_match_pascal_header():
    expected = {"WATER": 1, "GLUCOSE": 3, "SUCROSE": 21, "FRUCTOSE": 22, "MALTOSE": 23,
                "TRACER_MALTOSE": 24, "TRACER_GLUCOSE_1": 34, "TRACER_GLUCOSE_2": 35,
                "TRACER_SUCROSE": 36, "CONJUGATED_BILE_SALT": 51, "BILE_SALT_LCFA_COMPLEX": 54,
                "UNCONJUGATED_BILE_SALT": 59, "NACL": 78, "H_ION": 80}
    for name, number in expected.items():
        assert Component[name] == number
    assert max(Component) <= 110 and st.liquid.shape[1] > max(Component)
