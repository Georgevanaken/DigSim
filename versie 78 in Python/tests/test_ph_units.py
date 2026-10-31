import math

from digsim.calc_mucin_hydrogen import mucin_hydrogen
from digsim.calc_protein_hydrogen import calc_ph, protein_hydrogen


def test_ber_ph_returns_plausible_value():
    value = calc_ph(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert math.isfinite(value)
    assert 1.0 <= value <= 14.0


def test_protein_and_mucin_hydrogen_are_finite_and_decrease_in_expected_range():
    assert math.isfinite(protein_hydrogen(3.0, 1.0))
    assert protein_hydrogen(2.0, 1.0) > protein_hydrogen(3.0, 1.0)
    assert math.isfinite(mucin_hydrogen(3.0, 1.0))
    assert math.isfinite(mucin_hydrogen(4.0, 1.0))
    assert mucin_hydrogen(3.0, 1.0) > mucin_hydrogen(4.0, 1.0)
