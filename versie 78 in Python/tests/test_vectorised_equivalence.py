"""The vectorised volumina/decomposition must reproduce the original loop translation."""
import copy
from pathlib import Path

import numpy as np
import pytest
import ref_compartment_properties as ref_cp
import ref_main_part2 as ref_p2

from digsim import state as st
from digsim.compartment_properties import volumina
from digsim.diet_from_file import read_compounds, read_diet
from digsim.digsim_v78 import make_step
from digsim.enums import Component as C
from digsim.main_part1 import initialise
from digsim.main_part2 import decomposition
from digsim.model_parameters import set_model_parameters

DATA = Path(__file__).parent / "data"

# Loop counters that used to be module globals; they are local variables now, so their leftovers are not results.
LOOP_VARS = ("i", "j", "s", "compound_number", "dish_number")


GROUP_NAMES = (
    "hormones", "stomach", "glucose", "bile", "diet", "calories", "hunger", "beta_cell", "lipolysis",
    "secretion", "osmotic", "brush_enzymes", "clock", "solid_props",
)


def snapshot():
    snap = {k: copy.deepcopy(v) for k, v in vars(st).items()
            if not k.startswith("__") and (isinstance(v, (np.ndarray, int, float, bool)))}
    for g in GROUP_NAMES:
        for k, v in vars(getattr(st, g)).items():
            snap[f"{g}.{k}"] = copy.deepcopy(v)
    return snap


def restore(snap):
    for k, v in snap.items():
        group, _, field = k.partition(".")
        if field and group in GROUP_NAMES:
            setattr(getattr(st, group), field, copy.deepcopy(v))
        else:
            setattr(st, k, copy.deepcopy(v))


@pytest.fixture(scope="module")
def running_model():
    import os
    old = os.getcwd()
    os.chdir(DATA)
    try:
        np.seterr(all="ignore")
        st.comp = 110
        set_model_parameters()
        read_diet()
        initialise()
        st.clock.substeps = 60
        st.clock.dt = 1 / 60
        volumina()
        st.dish_number = 1
        st.clock.minutes = -10
        st.ph[1] = 7
        st.diet.liquid_viscosity_dish[1] = 1
        st.liquid_dish[1, 1:5] = [250, 1, 50, 10]
        st.eating_speed_dish[1] = 50
        st.diet.new_meal_number = 1
        read_compounds()
        for _ in range(30):
            st.clock.minutes += st.clock.dt
            make_step()
        return snapshot()
    finally:
        os.chdir(old)


def test_state_has_compounds(running_model):
    assert running_model["diet.number_of_compounds"] > 0
    assert np.any(running_model["compound"] != 0)


def test_decomposition_identical(running_model):
    restore(running_model)
    ref_p2.decomposition()
    expected = snapshot()
    restore(running_model)
    decomposition()
    got = snapshot()
    for k in expected:
        if k in LOOP_VARS:
            continue
        if isinstance(expected[k], np.ndarray):
            assert np.array_equal(expected[k], got[k], equal_nan=True), k
        else:
            assert expected[k] == got[k], k


def test_calc_all_ph_identical(running_model):
    import ref_ph

    from digsim.calc_protein_hydrogen import calc_all_ph

    rng = np.random.default_rng(7)

    def prepare():
        restore(running_model)
        for comp in (C.CASEIN, C.CASEINPEP, C.CASEINPAN, C.CASEIN_PEPTIDES, C.SMALL_PEPTIDES, C.BSA,
                     C.MUCUS_PROTEIN, C.BLG, C.PROTEIN, C.CITRIC_ACID, C.H_ION):
            st.liquid[0:17, comp] = rng.random(17) * 10 ** rng.uniform(-6, -1, 17)
        st.volume[0:17] = rng.uniform(0.01, 500, 17)

    for _ in range(25):
        prepare()
        ref_ph.calc_all_ph()
        expected = snapshot()
        restore(expected)
        calc_all_ph()
        got = snapshot()
        assert np.array_equal(expected["ph"], got["ph"]), (expected["ph"], got["ph"])
        assert expected["stomach.ph_fundus"] == got["stomach.ph_fundus"]


@pytest.mark.parametrize("cycle", [0.01, 0.15, 0.25])
def test_gastric_redistribution_identical(running_model, cycle):
    import ref_gastric_motility as ref_gm

    from digsim.gastric_motility import gastric_redistribution

    def prepare():
        restore(running_model)
        rng = np.random.default_rng(3)
        n = st.diet.number_of_compounds
        st.cycle = cycle
        st.liquid[1:17, 1:st.comp + 1] += rng.random((16, st.comp)) * 1e-3
        st.compound[1:17, 1:n + 1] += rng.random((16, n)) * 1e-3
        st.stomach.fundus[1:st.comp + 1] += rng.random(st.comp) * 1e-3
        st.stomach.solids_corpus[1:st.comp + 1] += rng.random(st.comp) * 1e-3
        st.stomach.solids_fundus[1:st.comp + 1] += rng.random(st.comp) * 1e-3

    prepare()
    ref_gm.gastric_redistribution()
    expected = snapshot()
    prepare()
    gastric_redistribution()
    got = snapshot()
    for k in expected:
        if k in LOOP_VARS:
            continue
        if isinstance(expected[k], np.ndarray):
            assert np.array_equal(expected[k], got[k], equal_nan=True), k
        else:
            assert expected[k] == got[k], k


def test_calories_identical(running_model):
    import ref_main_part1 as ref_p1

    from digsim.main_part1 import caloric_densities_segments, calories

    def prepare():
        restore(running_model)
        rng = np.random.default_rng(5)
        n = st.diet.number_of_compounds
        st.liquid[1:17, 1:st.comp + 1] = rng.random((16, st.comp)) * 1e-2
        st.compound[1:17, 1:n + 1] = rng.random((16, n)) * 1e-2
        st.concentration_brush[1:17, 1:st.comp + 1] = rng.random((16, st.comp)) * 1e-2
        st.volume[1:17] = rng.random(16)
        st.volume[7] = 1e-12

    for ref, new in ((ref_p1.caloric_densities_segments, caloric_densities_segments),
                     (ref_p1.calories, calories)):
        prepare()
        ref()
        expected = snapshot()
        prepare()
        new()
        got = snapshot()
        for k in expected:
            if isinstance(expected[k], np.ndarray):
                assert np.array_equal(expected[k], got[k], equal_nan=True), k
            else:
                assert expected[k] == got[k], k


def test_digestion_identical(running_model):
    from digsim.main_part2 import digestion

    def prepare(seed):
        rng = np.random.default_rng(seed)
        restore(running_model)
        n = st.diet.number_of_compounds
        st.liquid[1:17, 1:st.comp + 1] = rng.random((16, st.comp)) * 1e-3 * (rng.random((16, st.comp)) > 0.3)
        st.compound[1:17, 1:n + 1] = rng.random((16, n)) * 1e-3
        st.volume[1:17] = rng.uniform(0.01, 300, 16)
        st.ph[1:17] = rng.uniform(0.3, 9, 16)
        st.concentration_lumen[1:17, 1:st.comp + 1] = rng.random((16, st.comp)) * 1e-2

    for seed in range(10):
        prepare(seed)
        ref_p2.digestion()
        expected = snapshot()
        prepare(seed)
        digestion()
        got = snapshot()
        for k in expected:
            if k in LOOP_VARS:
                continue
            if isinstance(expected[k], np.ndarray):
                assert np.array_equal(expected[k], got[k], equal_nan=True), k
            else:
                assert expected[k] == got[k] and isinstance(expected[k], int) == isinstance(got[k], int), k  # 0 and 0.0 print differently


def test_gut_hormones_equivalent(running_model):
    import ref_main_part3 as ref_p3

    from digsim.main_part3 import gut_hormones

    def prepare():
        restore(running_model)
        rng = np.random.default_rng(1)
        st.concentration_brush[5:16, 1:111] = rng.random((11, 110)) * 1e-4
        st.concentration_cell[5:16, 1:111] = rng.random((11, 110)) * 1e-4
        st.liquid[5:16, 1:111] = rng.random((11, 110)) * 1e-2
        st.ph[5:16] = rng.uniform(4, 8, 11)
        st.hormones.pyy = 0.3

    prepare()
    ref_p3.gut_hormones()
    expected = snapshot()
    prepare()
    gut_hormones()
    got = snapshot()
    # Loop counters are leftovers, not results.
    for k in expected:
        if k in LOOP_VARS:
            continue
        if isinstance(expected[k], np.ndarray):
            np.testing.assert_allclose(got[k], expected[k], rtol=1e-11, atol=1e-300, err_msg=k)
        else:
            assert expected[k] == pytest.approx(got[k], rel=1e-11), k
    assert expected["hormones.total_rate_cck_release"] != 0


def _compare(expected, got, rtol=1e-9):
    for k in expected:
        if k in LOOP_VARS:
            continue
        if isinstance(expected[k], np.ndarray):
            np.testing.assert_allclose(got[k], expected[k], rtol=rtol, atol=1e-300, err_msg=k)
        else:
            assert expected[k] == pytest.approx(got[k], rel=rtol), k


def _randomise(seed):
    rng = np.random.default_rng(seed)
    shape = st.liquid[5:16, 1:111].shape
    st.concentration_brush[5:16, 1:111] = rng.random(shape) * rng.integers(0, 2, shape) * 1e-2
    st.concentration_lumen[2:16, 1:111] = rng.random((14, 110)) * 1e-2
    st.liquid[2:17, 1:111] = rng.random((15, 110)) * 5
    st.liquid[5:16, 22] = 1e-9
    st.mesenteric[5:16, 1:111] = rng.random(shape) * 1e-2
    st.absorption_rate[5:16, 1:111] = rng.random(shape) * 1e-3
    st.outflow[1:16] = rng.random(15)
    st.outflow_dish[1:6] = rng.random(5)
    st.concentration_dish[1:6, 1:111] = rng.random((5, 110))
    st.concentration_compound[2:16, 1:51] = rng.random((14, 50))
    st.concentration_compound_dish[1:6, 1:51] = rng.random((5, 50))
    st.compound[2:17, 1:51] = rng.random((15, 50))
    st.compound_dish[1:6, 1:51] = rng.random((5, 50))
    st.solid[1:51] = rng.integers(0, 20, 50)


def test_absorption_equivalent(running_model):
    import ref_main_part3 as ref_p3

    from digsim.main_part3 import absorption_in_small_intestine
    for seed in (1, 2):
        restore(running_model)
        _randomise(seed)
        pre = snapshot()
        ref_p3.absorption_in_small_intestine()
        expected = snapshot()
        restore(pre)
        absorption_in_small_intestine()
        _compare(expected, snapshot())
    assert np.any(expected["absorption_rate"] != 0)


def test_brush_border_enzymes_identical(running_model):
    import ref_brush_border_enzymes as ref_bbe

    from digsim import brush_border_enzymes as bbe

    def prepare(seed):
        rng = np.random.default_rng(seed)
        restore(running_model)
        st.concentration_brush[1:17, 1:111] = rng.random((16, 110)) * 1e-3
        st.concentration_lumen[1:17, 1:111] = rng.random((16, 110)) * 1e-2 * (rng.random((16, 110)) > 0.2)
        st.absorption_rate[1:17, 1:111] = rng.random((16, 110)) * 1e-6
        st.brush_enzymes.maltase_inhibition[:] = rng.random(len(st.brush_enzymes.maltase_inhibition))
        st.brush_enzymes.sucrase_inhibition[:] = rng.random(len(st.brush_enzymes.sucrase_inhibition))

    for seed in range(8):
        for ref, new in ((ref_bbe.sucrase, bbe.sucrase), (ref_bbe.maltase, bbe.maltase)):
            prepare(seed)
            ref()
            expected = snapshot()
            prepare(seed)
            new()
            _compare(expected, snapshot())


def test_mucus_layer_equivalent(running_model):
    import ref_main_part3 as ref_p3

    from digsim.main_part3 import mucus_layer
    for seed in (1, 2):
        restore(running_model)
        _randomise(seed)
        pre = snapshot()
        ref_p3.mucus_layer()
        expected = snapshot()
        restore(pre)
        mucus_layer()
        _compare(expected, snapshot())


def test_transport_equivalent(running_model):
    import ref_transport

    from digsim.digsim_v78 import transport
    restore(running_model)
    _randomise(3)
    pre = snapshot()
    ref_transport.transport()
    expected = snapshot()
    restore(pre)
    transport()
    _compare(expected, snapshot(), rtol=1e-14)


def test_volumina_equivalent(running_model):
    restore(running_model)
    ref_cp.volumina()
    expected = snapshot()
    restore(running_model)
    volumina()
    got = snapshot()
    for k in expected:
        if k in LOOP_VARS:
            continue
        if isinstance(expected[k], np.ndarray):
            np.testing.assert_allclose(got[k], expected[k], rtol=1e-12, atol=0, err_msg=k)
        else:
            assert expected[k] == pytest.approx(got[k], rel=1e-12), k
