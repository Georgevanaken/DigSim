"""Regression tests against a Free Pascal build of DigSim v78 (20 steps/min, tests/data/meal.txt).

The Pascal baselines are in tests/baselines. The simulations are slow (the 500-minute
initialisation alone takes about 7 minutes), so they only run with DIGSIM_SLOW=1:

    DIGSIM_SLOW=1 python -m pytest tests/test_regression.py -q

The tolerance is 5-6% of each column's range: the pH bisection (Tol 0.001) and x87 versus
double arithmetic make pulse-like and discrete columns differ locally, so only physiologically
meaningful columns are compared.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from digsim import state as st
from digsim.diet_from_file import read_compounds

HERE = Path(__file__).parent
BASELINES = HERE / "baselines"
NCOLS = 166

KEY_COLUMNS = [
    "corpus", "antrum", "duodenum 1", "jejunum 1", "ileum 1", "colon", "total volume SI",
    "CCK", "PYY", "gastric emptying", "calories proximal", "calories distal", "calories colon",
    "calories absorbed", "fullness", "hunger", "GIP", "GLP1", "glucose absorbed",
    "fatty acids absorbed", "pH upper gastric", "pH antrum", "serum glucose (mMol/l)",
    "serum insulin (pMol/l)", "Gastric Volume (mL)", "HbA1C (%)", "Ghrelin", "Motilin",
    "glycogen (g)",
]

slow = pytest.mark.skipif(not os.environ.get("DIGSIM_SLOW"), reason="set DIGSIM_SLOW=1 to run")


def load_results(path):
    lines = Path(path).read_text().splitlines()
    header = [h.strip() for h in lines[0].split(",")]
    rows = [[float(x) for x in l.replace(",", " ").split()[:NCOLS]] for l in lines[1:] if l.strip()]
    return header, np.array(rows)


def compare(baseline_file, result_rows, tol, row_slice=slice(None)):
    header, base = load_results(BASELINES / baseline_file)
    idx = [header.index(k) for k in KEY_COLUMNS]
    scale = np.maximum(np.abs(base).max(axis=0), 1e-6)
    err = np.abs(result_rows[:, idx] - base[:, idx]) / scale[idx]
    worst = [(KEY_COLUMNS[j], round(float(err[:, j].max()), 3)) for j in np.argsort(-err.max(axis=0))[:5]]
    assert err.max() < tol, f"max error {err.max():.3f} of column range; worst: {worst}"


def test_baseline_files_are_consistent():
    for name in ("pascal_n20_head.txt", "pascal_n20_sampled.txt"):
        header, rows = load_results(BASELINES / name)
        assert len(header) >= NCOLS
        assert rows.shape[1] == NCOLS
        assert all(k in header for k in KEY_COLUMNS)
    _, head = load_results(BASELINES / "pascal_n20_head.txt")
    assert head.shape[0] == 40 and head[0, 0] == 1.0


def test_compoundreader_out_of_range_compound_gives_zero_ph_and_viscosity():
    # meal.txt contains compound number 78 (> 50 array size); Pascal reads 0 beyond the array.
    d = st.diet
    saved = (d.new_meal_number, d.meal[1], d.dish[1], d.meal_amount[1], d.eating_speed[1])
    d.new_meal_number = 1
    d.meal[1], d.dish[1], d.meal_amount[1], d.eating_speed[1] = 78, 1, 1.0, 10.0
    d.ph_dish[1] = 7.0
    d.liquid_viscosity_dish[1] = 5.0
    before = st.compound_dish[1, 78]
    read_compounds()
    try:
        assert d.ph_dish[1] == 0.0
        assert d.liquid_viscosity_dish[1] == 0.0
        assert st.compound_dish[1, 78] == before + 1.0
    finally:
        st.compound_dish[1, 78] = before
        d.new_meal_number, d.meal[1], d.dish[1], d.meal_amount[1], d.eating_speed[1] = saved


def _run(tmp_path, minutes_max=None):
    # A fresh process per run: the simulation state lives in module globals.
    out = tmp_path / "Results.txt"
    cmd = [sys.executable, "-m", "digsim.digsim_v78", "-n", "20", "--results", str(out)]
    if minutes_max is not None:
        cmd += ["--minutes-max", str(minutes_max)]
    env = dict(os.environ, PYTHONPATH=str(HERE.parent))
    subprocess.run(cmd, cwd=HERE / "data", env=env, check=True, capture_output=True)
    return load_results(out)[1]


@slow
def test_first_40_minutes_match_pascal(tmp_path):
    rows = _run(tmp_path, minutes_max=40)
    assert rows.shape[0] >= 40
    compare("pascal_n20_head.txt", rows[:40], tol=0.05)


@slow
def test_full_run_matches_pascal_sampled(tmp_path):
    rows = _run(tmp_path)
    assert rows.shape[0] == 3502
    compare("pascal_n20_sampled.txt", rows[0::50], tol=0.06)
