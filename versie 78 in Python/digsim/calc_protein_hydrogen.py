"""Buffering of the luminal proteins: H+ binding of the mixture and the resulting pH per segment."""
import math

import numpy as np

from . import state as st
from .calc_blg_hydrogen import blg_hydrogen
from .calc_bsa_hydrogen import bsa_hydrogen
from .calc_casein_hydrogen import casein_hydrogen
from .calc_citric_acid_hydrogen import citric_acid_hydrogen
from .calc_mucin_hydrogen import mucin_hydrogen
from .enums import Component as C
from .enums import Segment
from .state import stomach


# Based on Curvale (2009), Buffer capacity of bovine serum albumin; estimates H+ added per gram of protein.
def protein_hydrogen(pH: float, Protein: float) -> float:
    """Return the H+ bound by generic protein (modelled on albumin) at the given pH."""
    correction_to_IEP = 0.0267  # IEP is 4.7.
    result = 0.0  # NOTE(port): uninitialised in Pascal when trunc(pH) is outside 1..12.
    case = int(pH)
    if case == 1:
        result = 0.14 - (pH - int(pH)) * 0.01
    elif case == 2:
        result = 0.13 - (pH - int(pH)) * 0.05
    elif case == 3:
        result = 0.08 - (pH - int(pH)) * 0.037
    elif case == 4:
        result = 0.043 - (pH - int(pH)) * 0.018
    elif case == 5:
        result = 0.025 - (pH - int(pH)) * 0.02
    elif case == 6:
        result = 0.005 - (pH - int(pH)) * 0.015
    elif case == 7:
        result = -0.01 - (pH - int(pH)) * 0.01
    elif case == 8:
        result = -0.02 - (pH - int(pH)) * 0.005
    elif case == 9:
        result = -0.025 - (pH - int(pH)) * 0.01
    elif case == 10:
        result = -0.035 - (pH - int(pH)) * 0.06
    elif case == 11:
        result = -0.095 - (pH - int(pH)) * 0.105
    elif case == 12:
        result = -0.20 - (pH - int(pH)) * 1
    result = Protein * (result - correction_to_IEP) / (15e-4 * 66463)  # Bound hydrogen ions relative to the protein IEP.
    return result


# Bound H+ for a protein mixture relative to the isoelectric points of the present proteins.
def calc_hydrogen(
    pH: float, casein: float, BSA: float, Mucin: float, BLG: float, Protein: float, citric_acid: float
) -> float:
    """Return the total H+ bound by the buffering components at the given pH."""
    # NOTE(port): Pascal var parameters are only read; no values need returning.
    result = (
        casein_hydrogen(pH, casein)
        + bsa_hydrogen(pH, BSA)
        + mucin_hydrogen(pH, Mucin)
        + blg_hydrogen(pH, BLG)
        + protein_hydrogen(pH, Protein)
        + citric_acid_hydrogen(pH, citric_acid)
    )
    return result


# Solve for pH at a given amount of hydrogen ions bound to a protein mixture.
def calc_ph(
    Volume: float, hydrogen: float, casein: float, BSA: float, Mucin: float, BLG: float, Protein: float,
    citric_acid: float,
) -> float:
    """Return the pH (bisection, tolerance 0.001) of a liquid of the given volume holding `hydrogen` H+."""
    # NOTE(port): Pascal var parameters are only read; no values need returning.
    Tol = 0.001
    Min = 0
    Max = 10 + Tol
    while Max - Min > Tol:
        pH = (Max + Min) / 2
        Guess_hydrogen = calc_hydrogen(
            pH, casein, BSA, Mucin, BLG, Protein, citric_acid
        )
        if Guess_hydrogen > hydrogen - Volume * 0.001 * (
            math.exp(-2.302585093 * pH) - 0.0000001
        ):
            Min = pH
        else:
            Max = pH
    result = (Max + Min) / 2
    return result


# Vectorised versions of the pH tables above. Every element goes through exactly the same floating-point
# operations, in the same order, as the scalar functions (tests/ref_ph.py checks this bit for bit).
def _table(a: list[float], b: list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Coefficients (index 0 unused: value 0, as the scalar code leaves `result` at 0.0 for trunc(pH) = 0)."""
    return np.array([0.0] + a), np.array([0.0] + b)


_PROTEIN_TABLE = _table(
    [0.14, 0.13, 0.08, 0.043, 0.025, 0.005, -0.01, -0.02, -0.025, -0.035, -0.095, -0.20],
    [0.01, 0.05, 0.037, 0.018, 0.02, 0.015, 0.01, 0.005, 0.01, 0.06, 0.105, 1],
)
_BLG_TABLE = _table(
    [0.002173913, 0.002081522, 0.00169837, 0.000434783, -0.000108696, -0.000380435, -0.000641304,
     -0.000869565, -0.000978261, -0.1032609, -0.1032609, -0.1032609],
    [9.23913e-05, 0.000383152, 0.001263587, 0.000543478, 0.000271739, 0.00026087, 0.000228261,
     0.000108696, 5.43478e-05, 0, 0, 0],
)
_MUCIN_TABLE = _table(
    [0.001333333, 0.000766667, 0.0003, 0.00014, 0.00001333, -0.0000333, -0.0000633, -0.0000733,
     -0.000103333, -0.000146667, -0.000333333, -0.001],
    [0.000566667, 0.000466667, 0.00016, 0.00012667, 4.663e-05, 0.00003, 0.00001, 0.000030033,
     (4.3334 - 5), 0.000186666, 0.000333333, 0.000666667],
)
_CITRIC_TABLE = _table(
    [3, 3, 2.5, 2, 1.4, 0.8, 0.3, 0, 0, 0, 0, 0],
    [0, 0.5, 0.5, 0.6, 0.6, 0.5, 0.3, 0, 0, 0, 0, 0],
)


# One row per trunc(pH): columns protein/BSA, BLG, mucin, citric acid. A single fancy-index lookup then serves
# all four tables.
_TABLES_A = np.stack([_PROTEIN_TABLE[0], _BLG_TABLE[0], _MUCIN_TABLE[0], _CITRIC_TABLE[0]], axis=1)
_TABLES_B = np.stack([_PROTEIN_TABLE[1], _BLG_TABLE[1], _MUCIN_TABLE[1], _CITRIC_TABLE[1]], axis=1)


_CASEIN_EDGES = np.array([2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9])
# Per pH range: a + (c - pH) * k. Doubling k is exact, so (c-pH)*d*2 == (c-pH)*(d*2) bit for bit.
_CASEIN_A = np.array([0.0016, 0.0009, 0.0004, 0.0001, 0.00005, 0.00002, 0, -0.00002, -0.0003, -0.0007, -0.0008, -0.0011])
_CASEIN_C = np.array([2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10])
_CASEIN_K = np.array([
    0.002 - 0.0016, (0.0016 - 0.0009) * 2, (0.0009 - 0.0004) * 2, (0.0004 - 0.0001) * 2,
    (0.0001 - 0.00005) * 2, (0.00005 - 0.00002) * 2, (0.00002 - 0) * 2, 0 + 0.00002,
    -0.00002 + 0.0003, -0.0003 + 0.0007, -0.0007 + 0.0008, -0.0008 + 0.0011,
])


# Casein coefficients per pH range, one row each: a, c, k.
_CASEIN_ROWS = np.stack([_CASEIN_A, _CASEIN_C, _CASEIN_K], axis=1)


def calc_hydrogen_array(
    pH: np.ndarray, casein: np.ndarray, BSA: np.ndarray, Mucin: np.ndarray, BLG: np.ndarray,
    Protein: np.ndarray, citric_acid: np.ndarray,
) -> np.ndarray:
    """Array version of calc_hydrogen (one element per segment); same operations and summation order."""
    idx = pH.astype(np.intp)
    frac = pH - idx
    t = _TABLES_A[idx] - frac[:, None] * _TABLES_B[idx]
    cc = _CASEIN_ROWS[np.searchsorted(_CASEIN_EDGES, pH, side="right")]
    bound = cc[:, 0] + (cc[:, 1] - pH) * cc[:, 2]
    protein_scaled = (t[:, 0] - 0.0267)
    casein_term = casein * (bound + -0.00002)
    bsa_term = BSA * protein_scaled / (15e-4 * 66463)
    mucin_term = Mucin * (t[:, 2] - 0.0005)
    blg_term = BLG * t[:, 1]
    protein_term = Protein * protein_scaled / (15e-4 * 66463)
    citric_term = citric_acid * -t[:, 3] / 192.124
    return casein_term + bsa_term + mucin_term + blg_term + protein_term + citric_term


def _count_bisection_steps(tol: float = 0.001) -> int:
    """Iterations of the interval halving; independent of the data, since the width always halves."""
    low, high, steps = 0.0, 10 + tol, 0
    while high - low > tol:
        high = (high + low) / 2
        steps += 1
    return steps


_BISECTION_STEPS = _count_bisection_steps()


def calc_ph_array(
    Volume: np.ndarray, hydrogen: np.ndarray, casein: np.ndarray, BSA: np.ndarray, Mucin: np.ndarray,
    BLG: np.ndarray, Protein: np.ndarray, citric_acid: np.ndarray,
) -> np.ndarray:
    """Array version of calc_ph: the same bisection for all elements at once.

    The interval halves every iteration whatever the comparison outcome, so all elements need the same
    number of iterations. math.exp is used (not np.exp) to get the same last bit as the scalar code.
    """
    Tol = 0.001
    n = len(Volume)
    low = np.zeros(n)
    high = np.full(n, 10 + Tol)
    volume_scaled = Volume * 0.001
    exp = math.exp
    for _ in range(_BISECTION_STEPS):
        pH = (high + low) / 2
        guess = calc_hydrogen_array(pH, casein, BSA, Mucin, BLG, Protein, citric_acid)
        decay = np.fromiter((exp(-2.302585093 * x) for x in pH.tolist()), dtype=float, count=n)
        raised = guess > hydrogen - volume_scaled * (decay - 0.0000001)
        low = np.where(raised, pH, low)
        high = np.where(raised, high, pH)
    return (high + low) / 2


def calc_all_ph() -> None:
    """Update st.ph for every segment from mouth to ileum 6 and the fundus (all solved together)."""
    rows = slice(Segment.MOUTH, Segment.ILEUM_6 + 1)
    # The fundus is appended as one extra element.
    liq = np.vstack((st.liquid[rows], stomach.fundus[None, :]))
    volume = np.append(st.volume[rows], stomach.v_fundus)
    casein = (
        liq[:, C.CASEIN] + liq[:, C.CASEINPEP] + liq[:, C.CASEINPAN] + liq[:, C.CASEIN_PEPTIDES]
        + liq[:, C.SMALL_PEPTIDES]
    )
    result = calc_ph_array(
        volume, liq[:, C.H_ION], casein, liq[:, C.BSA], liq[:, C.MUCUS_PROTEIN], liq[:, C.BLG],
        liq[:, C.PROTEIN], liq[:, C.CITRIC_ACID],
    )
    st.ph[rows] = result[:-1]
    stomach.ph_fundus = float(result[-1])



