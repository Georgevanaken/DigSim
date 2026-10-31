"""Main-program digestion and decomposition routines ported from Pascal."""

import numpy as np

from . import state as st
from .enums import Component as C
from .enums import Segment
from .main_part1 import _sequential_total, volumina
from .state import diet, stomach


def pepsin_units(pH: float, pepsin: float) -> float:
    """Return the pepsin activity (units) at the given pH."""
    if pH < 1:
        pepsin_act = 0 + 2 * (pH - 0.5) * (0.1 - 0)
    elif pH < 1.5:
        pepsin_act = 0.1 + 2 * (pH - 1) * (0.66 - 0.1)
    elif pH < 2:
        pepsin_act = 0.66 + 2 * (pH - 1.5) * (1 - 0.66)
    elif pH < 2.5:
        pepsin_act = 1 + 2 * (pH - 2) * (0.82 - 1)
    elif pH < 3:
        pepsin_act = 0.82 + 2 * (pH - 2.5) * (0.75 - 0.82)
    elif pH < 3.5:
        pepsin_act = 0.75 + 2 * (pH - 3) * (0.73 - 0.75)
    elif pH < 4:
        pepsin_act = 0.73 + 2 * (pH - 3.5) * (0.7 - 0.73)
    elif pH < 4.5:
        pepsin_act = 0.7 + 2 * (pH - 4) * (0.65 - 0.7)
    elif pH < 5:
        pepsin_act = 0.65 + 2 * (pH - 4.5) * (0.4 - 0.65)
    elif pH < 5.5:
        pepsin_act = 0.4 + 2 * (pH - 5) * (0.06 - 0.4)
    elif pH < 6:
        pepsin_act = 0.06 + 2 * (pH - 5.5) * (0.03 - 0.06)
    elif pH < 6.5:
        pepsin_act = 0.03 + 2 * (pH - 6) * (0.01 - 0.03)
    elif pH < 7:
        pepsin_act = 0.01 + 2 * (pH - 6.5) * (0.005 - 0.01)
    elif pH < 7.5:
        pepsin_act = 0.005 + 2 * (pH - 7) * (0 - 0.005)
    else:
        pepsin_act = 0
    return pepsin_act * pepsin * 4e6


def lingual_lipase_units(pH: float, lingual_lipase: float) -> float:
    """Return the lingual lipase activity (units) at the given pH."""
    if pH < 1:
        lingual_lipase_act = 0 + 2 * (pH - 0.5) * (0.1 - 0)
    elif pH < 1.5:
        lingual_lipase_act = 0.1 + 2 * (pH - 1) * (0.66 - 0.1)
    elif pH < 2:
        lingual_lipase_act = 0.66 + 2 * (pH - 1.5) * (1 - 0.66)
    elif pH < 2.5:
        lingual_lipase_act = 1 + 2 * (pH - 2) * (0.82 - 1)
    elif pH < 3:
        lingual_lipase_act = 0.82 + 2 * (pH - 2.5) * (0.75 - 0.82)
    elif pH < 3.5:
        lingual_lipase_act = 0.75 + 2 * (pH - 3) * (0.73 - 0.75)
    elif pH < 4:
        lingual_lipase_act = 0.73 + 2 * (pH - 3.5) * (0.7 - 0.73)
    elif pH < 4.5:
        lingual_lipase_act = 0.7 + 2 * (pH - 4) * (0.65 - 0.7)
    elif pH < 5:
        lingual_lipase_act = 0.65 + 2 * (pH - 4.5) * (0.4 - 0.65)
    elif pH < 5.5:
        lingual_lipase_act = 0.4 + 2 * (pH - 5) * (0.06 - 0.4)
    elif pH < 6:
        lingual_lipase_act = 0.06 + 2 * (pH - 5.5) * (0.03 - 0.06)
    elif pH < 6.5:
        lingual_lipase_act = 0.03 + 2 * (pH - 6) * (0.01 - 0.03)
    elif pH < 7:
        lingual_lipase_act = 0.01 + 2 * (pH - 6.5) * (0.005 - 0.01)
    elif pH < 7.5:
        lingual_lipase_act = 0.005 + 2 * (pH - 7) * (0 - 0.005)
    else:
        lingual_lipase_act = 0
    return lingual_lipase_act * lingual_lipase * 4e6


def gastric_lipase_units(pH: float, gastric_lipase: float) -> float:
    """Return the gastric lipase activity (units) at the given pH."""
    if pH < 1:
        gastric_lipase_act = 0 + 2 * (pH - 0.5) * (0.1 - 0)
    elif pH < 1.5:
        gastric_lipase_act = 0.1 + 2 * (pH - 1) * (0.66 - 0.1)
    elif pH < 2:
        gastric_lipase_act = 0.66 + 2 * (pH - 1.5) * (1 - 0.66)
    elif pH < 2.5:
        gastric_lipase_act = 1 + 2 * (pH - 2) * (0.82 - 1)
    elif pH < 3:
        gastric_lipase_act = 0.82 + 2 * (pH - 2.5) * (0.75 - 0.82)
    elif pH < 3.5:
        gastric_lipase_act = 0.75 + 2 * (pH - 3) * (0.73 - 0.75)
    elif pH < 4:
        gastric_lipase_act = 0.73 + 2 * (pH - 3.5) * (0.7 - 0.73)
    elif pH < 4.5:
        gastric_lipase_act = 0.7 + 2 * (pH - 4) * (0.65 - 0.7)
    elif pH < 5:
        gastric_lipase_act = 0.65 + 2 * (pH - 4.5) * (0.4 - 0.65)
    elif pH < 5.5:
        gastric_lipase_act = 0.4 + 2 * (pH - 5) * (0.06 - 0.4)
    elif pH < 6:
        gastric_lipase_act = 0.06 + 2 * (pH - 5.5) * (0.03 - 0.06)
    elif pH < 6.5:
        gastric_lipase_act = 0.03 + 2 * (pH - 6) * (0.01 - 0.03)
    elif pH < 7:
        gastric_lipase_act = 0.01 + 2 * (pH - 6.5) * (0.005 - 0.01)
    elif pH < 7.5:
        gastric_lipase_act = 0.005 + 2 * (pH - 7) * (0 - 0.005)
    else:
        gastric_lipase_act = 0
    return gastric_lipase_act * gastric_lipase * 4e6


_PEPSIN_EDGES = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5])
# Per pH range: a + 2 * (pH - c) * d, written with the same literal differences as pepsin_units.
_PEPSIN_ROWS = np.array([
    (0, 0.5, 0.1 - 0), (0.1, 1, 0.66 - 0.1), (0.66, 1.5, 1 - 0.66), (1, 2, 0.82 - 1), (0.82, 2.5, 0.75 - 0.82),
    (0.75, 3, 0.73 - 0.75), (0.73, 3.5, 0.7 - 0.73), (0.7, 4, 0.65 - 0.7), (0.65, 4.5, 0.4 - 0.65),
    (0.4, 5, 0.06 - 0.4), (0.06, 5.5, 0.03 - 0.06), (0.03, 6, 0.01 - 0.03), (0.01, 6.5, 0.005 - 0.01),
    (0.005, 7, 0 - 0.005), (0, 0, 0),
])


def _pepsin_units_array(pH: np.ndarray, pepsin: np.ndarray) -> np.ndarray:
    """Array version of pepsin_units (one element per segment)."""
    row = _PEPSIN_ROWS[np.searchsorted(_PEPSIN_EDGES, pH, side="right")]
    return (row[:, 0] + 2 * (pH - row[:, 1]) * row[:, 2]) * pepsin * 4e6


def _clamp(change: np.ndarray, available: np.ndarray) -> np.ndarray:
    """Limit the change to the available amount."""
    return np.where(change > available, available, change)


def _percent_loss(change: np.ndarray, fat: np.ndarray) -> float | int:
    """Sum over the segments of 100 * change / fat, counting only segments with fat above 1e-10."""
    has_fat = fat > 1e-10
    if not has_fat.any():
        return 0  # the original total stays an integer 0 (it is printed as "0")
    return _sequential_total(np.where(has_fat, 100 * change / np.where(has_fat, fat, 1.0), 0.0))


def digestion() -> None:
    """Enzymatic digestion (amylases, pepsin, proteases, lipases) in all segments for one time step.

    The segments are independent, so every step is applied to all segments at once; the steps themselves
    keep the order of the original per-segment code.
    """
    pepsin_susceptibility_casein = 0.01
    pepsin_susceptibility_protein = 0.001
    pancreatic_protease_susceptibility_casein = 1
    pancreatic_protease_susceptibility_protein = 1
    activity_pancreatic_protease = 1
    pepsin_susceptibility_casein_fragments = 0.001
    pancreatic_protease_susceptibility_casein_fragments = 1
    salivary_amylase_susceptibility_HDstarch = 1
    pancreatic_amylase_susceptibility_HDstarch = 1
    salivary_amylase_susceptibility_SDstarch = 0.2
    pancreatic_amylase_susceptibility_SDstarch = 0.2
    salivary_amylase_susceptibility_Rstarch = 0.0
    pancreatic_amylase_susceptibility_Rstarch = 0.0
    K_ms = 0.001
    K_cs = 50e-6
    K_as = 50e-6
    K_ncs = 50e-6
    K_mp = 50e-6
    K_cp = 50e-6
    K_ap = 50e-6
    K_ncp = 50e-6

    rows = slice(Segment.MOUTH, Segment.ILEUM_6 + 1)
    liq = st.liquid[rows]  # a view: assignments below update st.liquid
    volume = st.volume[rows]
    volume_squared = volume ** 2
    dt = st.clock.dt

    def lipolysis(constant: float, lipase: int) -> np.ndarray:
        """Triglyceride -> diglyceride + fatty acids; returns the amount converted per segment."""
        change = _clamp(constant * liq[:, lipase] * liq[:, C.FAT] * dt / volume_squared, liq[:, C.FAT])
        liq[:, C.DIGLYCERIDE] = liq[:, C.DIGLYCERIDE] + 0.666666 * change
        liq[:, C.FATTY_ACIDS] = liq[:, C.FATTY_ACIDS] + 0.333333 * change
        liq[:, C.FAT] = liq[:, C.FAT] - change
        return change

    change = lipolysis(st.lipolysis.lingual_cst, C.LINGUAL_LIPASE)
    st.lipolysis.total_lingual_rate = _sequential_total(change)
    st.lipolysis.percentual_loss_lingual = _percent_loss(change, liq[:, C.FAT])

    change = _clamp(st.lipolysis.lingual_cst * liq[:, C.LINGUAL_LIPASE] * liq[:, C.DIGLYCERIDE] * dt / volume_squared,
                    liq[:, C.DIGLYCERIDE])
    liq[:, C.FATTY_ACIDS] = liq[:, C.FATTY_ACIDS] + 0.5 * change
    liq[:, C.MONOGLYCERIDE] = liq[:, C.MONOGLYCERIDE] + 0.5 * change
    liq[:, C.DIGLYCERIDE] = liq[:, C.DIGLYCERIDE] - change

    change = _clamp(st.lipolysis.lingual_cst * liq[:, C.LINGUAL_LIPASE] * liq[:, C.MONOGLYCERIDE] * dt / volume_squared,
                    liq[:, C.MONOGLYCERIDE])
    liq[:, C.FATTY_ACIDS] = liq[:, C.FATTY_ACIDS] + change
    # NOTE(port): Preserve Pascal assignment to diglyceride component 7 during monoglyceride hydrolysis.
    liq[:, C.DIGLYCERIDE] = liq[:, C.MONOGLYCERIDE] - change

    change = lipolysis(st.lipolysis.gastric_cst, C.GASTRIC_LIPASE)
    st.lipolysis.total_gastric_rate = _sequential_total(change)
    st.lipolysis.percentual_loss_gastric = _percent_loss(change, liq[:, C.FAT])

    change = lipolysis(st.lipolysis.pancreatic_cst, C.PANCREATIC_LIPASE)
    st.lipolysis.total_pancreatic_rate = _sequential_total(change)
    st.lipolysis.percentual_loss_pancreatic = _percent_loss(change, liq[:, C.FAT])

    change = _clamp(st.lipolysis.pancreatic_cst * liq[:, C.PANCREATIC_LIPASE] * liq[:, C.DIGLYCERIDE] * dt / volume_squared,
                    liq[:, C.DIGLYCERIDE])
    liq[:, C.FATTY_ACIDS] = liq[:, C.FATTY_ACIDS] + 0.5 * change
    liq[:, C.MONOGLYCERIDE] = liq[:, C.MONOGLYCERIDE] + 0.5 * change
    liq[:, C.DIGLYCERIDE] = liq[:, C.DIGLYCERIDE] - change

    change = 100 * liq[:, C.FATTY_ACIDS] * liq[:, C.CONJUGATED_BILE_SALT] / (0.1 + liq[:, C.CONJUGATED_BILE_SALT]) * dt / volume_squared
    change = _clamp(_clamp(change, liq[:, C.FATTY_ACIDS]), liq[:, C.CONJUGATED_BILE_SALT])
    liq[:, C.BILE_SALT_LCFA_COMPLEX] = liq[:, C.BILE_SALT_LCFA_COMPLEX] + change
    liq[:, C.FATTY_ACIDS] = liq[:, C.FATTY_ACIDS] - change
    liq[:, C.CONJUGATED_BILE_SALT] = liq[:, C.CONJUGATED_BILE_SALT] - change

    change = 100 * liq[:, C.MONOGLYCERIDE] * liq[:, C.CONJUGATED_BILE_SALT] / (0.1 + liq[:, C.CONJUGATED_BILE_SALT]) * dt / volume_squared
    # NOTE(port): Pascal clamps micelle formation to fatty acid component 20, not monoglyceride 8.
    change = np.where(change > liq[:, C.MONOGLYCERIDE], liq[:, C.FATTY_ACIDS], change)
    change = _clamp(change, liq[:, C.CONJUGATED_BILE_SALT])
    liq[:, C.BILE_SALT_LCFA_COMPLEX] = liq[:, C.BILE_SALT_LCFA_COMPLEX] + change
    liq[:, C.MONOGLYCERIDE] = liq[:, C.MONOGLYCERIDE] - change
    liq[:, C.CONJUGATED_BILE_SALT] = liq[:, C.CONJUGATED_BILE_SALT] - change

    change = _clamp(0.1 * liq[:, C.BILE_PHOSPHOLIPIDS] * liq[:, C.PANCREATIC_PHOSPHOLIPASE] * dt / volume_squared,
                    liq[:, C.BILE_PHOSPHOLIPIDS])
    liq[:, C.BILE_PHOSPHOLIPIDS] = liq[:, C.BILE_PHOSPHOLIPIDS] - change
    liq[:, C.BILE_LYSOPHOSPHOLIPIDS] = liq[:, C.BILE_LYSOPHOSPHOLIPIDS] + change
    liq[:, C.FATTY_ACIDS] = liq[:, C.FATTY_ACIDS] + change

    # Left over from the Pascal global; holds the value of the last segment.
    st.segment_ph = st.ph[Segment.ILEUM_6]
    pepsin = _pepsin_units_array(st.ph[rows], liq[:, C.PEPSIN])

    def convert(source: int, target: int, rate: np.ndarray, clamp_to: np.ndarray | None = None) -> None:
        """Move `rate` of source to target, limited to the amount of source (or to `clamp_to`, as in Pascal)."""
        amount = _clamp(rate, liq[:, source]) if clamp_to is None else np.where(rate > liq[:, source], clamp_to, rate)
        liq[:, target] = liq[:, target] + amount
        liq[:, source] = liq[:, source] - amount

    convert(C.PROTEIN, C.SMALL_PEPTIDES, liq[:, C.PROTEIN] * pepsin_susceptibility_protein * pepsin * dt / volume)
    convert(C.BSA, C.SMALL_PEPTIDES, liq[:, C.BSA] * pepsin_susceptibility_protein * pepsin * dt / volume)
    convert(C.BLG, C.SMALL_PEPTIDES, liq[:, C.BLG] * pepsin_susceptibility_protein * pepsin * dt / volume)
    convert(C.CASEIN, C.CASEINPEP, liq[:, C.CASEIN] * pepsin_susceptibility_casein * pepsin * dt / volume)

    pancreatic = liq[:, C.PANCREATIC_PROTEASES]
    convert(C.PROTEIN, C.SMALL_PEPTIDES,
            liq[:, C.PROTEIN] * pancreatic * pancreatic_protease_susceptibility_protein * activity_pancreatic_protease * dt / volume)
    # NOTE(port): Pascal uses component 64 as the clamp for BSA/BLG digestion; preserve it.
    convert(C.BSA, C.SMALL_PEPTIDES,
            liq[:, C.BSA] * pancreatic * pancreatic_protease_susceptibility_protein * activity_pancreatic_protease * dt / volume,
            clamp_to=liq[:, C.CASEIN_PEPTIDES])
    convert(C.BLG, C.SMALL_PEPTIDES,
            liq[:, C.BLG] * pancreatic * pancreatic_protease_susceptibility_protein * activity_pancreatic_protease * dt / volume,
            clamp_to=liq[:, C.CASEIN_PEPTIDES])
    convert(C.CASEIN, C.CASEINPAN,
            liq[:, C.CASEIN] * pancreatic * pancreatic_protease_susceptibility_casein * activity_pancreatic_protease * dt / volume)
    convert(C.CASEINPEP, C.CASEIN_PEPTIDES,
            liq[:, C.CASEINPEP] * pepsin_susceptibility_casein_fragments * pepsin * dt / volume)
    convert(C.CASEINPAN, C.CASEIN_PEPTIDES,
            liq[:, C.CASEINPAN] * pancreatic * pancreatic_protease_susceptibility_casein_fragments * activity_pancreatic_protease * dt / volume)
    liq[:, C.SMALL_PEPTIDES] = liq[:, C.SMALL_PEPTIDES] + liq[:, C.CASEIN_PEPTIDES]
    liq[:, C.CASEIN_PEPTIDES] = 0

    conc = st.concentration_lumen[rows]
    ph = st.ph[rows]
    liq[:, C.SALIVARY_AMYLASE] = liq[:, C.SALIVARY_AMYLASE] - dt * liq[:, C.SALIVARY_AMYLASE] * np.exp(5 * np.log(5 / (5 + ph)))
    sugars = 1 + conc[:, C.I_A] / K_as + conc[:, C.GLUCOSE] / 0.1 + conc[:, C.MALTOSE] / 0.035
    st.salivary_amylase_activity[rows] = liq[:, C.SALIVARY_AMYLASE] / (
        (volume * K_ms * (1 + conc[:, C.I_C] / K_cs) + liq[:, C.HD_STARCH] * sugars) * (1 + conc[:, C.I_NC] / K_ncs))
    sugars = 1 + conc[:, C.I_A] / K_ap + conc[:, C.GLUCOSE] / 0.1 + conc[:, C.MALTOSE] / 0.035
    st.pancreatic_amylase_activity[rows] = liq[:, C.PANCREATIC_AMYLASE] / (
        (volume * K_mp * (1 + conc[:, C.I_C] / K_cp) + liq[:, C.HD_STARCH] * sugars) * (1 + conc[:, C.I_NC] / K_ncp))
    salivary_activity = st.salivary_amylase_activity[rows]
    pancreatic_activity = st.pancreatic_amylase_activity[rows]

    def starch(source: int, salivary_susceptibility: float, pancreatic_susceptibility: float) -> None:
        """Starch -> maltose (99 %) and tracer maltose (1 %)."""
        change = liq[:, source] * salivary_susceptibility * salivary_activity * dt
        change = change + liq[:, source] * pancreatic_susceptibility * pancreatic_activity * dt
        change = _clamp(change, liq[:, source])
        liq[:, C.MALTOSE] = liq[:, C.MALTOSE] + 0.99 * change
        liq[:, C.TRACER_MALTOSE] = liq[:, C.TRACER_MALTOSE] + 0.01 * change
        liq[:, source] = liq[:, source] - change

    starch(C.HD_STARCH, salivary_amylase_susceptibility_HDstarch, pancreatic_amylase_susceptibility_HDstarch)
    starch(C.SD_STARCH, salivary_amylase_susceptibility_SDstarch, pancreatic_amylase_susceptibility_SDstarch)
    starch(C.R_STARCH, salivary_amylase_susceptibility_Rstarch, pancreatic_amylase_susceptibility_Rstarch)

    volumina()


def _decompose_rows(rows: slice, n: int, dec_time: np.ndarray, dec_comp: np.ndarray, dec_frac: np.ndarray) -> None:
    """Decompose compound n in the given segment rows (rows are independent; same operations as the scalar loop)."""
    comp = st.compound
    column = comp[rows, n]
    if not column.any():
        return  # nothing of this compound present: increase would be 0 and all updates add 0
    comps = slice(1, st.comp + 1)
    increase = column * st.clock.dt / dec_time[n]
    aux = column - increase
    degcomp = dec_comp[n, 1]
    comp[rows, degcomp] = comp[rows, degcomp] + increase * dec_frac[n, 1] / 100
    degcomp = dec_comp[n, 2]
    comp[rows, degcomp] = comp[rows, degcomp] + increase * dec_frac[n, 2] / 100
    st.liquid[rows, comps] += (
        0.01 * (100 - dec_frac[n, 1] - dec_frac[n, 2]) * diet.compound_composition[n, comps] * increase[:, None]
    )
    comp[rows, n] = aux


def decomposition() -> None:
    """Decompose food compounds (oral, gastric, intestinal) into their components per the diet tables."""
    corpus_antrum = slice(Segment.CORPUS, Segment.ANTRUM + 1)
    mouth = slice(Segment.MOUTH, Segment.MOUTH + 1)
    intestine = slice(Segment.DUODENUM_1, Segment.ILEUM_6 + 1)
    for n in range(1, diet.number_of_compounds + 1):
        _decompose_rows(corpus_antrum, n, diet.gastric_dec_time, diet.gastric_dec_comp, diet.gastric_dec_frac)
        _decompose_rows(mouth, n, diet.oral_dec_time, diet.oral_dec_comp, diet.oral_dec_frac)

        if st.compound_fundus[n] != 0:
            increase = st.compound_fundus[n] * st.clock.dt / diet.gastric_dec_time[n]
            aux = st.compound_fundus[n] - increase
            degcomp = diet.gastric_dec_comp[n, 1]
            st.compound_fundus[degcomp] = st.compound_fundus[degcomp] + increase * diet.gastric_dec_frac[n, 1] / 100
            degcomp = diet.gastric_dec_comp[n, 2]
            st.compound_fundus[degcomp] = st.compound_fundus[degcomp] + increase * diet.gastric_dec_frac[n, 2] / 100
            stomach.fundus[1:st.comp + 1] += 0.01 * (100 - diet.gastric_dec_frac[n, 1] - diet.gastric_dec_frac[n, 2]) * diet.compound_composition[n, 1:st.comp + 1] * increase
            st.compound_fundus[n] = aux

        _decompose_rows(intestine, n, diet.intestinal_dec_time, diet.intestinal_dec_comp, diet.intestinal_dec_frac)



