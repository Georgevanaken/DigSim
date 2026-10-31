"""Brush-border enzymes: sucrase and maltase (all small-intestinal segments are solved together)."""
import numpy as np

from . import state as st
from .enums import Component as Comp
from .enums import Segment

drag = 1  # Drag coefficient.
Km_sucrase = 8.56e-3  # g/ml, based on plots in Seri et al. (1996).
Km_maltase = 0.01  # g/ml.

_ROWS = slice(Segment.DUODENUM_1, Segment.ILEUM_6 + 1)


def _hydrolysis(
    substrate: int, tracer: int, Km: float, presence: np.ndarray, activity: np.ndarray | float,
    inhibition_competitive: np.ndarray, inhibition_uncompetitive: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve the brush-border concentration of a disaccharide (quadratic) for every segment.

    Writes the brush concentrations and returns (substrate transformation, tracer transformation,
    inhibition, present). Where there is no substrate the transformations are 0 and `present` is False.
    The expressions keep the operation order of the original per-segment code.
    """
    surface = st.surface_factor[_ROWS]
    delta = st.delta[_ROWS]
    D = st.diffusion_coefficient[substrate]
    absorption_rate = st.absorption_rate[_ROWS, substrate] + st.absorption_rate[_ROWS, tracer]
    lumen_sub = st.concentration_lumen[_ROWS, substrate]
    lumen_tracer = st.concentration_lumen[_ROWS, tracer]
    lumen = lumen_sub + lumen_tracer

    sugar_flux = surface * absorption_rate
    water_flux = surface * st.absorption_rate[_ROWS, Comp.WATER]
    st.brush_enzymes.sugar_flux = float(sugar_flux[-1])  # leftovers of the last segment, as in the original
    st.brush_enzymes.water_flux = float(water_flux[-1])

    A = -inhibition_uncompetitive
    B = (
        -Km * inhibition_competitive
        + lumen * inhibition_uncompetitive
        + (delta / D)
        * (
            -presence * activity
            - sugar_flux * inhibition_uncompetitive
            + drag * lumen * water_flux * inhibition_uncompetitive
        )
    )
    C = (
        Km * inhibition_competitive * lumen
        + (delta / D)
        * (
            -sugar_flux * Km * inhibition_competitive
            + drag * lumen * Km * inhibition_competitive * water_flux
        )
    )
    brush = (-B - np.sqrt(B * B - 4 * A * C)) / (2 * A)
    overall = (
        (1 / surface)
        * presence
        * activity
        * brush
        / (Km * inhibition_competitive + brush * inhibition_uncompetitive)
    )
    present = (brush > 0) & (lumen > 0)
    safe_lumen = np.where(present, lumen, 1.0)
    safe_brush = np.where(present, brush, 1.0)
    brush_sub = np.where(present, (lumen_sub / safe_lumen) * brush, 0.0)
    brush_tracer = np.where(present, (lumen_tracer / safe_lumen) * brush, 0.0)
    st.concentration_brush[_ROWS, substrate] = brush_sub
    st.concentration_brush[_ROWS, tracer] = brush_tracer
    transformation_sub = np.where(present, (brush_sub / safe_brush) * overall, 0.0)
    transformation_tracer = np.where(present, (brush_tracer / safe_brush) * overall, 0.0)
    inhibition = (Km * inhibition_competitive + safe_brush * inhibition_uncompetitive) / (Km + safe_brush)
    return transformation_sub, transformation_tracer, inhibition, present


def sucrase() -> None:
    """Hydrolyse sucrose to glucose and fructose at the brush border of each small-intestinal segment."""
    st.brush_enzymes.sucrase_presence[1:16] = 0  # Reset presence distribution before segment assignments.
    st.brush_enzymes.sucrase_presence[Segment.DUODENUM_1] = 2
    st.brush_enzymes.sucrase_presence[Segment.DUODENUM_2] = 1.5
    st.brush_enzymes.sucrase_presence[Segment.JEJUNUM_1] = 1.2
    st.brush_enzymes.sucrase_presence[Segment.JEJUNUM_2] = 1
    st.brush_enzymes.sucrase_presence[Segment.JEJUNUM_3] = 1
    st.brush_enzymes.sucrase_presence[Segment.ILEUM_1:Segment.ILEUM_6 + 1] = 0.5

    activity = 0.15 * st.surface_factor[_ROWS]  # g/min per intestinal segment.
    conc = st.concentration_brush[_ROWS]
    # Uncompetitive inhibition by arabinose and xylose; competitive inhibition by acarbose.
    uncompetitive = (
        1
        + 1 * conc[:, Comp.ARABINOSE] / 6.85e-4
        + conc[:, Comp.XYLOSE] / 2e-5
        + 0 * conc[:, Comp.FRUCTOSE] / 6.85e-4
    )
    competitive = 1 + conc[:, Comp.ACARBOSE] / 1e-6
    st.brush_enzymes.inhibition_uncompetitive[_ROWS] = uncompetitive
    st.brush_enzymes.inhibition_competitive[_ROWS] = competitive
    st.brush_enzymes.sucrase_transformation[_ROWS, 1:st.comp + 1] = 0

    sucrose, tracer, inhibition, present = _hydrolysis(
        Comp.SUCROSE, Comp.TRACER_SUCROSE, Km_sucrase, st.brush_enzymes.sucrase_presence[_ROWS], activity,
        competitive, uncompetitive,
    )
    st.brush_enzymes.sucrase_transformation[_ROWS, Comp.SUCROSE] = sucrose
    st.brush_enzymes.sucrase_transformation[_ROWS, Comp.TRACER_SUCROSE] = tracer
    st.brush_enzymes.sucrase_inhibition[_ROWS] = np.where(present, inhibition, 1)


def maltase() -> None:
    """Hydrolyse maltose to glucose at the brush border of each small-intestinal segment."""
    # NOTE(port): Pascal reads the stale shared loop variable (left at 15 by sucrase's loop).
    stale_i = 15
    maltase_activity = 0.2 * st.surface_factor[stale_i]  # g/min per intestinal segment.
    st.brush_enzymes.inhibition_uncompetitive[stale_i] = 1 * (1 + st.concentration_brush[stale_i, Comp.XYLOSE] / 2e-5)
    st.brush_enzymes.inhibition_competitive[stale_i] = 1 + st.concentration_brush[stale_i, Comp.ACARBOSE] / 7e-7

    st.brush_enzymes.maltase_presence[1:16] = st.brush_enzymes.sucrase_presence[1:16]  # Maltase follows the sucrase distribution.
    st.brush_enzymes.maltase_transformation[_ROWS, 1:st.comp + 1] = 0

    maltose, tracer, inhibition, present = _hydrolysis(
        Comp.MALTOSE, Comp.TRACER_MALTOSE, Km_maltase, st.brush_enzymes.maltase_presence[_ROWS], maltase_activity,
        st.brush_enzymes.inhibition_competitive[_ROWS], st.brush_enzymes.inhibition_uncompetitive[_ROWS],
    )
    st.brush_enzymes.maltase_transformation[_ROWS, Comp.MALTOSE] = maltose
    st.brush_enzymes.maltase_transformation[_ROWS, Comp.TRACER_MALTOSE] = tracer
    # As in the original, the maltase inhibition is left unchanged where no maltose is present.
    st.brush_enzymes.maltase_inhibition[_ROWS] = np.where(present, inhibition, st.brush_enzymes.maltase_inhibition[_ROWS])
