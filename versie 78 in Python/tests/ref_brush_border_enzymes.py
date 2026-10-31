"""Brush-border enzymes: sucrase and maltase."""
import math

from digsim import state as st
from digsim.enums import Component as Comp
from digsim.enums import Segment

sucrase_activity = 0.0
maltase_activity = 0.0
A = 0.0
B = 0.0
C = 0.0
effective_concentration = 0.0
sucrose_absorption_rate = 0.0
sucrose_concentration_lumen = 0.0
sucrase_transformation_overall = 0.0
sucrose_concentration_brush = 0.0
maltose_absorption_rate = 0.0
maltose_concentration_lumen = 0.0
maltase_transformation_overall = 0.0
maltose_concentration_brush = 0.0

drag = 1  # Drag coefficient.
Km_sucrase = 8.56e-3  # g/ml, based on plots in Seri et al. (1996).
Km_maltase = 0.01  # g/ml.


def sucrase() -> None:
    """Hydrolyse sucrose to glucose and fructose at the brush border of each small-intestinal segment."""
    global sucrase_activity, A, B, C
    global sucrose_absorption_rate, sucrose_concentration_lumen
    global sucrase_transformation_overall, sucrose_concentration_brush

    for i in range(1, 16):  # Reset presence distribution before segment assignments.
        st.brush_enzymes.sucrase_presence[i] = 0
    st.brush_enzymes.sucrase_presence[Segment.DUODENUM_1] = 2
    st.brush_enzymes.sucrase_presence[Segment.DUODENUM_2] = 1.5
    st.brush_enzymes.sucrase_presence[Segment.JEJUNUM_1] = 1.2
    st.brush_enzymes.sucrase_presence[Segment.JEJUNUM_2] = 1
    st.brush_enzymes.sucrase_presence[Segment.JEJUNUM_3] = 1
    for seg in range(Segment.ILEUM_1, Segment.ILEUM_6 + 1):
        st.brush_enzymes.sucrase_presence[seg] = 0.5

    for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
        sucrase_activity = 0.15 * st.surface_factor[i]  # g/min per intestinal segment.
        st.brush_enzymes.inhibition_uncompetitive[i] = (  # Uncompetitive inhibition by arabinose and xylose.
            1
            + 1 * st.concentration_brush[i, Comp.ARABINOSE] / 6.85e-4
            + st.concentration_brush[i, Comp.XYLOSE] / 2e-5
            + 0 * st.concentration_brush[i, Comp.FRUCTOSE] / 6.85e-4
        )
        st.brush_enzymes.inhibition_competitive[i] = 1 + st.concentration_brush[i, Comp.ACARBOSE] / 1e-6  # Competitive inhibition by acarbose.
        for j in range(1, st.comp + 1):
            st.brush_enzymes.sucrase_transformation[i, j] = 0
        sucrose_absorption_rate = (
            st.absorption_rate[i, Comp.SUCROSE] + st.absorption_rate[i, Comp.TRACER_SUCROSE]
        )
        sucrose_concentration_lumen = (
            st.concentration_lumen[i, Comp.SUCROSE] + st.concentration_lumen[i, Comp.TRACER_SUCROSE]
        )

        st.brush_enzymes.sugar_flux = st.surface_factor[i] * sucrose_absorption_rate
        st.brush_enzymes.water_flux = st.surface_factor[i] * st.absorption_rate[i, Comp.WATER]
        A = -st.brush_enzymes.inhibition_uncompetitive[i]
        B = (
            -Km_sucrase * st.brush_enzymes.inhibition_competitive[i]
            + sucrose_concentration_lumen * st.brush_enzymes.inhibition_uncompetitive[i]
            + (st.delta[i] / st.diffusion_coefficient[Comp.SUCROSE])
            * (
                -st.brush_enzymes.sucrase_presence[i] * sucrase_activity
                - st.brush_enzymes.sugar_flux * st.brush_enzymes.inhibition_uncompetitive[i]
                + drag * sucrose_concentration_lumen * st.brush_enzymes.water_flux * st.brush_enzymes.inhibition_uncompetitive[i]
            )
        )
        C = (
            Km_sucrase * st.brush_enzymes.inhibition_competitive[i] * sucrose_concentration_lumen
            + (st.delta[i] / st.diffusion_coefficient[Comp.SUCROSE])
            * (
                -st.brush_enzymes.sugar_flux * Km_sucrase * st.brush_enzymes.inhibition_competitive[i]
                + drag * sucrose_concentration_lumen * Km_sucrase * st.brush_enzymes.inhibition_competitive[i] * st.brush_enzymes.water_flux
            )
        )
        sucrose_concentration_brush = (-B - math.sqrt(B * B - 4 * A * C)) / (2 * A)
        sucrase_transformation_overall = (
            (1 / st.surface_factor[i])
            * st.brush_enzymes.sucrase_presence[i]
            * sucrase_activity
            * sucrose_concentration_brush
            / (Km_sucrase * st.brush_enzymes.inhibition_competitive[i] + sucrose_concentration_brush * st.brush_enzymes.inhibition_uncompetitive[i])
        )
        if sucrose_concentration_brush > 0 and sucrose_concentration_lumen > 0:
            st.concentration_brush[i, Comp.SUCROSE] = (
                st.concentration_lumen[i, Comp.SUCROSE] / sucrose_concentration_lumen
            ) * sucrose_concentration_brush
            st.concentration_brush[i, Comp.TRACER_SUCROSE] = (
                st.concentration_lumen[i, Comp.TRACER_SUCROSE] / sucrose_concentration_lumen
            ) * sucrose_concentration_brush
            st.brush_enzymes.sucrase_transformation[i, Comp.SUCROSE] = (
                st.concentration_brush[i, Comp.SUCROSE] / sucrose_concentration_brush
            ) * sucrase_transformation_overall
            st.brush_enzymes.sucrase_transformation[i, Comp.TRACER_SUCROSE] = (
                st.concentration_brush[i, Comp.TRACER_SUCROSE] / sucrose_concentration_brush
            ) * sucrase_transformation_overall
            st.brush_enzymes.sucrase_inhibition[i] = (
                Km_sucrase * st.brush_enzymes.inhibition_competitive[i] + sucrose_concentration_brush * st.brush_enzymes.inhibition_uncompetitive[i]
            ) / (Km_sucrase + sucrose_concentration_brush)
        else:
            st.concentration_brush[i, Comp.SUCROSE] = 0
            st.concentration_brush[i, Comp.TRACER_SUCROSE] = 0
            st.brush_enzymes.sucrase_transformation[i, Comp.SUCROSE] = 0
            st.brush_enzymes.sucrase_transformation[i, Comp.TRACER_SUCROSE] = 0
            st.brush_enzymes.sucrase_inhibition[i] = 1


def maltase() -> None:
    """Hydrolyse maltose to glucose at the brush border of each small-intestinal segment."""
    global maltase_activity, A, B, C
    global maltose_absorption_rate, maltose_concentration_lumen
    global maltase_transformation_overall, maltose_concentration_brush

    # NOTE(port): Pascal reads the stale shared loop variable (left at 15 by sucrase's loop).
    stale_i = 15
    maltase_activity = 0.2 * st.surface_factor[stale_i]  # g/min per intestinal segment.
    st.brush_enzymes.inhibition_uncompetitive[stale_i] = 1 * (1 + st.concentration_brush[stale_i, Comp.XYLOSE] / 2e-5)
    st.brush_enzymes.inhibition_competitive[stale_i] = 1 + st.concentration_brush[stale_i, Comp.ACARBOSE] / 7e-7

    for i in range(1, 16):  # Reset presence distribution before segment assignments.
        st.brush_enzymes.maltase_presence[i] = st.brush_enzymes.sucrase_presence[i]  # Maltase follows the sucrase distribution.

    for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
        for j in range(1, st.comp + 1):
            st.brush_enzymes.maltase_transformation[i, j] = 0
        maltose_absorption_rate = (
            st.absorption_rate[i, Comp.MALTOSE] + st.absorption_rate[i, Comp.TRACER_MALTOSE]
        )
        maltose_concentration_lumen = (
            st.concentration_lumen[i, Comp.MALTOSE] + st.concentration_lumen[i, Comp.TRACER_MALTOSE]
        )

        st.brush_enzymes.sugar_flux = st.surface_factor[i] * maltose_absorption_rate
        st.brush_enzymes.water_flux = st.surface_factor[i] * st.absorption_rate[i, Comp.WATER]
        A = -st.brush_enzymes.inhibition_uncompetitive[i]
        B = (
            -Km_maltase * st.brush_enzymes.inhibition_competitive[i]
            + maltose_concentration_lumen * st.brush_enzymes.inhibition_uncompetitive[i]
            + (st.delta[i] / st.diffusion_coefficient[Comp.MALTOSE])
            * (
                -st.brush_enzymes.maltase_presence[i] * maltase_activity
                - st.brush_enzymes.sugar_flux * st.brush_enzymes.inhibition_uncompetitive[i]
                + drag * maltose_concentration_lumen * st.brush_enzymes.water_flux * st.brush_enzymes.inhibition_uncompetitive[i]
            )
        )
        C = (
            Km_maltase * st.brush_enzymes.inhibition_competitive[i] * maltose_concentration_lumen
            + (st.delta[i] / st.diffusion_coefficient[Comp.MALTOSE])
            * (
                -st.brush_enzymes.sugar_flux * Km_maltase * st.brush_enzymes.inhibition_competitive[i]
                + drag * maltose_concentration_lumen * Km_maltase * st.brush_enzymes.inhibition_competitive[i] * st.brush_enzymes.water_flux
            )
        )
        maltose_concentration_brush = (-B - math.sqrt(B * B - 4 * A * C)) / (2 * A)
        maltase_transformation_overall = (
            (1 / st.surface_factor[i])
            * st.brush_enzymes.maltase_presence[i]
            * maltase_activity
            * maltose_concentration_brush
            / (Km_maltase * st.brush_enzymes.inhibition_competitive[i] + maltose_concentration_brush * st.brush_enzymes.inhibition_uncompetitive[i])
        )
        if maltose_concentration_brush > 0 and maltose_concentration_lumen > 0:
            st.concentration_brush[i, Comp.MALTOSE] = (
                st.concentration_lumen[i, Comp.MALTOSE] / maltose_concentration_lumen
            ) * maltose_concentration_brush
            st.concentration_brush[i, Comp.TRACER_MALTOSE] = (
                st.concentration_lumen[i, Comp.TRACER_MALTOSE] / maltose_concentration_lumen
            ) * maltose_concentration_brush
            st.brush_enzymes.maltase_transformation[i, Comp.MALTOSE] = (
                st.concentration_brush[i, Comp.MALTOSE] / maltose_concentration_brush
            ) * maltase_transformation_overall
            st.brush_enzymes.maltase_transformation[i, Comp.TRACER_MALTOSE] = (
                st.concentration_brush[i, Comp.TRACER_MALTOSE] / maltose_concentration_brush
            ) * maltase_transformation_overall
            st.brush_enzymes.maltase_inhibition[i] = (
                Km_maltase * st.brush_enzymes.inhibition_competitive[i] + maltose_concentration_brush * st.brush_enzymes.inhibition_uncompetitive[i]
            ) / (Km_maltase + maltose_concentration_brush)
        else:
            st.concentration_brush[i, Comp.MALTOSE] = 0
            st.concentration_brush[i, Comp.TRACER_MALTOSE] = 0
            st.brush_enzymes.maltase_transformation[i, Comp.MALTOSE] = 0
            st.brush_enzymes.maltase_transformation[i, Comp.TRACER_MALTOSE] = 0


