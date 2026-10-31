"""Close translation of the opening procedures in the DigSim v78 main program."""

import math

import numpy as np

from . import state as st
from .calc_mucin_hydrogen import mucin_hydrogen
from .compartment_properties import volumina
from .enums import Component as C
from .enums import MMCStatus, Segment
from .gastric_motility import gastric_redistribution
from .glucose_homeostasis2 import initialise_glucose_dynamics
from .state import bile, diet, glucose, hormones, stomach


def initialise() -> None:
    """Define the starting conditions for the model calculations."""
    st.mmc_status = MMCStatus.FED
    st.mmc_time[0] = 0
    st.mmc_time[1] = 0
    st.mmc_time[2] = 0
    st.mmc_time[3] = 0
    st.mmc_time[4] = 0
    stomach.emptying_cycle = 0
    glucose.glycogen_store0 = 10
    glucose.glycogen_store1 = 0
    glucose.glycogen_store2 = 0
    glucose.glycogen_store = 100
    st.cycle = 0
    bile.blood_salt_pool = 0.008
    stomach.tone = 0  # Gastric pressure initially set to 0.
    st.osmotic.lumen_duodenum = 0.05
    st.secretion.bicarbonate = 0
    stomach.corpus_viscosity = 10
    stomach.antral_viscosity = 10

    for dish in range(1, 6):
        st.liquid[dish, C.H_ION] = 1e-7
    st.calories.colon = 0
    st.calories.absorbed = 0
    for i in range(Segment.MOUTH, Segment.COLON + 1):
        # Set non-water liquid components to zero; keep a small water amount to avoid division by zero.
        st.liquid[i, C.WATER] = 1
        st.liquid_next[i, C.WATER] = 1
        stomach.fundus[C.WATER] = 1
        st.absorption_rate[i, C.WATER] = 0
        st.concentration_brush[i, C.WATER] = 1
        for n in range(1, diet.number_of_compounds + 1):
            st.compound[i, n] = 0
            st.compound_fundus[n] = 0
        for j in range(2, st.comp + 1):
            st.liquid[i, j] = 0
            st.liquid_next[i, j] = 0
            stomach.fundus[j] = 0
            for s in range(1, 100):
                st.solids[s, i, j] = 0
                st.solids_next[s, i, j] = 0
            st.absorption_rate[i, j] = 0
            st.enterocyte[i, j] = 0
            st.concentration_brush[i, j] = 0

    st.liquid[Segment.CORPUS, C.WATER] = 25
    st.liquid[Segment.ANTRUM, C.WATER] = 10
    for j in range(1, st.comp + 1):
        stomach.solids_corpus[j] = 0
        stomach.solids_antrum[j] = 0
        stomach.solids_fundus[j] = 0
    stomach.solids_in_fundus = 0
    stomach.solids_in_corpus = 0
    stomach.solids_in_antrum = 0
    stomach.total_solids_in_stomach = 0
    for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
        st.liquid[i, C.PANCREATIC_LIPASE] = 0.01
    st.liquid[Segment.DUODENUM_2, C.CONJUGATED_BILE_SALT] = 0.1
    st.liquid[Segment.DUODENUM_2, C.PANCREATIC_LIPASE] = 0.012 / 5

    initialise_glucose_dynamics()

    for j in range(1, st.comp + 1):
        for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
            st.mesenteric[i, j] = 0
        st.blood[j] = 0
    for i in range(Segment.MOUTH, Segment.ILEUM_6 + 1):
        hormones.mesenteric_cck[i] = 0
        hormones.mesenteric_pyy[i] = 0
        hormones.mesenteric_gip[i] = 0
        hormones.mesenteric_glp1[i] = 0
        hormones.mesenteric_ghrelin[i] = 0
        hormones.mesenteric_motilin[i] = 0
    glucose.hba1c = 0


def _sequential_total(terms: np.ndarray) -> float:
    """Left-to-right sum of a 1-D array, starting from 0 (np.sum would round differently)."""
    return float(np.cumsum(np.concatenate(([0.0], terms)))[-1])


def caloric_densities_segments() -> None:
    """Calculate caloric values of luminal contents and brush borders."""
    segs = slice(Segment.DUODENUM_1, Segment.COLON + 1)
    comps = slice(1, st.comp + 1)
    comps_n = slice(1, diet.number_of_compounds + 1)
    terms = np.concatenate(
        (
            np.zeros((Segment.COLON - Segment.DUODENUM_1 + 1, 1)),
            st.liquid[segs, comps] * st.calories.value[comps],
            st.compound[segs, comps_n] * st.calories.value[comps_n],
        ),
        axis=1,
    )
    # cumsum along the rows adds left to right, as the original loops did.
    calories_segment = np.cumsum(terms, axis=1)[:, -1]
    volume = st.volume[segs]
    big = volume > 1e-10
    st.calories.density_local[segs] = np.where(big, calories_segment / np.where(big, volume, 1), 0)
    # The original recomputes this for every segment; only the value for the last one (colon) survives.
    st.calories.brush = _sequential_total(st.concentration_brush[Segment.COLON, comps] * st.calories.value[comps])


def _osmotic_value(concentration: np.ndarray, i: int) -> float:
    """Small sugars and salt as equivalent glucose concentration (terms added in the original order)."""
    return (
        (180.156 / 58.44) * concentration[i, C.NACL]
        + concentration[i, C.GLUCOSE]
        + concentration[i, C.TRACER_GLUCOSE_1]
        + concentration[i, C.TRACER_GLUCOSE_2]
        + 0.5 * concentration[i, C.SUCROSE]
        + 0.5 * concentration[i, C.TRACER_SUCROSE]
        + concentration[i, C.FRUCTOSE]
        + 0.5 * concentration[i, C.MALTOSE]
        + 0.5 * concentration[i, C.TRACER_MALTOSE]
        + concentration[i, C.ARABINOSE]
        + concentration[i, C.XYLOSE]
        + 0.2 * concentration[i, C.ACARBOSE]
    )


def osmotic_values_lumen() -> None:
    """Calculate luminal osmotic values as equivalent glucose concentrations."""
    for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
        st.osmotic.lumen[i] = _osmotic_value(st.concentration_lumen, i)

    ov = st.osmotic.lumen
    st.osmotic.lumen_duodenum = (
        ov[Segment.DUODENUM_1]
        + 1 * ov[Segment.DUODENUM_2]
        + 0 * ov[Segment.JEJUNUM_1]
        + 0 * ov[Segment.JEJUNUM_2]
    )
    st.osmotic.lumen_distal = (
        ov[Segment.ILEUM_1]
        + 1 * ov[Segment.ILEUM_2]
        + ov[Segment.ILEUM_3]
        + ov[Segment.ILEUM_4]
        + ov[Segment.ILEUM_5]
        + ov[Segment.ILEUM_6]
    )


def osmotic_values_brush() -> None:
    """Calculate brush-border osmotic values as equivalent glucose concentrations."""
    aux = 0.0
    decelerator = 0.0
    st.osmotic.brush_duodenum = 0
    for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
        st.osmotic.brush[i] = _osmotic_value(st.concentration_brush, i)
    decelerator = 0 * 60 * st.clock.dt
    ov = st.osmotic.brush
    aux = (
        1 * ov[Segment.DUODENUM_1]
        + 1 * ov[Segment.DUODENUM_2]
        + 0.5 * ov[Segment.JEJUNUM_1]
        + 0.2 * ov[Segment.JEJUNUM_2]
    )
    st.osmotic.brush_duodenum = (
        st.osmotic.brush_duodenum * decelerator + aux
    ) / (decelerator + 1)
    st.osmotic.brush_distal = (
        ov[Segment.ILEUM_1]
        + 1 * ov[Segment.ILEUM_2]
        + ov[Segment.ILEUM_3]
        + ov[Segment.ILEUM_4]
        + ov[Segment.ILEUM_5]
        + ov[Segment.ILEUM_6]
    )


def calories() -> None:
    """Calculate calories in the proximal small bowel, distal small bowel, and colon."""
    st.calories.proximal = 0
    st.calories.distal = 0
    st.calories.brush = 0
    st.calories.colon = 0
    proximal = slice(Segment.DUODENUM_1, Segment.JEJUNUM_2 + 1)
    distal = slice(Segment.JEJUNUM_3, Segment.ILEUM_6 + 1)
    comps = slice(1, st.comp + 1)
    comps_n = slice(1, diet.number_of_compounds + 1)

    def total(rows: slice) -> float:
        # Original order: component by component, segments inside; then compounds the same way.
        liquid_terms = (st.calories.value[comps] * st.liquid[rows, comps]).T.ravel()
        compound_terms = (st.compound[rows, comps_n] * st.calories.compound[comps_n]).T.ravel()
        return _sequential_total(np.concatenate((liquid_terms, compound_terms)))

    st.calories.proximal = total(proximal)
    st.calories.distal = total(distal)
    colon = slice(Segment.COLON, Segment.COLON + 1)
    st.calories.colon = total(colon)


def time_of_day(minutes: float, referencetime: float) -> float:
    """Return the fixed time used for averaged unstimulated saliva secretion."""
    return 11 * 60


def saliva_secretion() -> None:
    """Add saliva (water, mucin, amylase, lingual lipase) to the mouth for one time step."""
    mucus_content_in_saliva = 0.003
    amylase_content_in_saliva = 0.05
    lingual_lipase_content_in_saliva = 0.01

    aux1 = 0.1 + 0.4 * math.sin(
        6.283 * (time_of_day(st.clock.minutes, st.clock.reference_time) - 11 * 60) / (24 * 60)
    )
    mouth = Segment.MOUTH
    aux2 = 5 * st.volume[mouth] / (1 + st.volume[mouth])
    aux3 = 4 * abs(6 - st.ph[mouth]) / (1 + abs(6 - st.ph[mouth]))
    aux4 = 100 * st.concentration_lumen[mouth, C.NACL]
    st.secretion.saliva_addition = (aux1 + aux2 + aux3 + aux4) * st.clock.dt
    st.liquid[mouth, C.WATER] = st.liquid[mouth, C.WATER] + st.secretion.saliva_addition
    mucus_added = mucus_content_in_saliva * st.secretion.saliva_addition
    amylase_added = amylase_content_in_saliva * st.secretion.saliva_addition
    lingual_lipase_added = lingual_lipase_content_in_saliva * st.secretion.saliva_addition
    st.liquid[mouth, C.MUCUS_PROTEIN] = st.liquid[mouth, C.MUCUS_PROTEIN] + mucus_added
    st.liquid[mouth, C.SALIVARY_AMYLASE] = st.liquid[mouth, C.SALIVARY_AMYLASE] + amylase_added
    st.liquid[mouth, C.LINGUAL_LIPASE] = st.liquid[mouth, C.LINGUAL_LIPASE] + lingual_lipase_added
    saliva_pH = 6
    st.liquid[mouth, C.H_ION] = st.liquid[mouth, C.H_ION] + mucin_hydrogen(saliva_pH, mucus_added)
    st.secretion.total_saliva = st.secretion.total_saliva + st.secretion.saliva_addition


def pancreatobiliary_secretion() -> None:
    """Pancreatic and biliary secretion, gall-bladder filling/emptying and bile-salt circulation."""
    liver_secretion_rate = 0.3
    max_gall_bladder_volume = 50
    max_emptying_rate = 0.3
    bile_content = 0.2 * 0.007
    phospholipids_content = 0.007
    target_blood_bile_salt_pool = 0.008

    bile.blood_salt_pool = bile.blood_salt_pool + 0.2 * (
        target_blood_bile_salt_pool - bile.blood_salt_pool
    ) * st.clock.dt
    circulated_bile_salt_rate = 0.2 * bile.blood_salt_pool
    bile.blood_salt_pool = bile.blood_salt_pool + (
        st.total_absorption_rate[C.CONJUGATED_BILE_SALT] - circulated_bile_salt_rate
    ) * st.clock.dt
    water_rate = (1 - bile_content - phospholipids_content) * liver_secretion_rate

    bile_salt_rate = bile_content * water_rate + circulated_bile_salt_rate
    bile_phospholipids_rate = phospholipids_content * water_rate
    filling_rate_fraction = (
        max_gall_bladder_volume - bile.gall_bladder_volume
    ) / max_gall_bladder_volume
    gb = bile.gall_bladder
    gb[C.WATER] = (
        gb[C.WATER]
        + water_rate * filling_rate_fraction * st.clock.dt
        - bile.gall_bladder_volume * 0.01 * st.clock.dt
    )
    gb[C.CONJUGATED_BILE_SALT] = (
        gb[C.CONJUGATED_BILE_SALT] + bile_salt_rate * filling_rate_fraction * st.clock.dt
    )
    gb[C.BILE_PHOSPHOLIPIDS] = (
        gb[C.BILE_PHOSPHOLIPIDS]
        + bile_phospholipids_rate * filling_rate_fraction * st.clock.dt
    )
    bile.gall_bladder_volume = (
        gb[C.WATER] + gb[C.CONJUGATED_BILE_SALT] + gb[C.BILE_PHOSPHOLIPIDS]
    )
    bile.conjugated_salt_concentration_in_gall_bladder = (
        gb[C.CONJUGATED_BILE_SALT] / bile.gall_bladder_volume
    )
    if bile.gall_bladder_volume > max_gall_bladder_volume:
        bile.gall_bladder_emptying_rate = (
            0.1 * (bile.gall_bladder_volume - max_gall_bladder_volume) / st.clock.dt
        )
    bile.gall_bladder_emptying_rate = bile.gall_bladder_emptying_rate + 0 + (
        1 * hormones.total_rate_cck_release / 5
        + 2 * math.exp(-bile.duodenal_conjugated_salt_concentration / 0.0002)
    )
    if st.mmc_status == MMCStatus.PHASE_III:
        bile.gall_bladder_emptying_rate = bile.gall_bladder_emptying_rate + 2
    bile.gall_bladder_emptying_rate = (
        max_emptying_rate
        * bile.gall_bladder_emptying_rate
        / (max_emptying_rate + bile.gall_bladder_emptying_rate)
    )

    if bile.gall_bladder_emptying_rate * st.clock.dt > bile.gall_bladder_volume:
        bile.gall_bladder_emptying_rate = 0.9 * bile.gall_bladder_volume / st.clock.dt
    aux = st.clock.dt * bile.gall_bladder_emptying_rate / bile.gall_bladder_volume
    gb[C.WATER] = gb[C.WATER] - gb[C.WATER] * aux
    gb[C.CONJUGATED_BILE_SALT] = gb[C.CONJUGATED_BILE_SALT] - gb[C.CONJUGATED_BILE_SALT] * aux
    gb[C.BILE_PHOSPHOLIPIDS] = gb[C.BILE_PHOSPHOLIPIDS] - gb[C.BILE_PHOSPHOLIPIDS] * aux

    duo2 = Segment.DUODENUM_2
    st.liquid[duo2, C.WATER] = (
        st.liquid[duo2, C.WATER]
        + gb[C.WATER] * aux
        + (1 - filling_rate_fraction) * water_rate * st.clock.dt
    )
    st.liquid[duo2, C.CONJUGATED_BILE_SALT] = (
        st.liquid[duo2, C.CONJUGATED_BILE_SALT]
        + gb[C.CONJUGATED_BILE_SALT] * aux
        + (1 - filling_rate_fraction) * bile_salt_rate * st.clock.dt
    )
    st.liquid[duo2, C.BILE_PHOSPHOLIPIDS] = (
        st.liquid[duo2, C.BILE_PHOSPHOLIPIDS]
        + gb[C.BILE_PHOSPHOLIPIDS] * aux
        + (1 - filling_rate_fraction) * bile_phospholipids_rate * st.clock.dt
    )
    bile.conjugated_salt_secretion_rate = (
        gb[C.CONJUGATED_BILE_SALT] * aux / st.clock.dt
        + (1 - filling_rate_fraction) * bile_salt_rate
    )
    # NOTE(port): the following addition appears to count bile salt already added above again; retained from Pascal.
    st.liquid[duo2, C.CONJUGATED_BILE_SALT] = (
        st.liquid[duo2, C.CONJUGATED_BILE_SALT] + bile.conjugated_salt_secretion_rate * st.clock.dt
    )
    bile.total_conjugated_salt_secretion = (
        bile.total_conjugated_salt_secretion
        + bile.conjugated_salt_secretion_rate * st.clock.dt
    )
    bile.total_secretion_rate = bile.gall_bladder_emptying_rate + (
        1 - filling_rate_fraction
    ) * liver_secretion_rate
    bile.total_secreted = (
        bile.total_secreted + bile.total_secretion_rate * st.clock.dt
    )
    volumina()
    bile.conjugated_salt_in_secreted_bile = (
        bile.conjugated_salt_secretion_rate
        / (
            bile.gall_bladder_emptying_rate
            + (1 - filling_rate_fraction) * liver_secretion_rate
        )
    )
    bile.duodenal_conjugated_salt_concentration = (
        st.liquid[duo2, C.CONJUGATED_BILE_SALT] + st.liquid[duo2, C.BILE_SALT_LCFA_COMPLEX]
    ) / st.volume[duo2]

    pancreatic_bicarbonate_secretion_rate = 1e-8 * hormones.secretin
    pancreatic_enzyme_secretion_rate = 10 * 0.5 * (
        hormones.total_rate_cck_release / (1 + hormones.total_rate_cck_release)
        + (st.volume[Segment.DUODENUM_1] + st.volume[duo2]) / (100 + st.volume[Segment.DUODENUM_1] + st.volume[duo2])
    )
    pancreatic_fluid_secretion_rate = (
        0.02
        + 0.1 * pancreatic_enzyme_secretion_rate
        + 10 * pancreatic_bicarbonate_secretion_rate
    )

    if st.mmc_status == MMCStatus.PHASE_III:
        pancreatic_enzyme_secretion_rate = pancreatic_enzyme_secretion_rate + 1

    liq = st.liquid
    liq[duo2, C.H_ION] = liq[duo2, C.H_ION] - pancreatic_bicarbonate_secretion_rate * st.clock.dt
    liq[duo2, C.PANCREATIC_LIPASE] = (
        liq[duo2, C.PANCREATIC_LIPASE] + 0.1 * pancreatic_enzyme_secretion_rate * st.clock.dt
    )
    liq[duo2, C.PANCREATIC_PHOSPHOLIPASE] = (
        liq[duo2, C.PANCREATIC_PHOSPHOLIPASE] + 0.1 * pancreatic_enzyme_secretion_rate * st.clock.dt
    )
    liq[duo2, C.PANCREATIC_AMYLASE] = (
        liq[duo2, C.PANCREATIC_AMYLASE] + 0.5 * pancreatic_enzyme_secretion_rate * st.clock.dt
    )
    liq[duo2, C.PANCREATIC_PROTEASES] = (
        liq[duo2, C.PANCREATIC_PROTEASES] + 1e-1 * pancreatic_fluid_secretion_rate * st.clock.dt
    )
    liq[duo2, C.MUCUS_PROTEIN] = (
        liq[duo2, C.MUCUS_PROTEIN] + 0.01 * pancreatic_fluid_secretion_rate * st.clock.dt
    )
    liq[duo2, C.CCK_RELEASING_PROTEIN] = (
        liq[duo2, C.CCK_RELEASING_PROTEIN] + 1e-3 * pancreatic_fluid_secretion_rate * st.clock.dt
    )
    liq[duo2, C.WATER] = liq[duo2, C.WATER] + 1 * pancreatic_fluid_secretion_rate * st.clock.dt

    st.secretion.total_pancreatic_fluid = (
        st.secretion.total_pancreatic_fluid
        + pancreatic_fluid_secretion_rate * st.clock.dt
    )


def intestinal_secretion() -> None:
    """Bicarbonate and mucus secretion of the duodenum and the other small-intestinal segments."""
    decelerator = 0.0
    target_ph = 7
    bicarbonate_secretion = 0.0  # NOTE(port): uninitialised local in Pascal; used before assignment.
    bicarbonate_secretion = (
        st.clock.dt * (decelerator * bicarbonate_secretion + 1e-5 * hormones.secretin)
        / (decelerator + 1)
    )
    duo1 = Segment.DUODENUM_1
    st.liquid[duo1, C.H_ION] = st.liquid[duo1, C.H_ION] - bicarbonate_secretion
    st.liquid[duo1, C.WATER] = (
        st.liquid[duo1, C.WATER] + 1000 * bicarbonate_secretion + 0.01 * st.clock.dt
    )
    mucus_added = 0.01 * st.clock.dt
    # NOTE(port): Pascal uses the shared global loop variable i here, which still holds 15 (ileum 6)
    # from the previous step; kept as an explicit index.
    st.liquid[Segment.ILEUM_6, C.MUCUS_PROTEIN] = st.liquid[Segment.ILEUM_6, C.MUCUS_PROTEIN] + mucus_added
    st.liquid[duo1, C.H_ION] = (
        st.liquid[duo1, C.H_ION] + mucin_hydrogen(target_ph, mucus_added)
    )

    for i in range(Segment.DUODENUM_2, Segment.ILEUM_6 + 1):
        target_ph = 8
        bicarbonate_secretion = 1e-7 * math.exp(target_ph - st.ph[i]) * st.clock.dt
        st.liquid[i, C.H_ION] = st.liquid[i, C.H_ION] - bicarbonate_secretion
        st.liquid[i, C.WATER] = (
            st.liquid[i, C.WATER] + 100 * bicarbonate_secretion + 0.01 * st.clock.dt
        )
        mucus_added = 0.1 * bicarbonate_secretion
        mucus_added = mucus_added + 0.01 * st.clock.dt
        st.liquid[i, C.MUCUS_PROTEIN] = st.liquid[i, C.MUCUS_PROTEIN] + mucus_added
        st.liquid[i, C.H_ION] = (
            st.liquid[i, C.H_ION] + mucin_hydrogen(target_ph, mucus_added)
        )


def calc_outflows() -> None:
    """Compute the outflow rates (st.outflow, st.outflow_dish) of the dishes, mouth and intestinal segments."""
    pi = 3.14159
    max_mouth_volume_intake = 10
    max_mouth_volume_swallow = 1.5

    hormones.cck_trans = 1
    st.cck_transition = hormones.cck / hormones.cck_trans
    st.cck_transition = st.cck_transition / (1 + st.cck_transition)

    for dish in range(1, 6):
        if st.volume_dish[dish] > st.eating_speed_dish[dish] * st.clock.dt:
            st.outflow_dish[dish] = st.eating_speed_dish[dish]
        else:
            st.outflow_dish[dish] = 0.9 * st.volume_dish[dish] / st.clock.dt
        if st.volume_dish[dish] < 1e-9:
            st.outflow_dish[dish] = 0
        if st.volume[Segment.MOUTH] > max_mouth_volume_intake:
            st.outflow_dish[dish] = 0

    st.solid_volume = 0
    st.large_solids = 0
    st.surface_solids = 0
    for n in range(1, diet.number_of_compounds + 1):
        if st.solid[n] > 5:
            st.large_solids = (
                st.large_solids + st.compound[2, n]
            )
        if st.solid[n] > 0:
            st.solid_volume = (
                st.solid_volume + st.compound[2, n]
            )
        st.surface_solids = (
            st.surface_solids
            + st.compound[2, n]
            * pi
            * (st.solid[n] / 20) ** 2
        )
    mouth = Segment.MOUTH
    st.outflow[mouth] = 0
    if (
        st.volume[mouth]
        > st.solid_volume
        + 0.3 * st.surface_solids
        + max_mouth_volume_swallow
    ):
        st.outflow[mouth] = st.volume[mouth] / st.clock.dt

    gastric_redistribution()

    st.base_si_transit_rate = 0.05
    st.outflow[Segment.DUODENUM_1] = 5 * st.base_si_transit_rate * st.volume[Segment.DUODENUM_1]
    st.outflow[Segment.DUODENUM_2] = 3 * st.base_si_transit_rate * st.volume[Segment.DUODENUM_2]

    st.base_si_transit_rate = 0.05
    for i in range(Segment.JEJUNUM_1, Segment.JEJUNUM_3 + 1):
        Boost = 1 * st.volume[i] / (st.max_segment_volume / 2 + st.volume[i])
        st.outflow[i] = st.base_si_transit_rate * st.volume[i] * (1 + Boost)

    for i in range(Segment.ILEUM_1, Segment.ILEUM_6 + 1):
        brake = (
            (1 - 0.5 * hormones.cck_attenuation_si * hormones.cck / (1 + hormones.cck))
            * (1 - 0.5 * hormones.gip_attenuation_si * hormones.gip / (2 + hormones.gip))
            * (1 - 0.5 * hormones.pyy_attenuation_si * hormones.pyy / (10 + hormones.pyy))
            * (1 - 0.5 * hormones.glp1_attenuation_si * hormones.glp1 / (1 + hormones.glp1))
        )
        Boost = 1 * st.volume[i] / (st.max_segment_volume / 2 + st.volume[i])
        st.outflow[i] = st.base_si_transit_rate * st.volume[i] * (1 + Boost) * brake
    if st.mmc_status == MMCStatus.PHASE_III:
        for i in range(Segment.DUODENUM_1, Segment.ILEUM_6 + 1):
            st.outflow[i] = 1 + st.outflow[i] * 10




