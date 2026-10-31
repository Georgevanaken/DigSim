"""Gastric motility: stomach tone, solids, viscosities and emptying."""
import math

import numpy as np

from . import calc_protein_hydrogen as cph
from . import state as st
from .compartment_properties import volumina
from .enums import Component as C
from .enums import MMCStatus, Segment
from .state import diet, hormones, stomach

# Implementation-level Pascal variable.
target_corpus_comp = np.zeros((111,), dtype=float)


def pressures() -> None:
    """Update the stomach volume and the gastric tone (raised by CCK-driven relaxation)."""
    stomach.volume = stomach.v_fundus_corpus + st.volume[Segment.ANTRUM]
    stomach.tone = 3 + 0.0177 * (stomach.volume + stomach.total_solids_corpus) * (1 - 0.5 * st.cck_transition)
    if stomach.tone < 3:
        stomach.tone = 3


def solids_in_stomach() -> None:
    """Sum the solid mass in the antrum, corpus and fundus."""
    stomach.solids_in_antrum = 0
    stomach.solids_in_corpus = 0
    stomach.solids_in_fundus = 0
    for n in range(1, diet.number_of_compounds + 1):
        if diet.compound_viscosity[n] > 100 * stomach.viscosity_liquids_in_corpus and st.solid[n] == 0:
            st.solid[n] = 10
        if st.solid[n] != 0:
            stomach.solids_in_antrum += st.compound[Segment.ANTRUM, n]
            stomach.solids_in_corpus += st.compound[Segment.CORPUS, n]
            stomach.solids_in_fundus += st.compound_fundus[n]
    stomach.total_solids_in_stomach = stomach.solids_in_antrum + stomach.solids_in_corpus + stomach.solids_in_fundus


def gastric_compartment_viscosities() -> None:
    """Estimate the viscosity of the liquid in the fundus, corpus and antrum."""
    # NOTE(port): solids_in_corpus is an uninitialised local in Pascal and shadows the state variable.
    phi = 0.0  # NOTE(port): uninitialised in Pascal
    solids_in_corpus = 0.0  # NOTE(port): uninitialised in Pascal
    phimax = 0.55
    eta = 2.5

    corpus = Segment.CORPUS
    phi = solids_in_corpus / st.volume[corpus]
    stomach.viscosity_liquids_in_corpus = 5
    for n in range(1, diet.number_of_compounds + 1):
        if st.compound[corpus, n] > 1e-40:
            stomach.viscosity_liquids_in_corpus += diet.compound_viscosity[n] * (
                1 / (1 + st.volume[corpus] * diet.compound_density[n] / st.compound[corpus, n])
            )
    if phi > 0.5:
        stomach.corpus_viscosity = stomach.viscosity_liquids_in_corpus
        for n in range(1, diet.number_of_compounds + 1):
            stomach.corpus_viscosity += st.solid[n] * st.compound[corpus, n] * diet.compound_viscosity[n]
        stomach.corpus_viscosity /= solids_in_corpus
    else:
        stomach.corpus_viscosity = stomach.viscosity_liquids_in_corpus * math.exp(-eta * phimax * math.log(1 - phi / phimax))


def _sequential_sum(values: np.ndarray) -> float:
    """Left-to-right sum (np.sum would use pairwise summation and round differently)."""
    return float(np.cumsum(values)[-1]) + 0.0


def gastric_redistribution() -> None:
    """Redistribute liquid and solids between fundus, corpus and antrum and empty the antrum."""
    liquids_corpus_fundus = np.zeros((111,), dtype=float)
    target_corpus = np.zeros((111,), dtype=float)
    x = 0.0
    retropulsionvolume = 0.0
    flocculation_rate_corpus = 0.0
    flocculation_rate_antrum = 0.0
    floc_dissolution_rate_corpus = 0.0
    floc_ph = 0.0
    aux = 0.0
    bolus_surface = 0.0
    kineticfactor = 0.0
    change = 0.0
    discost = 0.0
    Corpus = 0.0
    up_corpus_to_fundus = 0.0
    down_fundus_to_corpus = 0.0
    target_antralvolume = 0.0
    transfer_efficiency_fundus_corpus = 0.0

    corpus = Segment.CORPUS
    antrum = Segment.ANTRUM
    duo1 = Segment.DUODENUM_1
    c = slice(1, st.comp + 1)
    cn = slice(1, diet.number_of_compounds + 1)

    volumina()
    st.cycle += st.clock.dt
    if st.cycle > 0.3333:
        st.cycle = 0
        st.count = 0
        st.outflow[antrum] = 0
    if st.cycle < 0.1111:
        # 1. Peristalsis from fundus to corpus.
        stomach.v_fundus_corpus = stomach.v_fundus + st.volume[corpus]
        transfer_efficiency_fundus_corpus = 1
        down_fundus_to_corpus = transfer_efficiency_fundus_corpus * stomach.v_fundus / 0.1111
        moved = down_fundus_to_corpus * st.clock.dt * stomach.fundus[c] / stomach.v_fundus
        st.liquid[corpus, c] += moved
        stomach.fundus[c] -= moved
        moved = down_fundus_to_corpus * st.clock.dt * st.compound_fundus[cn] / stomach.v_fundus
        st.compound[corpus, cn] += moved
        st.compound_fundus[cn] -= moved
        volumina()
    elif st.cycle < 0.2222:
        # Retropulsion from corpus to fundus and flow from corpus to antrum.
        Corpus = st.volume[corpus]
        up_corpus_to_fundus = 1 * st.volume[corpus] / 0.1111
        change = st.liquid[corpus, c] * st.clock.dt * up_corpus_to_fundus / st.volume[corpus]
        st.liquid[corpus, c] -= change
        stomach.fundus[c] += change
        change = st.compound[corpus, cn] * st.clock.dt * up_corpus_to_fundus / st.volume[corpus]
        st.compound[corpus, cn] -= change
        st.compound_fundus[cn] += change
        volumina()

        pressures()
        gastric_compartment_viscosities()
        st.max_antral_volume = 60
        target_antralvolume = 5 + st.max_antral_volume * stomach.tone / (1 + stomach.tone)
        st.outflow[corpus] = Corpus * ((1 - st.volume[antrum] / target_antralvolume) / 0.1111) * 0.1
        if st.mmc_status == MMCStatus.PHASE_III:
            st.outflow[corpus] = 2 * st.outflow[corpus]
        change = st.outflow[corpus] * st.liquid[corpus, c] * st.clock.dt / st.volume[corpus]
        st.liquid[corpus, c] -= change
        st.liquid[antrum, c] += change
        change = st.outflow[corpus] * st.clock.dt * st.compound[corpus, cn] / st.volume[corpus]
        st.compound[corpus, cn] -= change
        st.compound[antrum, cn] += change
        volumina()
    else:
        st.isotonic = 0.05
        st.osmotic.concentration = 2
        hormones.cck_gastric_emptying_stress = 2 * hormones.total_rate_cck_release + 1 * hormones.cck
        if st.mmc_status == MMCStatus.PHASE_II:
            hormones.cck_gastric_emptying_stress = 0.5 * hormones.cck_gastric_emptying_stress
        if st.mmc_status == MMCStatus.PHASE_III:
            hormones.cck_gastric_emptying_stress = 0 * hormones.cck_gastric_emptying_stress
        stomach.osmotic_emptying_stress = st.osmotic.concentration * st.osmotic.brush_duodenum / st.isotonic
        if st.mmc_status == MMCStatus.PHASE_II:
            stomach.osmotic_emptying_stress = 0.5 * stomach.osmotic_emptying_stress
        if st.mmc_status == MMCStatus.PHASE_III:
            stomach.osmotic_emptying_stress = 0 * stomach.osmotic_emptying_stress
        stomach.ph_emptying_stress = 5 * (abs(6 - st.ph[duo1]) / 6)
        if st.mmc_status == MMCStatus.PHASE_II:
            stomach.ph_emptying_stress = 0.5 * stomach.ph_emptying_stress
        if st.mmc_status == MMCStatus.PHASE_III:
            stomach.ph_emptying_stress = 0 * stomach.ph_emptying_stress
        hormones.glp1_gastric_emptying_stress = 10 * hormones.glp1
        if st.mmc_status == MMCStatus.PHASE_II:
            hormones.glp1_gastric_emptying_stress = 0.5 * hormones.glp1_gastric_emptying_stress
        if st.mmc_status == MMCStatus.PHASE_III:
            hormones.glp1_gastric_emptying_stress = 0 * hormones.glp1_gastric_emptying_stress
        stomach.viscous_emptying_stress = 0.21 * stomach.corpus_viscosity / stomach.tone
        if st.mmc_status == MMCStatus.PHASE_II:
            stomach.viscous_emptying_stress = 0.5 * stomach.viscous_emptying_stress
        if st.mmc_status == MMCStatus.PHASE_III:
            stomach.viscous_emptying_stress = 0 * stomach.viscous_emptying_stress
        if stomach.ph_emptying_stress < 0:
            stomach.ph_emptying_stress = 0
        stomach.emptying_stress = (
            hormones.cck_gastric_emptying_stress
            + stomach.osmotic_emptying_stress
            + stomach.ph_emptying_stress
            + hormones.glp1_gastric_emptying_stress
            + stomach.viscous_emptying_stress
        )
        stomach.emptying_stress = 0.2 * stomach.emptying_stress
        stomach.fraction_emptying = 1 / (1 + stomach.emptying_stress)
        if st.volume[antrum] > 0.1:
            st.outflow[antrum] = 1 * st.volume[antrum] * stomach.fraction_emptying / 0.1111
            retropulsionvolume = 1 * st.volume[antrum] * (1 - stomach.fraction_emptying) / 0.1111
        else:
            st.outflow[antrum] = 0
            retropulsionvolume = 0

        change = st.outflow[antrum] * st.liquid[antrum, c] * st.clock.dt / st.volume[antrum]
        retropulsed = retropulsionvolume * st.clock.dt * st.liquid[antrum, c] / st.volume[antrum]
        st.liquid[corpus, c] += retropulsed
        st.liquid[antrum, c] = st.liquid[antrum, c] - change - retropulsed
        st.liquid[duo1, c] += change
        stomach.gastric_emptying_water = st.outflow[antrum] * st.liquid[antrum, C.WATER] * st.clock.dt / st.volume[antrum]
        stomach.gastric_emptying_pharmaceutical = st.outflow[antrum] * st.liquid[antrum, C.PHARMACEUTICAL] * st.clock.dt / st.volume[antrum]
        for n in range(1, diet.number_of_compounds + 1):
            stomach.gastric_sieve[n] = math.exp(-st.solid[n])
        if st.mmc_status == MMCStatus.PHASE_III:
            stomach.gastric_sieve[cn] = 1
        change = st.outflow[antrum] * stomach.gastric_sieve[cn] * st.compound[antrum, cn] * st.clock.dt / st.volume[antrum]
        retropulsed = retropulsionvolume * st.clock.dt * st.compound[antrum, cn] / st.volume[antrum]
        st.compound[corpus, cn] += retropulsed
        st.compound[antrum, cn] = st.compound[antrum, cn] - change - retropulsed
        st.compound[duo1, cn] += change
        volumina()

    # 2a. Redistribution between fundus and corpus through slow sedimentation or creaming in the liquid medium.
    liquids_corpus_fundus[c] = st.liquid[corpus, c] + stomach.fundus[c]
    target_corpus[c] = st.liquid[corpus, c]
    target_corpus[C.FAT] = 0 * liquids_corpus_fundus[C.FAT]
    x = (target_corpus[c] - st.liquid[corpus, c]) * 10 * st.clock.dt / stomach.corpus_viscosity
    stomach.fundus[c] -= x
    st.liquid[corpus, c] += x
    # Sequential subtraction (cumsum) keeps the original left-to-right rounding.
    stomach.v_fundus = float(np.cumsum(np.concatenate(([stomach.v_fundus], -(x / st.density[c]))))[-1])
    st.volume[corpus] = float(np.cumsum(np.concatenate(([st.volume[corpus]], x / st.density[c])))[-1])

    # NOTE(port): Pascal evaluates the settling speed with the stale loop variable, which still holds
    # the last compound number (25); kept as an explicit index. The same applies to step 2b below.
    last_compound = diet.number_of_compounds
    x = diet.compound_density[last_compound] - 1
    x = x * 0.6 * 10 * st.solid[last_compound] * 9.81 * (st.solid[last_compound] * st.solid[last_compound]) * st.clock.dt / stomach.corpus_viscosity
    if x > 0:
        moved = x * st.compound_fundus[cn] / stomach.v_fundus
        st.compound[corpus, cn] += moved
        st.compound_fundus[cn] -= moved
    if x < 0:
        moved = x * st.compound[corpus, cn] / st.volume[corpus]
        st.compound_fundus[cn] -= moved
        st.compound[corpus, cn] += moved
    volumina()

    # 2b. Redistribution of compounds between corpus and antrum by slow sedimentation or creaming.
    x = diet.compound_density[last_compound] - 1
    x = x * 0.6 * 10 * st.solid[last_compound] * 9.81 * (st.solid[last_compound] * st.solid[last_compound]) * st.clock.dt / stomach.corpus_viscosity
    if x > 0:
        moved = x * st.compound[corpus, cn] / st.volume[corpus]
        st.compound[antrum, cn] += moved
        st.compound[corpus, cn] -= moved
    if x < 0:
        moved = x * st.compound[antrum, cn] / st.volume[antrum]
        st.compound[corpus, cn] -= moved
        st.compound[antrum, cn] += moved
    volumina()

    stomach.fundus_viscosity = 2
    stomach.corpus_viscosity = 2
    stomach.antral_viscosity = 2

    cph.calc_all_ph()
    st.flocculation_composition[c] = 0
    floc_ph = 4.5
    st.flocculation_composition[C.CASEIN] = 1
    st.flocculation_composition[C.CASEINPEP] = 1
    aux = st.flocculation_composition[C.CASEIN] + st.flocculation_composition[C.CASEINPEP]
    st.flocculation_composition[C.H_ION] = 1 * cph.calc_hydrogen(
        floc_ph,
        aux,
        st.flocculation_composition[C.MUCUS_PROTEIN],
        st.flocculation_composition[C.BSA],
        st.flocculation_composition[C.BLG],
        st.flocculation_composition[C.PROTEIN],
        st.flocculation_composition[C.CITRIC_ACID],
    )

    flocculation_rate_corpus = 0
    kineticfactor = 3000
    if st.volume[corpus] != 0:
        flocculation_rate_corpus = kineticfactor * math.exp(-8 * abs(st.ph[corpus] - floc_ph))
    change = flocculation_rate_corpus * st.flocculation_composition[c] * st.liquid[corpus, c] * st.clock.dt
    # Where the change would exceed the liquid amount nothing is transferred (as in the original).
    moved = np.where(change > st.liquid[corpus, c], 0.0, change)
    stomach.solids_corpus[c] += moved
    st.liquid[corpus, c] -= moved

    discost = 1
    aux = 10 * abs(st.ph[corpus] - floc_ph)
    floc_dissolution_rate_corpus = discost * aux / (1 + aux)
    stomach.total_solids_corpus = _sequential_sum(stomach.solids_corpus[c])
    if stomach.total_solids_corpus > 0.00001:
        bolus_surface = math.exp(0.666666666666 * math.log(stomach.total_solids_corpus))
    else:
        bolus_surface = 0
    if stomach.total_solids_corpus > 0:
        dissolved = floc_dissolution_rate_corpus * bolus_surface * st.flocculation_composition[c] * st.clock.dt / stomach.total_solids_corpus
        dissolved = np.where(dissolved > stomach.solids_corpus[c], stomach.solids_corpus[c], dissolved)
        stomach.solids_corpus[c] -= dissolved
        st.liquid[corpus, c] += dissolved

    stomach.solids_in_corpus = _sequential_sum(stomach.solids_corpus[c])

    flocculation_rate_antrum = 0
    kineticfactor = 1e3
    if st.ph[antrum] < 5.0 and st.ph[antrum] > 4 and st.volume[antrum] != 0:
        flocculation_rate_antrum = kineticfactor * math.exp(-abs(st.ph[antrum] - floc_ph))
    change = flocculation_rate_antrum * st.flocculation_composition[c] * (st.liquid[antrum, c] / st.volume[antrum]) * st.clock.dt
    stomach.solids_antrum[c] += change
    st.liquid[antrum, c] -= change

    stomach.total_solids_antrum = _sequential_sum(stomach.solids_antrum[c])
    aux = stomach.solids_antrum[C.H_ION]
    stomach.solids_antrum[C.H_ION] = st.flocculation_composition[C.H_ION] * stomach.total_solids_antrum
    st.liquid[antrum, C.H_ION] = st.liquid[antrum, C.H_ION] - stomach.solids_antrum[C.H_ION] + aux

    # Syneresis is switched off in the original (rate 0), so the antrum solids are unchanged.
    stomach.solids_corpus[c] += stomach.solids_antrum[c]
    stomach.solids_antrum[c] = 0
    stomach.total_solids_corpus = _sequential_sum(stomach.solids_corpus[c])
    stomach.total_solids_antrum = 0

    stomach.total_solids_corpus = _sequential_sum(stomach.solids_corpus[c] + stomach.solids_fundus[c])
    aux = stomach.solids_corpus[C.H_ION] + stomach.solids_fundus[C.H_ION]
    stomach.solids_corpus[C.H_ION] = st.flocculation_composition[C.H_ION] * stomach.total_solids_corpus
    st.liquid[corpus, C.H_ION] = st.liquid[corpus, C.H_ION] - stomach.solids_corpus[C.H_ION] - stomach.solids_fundus[C.H_ION] + aux
