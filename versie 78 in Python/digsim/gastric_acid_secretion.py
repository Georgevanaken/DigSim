"""Gastric acid and juice secretion."""
from . import state as st
from .calc_mucin_hydrogen import mucin_hydrogen
from .enums import Component as C
from .enums import MMCStatus, Segment
from .state import stomach


def gastric_secretion() -> None:
    """Add gastric juice (water, HCl, pepsin, gastric lipase, mucin) to the corpus for one time step."""
    pH_mucus = 2
    corpus = Segment.CORPUS
    antrum = Segment.ANTRUM
    basis = 0.1  # Basis secretion rate.
    # MMC3 is the cleaning phase; extra gastric juice is secreted (amount remains to be investigated).
    if st.mmc_status == MMCStatus.PHASE_III:
        basis = 1
    target_ph = 1  # Standard acid addition is proportional to antral pH minus 1 (Konturec, 1974).
    mole_HCl = 1 * (10 / 15000) * (st.ph[antrum] - target_ph) / (5.5 - target_ph)
    mole_HCl = mole_HCl * (1 + 1 * st.concentration_lumen[antrum, C.NACL])  # Salt stimulation estimate.
    if mole_HCl < 0:
        mole_HCl = 0
    if st.ph[Segment.DUODENUM_1] < 5:
        mole_HCl = 0.01 * (10 / 15000) * (st.ph[antrum] - target_ph) / (5.5 - target_ph)

    stomach.fluid_secretion = mole_HCl * 1000 + basis
    stomach.fluid_secretion = stomach.fluid_secretion + st.volume[antrum] / (st.volume[antrum] + 5)  # Antral distension effect.
    stomach.fluid_secretion = (
        stomach.fluid_secretion
        + st.liquid[antrum, C.CASEINPEP] / (st.liquid[antrum, C.CASEINPEP] + 1e-5)
        + st.liquid[antrum, C.CASEIN_PEPTIDES] / (st.liquid[antrum, C.CASEIN_PEPTIDES] + 1e-5)
        + st.liquid[antrum, C.SMALL_PEPTIDES] / (st.liquid[antrum, C.SMALL_PEPTIDES] + 1e-5)
    )

    st.liquid[corpus, C.WATER] = st.liquid[corpus, C.WATER] + stomach.fluid_secretion * st.clock.dt
    st.liquid[corpus, C.H_ION] = st.liquid[corpus, C.H_ION] + mole_HCl * st.clock.dt
    st.liquid[corpus, C.GASTRIC_LIPASE] = st.liquid[corpus, C.GASTRIC_LIPASE] + 0.001 * stomach.fluid_secretion * st.clock.dt
    mucin_secretion = 0.01 * stomach.fluid_secretion * st.clock.dt  # Estimated mucin protein fraction in gastric secretion.
    mucin_secretion = mucin_secretion + 0.01 * st.clock.dt  # Continuous mucin secretion by the gastric wall.
    st.liquid[corpus, C.MUCUS_PROTEIN] = st.liquid[corpus, C.MUCUS_PROTEIN] + mucin_secretion
    st.liquid[corpus, C.H_ION] = st.liquid[corpus, C.H_ION] + mucin_hydrogen(pH_mucus, mucin_secretion)
    st.liquid[corpus, C.PEPSIN] = st.liquid[corpus, C.PEPSIN] + 0.005 * stomach.fluid_secretion * st.clock.dt

    stomach.total_fluid_secretion = (
        stomach.total_fluid_secretion + stomach.fluid_secretion * st.clock.dt
    )


