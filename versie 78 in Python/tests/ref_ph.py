"""Original scalar implementation of calc_all_ph (reference for the vectorised pH bisection)."""
from digsim import state as st
from digsim.calc_protein_hydrogen import calc_ph
from digsim.enums import Component as C
from digsim.enums import Segment
from digsim.state import stomach


def calc_all_ph():
    for i in range(Segment.MOUTH, Segment.ILEUM_6 + 1):
        liq = st.liquid[i]
        aux1 = (
            liq[C.CASEIN]
            + liq[C.CASEINPEP]
            + liq[C.CASEINPAN]
            + liq[C.CASEIN_PEPTIDES]
            + liq[C.SMALL_PEPTIDES]
        )
        st.ph[i] = calc_ph(
            st.volume[i],
            liq[C.H_ION],
            aux1,
            liq[C.BSA],
            liq[C.MUCUS_PROTEIN],
            liq[C.BLG],
            liq[C.PROTEIN],
            liq[C.CITRIC_ACID],
        )
    f = stomach.fundus
    aux1 = f[C.CASEIN] + f[C.CASEINPEP] + f[C.CASEINPAN] + f[C.CASEIN_PEPTIDES] + f[C.SMALL_PEPTIDES]
    stomach.ph_fundus = calc_ph(
        stomach.v_fundus, f[C.H_ION], aux1, f[C.BSA], f[C.MUCUS_PROTEIN], f[C.BLG], f[C.PROTEIN], f[C.CITRIC_ACID]
    )
