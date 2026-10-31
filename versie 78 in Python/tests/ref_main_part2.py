"""Main-program digestion and decomposition routines ported from Pascal."""

import numpy as np

from digsim import state as st
from digsim.main_part1 import volumina
from digsim.state import diet, stomach


def pepsin_units(pH, pepsin):
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


def lingual_lipase_units(pH, lingual_lipase):
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


def gastric_lipase_units(pH, gastric_lipase):
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


def digestion():
    change = 0.0
    linguallipolysis_rate = np.zeros(111, dtype=float)
    gastriclipolysis_rate = np.zeros(111, dtype=float)
    pancreaticlipolysis_rate = np.zeros(111, dtype=float)
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

    st.lipolysis.total_lingual_rate = 0
    st.lipolysis.total_gastric_rate = 0
    st.lipolysis.total_pancreatic_rate = 0
    st.lipolysis.percentual_loss_lingual = 0
    st.lipolysis.percentual_loss_gastric = 0
    st.lipolysis.percentual_loss_pancreatic = 0
    for st.i in range(2, 16):
        change = st.lipolysis.lingual_cst * st.liquid[st.i, 48] * st.liquid[st.i, 2] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 2]:
            change = st.liquid[st.i, 2]
        st.liquid[st.i, 7] = st.liquid[st.i, 7] + 0.666666 * change
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + 0.333333 * change
        st.liquid[st.i, 2] = st.liquid[st.i, 2] - change
        linguallipolysis_rate[st.i] = change
        if st.liquid[st.i, 2] > 1e-10:
            st.lipolysis.percentual_loss_lingual = st.lipolysis.percentual_loss_lingual + 100 * change / st.liquid[st.i, 2]
        st.lipolysis.total_lingual_rate = st.lipolysis.total_lingual_rate + linguallipolysis_rate[st.i]

        change = st.lipolysis.lingual_cst * st.liquid[st.i, 48] * st.liquid[st.i, 7] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 7]:
            change = st.liquid[st.i, 7]
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + 0.5 * change
        st.liquid[st.i, 8] = st.liquid[st.i, 8] + 0.5 * change
        st.liquid[st.i, 7] = st.liquid[st.i, 7] - change

        change = st.lipolysis.lingual_cst * st.liquid[st.i, 48] * st.liquid[st.i, 8] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 8]:
            change = st.liquid[st.i, 8]
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + change
        # NOTE(port): Preserve Pascal assignment to diglyceride component 7 during monoglyceride hydrolysis.
        st.liquid[st.i, 7] = st.liquid[st.i, 8] - change

        change = st.lipolysis.gastric_cst * st.liquid[st.i, 49] * st.liquid[st.i, 2] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 2]:
            change = st.liquid[st.i, 2]
        st.liquid[st.i, 7] = st.liquid[st.i, 7] + 0.666666 * change
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + 0.333333 * change
        st.liquid[st.i, 2] = st.liquid[st.i, 2] - change
        gastriclipolysis_rate[st.i] = change
        if st.liquid[st.i, 2] > 1e-10:
            st.lipolysis.percentual_loss_gastric = st.lipolysis.percentual_loss_gastric + 100 * change / st.liquid[st.i, 2]
        st.lipolysis.total_gastric_rate = st.lipolysis.total_gastric_rate + gastriclipolysis_rate[st.i]

        change = st.lipolysis.pancreatic_cst * st.liquid[st.i, 50] * st.liquid[st.i, 2] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 2]:
            change = st.liquid[st.i, 2]
        st.liquid[st.i, 7] = st.liquid[st.i, 7] + 0.666666 * change
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + 0.333333 * change
        st.liquid[st.i, 2] = st.liquid[st.i, 2] - change
        pancreaticlipolysis_rate[st.i] = change
        if st.liquid[st.i, 2] > 1e-10:
            st.lipolysis.percentual_loss_pancreatic = st.lipolysis.percentual_loss_pancreatic + 100 * change / st.liquid[st.i, 2]
        st.lipolysis.total_pancreatic_rate = st.lipolysis.total_pancreatic_rate + pancreaticlipolysis_rate[st.i]

        change = st.lipolysis.pancreatic_cst * st.liquid[st.i, 50] * st.liquid[st.i, 7] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 7]:
            change = st.liquid[st.i, 7]
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + 0.5 * change
        st.liquid[st.i, 8] = st.liquid[st.i, 8] + 0.5 * change
        st.liquid[st.i, 7] = st.liquid[st.i, 7] - change

        change = 100 * st.liquid[st.i, 20] * st.liquid[st.i, 51] / (0.1 + st.liquid[st.i, 51]) * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 20]:
            change = st.liquid[st.i, 20]
        if change > st.liquid[st.i, 51]:
            change = st.liquid[st.i, 51]
        st.liquid[st.i, 54] = st.liquid[st.i, 54] + change
        st.liquid[st.i, 20] = st.liquid[st.i, 20] - change
        st.liquid[st.i, 51] = st.liquid[st.i, 51] - change

        change = 100 * st.liquid[st.i, 8] * st.liquid[st.i, 51] / (0.1 + st.liquid[st.i, 51]) * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 8]:
            # NOTE(port): Pascal clamps micelle formation to fatty acid component 20, not monoglyceride 8.
            change = st.liquid[st.i, 20]
        if change > st.liquid[st.i, 51]:
            change = st.liquid[st.i, 51]
        st.liquid[st.i, 54] = st.liquid[st.i, 54] + change
        st.liquid[st.i, 8] = st.liquid[st.i, 8] - change
        st.liquid[st.i, 51] = st.liquid[st.i, 51] - change

        change = 0.1 * st.liquid[st.i, 52] * st.liquid[st.i, 56] * st.clock.dt / (st.volume[st.i] ** 2)
        if change > st.liquid[st.i, 52]:
            change = st.liquid[st.i, 52]
        st.liquid[st.i, 52] = st.liquid[st.i, 52] - change
        st.liquid[st.i, 53] = st.liquid[st.i, 53] + change
        st.liquid[st.i, 20] = st.liquid[st.i, 20] + change

        st.segment_ph = st.ph[st.i]
        change = st.liquid[st.i, 4] * pepsin_susceptibility_protein * pepsin_units(st.segment_ph, st.liquid[st.i, 40]) * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 4]:
            change = st.liquid[st.i, 4]
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + change
        st.liquid[st.i, 4] = st.liquid[st.i, 4] - change

        change = st.liquid[st.i, 66] * pepsin_susceptibility_protein * pepsin_units(st.segment_ph, st.liquid[st.i, 40]) * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 66]:
            change = st.liquid[st.i, 66]
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + change
        st.liquid[st.i, 66] = st.liquid[st.i, 66] - change

        change = st.liquid[st.i, 67] * pepsin_susceptibility_protein * pepsin_units(st.segment_ph, st.liquid[st.i, 40]) * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 67]:
            change = st.liquid[st.i, 67]
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + change
        st.liquid[st.i, 67] = st.liquid[st.i, 67] - change

        change = st.liquid[st.i, 61] * pepsin_susceptibility_casein * pepsin_units(st.segment_ph, st.liquid[st.i, 40]) * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 61]:
            change = st.liquid[st.i, 61]
        st.liquid[st.i, 62] = st.liquid[st.i, 62] + change
        st.liquid[st.i, 61] = st.liquid[st.i, 61] - change

        change = st.liquid[st.i, 4] * st.liquid[st.i, 41] * pancreatic_protease_susceptibility_protein * activity_pancreatic_protease * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 4]:
            change = st.liquid[st.i, 4]
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + change
        st.liquid[st.i, 4] = st.liquid[st.i, 4] - change

        change = st.liquid[st.i, 66] * st.liquid[st.i, 41] * pancreatic_protease_susceptibility_protein * activity_pancreatic_protease * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 66]:
            # NOTE(port): Pascal uses component 64 as the clamp for BSA/BLG digestion; preserve it.
            change = st.liquid[st.i, 64]
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + change
        st.liquid[st.i, 66] = st.liquid[st.i, 66] - change

        change = st.liquid[st.i, 67] * st.liquid[st.i, 41] * pancreatic_protease_susceptibility_protein * activity_pancreatic_protease * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 67]:
            # NOTE(port): Pascal uses component 64 as the clamp for BSA/BLG digestion; preserve it.
            change = st.liquid[st.i, 64]
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + change
        st.liquid[st.i, 67] = st.liquid[st.i, 67] - change

        change = st.liquid[st.i, 61] * st.liquid[st.i, 41] * pancreatic_protease_susceptibility_casein * activity_pancreatic_protease * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 61]:
            change = st.liquid[st.i, 61]
        st.liquid[st.i, 63] = st.liquid[st.i, 63] + change
        st.liquid[st.i, 61] = st.liquid[st.i, 61] - change

        change = st.liquid[st.i, 62] * pepsin_susceptibility_casein_fragments * pepsin_units(st.segment_ph, st.liquid[st.i, 40]) * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 62]:
            change = st.liquid[st.i, 62]
        st.liquid[st.i, 64] = st.liquid[st.i, 64] + change
        st.liquid[st.i, 62] = st.liquid[st.i, 62] - change

        change = st.liquid[st.i, 63] * st.liquid[st.i, 41] * pancreatic_protease_susceptibility_casein_fragments * activity_pancreatic_protease * st.clock.dt / st.volume[st.i]
        if change > st.liquid[st.i, 63]:
            change = st.liquid[st.i, 63]
        st.liquid[st.i, 64] = st.liquid[st.i, 64] + change
        st.liquid[st.i, 63] = st.liquid[st.i, 63] - change
        st.liquid[st.i, 65] = st.liquid[st.i, 65] + st.liquid[st.i, 64]
        st.liquid[st.i, 64] = 0

        st.liquid[st.i, 57] = st.liquid[st.i, 57] - st.clock.dt * st.liquid[st.i, 57] * np.exp(5 * np.log(5 / (5 + st.ph[st.i])))
        st.salivary_amylase_activity[st.i] = st.liquid[st.i, 57] / ((st.volume[st.i] * K_ms * (1 + st.concentration_lumen[st.i, 72] / K_cs) + st.liquid[st.i, 10] * (1 + st.concentration_lumen[st.i, 73] / K_as + st.concentration_lumen[st.i, 3] / 0.1 + st.concentration_lumen[st.i, 23] / 0.035)) * (1 + st.concentration_lumen[st.i, 74] / K_ncs))
        st.pancreatic_amylase_activity[st.i] = st.liquid[st.i, 58] / ((st.volume[st.i] * K_mp * (1 + st.concentration_lumen[st.i, 72] / K_cp) + st.liquid[st.i, 10] * (1 + st.concentration_lumen[st.i, 73] / K_ap + st.concentration_lumen[st.i, 3] / 0.1 + st.concentration_lumen[st.i, 23] / 0.035)) * (1 + st.concentration_lumen[st.i, 74] / K_ncp))
        change = st.liquid[st.i, 10] * salivary_amylase_susceptibility_HDstarch * st.salivary_amylase_activity[st.i] * st.clock.dt
        change = change + st.liquid[st.i, 10] * pancreatic_amylase_susceptibility_HDstarch * st.pancreatic_amylase_activity[st.i] * st.clock.dt
        if change > st.liquid[st.i, 10]:
            change = st.liquid[st.i, 10]
        st.liquid[st.i, 23] = st.liquid[st.i, 23] + 0.99 * change
        st.liquid[st.i, 24] = st.liquid[st.i, 24] + 0.01 * change
        st.liquid[st.i, 10] = st.liquid[st.i, 10] - change

        change = st.liquid[st.i, 11] * salivary_amylase_susceptibility_SDstarch * st.salivary_amylase_activity[st.i] * st.clock.dt
        change = change + st.liquid[st.i, 11] * pancreatic_amylase_susceptibility_SDstarch * st.pancreatic_amylase_activity[st.i] * st.clock.dt
        if change > st.liquid[st.i, 11]:
            change = st.liquid[st.i, 11]
        st.liquid[st.i, 23] = st.liquid[st.i, 23] + 0.99 * change
        st.liquid[st.i, 24] = st.liquid[st.i, 24] + 0.01 * change
        st.liquid[st.i, 11] = st.liquid[st.i, 11] - change

        change = st.liquid[st.i, 12] * salivary_amylase_susceptibility_Rstarch * st.salivary_amylase_activity[st.i] * st.clock.dt
        change = change + st.liquid[st.i, 12] * pancreatic_amylase_susceptibility_Rstarch * st.pancreatic_amylase_activity[st.i] * st.clock.dt
        if change > st.liquid[st.i, 12]:
            change = st.liquid[st.i, 12]
        st.liquid[st.i, 23] = st.liquid[st.i, 23] + 0.99 * change
        st.liquid[st.i, 24] = st.liquid[st.i, 24] + 0.01 * change
        st.liquid[st.i, 12] = st.liquid[st.i, 12] - change

    volumina()


def decomposition():
    aux = 0.0
    increase = 0.0
    degcomp = 0
    for st.compound_number in range(1, diet.number_of_compounds + 1):
        for st.i in range(3, 5):
            increase = st.compound[st.i, st.compound_number] * st.clock.dt / diet.gastric_dec_time[st.compound_number]
            aux = st.compound[st.i, st.compound_number] - increase
            degcomp = diet.gastric_dec_comp[st.compound_number, 1]
            st.compound[st.i, degcomp] = st.compound[st.i, degcomp] + increase * diet.gastric_dec_frac[st.compound_number, 1] / 100
            degcomp = diet.gastric_dec_comp[st.compound_number, 2]
            st.compound[st.i, degcomp] = st.compound[st.i, degcomp] + increase * diet.gastric_dec_frac[st.compound_number, 2] / 100
            for st.j in range(1, st.comp + 1):
                st.liquid[st.i, st.j] = st.liquid[st.i, st.j] + 0.01 * (100 - diet.gastric_dec_frac[st.compound_number, 1] - diet.gastric_dec_frac[st.compound_number, 2]) * diet.compound_composition[st.compound_number, st.j] * increase
            st.compound[st.i, st.compound_number] = aux

        st.i = 2
        increase = st.compound[st.i, st.compound_number] * st.clock.dt / diet.oral_dec_time[st.compound_number]
        aux = st.compound[st.i, st.compound_number] - increase
        degcomp = diet.oral_dec_comp[st.compound_number, 1]
        st.compound[st.i, degcomp] = st.compound[st.i, degcomp] + increase * diet.oral_dec_frac[st.compound_number, 1] / 100
        degcomp = diet.oral_dec_comp[st.compound_number, 2]
        st.compound[st.i, degcomp] = st.compound[st.i, degcomp] + increase * diet.oral_dec_frac[st.compound_number, 2] / 100
        for st.j in range(1, st.comp + 1):
            st.liquid[st.i, st.j] = st.liquid[st.i, st.j] + 0.01 * (100 - diet.oral_dec_frac[st.compound_number, 1] - diet.oral_dec_frac[st.compound_number, 2]) * diet.compound_composition[st.compound_number, st.j] * increase
        st.compound[st.i, st.compound_number] = aux

        increase = st.compound_fundus[st.compound_number] * st.clock.dt / diet.gastric_dec_time[st.compound_number]
        aux = st.compound_fundus[st.compound_number] - increase
        degcomp = diet.gastric_dec_comp[st.compound_number, 1]
        st.compound_fundus[degcomp] = st.compound_fundus[degcomp] + increase * diet.gastric_dec_frac[st.compound_number, 1] / 100
        degcomp = diet.gastric_dec_comp[st.compound_number, 2]
        st.compound_fundus[degcomp] = st.compound_fundus[degcomp] + increase * diet.gastric_dec_frac[st.compound_number, 2] / 100
        for st.j in range(1, st.comp + 1):
            stomach.fundus[st.j] = stomach.fundus[st.j] + 0.01 * (100 - diet.gastric_dec_frac[st.compound_number, 1] - diet.gastric_dec_frac[st.compound_number, 2]) * diet.compound_composition[st.compound_number, st.j] * increase
        st.compound_fundus[st.compound_number] = aux

        for st.i in range(5, 16):
            increase = st.compound[st.i, st.compound_number] * st.clock.dt / diet.intestinal_dec_time[st.compound_number]
            aux = st.compound[st.i, st.compound_number] - increase
            degcomp = diet.intestinal_dec_comp[st.compound_number, 1]
            st.compound[st.i, degcomp] = st.compound[st.i, degcomp] + increase * diet.intestinal_dec_frac[st.compound_number, 1] / 100
            degcomp = diet.intestinal_dec_comp[st.compound_number, 2]
            st.compound[st.i, degcomp] = st.compound[st.i, degcomp] + increase * diet.intestinal_dec_frac[st.compound_number, 2] / 100
            for st.j in range(1, st.comp + 1):
                st.liquid[st.i, st.j] = st.liquid[st.i, st.j] + 0.01 * (100 - diet.intestinal_dec_frac[st.compound_number, 1] - diet.intestinal_dec_frac[st.compound_number, 2]) * diet.compound_composition[st.compound_number, st.j] * increase
            st.compound[st.i, st.compound_number] = aux



