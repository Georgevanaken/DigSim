"""Original scalar implementation of main_part3, kept as a reference for equivalence tests."""
import math

import numpy as np
from ref_brush_border_enzymes import drag, maltase, sucrase

from digsim import state as st
from digsim.compartment_properties import volumina
from digsim.state import bile, glucose, hormones, stomach


def absorption_in_small_intestine():
    totaltargetabsorption_rate = 0.0
    delay = 0.0
    change = 0.0
    Mestentric_transfomation_of_fructose_to_glucose = 0.0
    R = 0.8  # Recovery of glucose after absorption by enterocytes; some fuels enterocytes.

    st.calories.total_absorption_rate = 0
    glucose.g0_in = 0
    glucose.g1_in = 0
    glucose.g2_in = 0
    glucose.fructose_decay_rate = 0
    for st.i in range(5, 16):
        totaltargetabsorption_rate = 0
        for st.j in range(2, st.comp + 1):
            if (st.j != 3) and (st.j != 34) and (st.j != 35):
                if (st.concentration_brush[st.i, st.j] > 1e-40) and (st.j != 80):
                    st.absorption_rate[st.i, st.j] = (
                        st.vmax_absorption[st.i, st.j] * st.concentration_brush[st.i, st.j]
                        / (st.concentration_brush[st.i, st.j] + st.km_absorption[st.i, st.j])
                    )
                else:
                    st.absorption_rate[st.i, st.j] = 0
            else:
                st.scratch = (
                    st.concentration_brush[st.i, 3]
                    + st.concentration_brush[st.i, 34]
                    + st.concentration_brush[st.i, 35]
                )
                st.scratch = st.scratch * (1 + 100 * st.concentration_brush[st.i, 78])
                if st.scratch > 0:
                    st.absorption_rate[st.i, 3] = (
                        st.vmax_absorption[st.i, 3] * st.scratch / (st.scratch + st.km_absorption[st.i, 3])
                        + 0.5 * st.scratch / (st.scratch + 0.1 * 0.0225560)
                    )
                    st.absorption_rate[st.i, 34] = 0
                    st.absorption_rate[st.i, 35] = 0
                else:
                    st.absorption_rate[st.i, 3] = 0
                    st.absorption_rate[st.i, 34] = 0
                    st.absorption_rate[st.i, 35] = 0
            # Weber competition among macronutrients.
            totaltargetabsorption_rate = totaltargetabsorption_rate + st.absorption_rate[st.i, st.j]

        if st.max_total_absorption_rate[st.i] == 0:
            for st.j in range(2, st.comp + 1):
                st.absorption_rate[st.i, st.j] = 0
        else:
            for st.j in range(2, st.comp + 1):
                if st.calories.value[st.j] != 0:
                    st.absorption_rate[st.i, st.j] = (
                        st.absorption_rate[st.i, st.j] * st.max_total_absorption_rate[st.i]
                        / (totaltargetabsorption_rate + st.max_total_absorption_rate[st.i])
                    )

        if st.scratch != 0:
            glucose.absorption_rate[st.i] = st.absorption_rate[st.i, 3]
            st.absorption_rate[st.i, 3] = (
                st.concentration_brush[st.i, 3] / (st.scratch + 1e-3) * glucose.absorption_rate[st.i]
            )
            st.absorption_rate[st.i, 34] = (
                st.concentration_brush[st.i, 34] / (st.scratch + 1e-3) * glucose.absorption_rate[st.i]
            )
            st.absorption_rate[st.i, 35] = (
                st.concentration_brush[st.i, 35] / (st.scratch + 1e-3) * glucose.absorption_rate[st.i]
            )
        st.absorption_rate[st.i, 70] = 0.15 * st.absorption_rate[st.i, 64]
        st.absorption_rate[st.i, 71] = st.absorption_rate[st.i, 65]

        for st.j in range(1, st.comp + 1):
            if st.calories.value[st.j] != 0:
                st.absorption_rate[st.i, st.j] = st.absorption_rate[st.i, st.j] / st.calories.value[st.j]

        st.absorption_rate[st.i, 1] = (
            2 * (st.volume[st.i] / (st.volume[st.i] + 100)) * (0.2 - st.osmotic.brush[st.i])
            + (30 * 18 / 180.15588)
            * (st.absorption_rate[st.i, 3] + st.absorption_rate[st.i, 34] + st.absorption_rate[st.i, 35])
        )
        st.absorption_rate[st.i, 78] = (1 * 58.44 / 180.15588) * (
            st.absorption_rate[st.i, 3] + st.absorption_rate[st.i, 34] + st.absorption_rate[st.i, 35]
        )

    glucose.total_absorption_rate = 0
    bile.conjugated_salt_absorption_rate = 0
    for st.i in range(5, 16):
        glucose.total_absorption_rate = (
            glucose.total_absorption_rate + st.absorption_rate[st.i, 3]
            + st.absorption_rate[st.i, 34] + st.absorption_rate[st.i, 35]
        )
    glucose.total_absorption_rate = glucose.total_absorption_rate / 180.15588

    for st.j in range(1, st.comp + 1):
        st.total_absorption_rate[st.j] = 0
        st.absorbed[st.j] = 0
    for st.j in range(1, st.comp + 1):
        for st.i in range(5, 16):
            st.total_absorption_rate[st.j] = st.total_absorption_rate[st.j] + st.absorption_rate[st.i, st.j]
            st.calories.total_absorption_rate = (
                st.calories.total_absorption_rate + st.absorption_rate[st.i, st.j] * st.calories.value[st.j]
            )
            st.absorbed[st.j] = st.absorbed[st.j] + st.absorption_rate[st.i, st.j] * st.clock.dt
        st.calories.absorbed = st.calories.absorbed + st.calories.value[st.j] * 1 * st.absorbed[st.j]

    delay = 2
    delay = 1 / delay
    for st.i in range(5, 16):
        for st.j in range(1, st.comp + 1):
            st.mesenteric[st.i, st.j] = st.mesenteric[st.i, st.j] + st.absorption_rate[st.i, st.j] * st.clock.dt
            change = delay * st.mesenteric[st.i, st.j] * st.clock.dt
            if change > st.mesenteric[st.i, st.j]:
                change = st.mesenteric[st.i, st.j]
            st.blood[st.j] = st.blood[st.j] + change
            st.mesenteric[st.i, st.j] = st.mesenteric[st.i, st.j] - change
            Mestentric_transfomation_of_fructose_to_glucose = 0.5 * st.mesenteric[st.i, 22]
            if st.j == 3:
                glucose.g0_in = (
                    glucose.g0_in + delay * st.mesenteric[st.i, 3] * 1000
                    + Mestentric_transfomation_of_fructose_to_glucose * 1000
                )
            if st.j == 34:
                glucose.g1_in = glucose.g1_in + delay * st.mesenteric[st.i, 34] * 1000
            if st.j == 35:
                glucose.g2_in = glucose.g2_in + delay * st.mesenteric[st.i, 35] * 1000
            if st.j == 22:
                glucose.fructose_decay_rate = (
                    glucose.fructose_decay_rate + delay * st.mesenteric[st.i, 22] * 1000
                    - Mestentric_transfomation_of_fructose_to_glucose * 1000
                )

    glucose.g0_in = R * glucose.g0_in
    glucose.g1_in = R * glucose.g1_in
    glucose.g2_in = R * glucose.g2_in
    st.total_water_absorption_rate = st.total_absorption_rate[1]


def mucus_layer():
    divisor = 0.0
    aux = 0.0
    decelerator = 5
    divisor = 1 / (decelerator + 1)
    sucrase()
    maltase()
    volumina()
    for st.i in range(5, 16):
        aux = st.absorption_rate[st.i, 1] * st.clock.dt
        st.liquid_next[st.i, 1] = st.liquid[st.i, 1] - st.absorption_rate[st.i, 1] * st.clock.dt
        for st.j in range(2, st.comp + 1):
            if st.j != 80:
                if (st.j != 3) and (st.j != 21) and (st.j != 22) and (st.j != 23) and (st.j != 34) and (st.j != 35) and (st.j != 36):
                    aux = (
                        st.concentration_lumen[st.i, st.j]
                        * (1 + st.surface_factor[st.i] * drag * st.diffusion_time[st.i, st.j] * st.absorption_rate[st.i, 1])
                        - st.diffusion_time[st.i, st.j] * st.surface_factor[st.i] * st.absorption_rate[st.i, st.j]
                    )
                    st.concentration_brush[st.i, st.j] = divisor * (
                        decelerator * st.concentration_brush[st.i, st.j] + aux
                    )
                st.liquid_next[st.i, st.j] = st.liquid[st.i, st.j] - st.absorption_rate[st.i, st.j] * st.clock.dt
                if st.liquid_next[st.i, st.j] < 0:
                    st.liquid_next[st.i, st.j] = 0
        st.liquid_next[st.i, 80] = st.liquid[st.i, 80] - st.absorption_rate[st.i, 80] * st.clock.dt

        aux = (
            st.concentration_lumen[st.i, 3]
            * (1 + st.surface_factor[st.i] * st.diffusion_time[st.i, 3] * drag * st.absorption_rate[st.i, 1])
            + st.diffusion_time[st.i, 3] * st.surface_factor[st.i]
            * (-st.absorption_rate[st.i, 3] + st.brush_enzymes.sucrase_transformation[st.i, 21] * (180.15588 / 342.2965)
               + st.brush_enzymes.maltase_transformation[st.i, 23] * 2 * (180.15588 / 342.2965))
        )
        st.concentration_brush[st.i, 3] = divisor * (decelerator * st.concentration_brush[st.i, 3] + aux)
        st.liquid_next[st.i, 3] = st.liquid[st.i, 3] + st.clock.dt * (
            -st.absorption_rate[st.i, 3] + st.brush_enzymes.sucrase_transformation[st.i, 21] * (180.15588 / 342.2965)
            + st.brush_enzymes.maltase_transformation[st.i, 23] * 2 * (180.15588 / 342.2965)
        )

        aux = (
            st.concentration_lumen[st.i, 34]
            * (1 + st.surface_factor[st.i] * st.diffusion_time[st.i, 34] * drag * st.absorption_rate[st.i, 1])
            + st.diffusion_time[st.i, 34] * st.surface_factor[st.i]
            * (-st.absorption_rate[st.i, 34] + st.brush_enzymes.sucrase_transformation[st.i, 36] * (180.15588 / 342.2965)
               + st.brush_enzymes.maltase_transformation[st.i, 24] * 2 * (180.15588 / 342.2965))
        )
        st.concentration_brush[st.i, 34] = divisor * (decelerator * st.concentration_brush[st.i, 34] + aux)
        st.liquid_next[st.i, 34] = st.liquid[st.i, 34] + st.clock.dt * (
            -st.absorption_rate[st.i, 34] + st.brush_enzymes.sucrase_transformation[st.i, 36] * (180.15588 / 342.2965)
            + st.brush_enzymes.maltase_transformation[st.i, 24] * 2 * (180.15588 / 342.2965)
        )
        if st.liquid_next[st.i, 34] < 1e-10:
            st.liquid_next[st.i, 34] = 1e-10

        aux = (
            st.concentration_lumen[st.i, 35]
            * (1 + st.surface_factor[st.i] * st.diffusion_time[st.i, 35] * drag * st.absorption_rate[st.i, 1])
            + st.diffusion_time[st.i, 35] * st.surface_factor[st.i] * -st.absorption_rate[st.i, 35]
        )
        st.concentration_brush[st.i, 35] = divisor * (decelerator * st.concentration_brush[st.i, 35] + aux)

        st.liquid_next[st.i, 21] = st.liquid[st.i, 21] + st.clock.dt * (
            -st.absorption_rate[st.i, 21] - st.brush_enzymes.sucrase_transformation[st.i, 21]
        )
        st.liquid_next[st.i, 36] = st.liquid[st.i, 36] + st.clock.dt * (
            -st.absorption_rate[st.i, 36] - st.brush_enzymes.sucrase_transformation[st.i, 36]
        )
        aux = st.liquid[st.i, 22] + st.clock.dt * (
            -st.absorption_rate[st.i, 22]
            + st.brush_enzymes.sucrase_transformation[st.i, 21] * (180.15588 / 342.2965)
            + st.brush_enzymes.sucrase_transformation[st.i, 36] * (180.15588 / 342.2965)
        )
        if aux < 0:
            st.absorption_rate[st.i, 22] = (
                0.95 * st.liquid[st.i, 22] / st.clock.dt
                + st.brush_enzymes.sucrase_transformation[st.i, 21] * (180.15588 / 342.2965)
                + st.brush_enzymes.sucrase_transformation[st.i, 36] * (180.15588 / 342.2965)
            )
        st.liquid_next[st.i, 22] = st.liquid[st.i, 22] + st.clock.dt * (
            -st.absorption_rate[st.i, 22]
            + st.brush_enzymes.sucrase_transformation[st.i, 21] * (180.15588 / 342.2965)
            + st.brush_enzymes.sucrase_transformation[st.i, 36] * (180.15588 / 342.2965)
        )
        aux = (
            st.concentration_lumen[st.i, 22]
            * (1 + st.surface_factor[st.i] * st.diffusion_time[st.i, 22] * drag * st.absorption_rate[st.i, 1])
            + st.diffusion_time[st.i, 22] * st.surface_factor[st.i]
            * (-st.absorption_rate[st.i, 22]
               + st.brush_enzymes.sucrase_transformation[st.i, 21] * (180.15588 / 342.2965)
               + st.brush_enzymes.sucrase_transformation[st.i, 36] * (180.15588 / 342.2965))
        )
        st.concentration_brush[st.i, 22] = divisor * (decelerator * st.concentration_brush[st.i, 22] + aux)

        st.liquid_next[st.i, 23] = st.liquid[st.i, 23] + st.clock.dt * (
            -st.absorption_rate[st.i, 23] - st.brush_enzymes.maltase_transformation[st.i, 23]
        )
        st.liquid_next[st.i, 24] = st.liquid[st.i, 24] + st.clock.dt * (
            -st.absorption_rate[st.i, 24] - st.brush_enzymes.maltase_transformation[st.i, 24]
        )
        st.liquid_next[st.i, 51] = st.liquid[st.i, 51] + st.absorption_rate[st.i, 54] * st.clock.dt
        st.liquid_next[st.i, 54] = st.liquid[st.i, 54] - st.absorption_rate[st.i, 54] * st.clock.dt

        for st.j in range(1, st.comp + 1):
            st.liquid[st.i, st.j] = st.liquid_next[st.i, st.j]
    volumina()


def colon():
    st.liquid[16, 59] = st.liquid[16, 59] + st.liquid[16, 51]
    st.liquid[16, 51] = 0
    st.liquid[16, 59] = st.liquid[16, 59] - 0.01 * st.liquid[16, 59] * st.clock.dt


def gut_hormones():
    delay = 0.0
    change = 0.0
    aux = 0.0
    aux1 = 0.0
    delay_entry = 0.0
    rateCCKrelease = np.zeros(17)
    ratePYYrelease = np.zeros(17)
    rateGIPrelease = np.zeros(17)
    rateGLP1release = np.zeros(17)
    rateGhrelinrelease = np.zeros(17)
    rateMotilinrelease = np.zeros(17)
    rateSecretinrelease = np.zeros(17)
    rateCCKproduction = np.zeros(17)
    ratePYYproduction = np.zeros(17)
    rateGIPproduction = np.zeros(17)
    rateGLP1production = np.zeros(17)
    rateGhrelinproduction = np.zeros(17)
    rateMotilinproduction = np.zeros(17)
    rateSecretinproduction = np.zeros(17)

    delay = 0.5
    delay = 1 / delay
    delay_entry = 0.5
    delay_entry = 1 / delay_entry
    hormones.total_rate_cck_release = 0
    hormones.total_rate_gip_release = 0
    hormones.total_rate_pyy_release = 0
    hormones.total_rate_glp1_release = 0
    hormones.total_rate_ghrelin_release = 0
    hormones.total_rate_motilin_release = 0
    hormones.total_rate_secretin_release = 0

    for st.i in range(5, 16):
        rateCCKproduction[st.i] = 0
        ratePYYproduction[st.i] = 0
        rateGIPproduction[st.i] = 0
        rateGLP1production[st.i] = 0
        rateGhrelinproduction[st.i] = 0
        rateMotilinproduction[st.i] = 0
        rateSecretinproduction[st.i] = 0

    for st.i in range(5, 16):
        for st.j in range(1, st.comp + 1):
            st.concentration_cell[st.i, st.j] = st.concentration_cell[st.i, st.j] + (
                st.concentration_brush[st.i, st.j] - st.concentration_cell[st.i, st.j]
            ) * delay_entry * st.clock.dt

    for st.i in range(5, 16):
        aux1 = math.exp(4.5 - st.ph[st.i])
        rateSecretinproduction[st.i] = rateSecretinproduction[st.i] + 2 * st.s_cell[st.i] * aux1 / (1 + aux1)
        for st.j in range(1, st.comp + 1):
            rateCCKproduction[st.i] = (
                rateCCKproduction[st.i] + hormones.cck_response_factor[st.j] * st.i_cell[st.i] * st.concentration_cell[st.i, st.j]
                + 0.5 * st.proteolytic_susceptibility[st.j] * st.i_cell[st.i] * st.liquid[st.i, st.j]
            )
            ratePYYproduction[st.i] = (
                ratePYYproduction[st.i] + hormones.pyy_response_factor[st.j] * st.l_cell[st.i] * st.concentration_cell[st.i, st.j]
            )
            rateGIPproduction[st.i] = (
                rateGIPproduction[st.i] + hormones.gip_response_factor[st.j] * st.k_cell[st.i] * st.concentration_cell[st.i, st.j]
                / (1 + st.concentration_cell[st.i, st.j] / 1e-5)
            )
            aux = 1e2 * st.concentration_cell[st.i, st.j]
            rateGLP1production[st.i] = (
                rateGLP1production[st.i] + hormones.glp1_response_factor[st.j] * st.l_cell[st.i] * 5 * aux / (1 + aux)
            )
            rateGhrelinproduction[st.i] = (
                rateGhrelinproduction[st.i] + hormones.ghrelin_response_factor[st.j] * st.g_cell[st.i] * st.concentration_cell[st.i, st.j]
            )
            rateMotilinproduction[st.i] = (
                rateMotilinproduction[st.i] + hormones.motilin_response_factor[st.j] * st.mo_cell[st.i] * st.concentration_cell[st.i, st.j]
            )

        rateCCKproduction[st.i] = rateCCKproduction[st.i] * (1 - 0.2 * (hormones.pyy / (1 + hormones.pyy)))
        rateGIPproduction[st.i] = rateGIPproduction[st.i]
        ratePYYproduction[st.i] = ratePYYproduction[st.i] + 0 * hormones.total_rate_cck_release
        rateGLP1production[st.i] = rateGLP1production[st.i] + 0 * rateGIPproduction[st.i]

        hormones.mesenteric_cck[st.i] = hormones.mesenteric_cck[st.i] + rateCCKproduction[st.i] * st.clock.dt
        rateCCKrelease[st.i] = delay * hormones.mesenteric_cck[st.i]
        change = rateCCKrelease[st.i] * st.clock.dt
        hormones.mesenteric_cck[st.i] = hormones.mesenteric_cck[st.i] - change

        hormones.mesenteric_pyy[st.i] = hormones.mesenteric_pyy[st.i] + ratePYYproduction[st.i] * st.clock.dt
        ratePYYrelease[st.i] = delay * hormones.mesenteric_pyy[st.i]
        change = ratePYYrelease[st.i] * st.clock.dt
        hormones.mesenteric_pyy[st.i] = hormones.mesenteric_pyy[st.i] - change

        hormones.mesenteric_gip[st.i] = hormones.mesenteric_gip[st.i] + rateGIPproduction[st.i] * st.clock.dt
        rateGIPrelease[st.i] = delay * hormones.mesenteric_gip[st.i]
        change = rateGIPrelease[st.i] * st.clock.dt
        hormones.mesenteric_gip[st.i] = hormones.mesenteric_gip[st.i] - change

        hormones.mesenteric_glp1[st.i] = hormones.mesenteric_glp1[st.i] + rateGLP1production[st.i] * st.clock.dt
        rateGLP1release[st.i] = delay * hormones.mesenteric_glp1[st.i]
        change = rateGLP1release[st.i] * st.clock.dt
        hormones.mesenteric_glp1[st.i] = hormones.mesenteric_glp1[st.i] - change

        hormones.mesenteric_ghrelin[st.i] = hormones.mesenteric_ghrelin[st.i] + rateGhrelinproduction[st.i] * st.clock.dt
        rateGhrelinrelease[st.i] = delay * hormones.mesenteric_ghrelin[st.i]
        change = rateGhrelinrelease[st.i] * st.clock.dt
        hormones.mesenteric_ghrelin[st.i] = hormones.mesenteric_ghrelin[st.i] - change

        hormones.mesenteric_motilin[st.i] = hormones.mesenteric_motilin[st.i] + rateMotilinproduction[st.i] * st.clock.dt
        rateMotilinrelease[st.i] = delay * hormones.mesenteric_motilin[st.i]
        change = rateMotilinrelease[st.i] * st.clock.dt
        hormones.mesenteric_motilin[st.i] = hormones.mesenteric_motilin[st.i] - change

        hormones.mesenteric_secretin[st.i] = hormones.mesenteric_secretin[st.i] + rateSecretinproduction[st.i] * st.clock.dt
        rateSecretinrelease[st.i] = delay * hormones.mesenteric_secretin[st.i]
        change = rateSecretinrelease[st.i] * st.clock.dt
        hormones.mesenteric_secretin[st.i] = hormones.mesenteric_secretin[st.i] - change

        hormones.total_rate_cck_release = hormones.total_rate_cck_release + rateCCKrelease[st.i]
        hormones.total_rate_gip_release = hormones.total_rate_gip_release + rateGIPrelease[st.i]
        hormones.total_rate_pyy_release = hormones.total_rate_pyy_release + ratePYYrelease[st.i]
        hormones.total_rate_glp1_release = hormones.total_rate_glp1_release + rateGLP1release[st.i]
        hormones.total_rate_ghrelin_release = hormones.total_rate_ghrelin_release + rateGhrelinrelease[st.i]
        hormones.total_rate_motilin_release = hormones.total_rate_motilin_release + rateMotilinrelease[st.i]
        hormones.total_rate_secretin_release = hormones.total_rate_secretin_release + rateSecretinrelease[st.i]

    hormones.total_rate_ghrelin_release = hormones.total_rate_ghrelin_release + 10
    hormones.total_rate_motilin_release = hormones.total_rate_motilin_release + 10
    hormones.cck = hormones.cck + hormones.total_rate_cck_release * st.clock.dt
    hormones.pyy = hormones.pyy + hormones.total_rate_pyy_release * st.clock.dt
    hormones.gip = hormones.gip + hormones.total_rate_gip_release * st.clock.dt
    hormones.glp1 = hormones.glp1 + hormones.total_rate_glp1_release * st.clock.dt
    hormones.ghrelin = (
        hormones.ghrelin + hormones.total_rate_ghrelin_release * st.clock.dt * math.exp(-glucose.serum_insulin_mu_l / 30)
        * math.exp(-hormones.cck / 2) * math.exp(-hormones.glp1 / 2) * math.exp(-st.osmotic.lumen_duodenum / 0.02)
    )
    hormones.motilin = (
        hormones.motilin + hormones.total_rate_motilin_release * st.clock.dt * math.exp(-glucose.serum_insulin_mu_l / 20)
        * math.exp(-hormones.cck / 2) * math.exp(-hormones.glp1 / 0.2)
        * math.exp(-(abs(st.ph[6] - 6) / 6))
    )
    hormones.secretin = hormones.secretin + hormones.total_rate_secretin_release * st.clock.dt

    hormones.cck = hormones.cck - hormones.cck * hormones.cck_degradation_constant * st.clock.dt
    hormones.pyy = hormones.pyy - hormones.pyy * hormones.pyy_degradation_constant * st.clock.dt
    hormones.gip = hormones.gip - hormones.gip * hormones.gip_degradation_constant * st.clock.dt
    hormones.glp1 = hormones.glp1 - hormones.glp1 * hormones.glp1_degradation_constant * st.clock.dt
    hormones.ghrelin = hormones.ghrelin - hormones.ghrelin * hormones.ghrelin_degradation_constant * st.clock.dt
    hormones.motilin = hormones.motilin - hormones.motilin * hormones.motilin_degradation_constant * st.clock.dt
    hormones.secretin = hormones.secretin - hormones.secretin * hormones.secretin_degradation_constant * st.clock.dt

    hormones.real_glp1 = 0.5 * hormones.glp1
    hormones.real_gip = 1 * hormones.gip


def sensations():
    st.fullness = (10 / 15) * (stomach.tone - 2)
    st.hunger.total_caloric_absorption_rate = 10 - 25 * (st.calories.total_absorption_rate / st.calories.max_absorption_capacity)
    st.hunger.intestinal_volume = (
        st.volume[5] + st.volume[6] + st.volume[7] + st.volume[8] + st.volume[9]
        + st.volume[10] + st.volume[11] + st.volume[12] + st.volume[13] + st.volume[14]
    )
    st.hunger.intestinal_volume = 10 - st.hunger.intestinal_volume / 1000
    st.hunger.level = 0.8 * (st.hunger.total_caloric_absorption_rate + st.hunger.intestinal_volume) / 2
    if st.mmc_status == 3:
        st.hunger.level = 10


def intubation_in_antrum():
    if (st.clock.minutes > 250) and (st.clock.minutes < 251):
        st.liquid[4, 1] = st.liquid[4, 1] + 0 * st.clock.dt
        st.liquid[4, 3] = st.liquid[4, 3] + 0 * st.clock.dt


def intubation_in_duodenum():
    if (st.clock.minutes > 250) and (st.clock.minutes < 300):
        st.liquid[5, 1] = st.liquid[5, 1] + 0 * st.clock.dt
        st.liquid[5, 3] = st.liquid[5, 3] + 0 * st.clock.dt


def intravenous():
    glucose.g2_in = 1
    glucose.g2_in_mol_per_min = 0.001 / 180.15588



