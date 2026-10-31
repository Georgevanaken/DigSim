"""Default values of the model parameters."""
from . import state as st
from .state import hormones

i = 0
j = 0


def set_model_parameters() -> None:
    """Set all model parameters (transit, absorption, enzyme and cell-layer distributions) to their defaults."""
    global i, j

    st.base_si_transit_rate = 0.02  # 0.02; Basis Small Intestinal Transit rate in min-1
    hormones.cck_effect_on_gastric_emptying = 2
    st.osmotic.concentration = 200
    st.max_segment_volume = 50  # maximum volume of an SI segment

    for i in range(1, 17):
        st.hypertrophy[i] = 1
    st.hypertrophy[5] = 1
    st.hypertrophy[6] = 1

    for j in range(1, 81):
        hormones.cck_response_factor[j] = 0
        hormones.pyy_response_factor[j] = 0
        hormones.glp1_response_factor[j] = 0
        hormones.gip_response_factor[j] = 0
        hormones.ghrelin_response_factor[j] = 0
        hormones.motilin_response_factor[j] = 0
    hormones.cck_response_factor[20] = 2000
    hormones.cck_response_factor[54] = 2000
    hormones.cck_response_factor[64] = 200
    hormones.cck_response_factor[65] = 100

    hormones.pyy_response_factor[3] = 100
    hormones.pyy_response_factor[34] = 100
    hormones.pyy_response_factor[35] = 100
    hormones.pyy_response_factor[20] = 200
    hormones.pyy_response_factor[54] = 200
    hormones.pyy_response_factor[64] = 100
    hormones.pyy_response_factor[65] = 100

    hormones.glp1_response_factor[3] = 1
    hormones.glp1_response_factor[34] = hormones.glp1_response_factor[3]
    hormones.glp1_response_factor[35] = hormones.glp1_response_factor[3]
    hormones.glp1_response_factor[20] = 0
    hormones.glp1_response_factor[21] = 0
    hormones.glp1_response_factor[36] = hormones.glp1_response_factor[21]
    hormones.glp1_response_factor[22] = 0
    hormones.glp1_response_factor[54] = 0
    hormones.glp1_response_factor[64] = 0.1
    hormones.glp1_response_factor[65] = 0.1

    hormones.gip_response_factor[3] = 1000
    hormones.gip_response_factor[34] = hormones.gip_response_factor[3]
    hormones.gip_response_factor[35] = hormones.gip_response_factor[3]
    hormones.gip_response_factor[20] = 10
    hormones.gip_response_factor[21] = 0
    hormones.gip_response_factor[22] = 5
    hormones.gip_response_factor[25] = 5
    hormones.gip_response_factor[26] = 5
    hormones.gip_response_factor[54] = 20

    hormones.ghrelin_response_factor[3] = 0
    hormones.ghrelin_response_factor[34] = 0
    hormones.ghrelin_response_factor[35] = 0
    hormones.ghrelin_response_factor[20] = 0
    hormones.ghrelin_response_factor[64] = 0
    hormones.ghrelin_response_factor[65] = 0
    hormones.ghrelin_response_factor[25] = 0
    hormones.ghrelin_response_factor[26] = 0
    hormones.ghrelin_response_factor[54] = 0

    hormones.motilin_response_factor[3] = 0
    hormones.motilin_response_factor[34] = 0
    hormones.motilin_response_factor[35] = 0
    hormones.motilin_response_factor[20] = 0
    hormones.motilin_response_factor[64] = 0
    hormones.motilin_response_factor[65] = 0
    hormones.motilin_response_factor[25] = 0
    hormones.motilin_response_factor[26] = 0
    hormones.motilin_response_factor[51] = 0
    hormones.motilin_response_factor[54] = 0

    hormones.secretin_response_factor[3] = 0
    hormones.secretin_response_factor[34] = 0
    hormones.secretin_response_factor[35] = 0
    hormones.secretin_response_factor[20] = 0
    hormones.secretin_response_factor[64] = 0
    hormones.secretin_response_factor[65] = 0
    hormones.secretin_response_factor[25] = 0
    hormones.secretin_response_factor[26] = 0
    hormones.secretin_response_factor[51] = 0
    hormones.secretin_response_factor[54] = 0

    for j in range(1, st.comp + 1):
        st.proteolytic_susceptibility[j] = 0
    st.proteolytic_susceptibility[4] = 1
    st.proteolytic_susceptibility[61] = 1
    st.proteolytic_susceptibility[62] = 1
    st.proteolytic_susceptibility[63] = 1
    st.proteolytic_susceptibility[64] = 0
    st.proteolytic_susceptibility[65] = 0
    st.proteolytic_susceptibility[66] = 10
    st.proteolytic_susceptibility[67] = 10

    for i in range(1, 16):
        st.i_cell[i] = 0
    st.i_cell[5] = 1 * st.hypertrophy[5]
    st.i_cell[6] = 1 * st.hypertrophy[6]
    st.i_cell[7] = 1 * st.hypertrophy[7]
    st.i_cell[8] = 0.2 * st.hypertrophy[8]
    st.i_cell[9] = 0.1 * st.hypertrophy[9]

    for i in range(1, 16):
        st.l_cell[i] = 0
    st.l_cell[5] = 0 * st.hypertrophy[5]
    st.l_cell[6] = 0.1 * st.hypertrophy[6]
    st.l_cell[7] = 0.2 * st.hypertrophy[7]
    st.l_cell[8] = 0.3 * st.hypertrophy[8]
    st.l_cell[9] = 0.4 * st.hypertrophy[9]
    st.l_cell[10] = 0.5 * st.hypertrophy[10]
    st.l_cell[11] = 0.6 * st.hypertrophy[11]
    st.l_cell[12] = 0.8 * st.hypertrophy[12]
    st.l_cell[13] = 1.0 * st.hypertrophy[13]
    st.l_cell[14] = 2 * st.hypertrophy[14]
    st.l_cell[15] = 5 * st.hypertrophy[15]

    for i in range(1, 16):
        st.k_cell[i] = 0
    st.k_cell[5] = 2 * st.hypertrophy[5]
    st.k_cell[6] = 1.5 * st.hypertrophy[6]
    st.k_cell[7] = 1 * st.hypertrophy[7]
    st.k_cell[8] = 0.5 * st.hypertrophy[8]
    st.k_cell[9] = 0.2 * st.hypertrophy[9]
    st.k_cell[10] = 0.1 * st.hypertrophy[10]

    for i in range(1, 16):
        st.g_cell[i] = 0
    st.g_cell[3] = 0 * st.hypertrophy[3]
    st.g_cell[4] = 0 * st.hypertrophy[4]
    st.g_cell[5] = 0 * st.hypertrophy[5]
    st.g_cell[6] = 0 * st.hypertrophy[6]

    for i in range(1, 16):
        st.mo_cell[i] = 0
    st.mo_cell[5] = 2 * st.hypertrophy[5]
    st.mo_cell[6] = 5 * st.hypertrophy[6]
    st.mo_cell[7] = 1 * st.hypertrophy[7]
    st.mo_cell[8] = 1 * st.hypertrophy[8]
    st.mo_cell[9] = 1 * st.hypertrophy[9]
    st.mo_cell[10] = 1 * st.hypertrophy[10]
    st.mo_cell[11] = 1 * st.hypertrophy[11]
    st.mo_cell[12] = 1 * st.hypertrophy[12]
    st.mo_cell[13] = 1 * st.hypertrophy[13]
    st.mo_cell[14] = 1 * st.hypertrophy[14]
    st.mo_cell[15] = 1 * st.hypertrophy[15]

    for i in range(1, 16):
        st.s_cell[i] = 0
    st.s_cell[5] = 1 * st.hypertrophy[5]
    st.s_cell[6] = 1 * st.hypertrophy[6]
    st.s_cell[7] = 0 * st.hypertrophy[7]
    st.s_cell[8] = 0 * st.hypertrophy[8]
    st.s_cell[9] = 0 * st.hypertrophy[9]
    st.s_cell[10] = 0 * st.hypertrophy[10]
    st.s_cell[11] = 0 * st.hypertrophy[11]
    st.s_cell[12] = 0 * st.hypertrophy[12]
    st.s_cell[13] = 0 * st.hypertrophy[13]
    st.s_cell[14] = 0 * st.hypertrophy[14]
    st.s_cell[15] = 0 * st.hypertrophy[15]

    hormones.cck_degradation_constant = 0.5
    hormones.pyy_degradation_constant = 0.1
    hormones.gip_degradation_constant = 1 / 7.2
    hormones.glp1_degradation_constant = 1 / 2.3
    hormones.ghrelin_degradation_constant = 1 / 2
    hormones.motilin_degradation_constant = 1 / 15
    hormones.secretin_degradation_constant = 0.4

    hormones.cck_attenuation_si = 0
    hormones.gip_attenuation_si = 0
    hormones.pyy_attenuation_si = 0.9
    hormones.glp1_attenuation_si = 0.01
    hormones.ghrelin_attenuation_si = 0
    hormones.motilin_attenuation_si = 2
    hormones.secretin_attenuation_si = 0.5

    st.lipolysis.lingual_cst = 0
    st.lipolysis.gastric_cst = 6
    st.lipolysis.pancreatic_cst = 6

    for j in range(1, st.comp + 1):
        st.density[j] = 1e20
    st.density[1] = 1
    st.density[2] = 0.9
    st.density[3] = 1.54
    st.density[34] = 1.54
    st.density[35] = 1.54
    st.density[4] = 1.4
    st.density[5] = 2
    st.density[6] = 2
    st.density[10] = 1.54
    st.density[11] = 1.54
    st.density[12] = 1.54
    st.density[20] = 0.9
    st.density[21] = 1.59
    st.density[36] = 1.59
    st.density[22] = 1.69
    st.density[23] = 1.69
    st.density[25] = 1.69
    st.density[26] = 1.69
    st.density[27] = 1.69
    st.density[54] = 1.69
    st.density[61] = 1.2
    st.density[62] = 1.2
    st.density[63] = 1.2
    st.density[64] = 1.2
    st.density[65] = 1.2
    st.density[66] = 1.2
    st.density[67] = 1.2
    st.density[68] = 1.66
    st.density[78] = 2.16
    st.density[80] = 1e6

    for j in range(1, st.comp + 1):
        st.calories.value[j] = 0
    st.calories.value[2] = 9.54
    st.calories.value[3] = 3.86
    st.calories.value[34] = 3.86
    st.calories.value[35] = 3.86
    st.calories.value[4] = 4.76
    st.calories.value[10] = 3.86
    st.calories.value[11] = 3.86
    st.calories.value[12] = 3.86
    st.calories.value[20] = 9.54
    st.calories.value[21] = 3.86
    st.calories.value[36] = 3.86
    st.calories.value[22] = 3.86
    st.calories.value[23] = 3.86
    st.calories.value[25] = 0
    st.calories.value[26] = 0
    st.calories.value[27] = 0
    st.calories.value[54] = st.calories.value[20]
    st.calories.value[61] = st.calories.value[4]
    st.calories.value[62] = st.calories.value[4]
    st.calories.value[63] = st.calories.value[4]
    st.calories.value[64] = st.calories.value[4]
    st.calories.value[65] = st.calories.value[4]
    st.calories.value[66] = st.calories.value[4]
    st.calories.value[67] = st.calories.value[4]

    for i in range(1, 17):
        for j in range(1, st.comp + 1):
            st.km_absorption[i, j] = 0.0
            st.vmax_absorption[i, j] = 0
    st.vmax_absorption[16, 1] = 10
    st.km_absorption[16, 1] = 1
    st.water_exchange_rate = 1

    for i in range(5, 16):
        st.km_absorption[i, 1] = 1
        st.vmax_absorption[i, 1] = 100
        st.km_absorption[i, 20] = 0.1
        st.vmax_absorption[i, 20] = 5
        st.km_absorption[i, 21] = 0
        st.vmax_absorption[i, 21] = 0
        st.km_absorption[i, 54] = 0.1
        st.vmax_absorption[i, 54] = 5
        st.vmax_absorption[i, 3] = 3.48 / 0.2
        st.km_absorption[i, 3] = 1 * 0.022556
        st.vmax_absorption[i, 34] = st.vmax_absorption[i, 3]
        st.vmax_absorption[i, 35] = st.vmax_absorption[i, 3]
        st.km_absorption[i, 34] = st.km_absorption[i, 3]
        st.km_absorption[i, 35] = st.km_absorption[i, 3]
        st.vmax_absorption[i, 22] = 0.5 * 3.48 / 1
        st.km_absorption[i, 22] = 1 * 0.022556
        st.km_absorption[i, 25] = 1
        st.vmax_absorption[i, 25] = 12 * st.km_absorption[i, 25]
        st.km_absorption[i, 26] = 1
        st.vmax_absorption[i, 26] = 0.02 * st.km_absorption[i, 26]
        st.vmax_absorption[i, 27] = 0.0
        st.km_absorption[i, 27] = 0.007
        st.vmax_absorption[i, 4] = 0
        st.km_absorption[i, 4] = 0.02
        st.vmax_absorption[i, 64] = 0.3
        st.vmax_absorption[i, 65] = 0.3
        st.km_absorption[i, 64] = 0.02
        st.km_absorption[i, 65] = 0.02
        st.vmax_absorption[i, 78] = 0
        st.km_absorption[i, 78] = 1

    for i in range(5, 16):
        st.vmax_absorption[i, 5] = 1
        st.km_absorption[i, 5] = 1

    for i in range(10, 16):
        st.vmax_absorption[i, 25] = 0 * st.vmax_absorption[i, 25]
        st.vmax_absorption[i, 26] = 0 * st.vmax_absorption[i, 26]
        st.km_absorption[i, 51] = 0.1
        st.vmax_absorption[i, 51] = 0.02

    for j in range(1, st.comp + 1):
        st.vmax_absorption[5, j] = 1 * st.vmax_absorption[5, j]
        st.vmax_absorption[6, j] = 1 * st.vmax_absorption[6, j]
        st.vmax_absorption[7, j] = 1 * st.vmax_absorption[7, j]
        st.vmax_absorption[8, j] = 1 * st.vmax_absorption[8, j]
        st.vmax_absorption[9, j] = 1 * st.vmax_absorption[9, j]
        st.vmax_absorption[10, j] = 0.5 * st.vmax_absorption[10, j]
        st.vmax_absorption[11, j] = 0.5 * st.vmax_absorption[11, j]
        st.vmax_absorption[12, j] = 0.5 * st.vmax_absorption[12, j]
        st.vmax_absorption[13, j] = 0.5 * st.vmax_absorption[13, j]
        st.vmax_absorption[14, j] = 0.5 * st.vmax_absorption[14, j]
        st.vmax_absorption[15, j] = 0.5 * st.vmax_absorption[15, j]
    for i in range(5, 16):
        for j in range(1, st.comp + 1):
            st.vmax_absorption[i, j] = st.vmax_absorption[i, j] * st.hypertrophy[i]

    for i in range(5, 17):
        st.delta[i] = 1 * 5e-2
    for j in range(1, st.comp + 1):
        st.diffusion_coefficient[j] = 1e-20
    st.diffusion_coefficient[1] = 2.27e-5 * 60
    st.diffusion_coefficient[3] = 9.4e-6 * 60
    st.diffusion_coefficient[34] = st.diffusion_coefficient[3]
    st.diffusion_coefficient[35] = st.diffusion_coefficient[3]
    st.diffusion_coefficient[20] = 6e-6 * 60
    st.diffusion_coefficient[21] = 7.056e-6 * 60
    st.diffusion_coefficient[36] = st.diffusion_coefficient[21]
    st.diffusion_coefficient[22] = 10e-6 * 60
    st.diffusion_coefficient[23] = 4e-6 * 60
    st.diffusion_coefficient[24] = st.diffusion_coefficient[23]
    st.diffusion_coefficient[25] = 10e-6 * 60
    st.diffusion_coefficient[26] = 10e-6 * 60
    st.diffusion_coefficient[27] = 10e-6 * 60
    st.diffusion_coefficient[5] = 6e-6 * 60
    st.diffusion_coefficient[51] = 6e-6 * 60
    st.diffusion_coefficient[52] = 4e-6 * 60
    st.diffusion_coefficient[53] = 4e-6 * 60
    st.diffusion_coefficient[54] = 4e-6 * 60
    st.diffusion_coefficient[64] = 6.7e-6 * 60
    st.diffusion_coefficient[65] = 6.7e-6 * 60
    st.diffusion_coefficient[78] = 1.33e-5 * 60
    st.diffusion_coefficient[80] = 20e-5 * 60

    for i in range(5, 16):
        for j in range(1, st.comp + 1):
            st.surface_factor[i] = 1 / (3 * 100 * 100 * st.hypertrophy[i])
            st.diffusion_time[i, j] = st.delta[i] / st.diffusion_coefficient[j]
    for i in range(1, 5):
        st.max_total_absorption_rate[i] = 0
    for i in range(5, 10):
        st.max_total_absorption_rate[i] = 1
    for i in range(10, 16):
        st.max_total_absorption_rate[i] = 0.5

    st.calories.max_absorption_capacity = 0
    for i in range(1, 16):
        st.calories.max_absorption_capacity = (
            st.calories.max_absorption_capacity + st.max_total_absorption_rate[i]
        )
