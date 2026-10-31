"""Absorption, brush-border layer, colon, gut hormones and sensations (third part of the main program)."""
import math

import numpy as np

from . import state as st
from .brush_border_enzymes import drag, maltase, sucrase
from .compartment_properties import volumina
from .enums import SMALL_INTESTINE, MMCStatus, Segment
from .enums import Component as C
from .state import bile, glucose, hormones, stomach


def _seq_sum(a: np.ndarray) -> np.ndarray:
    """Sum along axis 0 strictly in order (same rounding as the original scalar loops)."""
    return np.cumsum(a, axis=0)[-1]


def absorption_in_small_intestine() -> None:
    """Transport between lumen, brush border and enterocytes: absorption, enzymes and cell metabolism."""
    R = 0.8  # Recovery of glucose after absorption by enterocytes; some fuels enterocytes.
    seg = SMALL_INTESTINE
    cols = slice(2, st.comp + 1)
    allc = slice(1, st.comp + 1)
    cb = st.concentration_brush[seg, :]
    rate = np.zeros((11, st.absorption_rate.shape[1]))
    glucoses = (C.GLUCOSE, C.TRACER_GLUCOSE_1, C.TRACER_GLUCOSE_2)

    with np.errstate(all="ignore"):
        vmax = st.vmax_absorption[seg, :]
        km = st.km_absorption[seg, :]
        active = cb > 1e-40
        active[:, C.H_ION] = False
        generic = np.where(active, vmax * cb / (cb + km), 0.0)
        rate[:, cols] = generic[:, cols]

        # The three glucose species share one transporter.
        for j in glucoses:
            rate[:, j] = 0
        aux = cb[:, C.GLUCOSE] + cb[:, C.TRACER_GLUCOSE_1] + cb[:, C.TRACER_GLUCOSE_2]
        aux = aux * (1 + 100 * cb[:, C.NACL])
        rate[:, C.GLUCOSE] = np.where(
            aux > 0,
            vmax[:, C.GLUCOSE] * aux / (aux + km[:, C.GLUCOSE]) + 0.5 * aux / (aux + 0.1 * 0.0225560),
            0.0,
        )
        st.scratch = float(aux[-1])

        # Weber competition among macronutrients.
        total_target = _seq_sum(rate[:, cols].T)
        maxrate = st.max_total_absorption_rate[seg]
        caloric = st.calories.value[cols] != 0
        factor = maxrate / (total_target + maxrate)
        scaled = np.where(caloric, rate[:, cols] * factor[:, None], rate[:, cols])
        rate[:, cols] = np.where((maxrate == 0)[:, None], 0.0, scaled)

        nonzero = aux != 0
        gar = np.where(nonzero, rate[:, C.GLUCOSE], glucose.absorption_rate[seg])
        glucose.absorption_rate[seg] = gar
        for j in glucoses:
            rate[:, j] = np.where(nonzero, cb[:, j] / (aux + 1e-3) * gar, rate[:, j])
        rate[:, C.LEUCINE] = 0.15 * rate[:, C.CASEIN_PEPTIDES]
        rate[:, C.LEUCINE_TRACER] = rate[:, C.SMALL_PEPTIDES]

        calv = st.calories.value[allc]
        rate[:, allc] = np.where(calv != 0, rate[:, allc] / np.where(calv != 0, calv, 1.0), rate[:, allc])

        sugars = rate[:, C.GLUCOSE] + rate[:, C.TRACER_GLUCOSE_1] + rate[:, C.TRACER_GLUCOSE_2]
        rate[:, C.WATER] = (
            2 * (st.volume[seg] / (st.volume[seg] + 100)) * (0.2 - st.osmotic.brush[seg])
            + (30 * 18 / 180.15588) * sugars
        )
        rate[:, C.NACL] = (1 * 58.44 / 180.15588) * sugars

    st.absorption_rate[seg, allc] = rate[:, allc]

    total = 0.0
    for i in range(rate.shape[0]):
        total = total + rate[i, C.GLUCOSE] + rate[i, C.TRACER_GLUCOSE_1] + rate[i, C.TRACER_GLUCOSE_2]
    glucose.total_absorption_rate = total / 180.15588
    bile.conjugated_salt_absorption_rate = 0

    rows = rate[:, allc]
    st.total_absorption_rate[allc] = _seq_sum(rows)
    st.absorbed[allc] = _seq_sum(rows * st.clock.dt)
    caloric_terms = rows * st.calories.value[allc]
    st.calories.total_absorption_rate = float(np.cumsum(caloric_terms.T.ravel())[-1])
    calorie_terms = st.calories.value[allc] * 1 * st.absorbed[allc]
    st.calories.absorbed = float(np.cumsum(np.concatenate(([st.calories.absorbed], calorie_terms)))[-1])

    delay = 1 / 2
    mes = st.mesenteric[seg, :]
    mes_old22 = mes[:, C.FRUCTOSE].copy()
    mes[:, allc] = mes[:, allc] + rate[:, allc] * st.clock.dt
    change = delay * mes[:, allc] * st.clock.dt
    change = np.minimum(change, mes[:, allc])
    st.blood[allc] = np.cumsum(np.vstack([st.blood[allc], change]), axis=0)[-1]
    mes[:, allc] = mes[:, allc] - change

    G0in = G1in = G2in = Fruc = 0.0
    for i in range(mes.shape[0]):
        # For glucose the fructose pool is still the one from before this step.
        to_glucose_before = 0.5 * mes_old22[i]
        G0in = G0in + delay * mes[i, C.GLUCOSE] * 1000 + to_glucose_before * 1000
        G1in = G1in + delay * mes[i, C.TRACER_GLUCOSE_1] * 1000
        G2in = G2in + delay * mes[i, C.TRACER_GLUCOSE_2] * 1000
        to_glucose_after = 0.5 * mes[i, C.FRUCTOSE]
        Fruc = Fruc + delay * mes[i, C.FRUCTOSE] * 1000 - to_glucose_after * 1000
    glucose.g0_in = R * G0in
    glucose.g1_in = R * G1in
    glucose.g2_in = R * G2in
    glucose.fructose_decay_rate = Fruc
    st.total_water_absorption_rate = st.total_absorption_rate[C.WATER]


def mucus_layer() -> None:
    """Run the brush-border enzymes and compute the brush-border concentrations from the lumen concentrations."""
    decelerator = 5
    divisor = 1 / (decelerator + 1)
    sucrase()
    maltase()
    volumina()

    seg = SMALL_INTESTINE
    allc = slice(1, st.comp + 1)
    dt = st.clock.dt
    c = 180.15588 / 342.2965
    ar = st.absorption_rate[seg, :]
    ar1 = ar[:, C.WATER]
    sf = st.surface_factor[seg]
    T = st.diffusion_time[seg, :]
    cl = st.concentration_lumen[seg, :]
    cb = st.concentration_brush[seg, :]
    liq = st.liquid[seg, :]
    suc = st.brush_enzymes.sucrase_transformation[seg, :]
    mal = st.brush_enzymes.maltase_transformation[seg, :]

    # Brush-border concentrations for the ordinary components. The sugars below and
    # the H+ ion are handled separately.
    special = [C.GLUCOSE, C.SUCROSE, C.FRUCTOSE, C.MALTOSE, C.TRACER_GLUCOSE_1,
               C.TRACER_GLUCOSE_2, C.TRACER_SUCROSE, C.H_ION]
    ordinary = np.zeros(st.concentration_brush.shape[1], dtype=bool)
    ordinary[2:st.comp + 1] = True
    ordinary[special] = False
    aux = (
        cl * (1 + sf[:, None] * drag * T * ar1[:, None])
        - T * sf[:, None] * ar
    )
    cb_new = divisor * (decelerator * cb + aux)
    cb[:, ordinary] = cb_new[:, ordinary]

    ln = liq - ar * dt
    ln[:, 2:st.comp + 1] = np.maximum(ln[:, 2:st.comp + 1], 0)
    ln[:, C.H_ION] = liq[:, C.H_ION] - ar[:, C.H_ION] * dt

    def brush(j: int, flux: np.ndarray) -> None:
        # flux is the net amount added to the brush-border layer by transformation and uptake
        h = cl[:, j] * (1 + sf * T[:, j] * drag * ar1) + T[:, j] * sf * flux
        cb[:, j] = divisor * (decelerator * cb[:, j] + h)

    # Glucose is released from sucrose (one half) and maltose (two halves) by the brush-border enzymes.
    glucose_in = -ar[:, C.GLUCOSE] + suc[:, C.SUCROSE] * c + mal[:, C.MALTOSE] * 2 * c
    brush(C.GLUCOSE, glucose_in)
    ln[:, C.GLUCOSE] = liq[:, C.GLUCOSE] + dt * glucose_in

    tracer_in = -ar[:, C.TRACER_GLUCOSE_1] + suc[:, C.TRACER_SUCROSE] * c + mal[:, C.TRACER_MALTOSE] * 2 * c
    brush(C.TRACER_GLUCOSE_1, tracer_in)
    ln[:, C.TRACER_GLUCOSE_1] = np.maximum(liq[:, C.TRACER_GLUCOSE_1] + dt * tracer_in, 1e-10)

    brush(C.TRACER_GLUCOSE_2, -ar[:, C.TRACER_GLUCOSE_2])

    ln[:, C.SUCROSE] = liq[:, C.SUCROSE] + dt * (-ar[:, C.SUCROSE] - suc[:, C.SUCROSE])
    ln[:, C.TRACER_SUCROSE] = liq[:, C.TRACER_SUCROSE] + dt * (-ar[:, C.TRACER_SUCROSE] - suc[:, C.TRACER_SUCROSE])
    sucrose_fructose = suc[:, C.SUCROSE] * c
    tracer_fructose = suc[:, C.TRACER_SUCROSE] * c
    projected = liq[:, C.FRUCTOSE] + dt * (-ar[:, C.FRUCTOSE] + sucrose_fructose + tracer_fructose)
    ar[:, C.FRUCTOSE] = np.where(
        projected < 0, 0.95 * liq[:, C.FRUCTOSE] / dt + sucrose_fructose + tracer_fructose, ar[:, C.FRUCTOSE])
    ln[:, C.FRUCTOSE] = liq[:, C.FRUCTOSE] + dt * (-ar[:, C.FRUCTOSE] + sucrose_fructose + tracer_fructose)
    brush(C.FRUCTOSE, -ar[:, C.FRUCTOSE] + sucrose_fructose + tracer_fructose)

    ln[:, C.MALTOSE] = liq[:, C.MALTOSE] + dt * (-ar[:, C.MALTOSE] - mal[:, C.MALTOSE])
    ln[:, C.TRACER_MALTOSE] = liq[:, C.TRACER_MALTOSE] + dt * (-ar[:, C.TRACER_MALTOSE] - mal[:, C.TRACER_MALTOSE])
    ln[:, C.CONJUGATED_BILE_SALT] = liq[:, C.CONJUGATED_BILE_SALT] + ar[:, C.BILE_SALT_LCFA_COMPLEX] * dt
    ln[:, C.BILE_SALT_LCFA_COMPLEX] = liq[:, C.BILE_SALT_LCFA_COMPLEX] - ar[:, C.BILE_SALT_LCFA_COMPLEX] * dt

    st.liquid_next[seg, allc] = ln[:, allc]
    st.liquid[seg, allc] = ln[:, allc]
    volumina()


def colon() -> None:
    """Deconjugate the bile salt in the colon and let part of it disappear."""
    colon_row = Segment.COLON
    st.liquid[colon_row, C.UNCONJUGATED_BILE_SALT] = (
        st.liquid[colon_row, C.UNCONJUGATED_BILE_SALT] + st.liquid[colon_row, C.CONJUGATED_BILE_SALT]
    )
    st.liquid[colon_row, C.CONJUGATED_BILE_SALT] = 0
    st.liquid[colon_row, C.UNCONJUGATED_BILE_SALT] = (
        st.liquid[colon_row, C.UNCONJUGATED_BILE_SALT]
        - 0.01 * st.liquid[colon_row, C.UNCONJUGATED_BILE_SALT] * st.clock.dt
    )


def gut_hormones() -> None:
    """Production of the gut hormones (CCK, PYY, GIP, GLP-1, ghrelin, motilin, secretin) and their release."""
    delay = 1 / 0.5
    delay_entry = 1 / 0.5
    seg = SMALL_INTESTINE
    comps = slice(1, st.comp + 1)

    cell = st.concentration_cell[seg, comps]
    cell += (st.concentration_brush[seg, comps] - cell) * delay_entry * st.clock.dt

    liquid = st.liquid[seg, comps]
    I_cell = st.i_cell[seg, None]
    L_cell = st.l_cell[seg, None]
    K_cell = st.k_cell[seg, None]
    G_cell = st.g_cell[seg, None]
    Mo_cell = st.mo_cell[seg, None]

    aux1 = np.exp(4.5 - st.ph[seg])
    rateSecretinproduction = 2 * st.s_cell[seg] * aux1 / (1 + aux1)
    rateCCKproduction = (
        hormones.cck_response_factor[comps] * I_cell * cell
        + 0.5 * st.proteolytic_susceptibility[comps] * I_cell * liquid
    ).sum(axis=1)
    ratePYYproduction = (hormones.pyy_response_factor[comps] * L_cell * cell).sum(axis=1)
    rateGIPproduction = (hormones.gip_response_factor[comps] * K_cell * cell / (1 + cell / 1e-5)).sum(axis=1)
    receptor = 1e2 * cell
    rateGLP1production = (hormones.glp1_response_factor[comps] * L_cell * 5 * receptor / (1 + receptor)).sum(axis=1)
    rateGhrelinproduction = (hormones.ghrelin_response_factor[comps] * G_cell * cell).sum(axis=1)
    rateMotilinproduction = (hormones.motilin_response_factor[comps] * Mo_cell * cell).sum(axis=1)

    rateCCKproduction = rateCCKproduction * (1 - 0.2 * (hormones.pyy / (1 + hormones.pyy)))

    def mesenteric_release(store: np.ndarray, production: float) -> float:
        # The mesentric cell layer accumulates the hormone and releases it with a delay.
        store[seg] = store[seg] + production * st.clock.dt
        release = delay * store[seg]
        store[seg] = store[seg] - release * st.clock.dt
        return release

    rateCCKrelease = mesenteric_release(hormones.mesenteric_cck, rateCCKproduction)
    ratePYYrelease = mesenteric_release(hormones.mesenteric_pyy, ratePYYproduction)
    rateGIPrelease = mesenteric_release(hormones.mesenteric_gip, rateGIPproduction)
    rateGLP1release = mesenteric_release(hormones.mesenteric_glp1, rateGLP1production)
    rateGhrelinrelease = mesenteric_release(hormones.mesenteric_ghrelin, rateGhrelinproduction)
    rateMotilinrelease = mesenteric_release(hormones.mesenteric_motilin, rateMotilinproduction)
    rateSecretinrelease = mesenteric_release(hormones.mesenteric_secretin, rateSecretinproduction)

    hormones.total_rate_cck_release = float(rateCCKrelease.sum())
    hormones.total_rate_gip_release = float(rateGIPrelease.sum())
    hormones.total_rate_pyy_release = float(ratePYYrelease.sum())
    hormones.total_rate_glp1_release = float(rateGLP1release.sum())
    hormones.total_rate_ghrelin_release = float(rateGhrelinrelease.sum())
    hormones.total_rate_motilin_release = float(rateMotilinrelease.sum())
    hormones.total_rate_secretin_release = float(rateSecretinrelease.sum())
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
        * math.exp(-(abs(st.ph[Segment.DUODENUM_2] - 6) / 6))
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


def sensations() -> None:
    """Compute fullness and hunger from the stomach tone and the intestinal content and absorption."""
    st.fullness = (10 / 15) * (stomach.tone - 2)
    st.hunger.total_caloric_absorption_rate = 10 - 25 * (st.calories.total_absorption_rate / st.calories.max_absorption_capacity)
    # Ileum 6 is deliberately excluded (as in the Pascal original); sum() adds left to right.
    st.hunger.intestinal_volume = sum(st.volume[Segment.DUODENUM_1:Segment.ILEUM_5 + 1])
    st.hunger.intestinal_volume = 10 - st.hunger.intestinal_volume / 1000
    st.hunger.level = 0.8 * (st.hunger.total_caloric_absorption_rate + st.hunger.intestinal_volume) / 2
    if st.mmc_status == MMCStatus.PHASE_III:
        st.hunger.level = 10


def intubation_in_antrum() -> None:
    """Optional infusion of water and glucose into the antrum (rates currently 0)."""
    if (st.clock.minutes > 250) and (st.clock.minutes < 251):
        st.liquid[Segment.ANTRUM, C.WATER] = st.liquid[Segment.ANTRUM, C.WATER] + 0 * st.clock.dt
        st.liquid[Segment.ANTRUM, C.GLUCOSE] = st.liquid[Segment.ANTRUM, C.GLUCOSE] + 0 * st.clock.dt


def intubation_in_duodenum() -> None:
    """Optional infusion of water and glucose into duodenum 1 (rates currently 0)."""
    if (st.clock.minutes > 250) and (st.clock.minutes < 300):
        st.liquid[Segment.DUODENUM_1, C.WATER] = st.liquid[Segment.DUODENUM_1, C.WATER] + 0 * st.clock.dt
        st.liquid[Segment.DUODENUM_1, C.GLUCOSE] = st.liquid[Segment.DUODENUM_1, C.GLUCOSE] + 0 * st.clock.dt


def intravenous() -> None:
    """Set the intravenous glucose infusion rates (glucose.g2_in and glucose.g2_in_mol_per_min)."""
    glucose.g2_in = 1
    glucose.g2_in_mol_per_min = 0.001 / 180.15588



