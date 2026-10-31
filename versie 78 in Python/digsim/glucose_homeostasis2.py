"""Glucose-insulin-glucagon homeostasis (Tolic et al., 2000)."""
import math

from . import state as st
from .state import glucose, hormones
from .subcellularinsulin_pederson import beta_cell_distribution

beta_cells_percentage = 1
alpha_cells_percentage = 1
Vp = 3
Vi = 11
Vg = 10
Vs = 1
E = 0.2
ti = 100
td1 = 3
td2 = 3
td3 = 4
Rm = 210
a1 = 300
C1 = 2000
Ub = 72
C2 = 144
C3 = 1000
U0 = 40
Um = 940
beta = 1.77
C4 = 80
Rg = 180
alpha = 0.29
C5 = 26
gamma = 5.0
C6 = 2.0
Sb = 20
Sm = 140
ddelta = -2.4
C7 = 2.0
insulin_sensitivity = 1.40
glucose_sensitivity = 1.20
ratio = 0.1
k1 = 0.1
k2 = 0.1
k3 = 0


def f2(Gluc: float) -> float:
    """Tolic f2: insulin-independent glucose utilisation."""
    return Ub * (1 - math.exp(-glucose_sensitivity * Gluc / (C2 * Vg)))


def f3(Gluc: float) -> float:
    """Tolic f3: glucose-dependent factor of the insulin-dependent utilisation (multiplied with f4)."""
    return glucose_sensitivity * Gluc / (C3 * Vg)


def f4(Ii: float) -> float:
    """Tolic f4: insulin-dependent factor of the glucose utilisation (multiplied with f3), from interstitial insulin Ii."""
    return U0 + (Um - U0) / (
        1 + math.exp(-beta * math.log((insulin_sensitivity * Ii / C4) * (1 / Vi + 1 / (E * ti))))
    )


def f5(x3: float) -> float:
    """Tolic f5: insulin-dependent suppression factor of endogenous glucose production (multiplied with f6)."""
    return Rg / (1 + math.exp(alpha * (x3 * insulin_sensitivity / Vp - C5)))


def f6(Gluc: float, Gin: float) -> float:
    """Tolic f6: glucose-dependent factor of endogenous glucose production (multiplied with f5)."""
    return 1 / (
        1
        + math.exp(
            gamma
            * (glucose_sensitivity * (Gluc / Vg + ratio * (Gin / 0.57 - Gluc / Vg)) / C3 - C6)
        )
    )


def f7(Gluc: float, Gin: float) -> float:
    """Additional glucose-dependent term (parameters Sb, Sm, ddelta, C7) that is added to the utilisation."""
    return Sb + (Sm - Sb) * (
        1
        / (
            1
            + math.exp(
                ddelta
                * (glucose_sensitivity * (Gluc / Vg + ratio * (Gin / 0.57 - Gluc / Vg)) / C3 - C7)
            )
        )
        - 1 / (1 + math.exp(-ddelta * C7))
    )


def initialise_glucose_dynamics() -> None:
    """Set the initial state (insulin, glucagon, glycogen, beta-cell distribution) of the glucose model."""
    beta_cell_distribution()
    glucose.serum_glucagon = 10

    glucose.x1 = 10
    glucose.x2 = 10
    glucose.x3 = 10
    glucose.gluc0 = 94 * 10 * Vg
    glucose.gluc1 = 0
    glucose.gluc2 = 0
    glucose.gluc = glucose.gluc0 + glucose.gluc1 + glucose.gluc2
    glucose.ip = 10 * Vp
    glucose.ii = 10 * Vi
    glucose.gluci = 10
    glucose.gluciskin = 1


def insulin_dynamics() -> None:
    """Advance the insulin and glucose compartments one time step (Tolic et al., 2000)."""
    # Tolic et al., J. Theor. Biol. (2000)207, 361-375
    Ip_increment = 0.0
    Ii_increment = 0.0
    x1_increment = 0.0
    x2_increment = 0.0
    x3_increment = 0.0
    f1_Tolic = 0.0
    incretin_factor = 0.0
    gluci_increment = 0.0
    Gluciskin_increment = 0.0
    Gi = 0.0
    glucose.tpf = 4
    glucose.incretin = 5 * (0.02 * hormones.real_gip + 0.05 * hormones.real_glp1)
    incretin_factor = 1 + glucose.incretin
    glucose.tp = glucose.tpf * (1 + 0.5 * glucose.incretin / (1 + glucose.incretin))

    glucose.g_in = glucose.g0_in + glucose.g1_in + glucose.g2_in
    glucose.gluc = glucose.gluc0 + glucose.gluc1 + glucose.gluc2

    f1_Tolic = incretin_factor * Rm / (
        1
        + math.exp(
            (C1 - glucose_sensitivity * (glucose.gluc / Vg + ratio * (glucose.g_in / 0.57))) / a1
        )
    )
    st.f1 = f1_Tolic

    Ip_increment = st.f1 - E * (glucose.ip / Vp - glucose.ii / Vi) - glucose.ip / glucose.tp
    Ii_increment = E * (glucose.ip / Vp - glucose.ii / Vi) - glucose.ii / ti
    gluci_increment = k1 * Vi * (glucose.gluc / Vg - glucose.gluci / Vi) - glucose.utilisation
    Gluciskin_increment = k2 * Vs * (glucose.gluc / Vg - glucose.gluciskin / Vs) - k3 * glucose.gluciskin
    Gi = glucose.gluc

    glucose.endogenous_production = f5(glucose.x3) * f6(Gi, glucose.g_in)
    glucose.utilisation = f3(Gi) * f4(glucose.ii) + f2(Gi) + f7(Gi, glucose.g_in)

    glucose.g0_increment = glucose.g0_in - (glucose.gluc0 / glucose.gluc) * glucose.utilisation + glucose.endogenous_production
    glucose.g1_increment = glucose.g1_in - (glucose.gluc1 / glucose.gluc) * glucose.utilisation
    glucose.g2_increment = glucose.g2_in - (glucose.gluc2 / glucose.gluc) * glucose.utilisation

    x1_increment = 3 * (glucose.ip - glucose.x1) / td1
    x2_increment = 3 * (glucose.x1 - glucose.x2) / td2
    x3_increment = 3 * (glucose.x2 - glucose.x3) / td3

    glucose.ip = glucose.ip + Ip_increment * st.clock.dt
    glucose.ii = glucose.ii + Ii_increment * st.clock.dt

    glucose.gluc0 = glucose.gluc0 + glucose.g0_increment * st.clock.dt
    glucose.gluc1 = glucose.gluc1 + glucose.g1_increment * st.clock.dt
    glucose.gluc2 = glucose.gluc2 + glucose.g2_increment * st.clock.dt
    glucose.gluc = glucose.gluc0 + glucose.gluc1 + glucose.gluc2

    glucose.gluci = glucose.gluci + gluci_increment * st.clock.dt
    glucose.gluciskin = glucose.gluciskin + Gluciskin_increment * st.clock.dt

    glucose.x1 = glucose.x1 + x1_increment * st.clock.dt
    glucose.x2 = glucose.x2 + x2_increment * st.clock.dt
    glucose.x3 = glucose.x3 + x3_increment * st.clock.dt

    glucose.serum_insulin_mu_l = glucose.ip / Vp
    glucose.serum_insulin = glucose.serum_insulin_mu_l * 6.67
    glucagon_dynamics()
    serum_glucose_dynamics()

    glucose.serum0_mg_dl = 0.1 * glucose.gluc0 / Vg
    glucose.serum1_mg_dl = 0.1 * glucose.gluc1 / Vg
    glucose.serum2_mg_dl = 0.1 * glucose.gluc2 / Vg
    glucose.serum_mg_dl = (
        glucose.serum0_mg_dl + glucose.serum1_mg_dl + glucose.serum2_mg_dl
    )
    glucose.intercellular_muscle_adipose_mg_dl = 0.1 * glucose.gluci / Vi
    glucose.intercellular_skin_mg_dl = 0.1 * glucose.gluciskin / Vs

    glucose.serum0 = glucose.serum0_mg_dl * 10 / 180.15588
    glucose.serum1 = glucose.serum1_mg_dl * 10 / 180.15588
    glucose.serum2 = glucose.serum2_mg_dl * 10 / 180.15588
    glucose.serum = glucose.serum0 + glucose.serum1 + glucose.serum2
    glucose.intercellular_muscle_adipose = (
        glucose.intercellular_muscle_adipose_mg_dl * 10 / 180.15588
    )
    glucose.intercellular_skin = glucose.intercellular_skin_mg_dl * 10 / 180.15588

    st.r_c13_c12 = 1000000 * glucose.serum1 / glucose.serum0


def glucagon_dynamics() -> None:
    """Advance serum glucagon one time step: production depends on the serum glucose deviation."""
    difference = 0.0
    glucagon_production = 0.0
    glucagon_decay = 0.0
    serum_glucose_criterion = 4.5
    difference = glucose.serum - serum_glucose_criterion
    if difference < 0:
        glucagon_production = 1 * math.exp(-difference) * alpha_cells_percentage
    else:
        glucagon_production = 0
    glucagon_decay = glucose.serum_glucagon / 5
    glucose.serum_glucagon = glucose.serum_glucagon + (glucagon_production - glucagon_decay) * st.clock.dt


def serum_glucose_dynamics() -> None:
    """Advance the serum glucose pools and the glycogen stores one time step."""
    if glucose.glycogen_store < 100:
        glucose.glycogen_store0 = glucose.glycogen_store0 + (
            0.01 * glucose.serum0 * glucose.serum_insulin_mu_l / 50
            - 0.01 * (glucose.glycogen_store0 / glucose.glycogen_store) * glucose.serum_glucagon
        ) * st.clock.dt
        glucose.glycogen_store1 = glucose.glycogen_store1 + (
            0.01 * glucose.serum1 * glucose.serum_insulin_mu_l / 50
            - 0.01 * (glucose.glycogen_store1 / glucose.glycogen_store) * glucose.serum_glucagon
        ) * st.clock.dt
        glucose.glycogen_store2 = glucose.glycogen_store2 + (
            0.01 * glucose.serum2 * glucose.serum_insulin_mu_l / 50
            - 0.01 * (glucose.glycogen_store2 / glucose.glycogen_store) * glucose.serum_glucagon
        ) * st.clock.dt
        glucose.glycogen_store = glucose.glycogen_store0 + glucose.glycogen_store1 + glucose.glycogen_store2
    glucose.serum0 = glucose.serum0 + 0.01 * (glucose.glycogen_store0 / glucose.glycogen_store) * glucose.serum_glucagon * st.clock.dt
    glucose.serum1 = glucose.serum1 + 0.01 * (glucose.glycogen_store1 / glucose.glycogen_store) * glucose.serum_glucagon * st.clock.dt
    glucose.serum2 = glucose.serum2 + 0.01 * (glucose.glycogen_store2 / glucose.glycogen_store) * glucose.serum_glucagon * st.clock.dt
    glucose.serum = glucose.serum0 + glucose.serum1 + glucose.serum2
    glucose.gluc = glucose.serum * 10 * Vg



