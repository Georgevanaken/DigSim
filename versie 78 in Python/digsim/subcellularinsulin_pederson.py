"""Sub-cellular beta-cell model of insulin secretion (Pedersen)."""
import math

from . import state as st
from .state import glucose

m_secretion = 0.62
fplus = 6.2
Pplus = 0.021
Pmin = 0.1695
k = 0.47
r = 3 * 0.0023
nM = 12
nH = 3.3
c = 1.25
M0 = 3 * 0.020
KmM = 250
G50 = 150
tau = 0.5
imax = 500
f1min = 0.01


def sigmoid(BASIS: float, MAX: float, POWER: float, x: float, xhalf: float) -> float:
    """Return a Hill-type sigmoid BASIS + MAX * x^POWER / (...) as written in the Pascal original."""
    result = BASIS + MAX * math.exp(POWER * math.log(x)) / math.exp(
        POWER * math.log(x) + math.exp(POWER * math.log(xhalf))
    )
    return result


def beta_cell_distribution() -> None:
    """Initialise the distribution (phi, dphi) of the beta-cell insulin granule pools."""
    for i in range(1, imax + 1):
        st.beta_cell.phi[i] = math.exp(nH * math.log(i)) / (
            math.exp(nH * math.log(imax / 2)) + math.exp(nH * math.log(i))
        )
    st.beta_cell.dphi[1] = 0
    for i in range(2, imax + 1):
        st.beta_cell.dphi[i] = st.beta_cell.phi[i] - st.beta_cell.phi[i - 1]
    for i in range(1, imax + 1):
        st.beta_cell.intermediate[i] = 1 * st.beta_cell.dphi[i]
        st.beta_cell.rrp[i] = 1.62 * st.beta_cell.dphi[i]
        st.beta_cell.fused[i] = 0


def calc_f1_pederson() -> None:
    """Advance the beta-cell granule model (Pedersen) and compute the insulin secretion rate st.beta_cell.f1."""
    total = 0
    if st.clock.minutes == st.clock.start_time:
        glucose.level = 40
        glucose.mobilization_level = M0 + c * math.exp(nM * math.log(glucose.level)) / (
            math.exp(nM * math.log(KmM)) + math.exp(nM * math.log(glucose.level))
        )
    else:
        mobilization = M0 + c * math.exp(nM * math.log(glucose.level)) / (
            math.exp(nM * math.log(KmM)) + math.exp(nM * math.log(glucose.level))
        )
        M_increment = (mobilization - glucose.mobilization_level) / tau
        glucose.mobilization_level = glucose.mobilization_level + M_increment * st.clock.dt

    for i in range(1, imax + 1):
        st.beta_cell.intermediate_increment[i] = (
            glucose.mobilization_level * st.beta_cell.dphi[i]
            - (r + Pplus) * st.beta_cell.intermediate[i]
            + Pmin * st.beta_cell.rrp[i]
        )
        st.beta_cell.rrp_increment[i] = (
            Pplus * st.beta_cell.intermediate[i] - Pmin * st.beta_cell.rrp[i]
        )
        # NOTE(port): Fused_increment remains unchanged when this condition is false, matching Pascal.
        if 2 * i < imax * glucose.level / G50:
            st.beta_cell.intermediate_increment[i] = (
                st.beta_cell.intermediate_increment[i] + k * st.beta_cell.fused[i]
            )
            st.beta_cell.rrp_increment[i] = st.beta_cell.rrp_increment[i] - fplus * st.beta_cell.rrp[i]
            st.beta_cell.fused_increment[i] = (
                fplus * st.beta_cell.rrp[i] - (m_secretion + k) * st.beta_cell.fused[i]
            )
            st.beta_cell.released[i] = m_secretion * st.beta_cell.fused[i]
            total = total + st.beta_cell.released[i]

    st.beta_cell.f1 = total / 0.0387
    st.beta_cell.f1 = st.beta_cell.f1
    if st.beta_cell.f1 < f1min:
        st.beta_cell.f1 = f1min

    for i in range(1, imax + 1):
        st.beta_cell.intermediate[i] = (
            st.beta_cell.intermediate[i] + st.beta_cell.intermediate_increment[i] * st.clock.dt
        )
        st.beta_cell.rrp[i] = st.beta_cell.rrp[i] + st.beta_cell.rrp_increment[i] * st.clock.dt
        st.beta_cell.fused[i] = st.beta_cell.fused[i] + st.beta_cell.fused_increment[i] * st.clock.dt


