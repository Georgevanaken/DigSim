"""Glycosylation of haemoglobin (HbA1c)."""
from . import state as st
from .state import glucose


def calc_hba1c() -> None:
    """Advance HbA1c (% of total Hb) by one time step: glycosylation by glucose minus decay."""
    glucose.hba1c_half_life = 3 * 30 * 24 * 60  # 3 months in minutes
    glucose.glycosylation_rate_constant = 1.23 * 1e-8
    glucose.hba1c = glucose.hba1c + (
        glucose.gluc * glucose.glycosylation_rate_constant
        - glucose.hba1c / (2 * glucose.hba1c_half_life)
    ) * st.clock.dt  # HbA1c in % of total Hb
