"""Serum concentrations and liver clearance."""
from . import state as st
from .enums import Component
from .state import glucose


def initialize_blood_serum_conc() -> None:
    """Reset the serum concentrations and set the decay rate of the pharmaceutical."""
    for j in range(1, 81):
        st.blood_serum_conc[j] = 0
        st.inverse_serum_half_time[j] = 0
    st.inverse_serum_half_time[Component.PHARMACEUTICAL] = 1 / 60  # inverse of the decaytime in minutes


def decay_in_blood() -> None:
    """Update the serum concentration (mg/l) of each component: absorption in, decay out."""
    # general
    for j in range(1, 81):
        st.blood_serum_conc[j] = st.blood_serum_conc[j] + (
            1000 * st.total_absorption_rate[j] / 5
            - st.blood_serum_conc[j] * st.inverse_serum_half_time[j]
        ) * st.clock.dt  # serum conc in mg/l


def fructose_to_glucose() -> None:
    """Let the liver remove fructose from the blood and convert part of it into glucose."""
    decay_fructose = 0.5 * glucose.fructose_decay_rate * st.clock.dt
    glucose.fructose_decay_rate = glucose.fructose_decay_rate - decay_fructose  # assumption is that fructose is removed from the blood by the liver in 2 minutes, fruc in mg
    glucose.gluc = glucose.gluc + 1 * decay_fructose  # assumption is that a part of the fructose is transformed by the liver into glucose; gluc in mg
