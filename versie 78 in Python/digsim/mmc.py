"""Migrating motor complex (MMC) state machine."""
from . import state as st
from .enums import MMCStatus
from .state import hormones


def update_mmc_status() -> None:
    """Advance the migrating motor complex (fed / phases I-IV) state machine, driven by motilin."""
    if st.clock.minutes > 0:
        if hormones.motilin > 40:
            Criterion_start_fasted = True
        else:
            Criterion_start_fasted = False
        if hormones.motilin < 40:
            criterion_start_fed = True
        else:
            criterion_start_fed = False

        if Criterion_start_fasted and st.mmc_status == MMCStatus.FED:
            st.mmc_status = MMCStatus.PHASE_I
            st.mmc_time[1] = 0
        if criterion_start_fed and st.mmc_status != MMCStatus.FED:
            st.mmc_status = MMCStatus.FED
            st.mmc_time[0] = 0

        # MMC_time is indexed by the numeric value of the phase.
        if st.mmc_status == MMCStatus.FED:
            st.mmc_time[0] = st.mmc_time[0] + st.clock.dt
        elif st.mmc_status == MMCStatus.PHASE_I:
            if st.mmc_time[1] > 50:
                st.mmc_status = MMCStatus.PHASE_II
                st.mmc_time[1] = 0
            else:
                st.mmc_time[1] = st.mmc_time[1] + st.clock.dt
        elif st.mmc_status == MMCStatus.PHASE_II:
            if st.mmc_time[2] > 20:
                st.mmc_status = MMCStatus.PHASE_III
                st.mmc_time[2] = 0
            else:
                st.mmc_time[2] = st.mmc_time[2] + st.clock.dt
        elif st.mmc_status == MMCStatus.PHASE_III:
            if st.mmc_time[3] > 7:
                st.mmc_status = MMCStatus.PHASE_IV
                st.mmc_time[3] = 0
            else:
                st.mmc_time[3] = st.mmc_time[3] + st.clock.dt
        elif st.mmc_status == MMCStatus.PHASE_IV:
            if st.mmc_time[4] > 10:
                st.mmc_status = MMCStatus.PHASE_I
                st.mmc_time[4] = 0
            else:
                st.mmc_time[4] = st.mmc_time[4] + st.clock.dt
