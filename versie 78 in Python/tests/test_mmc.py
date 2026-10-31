from digsim import state as st
from digsim.enums import MMCStatus
from digsim.mmc import update_mmc_status


def test_mmc_status_values_match_pascal_integers():
    assert [int(s) for s in MMCStatus] == [0, 1, 2, 3, 4]
    assert f"{MMCStatus.PHASE_III:5d}" == "    3"


def test_mmc_cycle_runs_through_all_phases():
    saved = (st.clock.minutes, st.hormones.motilin, st.mmc_status, st.clock.dt, st.mmc_time.copy())
    try:
        st.clock.minutes, st.hormones.motilin, st.clock.dt = 1.0, 100.0, 1.0  # fasted: motilin above 40
        st.mmc_status = MMCStatus.FED
        st.mmc_time[:] = 0
        seen = []
        for _ in range(200):
            update_mmc_status()
            if not seen or seen[-1] != st.mmc_status:
                seen.append(st.mmc_status)
        assert seen[:5] == [MMCStatus.PHASE_I, MMCStatus.PHASE_II, MMCStatus.PHASE_III,
                            MMCStatus.PHASE_IV, MMCStatus.PHASE_I]
        st.hormones.motilin = 10.0  # fed again: back to FED at once
        update_mmc_status()
        assert st.mmc_status == MMCStatus.FED
    finally:
        st.clock.minutes, st.hormones.motilin, st.mmc_status, st.clock.dt = saved[:4]
        st.mmc_time[:] = saved[4]
