"""Run the first step function by function and print key values after each one."""
import sys

import digsim.calc_protein_hydrogen as ph
import digsim.compartment_properties as cp
import digsim.gastric_acid_secretion as ga
import digsim.gastric_motility as gm
import digsim.glucose_homeostasis2 as gh
import digsim.hb_glycosylation as hb
import digsim.liver_clearance as lc
import digsim.main_part1 as p1
import digsim.main_part2 as p2
import digsim.main_part3 as p3
import digsim.mmc
from digsim import state as st

sys.argv = ["x", "6", "0"]
# Reuse the set-up part of trace_init.py (everything before its stepping loop).
exec(open("tools/trace_init.py").read().split("t0 = time.time()")[0])

print("start liquid[3,1]=", st.liquid[3, 1], "vol", st.volume[2:5], "density[1..4]", st.density[1:5])
st.clock.minutes += st.clock.dt
seq = [("update_mmc_status", digsim.mmc.update_mmc_status), ("volumina", cp.volumina), ("calories", p1.calories),
       ("decomposition", p2.decomposition), ("digestion", p2.digestion), ("osm_l", p1.osmotic_values_lumen),
       ("osm_b", p1.osmotic_values_brush), ("absorption", p3.absorption_in_small_intestine), ("mucus", p3.mucus_layer),
       ("intub_a", p3.intubation_in_antrum), ("intub_d", p3.intubation_in_duodenum), ("intrav", p3.intravenous),
       ("fructose", lc.fructose_to_glucose), ("insulin", gh.insulin_dynamics), ("hba1c", hb.calc_hba1c),
       ("hormones", p3.gut_hormones), ("decay", lc.decay_in_blood), ("pHs", ph.calc_all_ph),
       ("saliva", p1.saliva_secretion), ("gastric_sec", ga.gastric_secretion), ("pancreabil", p1.pancreatobiliary_secretion),
       ("intestinal", p1.intestinal_secretion), ("caloric_dens", p1.caloric_densities_segments),
       ("solids_in_stomach", gm.solids_in_stomach), ("colon", p3.colon), ("volumina2", cp.volumina), ("calc_outflows", p1.calc_outflows)]
for name, f in seq:
    f()
    print("%-18s liquid[3,1]=%-12.5g liquid[4,1]=%-10.5g vol3=%-10.5g up3=%-10.5g up2=%-10.5g pH3=%.3f" % (
        name, st.liquid[3, 1], st.liquid[4, 1], st.volume[3], st.outflow[3], st.outflow[2], st.ph[3]))
