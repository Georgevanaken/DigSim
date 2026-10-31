"""Step the initialisation phase and report where values first become suspicious."""
import sys
import time

import numpy as np

from digsim import state as st
from digsim.compartment_properties import volumina
from digsim.diet_from_file import read_diet
from digsim.digsim_v78 import make_step
from digsim.main_part1 import initialise
from digsim.model_parameters import set_model_parameters

n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
nsteps = int(sys.argv[2]) if len(sys.argv) > 2 else 200
np.seterr(all="ignore")
st.comp = 110
set_model_parameters()
read_diet()
initialise()
st.clock.substeps = n
st.clock.dt = 1 / n
volumina()
st.clock.minutes = -500
st.dish_number = 1
st.swallow_speed = 50
st.ph[1] = 7
st.diet.liquid_viscosity_dish[1] = 1
st.liquid_dish[1, 1] = 250
st.liquid_dish[1, 2] = 1
st.liquid_dish[1, 3] = 50
st.liquid_dish[1, 4] = 10
st.eating_speed_dish[1] = 50

t0 = time.time()
for k in range(nsteps):
    st.clock.minutes += st.clock.dt
    make_step()
    bad = [(name, np.argwhere(~np.isfinite(getattr(st, name))).tolist()[:3])
           for name in ("liquid", "volume", "outflow", "ph", "compound")
           if not np.all(np.isfinite(getattr(st, name)))]
    big = float(np.max(np.abs(st.liquid)))
    if bad or k % 20 == 0:
        print(k, round(st.clock.minutes, 3), "max|liquid|=%.3g" % big,
              "vol[3]=%.4g up[3]=%.4g" % (st.volume[3], st.outflow[3]), bad)
    if bad:
        break
print("time per step: %.4f s" % ((time.time() - t0) / (k + 1)))
