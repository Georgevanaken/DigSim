"""Segment volumes, densities and related properties derived from the liquid and solid contents."""
import numpy as np

from . import state as st
from .enums import Segment
from .state import diet, stomach


def volumina() -> None:
    """Compute the volume, density and related properties of each segment from liquid and solid contents."""

    # The inner loops over components/compounds/solid classes are vectorised (numpy slices).
    comps = slice(1, st.comp + 1)
    comps_n = slice(1, diet.number_of_compounds + 1)
    dens = st.density[comps]

    dishes = slice(1, 6)
    # Row sums over the component axis; each row is added as before: (1e-10 + liquid part) + compound part.
    dish_volume = (1e-10 + np.sum(st.liquid_dish[dishes, comps] / dens, axis=1)
                   + np.sum(st.compound_dish[dishes, comps_n], axis=1))
    present = dish_volume > 1e-20
    vol = np.where(present, dish_volume, 1.0)[:, None]
    st.volume_dish[dishes] = np.where(present, dish_volume, 1e-20)
    st.concentration_dish[dishes, comps] = np.where(present[:, None], st.liquid_dish[dishes, comps] / vol, 0)
    st.concentration_compound_dish[dishes, comps_n] = np.where(
        present[:, None], st.compound_dish[dishes, comps_n] / vol, 0)

    segs = slice(Segment.MOUTH, Segment.COLON + 1)
    n_seg = Segment.COLON - Segment.MOUTH + 1
    # Row-wise sums over contiguous copies use the same pairwise order as the former per-segment 1-D sums.
    total = 1e-5 + np.sum(st.liquid[segs, comps] / dens, axis=1)
    if st.s_max >= 1:
        solids = st.solids[1:st.s_max + 1, segs, comps] / dens
        total += np.sum(solids.transpose(1, 0, 2).reshape(n_seg, -1), axis=1)
    total += np.sum(st.compound[segs, comps_n] / diet.compound_density[comps_n], axis=1)
    st.volume[segs] = total

    present = total > 1e-5
    vol = total[:, None]
    st.concentration_lumen[segs, comps] = np.where(
        present[:, None], st.liquid[segs, comps] / vol, 1e-10)
    st.concentration_compound[segs, comps_n] = np.where(
        present[:, None], st.compound[segs, comps_n] / vol, 1e-10)
    # NOTE(port): Pascal writes Volume[j] (j is a stale loop counter) when the volume is too small; intended
    # is Volume[i]. This is a no-op because the volume starts at 1e-5 and only grows.
    st.volume[segs] = np.where(present, total, 1e-5)

    st.in_enterocytes = 0
    stomach.v_fundus = 0
    stomach.v_fundus += np.sum(stomach.fundus[comps] / dens)
    stomach.v_fundus += np.sum(st.compound_fundus[comps_n])

    if stomach.v_fundus > 1e-5:
        stomach.concentration_fundus[comps] = stomach.fundus[comps] / stomach.v_fundus
        stomach.concentration_compound_fundus[comps_n] = stomach.compound_fundus[comps_n] / stomach.v_fundus
    else:
        stomach.v_fundus = 1e-5
        stomach.concentration_fundus[comps] = 0
        stomach.concentration_compound_fundus[comps_n] = 0

    st.volume_total = float(np.cumsum(np.concatenate(([0.0], st.volume[Segment.DUODENUM_1:Segment.ILEUM_6 + 1])))[-1])
