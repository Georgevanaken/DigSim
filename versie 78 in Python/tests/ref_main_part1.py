"""Original scalar loops of caloric_densities_segments and calories (reference for the vectorised versions)."""
from digsim import state as st
from digsim.enums import Segment
from digsim.state import diet


def caloric_densities_segments():
    for i in range(Segment.DUODENUM_1, Segment.COLON + 1):
        calories_segment = 0
        for j in range(1, st.comp + 1):
            calories_segment = calories_segment + st.liquid[i, j] * st.calories.value[j]
        for n in range(1, diet.number_of_compounds + 1):
            calories_segment = calories_segment + st.compound[i, n] * st.calories.value[n]
        if st.volume[i] > 1e-10:
            st.calories.density_local[i] = calories_segment / st.volume[i]
        else:
            st.calories.density_local[i] = 0
        st.calories.brush = 0
        for j in range(1, st.comp + 1):
            st.calories.brush = st.calories.brush + st.concentration_brush[i, j] * st.calories.value[j]


def calories():
    st.calories.proximal = 0
    st.calories.distal = 0
    st.calories.brush = 0
    st.calories.colon = 0
    proximal = range(Segment.DUODENUM_1, Segment.JEJUNUM_2 + 1)
    distal = range(Segment.JEJUNUM_3, Segment.ILEUM_6 + 1)
    for j in range(1, st.comp + 1):
        for i in proximal:
            st.calories.proximal = st.calories.proximal + st.calories.value[j] * st.liquid[i, j]
        for i in distal:
            st.calories.distal = st.calories.distal + st.calories.value[j] * st.liquid[i, j]
        st.calories.colon = st.calories.colon + st.calories.value[j] * st.liquid[Segment.COLON, j]
    for n in range(1, diet.number_of_compounds + 1):
        for i in proximal:
            st.calories.proximal = st.calories.proximal + st.compound[i, n] * st.calories.compound[n]
        for i in distal:
            st.calories.distal = st.calories.distal + st.compound[i, n] * st.calories.compound[n]
        st.calories.colon = st.calories.colon + st.compound[Segment.COLON, n] * st.calories.compound[n]
