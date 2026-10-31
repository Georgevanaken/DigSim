from digsim import state as st
from digsim.state import stomach

# Implementation-level Pascal variable, distinct from st.dish_number.
dishnumber = 0


def volumina():
    global dishnumber

    for dishnumber in range(1, 6):
        st.volume_dish[dishnumber] = 1e-10
        for st.j in range(1, st.comp + 1):
            st.volume_dish[dishnumber] += st.liquid_dish[dishnumber, st.j] / st.density[st.j]
        for st.compound_number in range(1, st.diet.number_of_compounds + 1):
            st.volume_dish[dishnumber] += st.compound_dish[dishnumber, st.compound_number] / 1  # Compound density not yet introduced (currently set to 1 g/ml).

        if st.volume_dish[dishnumber] > 1e-20:
            for st.j in range(1, st.comp + 1):
                st.concentration_dish[dishnumber, st.j] = st.liquid_dish[dishnumber, st.j] / st.volume_dish[dishnumber]
            for st.compound_number in range(1, st.diet.number_of_compounds + 1):
                st.concentration_compound_dish[dishnumber, st.compound_number] = st.compound_dish[dishnumber, st.compound_number] / st.volume_dish[dishnumber]
        else:
            st.volume_dish[dishnumber] = 1e-20
            for st.j in range(1, st.comp + 1):
                st.concentration_dish[dishnumber, st.j] = 0
            for st.compound_number in range(1, st.diet.number_of_compounds + 1):
                st.concentration_compound_dish[dishnumber, st.compound_number] = 0

    for st.i in range(2, 17):
        st.volume[st.i] = 1e-5
        for st.j in range(1, st.comp + 1):
            st.volume[st.i] += st.liquid[st.i, st.j] / st.density[st.j]
            for st.s in range(1, st.s_max + 1):
                st.volume[st.i] += st.solids[st.s, st.i, st.j] / st.density[st.j]
        for st.compound_number in range(1, st.diet.number_of_compounds + 1):
            st.volume[st.i] += st.compound[st.i, st.compound_number] / st.diet.compound_density[st.compound_number]

        if st.volume[st.i] > 1e-5:
            for st.j in range(1, st.comp + 1):
                st.concentration_lumen[st.i, st.j] = st.liquid[st.i, st.j] / st.volume[st.i]
            for st.compound_number in range(1, st.diet.number_of_compounds + 1):
                st.concentration_compound[st.i, st.compound_number] = st.compound[st.i, st.compound_number] / st.volume[st.i]
        else:
            # NOTE(port): Pascal writes Volume[j] (j is a stale loop counter, up to 110, i.e. an
            # out-of-bounds write on array[1..16] without range checking). Intended is Volume[i];
            # this is a no-op anyway because volume[i] starts at 1e-5 and only grows.
            st.volume[st.i] = 1e-5
            for st.j in range(1, st.comp + 1):
                st.concentration_lumen[st.i, st.j] = 1e-10
            for st.compound_number in range(1, st.diet.number_of_compounds + 1):
                st.concentration_compound[st.i, st.compound_number] = 1e-10

    st.in_enterocytes = 0
    st.volume_total = 0
    stomach.v_fundus = 0
    for st.j in range(1, st.comp + 1):
        stomach.v_fundus += stomach.fundus[st.j] / st.density[st.j]
    for st.compound_number in range(1, st.diet.number_of_compounds + 1):
        stomach.v_fundus += st.compound_fundus[st.compound_number] / 1

    if stomach.v_fundus > 1e-5:
        for st.j in range(1, st.comp + 1):
            stomach.concentration_fundus[st.j] = stomach.fundus[st.j] / stomach.v_fundus
        for st.compound_number in range(1, st.diet.number_of_compounds + 1):
            stomach.concentration_compound_fundus[st.compound_number] = stomach.compound_fundus[st.compound_number] / stomach.v_fundus
    else:
        stomach.v_fundus = 1e-5
        for st.j in range(1, st.comp + 1):
            stomach.concentration_fundus[st.j] = 0
        for st.compound_number in range(1, st.diet.number_of_compounds + 1):
            stomach.concentration_compound_fundus[st.compound_number] = 0

    for st.i in range(5, 16):
        st.volume_total += st.volume[st.i]
