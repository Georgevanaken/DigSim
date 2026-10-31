"""Original scalar implementation of the transport part of make_step (reference for equivalence tests)."""
from digsim import state as st


def transport():
    for st.j in range(1, st.comp + 1):
        st.liquid_next[2, st.j] = 0
        for st.dish_number in range(1, 6):
            st.liquid_dish_next[st.dish_number, st.j] = st.liquid_dish[st.dish_number, st.j] - st.outflow_dish[st.dish_number] * st.concentration_dish[st.dish_number, st.j] * st.clock.dt
            st.liquid_next[2, st.j] = st.liquid_next[2, st.j] + st.outflow_dish[st.dish_number] * st.concentration_dish[st.dish_number, st.j] * st.clock.dt
        st.liquid_next[2, st.j] = st.liquid[2, st.j] + st.liquid_next[2, st.j] - st.outflow[2] * st.concentration_lumen[2, st.j] * st.clock.dt

        st.liquid_next[3, st.j] = st.liquid[3, st.j] + st.outflow[2] * st.concentration_lumen[2, st.j] * st.clock.dt
        st.liquid_next[4, st.j] = st.liquid[4, st.j]
        st.liquid_next[5, st.j] = st.liquid[5, st.j] - st.outflow[5] * st.concentration_lumen[5, st.j] * st.clock.dt

        for st.i in range(6, 16):
            st.liquid_next[st.i, st.j] = st.liquid[st.i, st.j] + st.outflow[st.i - 1] * st.concentration_lumen[st.i - 1, st.j] * st.clock.dt
            st.liquid_next[st.i, st.j] = st.liquid_next[st.i, st.j] - st.outflow[st.i] * st.concentration_lumen[st.i, st.j] * st.clock.dt
        st.liquid_next[16, st.j] = st.liquid[16, st.j] + st.outflow[15] * st.concentration_lumen[15, st.j] * st.clock.dt

    for compoundnumber in range(1, st.diet.number_of_compounds + 1):
        st.compound_next[2, compoundnumber] = 0
        for st.dish_number in range(1, 6):
            st.compound_dish_next[st.dish_number, compoundnumber] = st.compound_dish[st.dish_number, compoundnumber] - st.outflow_dish[st.dish_number] * st.concentration_compound_dish[st.dish_number, compoundnumber] * st.clock.dt
            st.compound_next[2, compoundnumber] = st.compound_next[2, compoundnumber] + st.outflow_dish[st.dish_number] * st.concentration_compound_dish[st.dish_number, compoundnumber] * st.clock.dt
        if st.solid[compoundnumber] > 10:
            st.swallow_solid = 0
        else:
            st.swallow_solid = 1
        st.compound_next[2, compoundnumber] = st.compound[2, compoundnumber] + st.compound_next[2, compoundnumber] - st.outflow[2] * st.swallow_solid * st.concentration_compound[2, compoundnumber] * st.clock.dt

        st.compound_next[3, compoundnumber] = st.compound[3, compoundnumber] + st.outflow[2] * st.swallow_solid * st.concentration_compound[2, compoundnumber] * st.clock.dt
        st.compound_next[4, compoundnumber] = st.compound[4, compoundnumber]
        st.compound_next[5, compoundnumber] = st.compound[5, compoundnumber] - st.outflow[5] * st.concentration_compound[5, compoundnumber] * st.clock.dt

        for st.i in range(6, 16):
            st.compound_next[st.i, compoundnumber] = st.compound[st.i, compoundnumber] + (st.outflow[st.i - 1] * st.concentration_compound[st.i - 1, compoundnumber] - st.outflow[st.i] * st.concentration_compound[st.i, compoundnumber]) * st.clock.dt
        st.compound_next[16, compoundnumber] = st.compound[16, compoundnumber] + st.outflow[15] * st.concentration_compound[15, compoundnumber] * st.clock.dt

    for st.dish_number in range(1, 6):
        for compoundnumber in range(1, st.diet.number_of_compounds + 1):
            st.compound_dish[st.dish_number, compoundnumber] = st.compound_dish_next[st.dish_number, compoundnumber]
        for st.j in range(1, st.comp + 1):
            st.liquid_dish[st.dish_number, st.j] = st.liquid_dish_next[st.dish_number, st.j]

    for st.i in range(2, 17):
        for compoundnumber in range(1, st.diet.number_of_compounds + 1):
            st.compound[st.i, compoundnumber] = st.compound_next[st.i, compoundnumber]
        for st.j in range(1, st.comp + 1):
            st.liquid[st.i, st.j] = st.liquid_next[st.i, st.j]
