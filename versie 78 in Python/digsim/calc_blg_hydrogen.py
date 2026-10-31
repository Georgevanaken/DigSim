"""H+ binding of beta-lactoglobulin as a function of pH."""
# Based on Tanford and Nozaki, Physico-chemical comparison of b-lactoglobulin A and B, J. Biol. Chem. 1959, 234(11), 2874-2877.
def blg_hydrogen(pH: float, BLG: float) -> float:
    """Return the H+ bound by beta-lactoglobulin (amount BLG) at the given pH, relative to its IEP."""
    correction_to_IEP = 0  # Source comment gives IEP as 000.
    result = 0.0  # NOTE(port): uninitialised in Pascal when trunc(pH) is outside 1..12.
    case = int(pH)
    if case == 1:
        result = 0.002173913 - (pH - int(pH)) * 9.23913e-05
    elif case == 2:
        result = 0.002081522 - (pH - int(pH)) * 0.000383152
    elif case == 3:
        result = 0.00169837 - (pH - int(pH)) * 0.001263587
    elif case == 4:
        result = 0.000434783 - (pH - int(pH)) * 0.000543478
    elif case == 5:
        result = -0.000108696 - (pH - int(pH)) * 0.000271739
    elif case == 6:
        result = -0.000380435 - (pH - int(pH)) * 0.00026087
    elif case == 7:
        result = -0.000641304 - (pH - int(pH)) * 0.000228261
    elif case == 8:
        result = -0.000869565 - (pH - int(pH)) * 0.000108696
    elif case == 9:
        result = -0.000978261 - (pH - int(pH)) * 5.43478e-05
    elif case == 10 or case == 11 or case == 12:
        result = -0.1032609 - (pH - int(pH)) * 0
    result = BLG * (result - correction_to_IEP)  # Hydrogen ions for deviation from the IEP.
    return result


