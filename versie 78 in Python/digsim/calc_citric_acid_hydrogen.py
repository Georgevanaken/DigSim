"""H+ binding of citric acid as a function of pH."""
# Based on Beltry (2019), considerations on the derivative of the titration curve of weak diacids.
def citric_acid_hydrogen(pH: float, citric_acid: float) -> float:
    """Return the H+ bound by citric acid (amount citric_acid) at the given pH."""
    correction_to_IEP = 0  # Source comment gives IEP as 000.
    result = 0.0  # NOTE(port): uninitialised in Pascal when trunc(pH) is outside 1..12.
    case = int(pH)  # First calculate on a mole-to-mole basis.
    if case == 1:
        result = 3 - (pH - int(pH)) * 0
    elif case == 2:
        result = 3 - (pH - int(pH)) * 0.5
    elif case == 3:
        result = 2.5 - (pH - int(pH)) * 0.5
    elif case == 4:
        result = 2 - (pH - int(pH)) * 0.6
    elif case == 5:
        result = 1.4 - (pH - int(pH)) * 0.6
    elif case == 6:
        result = 0.8 - (pH - int(pH)) * 0.5
    elif case == 7:
        result = 0.3 - (pH - int(pH)) * 0.3
    elif case in (8, 9, 10, 11, 12):
        result = 0 - (pH - int(pH)) * 0
    result = citric_acid * (-result - correction_to_IEP) / 192.124  # 192.124 is the molecular mass.
    return result



