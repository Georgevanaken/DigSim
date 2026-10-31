"""H+ binding of bovine serum albumin as a function of pH."""
# Based on Curvale (2009), Buffer capacity of bovine serum albumin.
def bsa_hydrogen(pH: float, BSA: float) -> float:
    """Return the H+ bound by bovine serum albumin (amount BSA) at the given pH, relative to its IEP."""
    correction_to_IEP = 0.0267  # IEP is 4.7.
    result = 0.0  # NOTE(port): uninitialised in Pascal when trunc(pH) is outside 1..12.
    case = int(pH)
    if case == 1:
        result = 0.14 - (pH - int(pH)) * 0.01
    elif case == 2:
        result = 0.13 - (pH - int(pH)) * 0.05
    elif case == 3:
        result = 0.08 - (pH - int(pH)) * 0.037
    elif case == 4:
        result = 0.043 - (pH - int(pH)) * 0.018
    elif case == 5:
        result = 0.025 - (pH - int(pH)) * 0.02
    elif case == 6:
        result = 0.005 - (pH - int(pH)) * 0.015
    elif case == 7:
        result = -0.01 - (pH - int(pH)) * 0.01
    elif case == 8:
        result = -0.02 - (pH - int(pH)) * 0.005
    elif case == 9:
        result = -0.025 - (pH - int(pH)) * 0.01
    elif case == 10:
        result = -0.035 - (pH - int(pH)) * 0.06
    elif case == 11:
        result = -0.095 - (pH - int(pH)) * 0.105
    elif case == 12:
        result = -0.20 - (pH - int(pH)) * 1
    result = BSA * (result - correction_to_IEP) / (15e-4 * 66463)  # Hydrogen ions for deviation from the IEP.
    return result


