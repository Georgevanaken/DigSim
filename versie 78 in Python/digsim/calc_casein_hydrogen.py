"""H+ binding of casein as a function of pH."""
# Based on Pertzoff and Carpenter (1932) and Erdos, effect of salts on the isoelectric point of casein.
def casein_hydrogen(pH: float, casein: float) -> float:
    """Return the H+ bound by casein (amount casein) at the given pH, relative to its IEP."""
    correction_to_IEP = -0.00002  # IEP set to 4.5; corrects the Pertzoff and Carpenter titration curve.
    if pH < 2:
        bound_hydrogen = 0.0016 + (2 - pH) * (0.002 - 0.0016)
    elif pH < 2.5:
        bound_hydrogen = 0.0009 + (2.5 - pH) * (0.0016 - 0.0009) * 2
    elif pH < 3:
        bound_hydrogen = 0.0004 + (3 - pH) * (0.0009 - 0.0004) * 2
    elif pH < 3.5:
        bound_hydrogen = 0.0001 + (3.5 - pH) * (0.0004 - 0.0001) * 2
    elif pH < 4:
        bound_hydrogen = 0.00005 + (4 - pH) * (0.0001 - 0.00005) * 2
    elif pH < 4.5:
        bound_hydrogen = 0.00002 + (4.5 - pH) * (0.00005 - 0.00002) * 2
    elif pH < 5:
        bound_hydrogen = 0 + (5 - pH) * (0.00002 - 0) * 2
    elif pH < 6:
        bound_hydrogen = -0.00002 + (6 - pH) * (0 + 0.00002)
    elif pH < 7:
        bound_hydrogen = -0.0003 + (7 - pH) * (-0.00002 + 0.0003)
    elif pH < 8:
        bound_hydrogen = -0.0007 + (8 - pH) * (-0.0003 + 0.0007)
    elif pH < 9:
        bound_hydrogen = -0.0008 + (9 - pH) * (-0.0007 + 0.0008)
    else:
        bound_hydrogen = -0.0011 + (10 - pH) * (-0.0008 + 0.0011)
    result = casein * (bound_hydrogen + correction_to_IEP)  # Hydrogen ions for deviation from the IEP.
    return result


