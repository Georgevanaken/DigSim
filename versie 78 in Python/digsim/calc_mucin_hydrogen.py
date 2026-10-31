"""H+ binding of mucin as a function of pH."""
# Based on Leon Li et al., Microfluidic in vitro model for quantitative study of stomachmucinacid barrier function, 2011, p. 771, Fig. C.
def mucin_hydrogen(pH: float, Mucin: float) -> float:
    """Return the H+ bound by mucin (amount Mucin) at the given pH, relative to its IEP."""
    correction_to_IEP = 0.0005  # IEP is between 2 and 3.
    result = 0.0  # NOTE(port): uninitialised in Pascal when trunc(pH) is outside 1..12.
    case = int(pH)
    if case == 1:
        result = 0.001333333 - (pH - int(pH)) * 0.000566667
    elif case == 2:
        result = 0.000766667 - (pH - int(pH)) * 0.000466667
    elif case == 3:
        result = 0.0003 - (pH - int(pH)) * 0.00016
    elif case == 4:
        result = 0.00014 - (pH - int(pH)) * 0.00012667
    elif case == 5:
        result = 0.00001333 - (pH - int(pH)) * 4.663e-05
    elif case == 6:
        result = -0.0000333 - (pH - int(pH)) * 0.00003
    elif case == 7:
        result = -0.0000633 - (pH - int(pH)) * 0.00001
    elif case == 8:
        result = -0.0000733 - (pH - int(pH)) * 0.000030033
    elif case == 9:
        # NOTE(port): source uses 4.3334-05 (subtraction), not scientific notation; preserve it.
        result = -0.000103333 - (pH - int(pH)) * (4.3334 - 5)
    elif case == 10:
        result = -0.000146667 - (pH - int(pH)) * 0.000186666
    elif case == 11:
        result = -0.000333333 - (pH - int(pH)) * 0.000333333
    elif case == 12:
        result = -0.001 - (pH - int(pH)) * 0.000666667
    result = Mucin * (result - correction_to_IEP)  # Hydrogen ions for deviation from the IEP.
    return result



