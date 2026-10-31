"""Enumerations for the segment, component and MMC-phase indices."""
from enum import IntEnum


class MMCStatus(IntEnum):
    """State of the migrating motor complex. Integer values match the Pascal original."""
    FED = 0
    PHASE_I = 1    # quiescence
    PHASE_II = 2   # irregular contractions
    PHASE_III = 3  # activity front
    PHASE_IV = 4   # transition back to quiescence


class Segment(IntEnum):
    """Row index of the location arrays (names and numbers from the header of digsim_v78.pas)."""
    MOUTH = 2
    CORPUS = 3
    ANTRUM = 4
    DUODENUM_1 = 5
    DUODENUM_2 = 6
    JEJUNUM_1 = 7
    JEJUNUM_2 = 8
    JEJUNUM_3 = 9
    ILEUM_1 = 10
    ILEUM_2 = 11
    ILEUM_3 = 12
    ILEUM_4 = 13
    ILEUM_5 = 14
    ILEUM_6 = 15
    COLON = 16  # the original distinguishes four colon segments, the model uses one


# Rows of the small intestine, as a slice for array indexing.
SMALL_INTESTINE = slice(int(Segment.DUODENUM_1), int(Segment.ILEUM_6) + 1)


class Component(IntEnum):
    """Column index of the component arrays (names and numbers from the header of digsim_v78.pas)."""
    WATER = 1
    FAT = 2
    GLUCOSE = 3
    PROTEIN = 4
    PHARMACEUTICAL = 5
    DIGLYCERIDE = 7
    MONOGLYCERIDE = 8
    HD_STARCH = 10
    SD_STARCH = 11
    R_STARCH = 12
    FATTY_ACIDS = 20
    PEPSIN = 40
    GASTRIC_LIPASE = 49
    BILE_LYSOPHOSPHOLIPIDS = 53
    CASEINPAN = 63
    I_C = 72
    I_A = 73
    I_NC = 74
    PANCREATIC_PROTEASES = 41
    LINGUAL_LIPASE = 48
    PANCREATIC_LIPASE = 50
    SUCROSE = 21
    FRUCTOSE = 22
    MALTOSE = 23
    TRACER_MALTOSE = 24      # 13C
    ARABINOSE = 25
    XYLOSE = 26
    ACARBOSE = 27
    TRACER_GLUCOSE_1 = 34    # 13C
    TRACER_GLUCOSE_2 = 35    # 2H
    TRACER_SUCROSE = 36      # 13C
    CONJUGATED_BILE_SALT = 51
    BILE_PHOSPHOLIPIDS = 52
    BILE_SALT_LCFA_COMPLEX = 54
    MUCUS_PROTEIN = 55
    PANCREATIC_PHOSPHOLIPASE = 56
    SALIVARY_AMYLASE = 57
    PANCREATIC_AMYLASE = 58
    UNCONJUGATED_BILE_SALT = 59
    CCK_RELEASING_PROTEIN = 60
    CASEIN = 61
    CASEINPEP = 62
    BSA = 66
    BLG = 67
    CITRIC_ACID = 68
    CASEIN_PEPTIDES = 64
    SMALL_PEPTIDES = 65
    LEUCINE = 70
    LEUCINE_TRACER = 71
    NACL = 78
    H_ION = 80               # in moles
