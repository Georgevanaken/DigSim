"""One-off helper: rename Dutch/Pascal function and local names to English (word-boundary rewrite).

Usage (from DigSim_py): python tools\\rename_identifiers.py [--apply]
Refuses to run if a new name already occurs as an identifier in the sources.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

RENAMES = {
    "initialiseer": "initialise",
    "parameter_invoer": "set_model_parameters",
    "betacelverdeling": "beta_cell_distribution",
    "compoundreader": "read_compounds",
    "ber_hydrogen": "calc_hydrogen",
    "ber_ph": "calc_ph",
    "ber_alle_phs": "calc_all_ph",
    "pancreabiliar_secretion": "pancreatobiliary_secretion",
    "ups": "calc_outflows",
    "mesentric_release": "mesenteric_release",
    "hulp": "aux",
    "hulp1": "aux1",
    "hulp2": "aux2",
    "hulp3": "aux3",
    "hulp4": "aux4",
    "hulp5": "aux5",
    "hulp6": "aux6",
    "hulp7": "aux7",
    "hulpreal": "time_value",
    "deler": "divisor",
    "vertrager": "decelerator",
    "richt_pH": "target_ph",
    "syneresesnelheid": "syneresis_rate",
    "vlokkingssnelheid_antrum": "flocculation_rate_antrum",
    "vlokkingssnelheid_corpus": "flocculation_rate_corpus",
    "oplossnelheid_vlokken_corpus": "floc_dissolution_rate_corpus",
    "pHvlokken": "floc_ph",
    "bolus_oppervlak": "bolus_surface",
    "max_gallbladdervolume": "max_gall_bladder_volume",
    "targetcorpus": "target_corpus",
    "targetcorpuscomp": "target_corpus_comp",
    "liquidscorpusfundus": "liquids_corpus_fundus",
    "alfa_cells_percentage": "alpha_cells_percentage",
}

files = [p for d in ("digsim", "tests", "tools") for p in sorted((ROOT / d).glob("*.py"))
         if p.name not in ("rename_identifiers.py", "rename_state.py", "state_renames.py")]
texts = {p: p.read_bytes().decode("utf-8") for p in files}

alltext = "\n".join(texts.values())
clash = [n for n in RENAMES.values() if re.search(rf"\b{re.escape(n)}\b", alltext)]
if clash:
    raise SystemExit(f"new names already in use: {clash}")

pattern = re.compile(r"\b(" + "|".join(sorted(RENAMES, key=len, reverse=True)) + r")\b")
apply = "--apply" in sys.argv
total = 0
for p, text in texts.items():
    new, n = pattern.subn(lambda m: RENAMES[m.group(1)], text)
    if n:
        total += n
        print(f"{p.relative_to(ROOT)}: {n}")
        if apply:
            p.write_bytes(new.encode("utf-8"))
print("total", total, "(applied)" if apply else "(dry run)")
