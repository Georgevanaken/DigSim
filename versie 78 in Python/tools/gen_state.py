"""Generate digsim/state.py from the interface-section variables of the Pascal units.

Pascal arrays are mapped to zero-initialised numpy arrays that are allocated with
size hi+1 per dimension, so the Pascal index can be used unchanged (index 0 of a
1-based array is simply unused).
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from state_renames import RENAMES

SRC = Path(sys.argv[1])
OUT = Path(sys.argv[2])

# units whose interface variables are shared, in Pascal 'uses' order (later wins on name clash)
UNITS = [
    "globalvariables",
    "diet_from_file",
    "glucose_homeostasis2",
    "subcellularinsulin_pederson",
    "brush_border_enzymes",
]


def strip_comments(t):
    t = re.sub(r"\{[^}]*\}", " ", t)
    t = re.sub(r"\(\*.*?\*\)", " ", t, flags=re.DOTALL)
    t = re.sub(r"//[^\n]*", " ", t)
    return t


def interface_text(name):
    t = (SRC / f"{name}.pas").read_text(encoding="utf-8", errors="replace")
    t = strip_comments(t)
    m = re.search(r"\binterface\b(.*?)\bimplementation\b", t, flags=re.DOTALL | re.IGNORECASE)
    return m.group(1)


types = {}
variables = {}  # name -> (type spec, unit); later units override earlier ones


def parse_array(spec):
    m = re.match(r"array\s*\[(.*?)\]\s*of\s+(\w+)", spec, flags=re.IGNORECASE)
    if not m:
        return None
    dims = []
    for part in m.group(1).split(","):
        lo, hi = [int(x) for x in part.strip().split("..")]
        dims.append(hi)
    return dims, m.group(2).lower()


def resolve(spec):
    spec = spec.strip()
    low = spec.lower()
    arr = parse_array(spec)
    if arr:
        return ("array", arr[0], arr[1])
    if low in types:
        return resolve(types[low])
    if low.startswith("string"):
        return ("scalar", "str")
    if low in ("real", "double", "extended", "single"):
        return ("scalar", "float")
    if low in ("integer", "word", "longint", "byte", "cardinal"):
        return ("scalar", "int")
    if low == "boolean":
        return ("scalar", "bool")
    if low == "text":
        return ("scalar", "none")
    raise ValueError(f"unknown type {spec}")


for unit in UNITS:
    text = interface_text(unit)
    mode = None
    for stmt in text.split(";"):
        s = " ".join(stmt.split())
        if not s:
            continue
        low = s.lower()
        m = re.match(r"(var|type|const)\b\s*(.*)", s, flags=re.IGNORECASE)
        if m:
            mode = m.group(1).lower()
            s = m.group(2)
            low = s.lower()
            if not s:
                continue
        elif re.match(r"(uses|procedure|function)\b", low):
            mode = None
            continue
        if mode == "type":
            n, spec = s.split("=", 1)
            types[n.strip().lower()] = spec.strip()
        elif mode == "var":
            names, spec = s.split(":", 1)
            for n in names.split(","):
                variables[n.strip()] = (spec.strip(), unit)


def init_expr(kind):
    if kind[0] == "scalar":
        return {"float": "0.0", "int": "0", "bool": "False", "str": "''", "none": "None"}[kind[1]]
    dims, base = kind[1], kind[2]
    shape = ", ".join(str(d + 1) for d in dims)
    dtype = {"real": "float", "double": "float", "integer": "int", "boolean": "bool", "word": "int"}.get(base)
    if dtype is None:
        raise ValueError(base)
    return f"np.zeros(({shape},), dtype={dtype})"


# Variables that live in a dataclass instead of a module global:
#   group name -> (class name, instance name, {Pascal variable: field name})
GROUPS = {
    "hormones": ("Hormones", "hormones", {
        "CCK": "cck",
        "CCKattenuationSI": "cck_attenuation_si",
        "CCKdegradationconstant": "cck_degradation_constant",
        "CCKEffectOnGastricEmptying": "cck_effect_on_gastric_emptying",
        "CCKenterocytes": "cck_enterocytes",
        "CCKgastricEmptying": "cck_gastric_emptying",
        "CCKgastricEmptying_stress": "cck_gastric_emptying_stress",
        "CCKresponsefactor": "cck_response_factor",
        "CCKtrans": "cck_trans",
        "Mesentric_CCK": "mesenteric_cck",
        "rateCCKproductionref": "rate_cck_production_ref",
        "totalrateCCKrelease": "total_rate_cck_release",
        "GLP1": "glp1",
        "GLP1attenuationSI": "glp1_attenuation_si",
        "GLP1degradationconstant": "glp1_degradation_constant",
        "GLP1enterocytes": "glp1_enterocytes",
        "GLP1gastricemptying_stress": "glp1_gastric_emptying_stress",
        "GLP1responsefactor": "glp1_response_factor",
        "Mesentric_GLP1": "mesenteric_glp1",
        "realGLP1": "real_glp1",
        "totalrateGLP1release": "total_rate_glp1_release",
        "PYY": "pyy",
        "PYYattenuationSI": "pyy_attenuation_si",
        "PYYdegradationconstant": "pyy_degradation_constant",
        "PYYenterocytes": "pyy_enterocytes",
        "PYYresponsefactor": "pyy_response_factor",
        "Mesentric_PYY": "mesenteric_pyy",
        "totalratePYYrelease": "total_rate_pyy_release",
        "Secretin": "secretin",
        "SecretinattenuationSI": "secretin_attenuation_si",
        "Secretindegradationconstant": "secretin_degradation_constant",
        "Secretinenterocytes": "secretin_enterocytes",
        "Secretinresponsefactor": "secretin_response_factor",
        "Mesentric_Secretin": "mesenteric_secretin",
        "totalrateSecretinrelease": "total_rate_secretin_release",
        "Gastrinenterocytes": "gastrin_enterocytes",
        "GIP": "gip",
        "GIPattenuationSI": "gip_attenuation_si",
        "GIPdegradationconstant": "gip_degradation_constant",
        "GIPenterocytes": "gip_enterocytes",
        "GIPresponsefactor": "gip_response_factor",
        "Mesentric_GIP": "mesenteric_gip",
        "RealGIP": "real_gip",
        "totalrateGIPrelease": "total_rate_gip_release",
        "Ghrelin": "ghrelin",
        "GhrelinattenuationSI": "ghrelin_attenuation_si",
        "Ghrelindegradationconstant": "ghrelin_degradation_constant",
        "Ghrelinresponsefactor": "ghrelin_response_factor",
        "Mesentric_Ghrelin": "mesenteric_ghrelin",
        "totalrateGhrelinrelease": "total_rate_ghrelin_release",
        "Motilin": "motilin",
        "MotilinattenuationSI": "motilin_attenuation_si",
        "Motilindegradationconstant": "motilin_degradation_constant",
        "Motilinenterocytes": "motilin_enterocytes",
        "Motilinresponsefactor": "motilin_response_factor",
        "Mesentric_Motilin": "mesenteric_motilin",
        "totalrateMotilinrelease": "total_rate_motilin_release",
    }),
    "glucose": ("Glucose", "glucose", {
        "Gluc": "gluc",
        "Gluc0": "gluc0",
        "Gluc1": "gluc1",
        "Gluc2": "gluc2",
        "glucose_absorption_rate": "absorption_rate",
        "glucose_in_colon": "in_colon",
        "total_glucose_absorption_rate": "total_absorption_rate",
        "intercellulair_glucose_muscle_adipose": "intercellular_muscle_adipose",
        "intercellulair_glucose_muscle_adipose_mg_dl": "intercellular_muscle_adipose_mg_dl",
        "intercellulair_glucose_skin": "intercellular_skin",
        "intercellulair_glucose_skin_mg_dl": "intercellular_skin_mg_dl",
        "serum_glucose": "serum",
        "serum_glucose0": "serum0",
        "serum_glucose0_mg_dl": "serum0_mg_dl",
        "serum_glucose1": "serum1",
        "serum_glucose1_mg_dl": "serum1_mg_dl",
        "serum_glucose2": "serum2",
        "serum_glucose2_eq": "serum2_eq",
        "serum_glucose2_mg_dl": "serum2_mg_dl",
        "serum_glucose_mg_dl": "serum_mg_dl",
        "endoglucprod": "endogenous_production",
        "Gin": "g_in",
        "gluci": "gluci",
        "gluciskin": "gluciskin",
        "glucutili": "utilisation",
        "glycogene_store": "glycogen_store",
        "Ii": "ii",
        "Ip": "ip",
        "serum_FA": "serum_fa",
        "serum_glucagon": "serum_glucagon",
        "serum_insulin": "serum_insulin",
        "serum_insulin_mU_l": "serum_insulin_mu_l",
        "serum_ketons": "serum_ketons",
        "serum_peptides": "serum_peptides",
        "tp": "tp",
        "tpf": "tpf",
        "x1": "x1",
        "x2": "x2",
        "x3": "x3",
    }),
    "bile": ("Bile", "bile", {
        "blood_bile_salt_pool": "blood_salt_pool",
        "conjugated_bile_salt_absorption_rate": "conjugated_salt_absorption_rate",
        "conjugated_bile_salt_concentration_in_gall_bladder": "conjugated_salt_concentration_in_gall_bladder",
        "conjugated_Bile_salt_in_secreted_bile": "conjugated_salt_in_secreted_bile",
        "conjugated_bile_salt_secretion_rate": "conjugated_salt_secretion_rate",
        "duodenal_conjugated_bile_salt_concentration": "duodenal_conjugated_salt_concentration",
        "gall_bladder": "gall_bladder",
        "gall_bladder_emptying_rate": "gall_bladder_emptying_rate",
        "gall_bladder_volume": "gall_bladder_volume",
        "total_bile_secreted": "total_secreted",
        "total_bile_secretion_rate": "total_secretion_rate",
        "total_conjugated_bile_salt_secretion": "total_conjugated_salt_secretion",
    }),
    "diet": ("Diet", "diet", {
        "active": "active",
        "compoundcomposition": "compound_composition",
        "compounddensity": "compound_density",
        "compoundGprime": "compound_g_prime",
        "compoundpH": "compound_ph",
        "compoundstrain": "compound_strain",
        "compoundviscosity": "compound_viscosity",
        "dish": "dish",
        "eating_speed": "eating_speed",
        "swallow_speed_meal": "swallow_speed_meal",
        "pHdish": "ph_dish",
        "liquidviscositydish": "liquid_viscosity_dish",
        "meal": "meal",
        "mealamount": "meal_amount",
        "mealnumber": "meal_number",
        "mealtime": "meal_time",
        "newmealnumber": "new_meal_number",
        "number_of_compounds": "number_of_compounds",
        "number_of_meals": "number_of_meals",
        "oral_dec_comp": "oral_dec_comp",
        "oral_dec_frac": "oral_dec_frac",
        "oral_dec_time": "oral_dec_time",
        "gastric_dec_comp": "gastric_dec_comp",
        "gastric_dec_frac": "gastric_dec_frac",
        "gastric_dec_time": "gastric_dec_time",
        "intestinal_dec_comp": "intestinal_dec_comp",
        "intestinal_dec_frac": "intestinal_dec_frac",
        "intestinal_dec_time": "intestinal_dec_time",
    }),
    "stomach": ("Stomach", "stomach", {
        "compound_fundus": "compound_fundus",
        "concentration_compound_fundus": "concentration_compound_fundus",
        "concentration_fundus": "concentration_fundus",
        "fundus": "fundus",
        "fundus_viscosity": "fundus_viscosity",
        "pHfundus": "ph_fundus",
        "Vfundus": "v_fundus",
        "Vcorpus": "v_corpus",
        "Vfunduscorpus": "v_fundus_corpus",
        "corpus_viscosity": "corpus_viscosity",
        "viscosity_liquids_in_corpus": "viscosity_liquids_in_corpus",
        "antral_viscosity": "antral_viscosity",
        "emptying_cycle": "emptying_cycle",
        "emptying_stress": "emptying_stress",
        "fraction_emptying": "fraction_emptying",
        "gastric_emptying": "gastric_emptying",
        "gastric_emptying_pharmaceutical": "gastric_emptying_pharmaceutical",
        "gastric_emptying_water": "gastric_emptying_water",
        "gastric_sieve": "gastric_sieve",
        "Gastric_tone": "tone",
        "Gastric_Volume": "volume",
        "BasisGastricEmptyingRate": "basis_emptying_rate",
        "OsmoticgastricEmptying": "osmotic_emptying",
        "OsmoticgastricEmptying_stress": "osmotic_emptying_stress",
        "pHgastricEmptying": "ph_emptying",
        "pHgastricEmptying_stress": "ph_emptying_stress",
        "ViscousgastricEmptying": "viscous_emptying",
        "Viscousgastricemptying_stress": "viscous_emptying_stress",
        "solids_antrum": "solids_antrum",
        "solids_corpus": "solids_corpus",
        "solids_fundus": "solids_fundus",
        "solids_in_antrum": "solids_in_antrum",
        "solids_in_corpus": "solids_in_corpus",
        "solids_in_fundus": "solids_in_fundus",
        "total_solids_antrum": "total_solids_antrum",
        "total_solids_corpus": "total_solids_corpus",
        "total_solids_in_stomach": "total_solids_in_stomach",
        "gastric_fluid_secretion": "fluid_secretion",
        "total_gastric_fluid_secretion": "total_fluid_secretion",
    }),
}
from group_state import NEW_GROUPS

_english_to_pascal = {RENAMES.get(p, p): p for p in variables}
for _inst, (_cls, _members) in NEW_GROUPS.items():
    _unknown = [e for e in _members if e not in _english_to_pascal]
    if _unknown:
        raise SystemExit(f"NEW_GROUPS names unknown variables: {_unknown}")
    _fields = {_english_to_pascal[e]: f for e, f in _members.items()}
    if _inst in GROUPS:
        GROUPS[_inst][2].update(_fields)
    else:
        GROUPS[_inst] = (_cls, _inst, _fields)
grouped =  {pv for _, _, fields in GROUPS.values() for pv in fields}
missing = grouped - set(variables)
if missing:
    raise SystemExit(f"GROUPS names unknown Pascal variables: {sorted(missing)}")


def field_line(field_name, spec, unit, pascal_name):
    kind = resolve(spec)
    comment = f"  # {unit}: {pascal_name}: {spec}"
    if kind[0] == "scalar":
        py_type = {"float": "float", "int": "int", "bool": "bool", "str": "str", "none": "None"}[kind[1]]
        return f"    {field_name}: {py_type} = {init_expr(kind)}{comment}"
    return f"    {field_name}: np.ndarray = field(default_factory=lambda: {init_expr(kind)}){comment}"


lines = [
    '"""Shared state of the DigSim model (generated by tools/gen_state.py -- do not edit by hand).',
    "",
    "Contains the interface variables of: " + ", ".join(UNITS) + ".",
    "Arrays keep their Pascal indices (array[1..n] is allocated with n+1 elements).",
    "Related variables are grouped in dataclass instances (see GROUPS in the generator).",
    '"""',
    "from dataclasses import dataclass, field",
    "",
    "import numpy as np",
    "",
]
for cls, inst, fields in GROUPS.values():
    lines += ["", "@dataclass", f"class {cls}:"]
    for pascal_name, field_name in fields.items():
        spec, unit = variables[pascal_name]
        lines.append(field_line(field_name, spec, unit, pascal_name))
    lines += ["", f"{inst} = {cls}()", ""]
for name, (spec, unit) in sorted(variables.items(), key=lambda kv: kv[0].lower()):
    if name not in grouped:
        lines.append(f"{RENAMES.get(name, name)} = {init_expr(resolve(spec))}  # {unit}: {name}: {spec}")
OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(len(variables), "variables written to", OUT)
