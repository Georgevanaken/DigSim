"""Definition of the second batch of dataclass groups, plus the one-off rewrite helper.

NEW_GROUPS (english global name -> field) is imported by tools/gen_state.py, which emits the dataclasses.
Run as a script (from DigSim_py) to rewrite the code that used the old globals:
  python tools\\group_state.py            dry run of the 'st.<name>' -> 'st.<group>.<field>' rewrite
  python tools\\group_state.py --apply    apply the rewrite to digsim/, tests/ and tools/
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# instance name -> (class name, {current global name: field name})
NEW_GROUPS = {
    "calories": ("Calories", {
        "calories_absorbed": "absorbed", "calories_brush": "brush", "calories_colon": "colon",
        "calories_distal": "distal", "calories_proximal": "proximal", "calories_compound": "compound",
        "caloric_density_local": "density_local", "caloric_density_mucus": "density_mucus",
        "caloric_value": "value", "total_caloric_absorption_rate": "total_absorption_rate",
        "max_caloric_absorption_capacity": "max_absorption_capacity",
    }),
    "hunger": ("Hunger", {
        "hunger": "level", "hunger_detectable_nutrients": "detectable_nutrients",
        "hunger_intestinal_volume": "intestinal_volume",
        "hunger_si_absorbable_calories": "si_absorbable_calories", "hunger_si_calories": "si_calories",
        "hunger_total_caloric_absorption_rate": "total_caloric_absorption_rate",
    }),
    "beta_cell": ("BetaCell", {
        "dphi": "dphi", "f1_pederson": "f1", "fused": "fused", "fused_increment": "fused_increment",
        "intermediate": "intermediate", "intermediate_increment": "intermediate_increment",
        "phi": "phi", "released": "released", "rrp": "rrp", "rrp_increment": "rrp_increment",
    }),
    "lipolysis": ("Lipolysis", {
        "gastric_lipolysis_cst": "gastric_cst", "lingual_lipolysis_cst": "lingual_cst",
        "pancreatic_lipolysis_cst": "pancreatic_cst",
        "percentual_loss_of_triglycerides_gastric": "percentual_loss_gastric",
        "percentual_loss_of_triglycerides_lingual": "percentual_loss_lingual",
        "percentual_loss_of_triglycerides_pancreatic": "percentual_loss_pancreatic",
        "total_gastric_lipase_lipolysis_rate": "total_gastric_rate",
        "total_lingual_lipase_lipolysis_rate": "total_lingual_rate",
        "total_pancreatic_lipase_lipolysis_rate": "total_pancreatic_rate",
    }),
    "secretion": ("Secretion", {
        "saliva_addition": "saliva_addition", "total_saliva_secretion": "total_saliva",
        "pancreatic_fluid_secretion": "pancreatic_fluid",
        "total_pancreatic_fluid_secretion": "total_pancreatic_fluid",
        "bicarbonate_secretion": "bicarbonate",
    }),
    "osmotic": ("Osmotic", {
        "osmotic_concentration": "concentration", "osmotic_value_brush": "brush",
        "osmotic_value_brush_distal": "brush_distal", "osmotic_value_brush_duodenum": "brush_duodenum",
        "osmotic_value_lumen": "lumen", "osmotic_value_lumen_distal": "lumen_distal",
        "osmotic_value_lumen_duodenum": "lumen_duodenum",
    }),
    "brush_enzymes": ("BrushEnzymes", {
        "sucrase_inhibition": "sucrase_inhibition", "sucrase_presence": "sucrase_presence",
        "sucrase_transformation": "sucrase_transformation",
        "maltase_inhibition": "maltase_inhibition", "maltase_presence": "maltase_presence",
        "maltase_transformation": "maltase_transformation",
        "inhibition_competitive": "inhibition_competitive",
        "inhibition_uncompetitive": "inhibition_uncompetitive",
        "sugar_flux": "sugar_flux", "water_flux": "water_flux",
    }),
    "clock": ("Clock", {
        "minutes": "minutes", "minutes_max": "minutes_max", "start_time": "start_time",
        "reference_time": "reference_time", "time_of_day": "time_of_day", "dt": "dt", "substeps": "substeps",
    }),
    "solid_props": ("SolidProperties", {
        "solids_g_prime": "g_prime", "solids_ph": "ph", "solids_strain": "strain",
        "solids_viscosity": "viscosity",
    }),
    # existing group: extra members
    "glucose": ("Glucose", {
        "glucose_level": "level", "mobilization_level": "mobilization_level",
        "fructose_decay_rate": "fructose_decay_rate", "incretin": "incretin",
        "g0_in": "g0_in", "g0_increment": "g0_increment", "g1_in": "g1_in", "g1_increment": "g1_increment",
        "g2_in": "g2_in", "g2_increment": "g2_increment", "g2_in_mol_per_min": "g2_in_mol_per_min",
        "glycogen_store0": "glycogen_store0", "glycogen_store1": "glycogen_store1",
        "glycogen_store2": "glycogen_store2",
        "hba1c": "hba1c", "hba1c_half_life": "hba1c_half_life",
        "glycosylation_rate_constant": "glycosylation_rate_constant",
    }),
}

MAPPING = {old: (inst, field) for inst, (_, fields) in NEW_GROUPS.items() for old, field in fields.items()}

def new_ref(inst, field):
    # The existing groups are imported by name (from .state import glucose); the new ones go through st.
    return f"{inst}.{field}" if inst == "glucose" else f"st.{inst}.{field}"


def main():
    pattern = re.compile(r"\bst\.(" + "|".join(sorted(MAPPING, key=len, reverse=True)) + r")\b")
    apply = "--apply" in sys.argv
    total = 0
    for folder in ("digsim", "tests", "tools"):
        for path in sorted((ROOT / folder).glob("*.py")):
            if path.name in ("state.py", "group_state.py", "gen_state.py"):
                continue
            text = path.read_bytes().decode("utf-8")
            new, n = pattern.subn(lambda m: new_ref(*MAPPING[m.group(1)]), text)
            if n:
                total += n
                print(f"{path.relative_to(ROOT)}: {n}")
                if apply:
                    path.write_bytes(new.encode("utf-8"))
    print("total", total, "(applied)" if apply else "(dry run)")


if __name__ == "__main__":
    main()
