"""One-off helper: rewrite 'st.<old>' references to 'st.<new>' using tools/state_renames.py.

Usage (from DigSim_py): python tools\\rename_state.py [--apply]
Without --apply it only reports what would change.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from state_renames import RENAMES

ROOT = Path(__file__).resolve().parent.parent
pattern = re.compile(r"\bst\.(" + "|".join(sorted(RENAMES, key=len, reverse=True)) + r")\b")
apply = "--apply" in sys.argv

total = 0
for folder in ("digsim", "tests", "tools"):
    for path in sorted((ROOT / folder).glob("*.py")):
        if path.name in ("state.py", "state_renames.py", "rename_state.py"):
            continue
        raw = path.read_bytes()
        text = raw.decode("utf-8")
        new, n = pattern.subn(lambda m: "st." + RENAMES[m.group(1)], text)
        if n:
            total += n
            print(f"{path.relative_to(ROOT)}: {n}")
            if apply:
                path.write_bytes(new.encode("utf-8"))
print("total", total, "(applied)" if apply else "(dry run)")
