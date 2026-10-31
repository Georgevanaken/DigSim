# digsim-py

A Python port of **DigSim v78** ("effect of salt added"), the literature-based simulator of human digestion and
absorption written in Free Pascal by George van Aken (www.insightfoodinside.com).

Licence: GNU Affero General Public License v3 or later (same as the original). See [LICENSE](LICENSE).
The scientific model, parameters and meal format are the original author's work; this repository only translates
the code. It is not an official release.

## Install and run

```
pip install -e .
digsim -n 20 --meal meal.txt --results Results.txt
```

- `-n` steps per minute. Use at least 20; 6 is unstable (the stomach cycle is 0.111 min).
- `--meal` meal definition file (same format as the Pascal `meal.txt`), default `meal.txt`.
- `--results` output table (166 columns, one row per simulated minute), default `Results.txt`.
- `--minutes-max` stop early.

Without `-n` the program asks for the step rate, like the Pascal original. A full 3500-minute run at 20 steps/min
takes about 5 minutes on a normal PC (the 500-minute fasted initialisation is part of that); the Pascal original needs
well under a minute.

## Validation

The port was checked against the Pascal source compiled with Free Pascal 3.2.2 at the same step rate. Results agree
on 29 key columns to within about 5% of each column's range (serum glucose and gastric volume within about 0.5%).
Exact agreement is not possible: the pH bisection is discontinuous and Pascal uses x87 extended precision.
The supplied `Results.txt` of the Pascal distribution is not reproduced by the Pascal code itself at 20 or 60
steps/min, so it was not used as the reference. Details are in [TODO.md](TODO.md).

## Tests

```
pip install -e .[test]
python -m pytest -q                      # fast tests, seconds
DIGSIM_SLOW=1 python -m pytest -q        # adds regression runs against Pascal baselines (about 6 min in total)
```

`tools/build_pascal_reference.ps1` rebuilds the Pascal reference output (needs Free Pascal).
