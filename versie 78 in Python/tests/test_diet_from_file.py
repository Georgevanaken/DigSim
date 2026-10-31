from pathlib import Path

from digsim import diet_from_file
from digsim import state as st


def test_read_diet_parses_meal_file_and_compounds(monkeypatch, capsys):
    meal_path = Path(__file__).parent / "data"
    monkeypatch.chdir(meal_path)

    diet_from_file.read_diet()

    d = st.diet
    assert d.number_of_meals == 16
    assert list(d.meal_time[1:17]) == [0, 1, 500, 501, 502, 1000, 1001, 1002, 2000, 2001, 2500, 2501, 2502, 3000, 3001, 3002]
    assert d.meal[1] == 1
    assert d.meal[2] == 2
    assert d.meal_amount[1] == 250
    assert d.dish[1] == 1
    assert d.number_of_compounds == 25
    assert d.compound_composition[1, 1] == 1.0
    assert d.compound_composition[2, 3] == 0.99
    assert d.compound_composition[4, 2] == 0.3273
    assert d.active[1]
    assert d.compound_density[1] == 1.0

    d.new_meal_number = 1
    diet_from_file.read_compounds()
    assert st.compound_dish[1, 1] == 250
    assert st.eating_speed_dish[1] == 10
    assert "mealnr = 1" in capsys.readouterr().out
