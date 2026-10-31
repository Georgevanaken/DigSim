"""Reader for the meal file (compound table and meal schedule)."""
from typing import IO, Any

from . import state as st
from .state import diet


class PascalTextReader:
    """Minimal reader for the whitespace and line semantics used by Pascal Text files."""

    def __init__(self, file: IO[str]) -> None:
        self.file = file
        self.text = file.read()
        self.position = 0

    def _skip_whitespace(self) -> None:
        """Advance past whitespace."""
        while self.position < len(self.text) and self.text[self.position].isspace():
            self.position += 1

    def read(self, value_type: type = str) -> Any:
        """Read one character (chr) or one whitespace-delimited token converted to int, float or str."""
        if value_type is chr:
            if self.position >= len(self.text):
                raise EOFError("read past end of file")
            value = self.text[self.position]
            self.position += 1
            return value
        self._skip_whitespace()
        start = self.position
        while self.position < len(self.text) and not self.text[self.position].isspace():
            self.position += 1
        if start == self.position:
            raise EOFError("read past end of file")
        token = self.text[start:self.position]
        if value_type is int:
            return int(float(token))
        if value_type is float:
            return float(token)
        return token

    def readln(self, value_type: type | None = str) -> Any:
        """Read a value (or the rest of the line for str) and skip to the start of the next line."""
        if value_type is str:
            value = ""
            while self.position < len(self.text) and self.text[self.position] not in "\r\n":
                value += self.text[self.position]
                self.position += 1
        elif value_type is None:
            value = ""
        else:
            value = self.read(value_type)
        while self.position < len(self.text) and self.text[self.position] not in "\r\n":
            self.position += 1
        if self.position < len(self.text) and self.text[self.position] == "\r":
            self.position += 1
        if self.position < len(self.text) and self.text[self.position] == "\n":
            self.position += 1
        return value
    def close(self) -> None:
        """Close the underlying file."""
        self.file.close()


def _read_until_closing_paren(reader: PascalTextReader) -> None:
    """Skip characters up to and including the next ')' (used to skip comments in the meal file)."""
    char = " "
    while char != ")":
        char = reader.read(chr)


def read_diet(meal_file: str = "meal.txt") -> None:
    """Read the compound table and the meal schedule from the meal file into the diet group."""
    reader = PascalTextReader(open(meal_file, encoding="utf-8"))
    st.data_file = reader
    for compoundnumber in range(51):
        for j in range(1, 111):
            diet.compound_composition[compoundnumber, j] = 0
            diet.oral_dec_frac[compoundnumber, 1] = 0
            diet.oral_dec_frac[compoundnumber, 2] = 0
            diet.gastric_dec_frac[compoundnumber, 1] = 0
            diet.gastric_dec_frac[compoundnumber, 2] = 0
            diet.intestinal_dec_frac[compoundnumber, 1] = 0
            diet.intestinal_dec_frac[compoundnumber, 2] = 0
    reader.readln(str)
    _read_until_closing_paren(reader)
    st.clock.reference_time = reader.readln(float)
    reader.readln(str)
    reader.readln(str)
    reader.readln(str)

    diet.meal_number = 0
    st.clock.minutes_max = -1
    while True:
        time_value = reader.read(float)
        if time_value > st.clock.minutes_max and time_value != 9999:
            st.clock.minutes_max = time_value
            diet.meal_number += 1
            diet.meal[diet.meal_number] = reader.read(int)
            diet.meal_time[diet.meal_number] = time_value
            diet.meal_amount[diet.meal_number] = reader.read(float)
            diet.dish[diet.meal_number] = reader.read(int)
            diet.eating_speed[diet.meal_number] = reader.read(float)
        if time_value == 9999:
            break
    diet.number_of_meals = diet.meal_number
    st.clock.minutes_max = 500 + diet.meal_time[diet.meal_number]

    reader.readln(str)
    diet.number_of_compounds = 0
    _read_until_closing_paren(reader)
    aux = reader.read(int)
    while aux != 9999:
        diet.number_of_compounds = aux
        aux = reader.read(int)

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.compound_ph[compoundnumber] = reader.read(float)
    reader.readln(str)

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        st.solid[compoundnumber] = reader.read(float)
    reader.readln(str)

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.compound_viscosity[compoundnumber] = reader.read(float)
    reader.readln(str)

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.oral_dec_time[compoundnumber] = reader.read(float)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.oral_dec_comp[compoundnumber, 1] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.oral_dec_frac[compoundnumber, 1] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.oral_dec_comp[compoundnumber, 2] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.oral_dec_frac[compoundnumber, 2] = reader.read(int)
    reader.readln(str)

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.gastric_dec_time[compoundnumber] = reader.read(float)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.gastric_dec_comp[compoundnumber, 1] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.gastric_dec_frac[compoundnumber, 1] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.gastric_dec_comp[compoundnumber, 2] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.gastric_dec_frac[compoundnumber, 2] = reader.read(int)
    reader.readln(str)

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.intestinal_dec_time[compoundnumber] = reader.read(float)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.intestinal_dec_comp[compoundnumber, 1] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.intestinal_dec_frac[compoundnumber, 1] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.intestinal_dec_comp[compoundnumber, 2] = reader.read(int)
    reader.readln(str)
    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.intestinal_dec_frac[compoundnumber, 2] = reader.read(int)
    reader.readln(str)

    reader.readln(str)
    reader.readln(str)
    reader.readln(str)

    for j in range(1, st.comp + 1):
        for compoundnumber in range(1, diet.number_of_compounds + 1):
            diet.compound_composition[compoundnumber, j] = 0
            diet.active[j] = False
    while True:
        j = reader.read(int)
        if j != 0:
            _read_until_closing_paren(reader)
            for compoundnumber in range(1, diet.number_of_compounds + 1):
                diet.compound_composition[compoundnumber, j] = reader.read(float)
                diet.compound_composition[compoundnumber, j] *= 0.01
                if diet.compound_composition[compoundnumber, j] != 0:
                    diet.active[j] = True
            reader.readln(str)
        else:
            break

    _read_until_closing_paren(reader)
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        diet.compound_density[compoundnumber] = reader.read(float)
    reader.close()

    for compoundnumber in range(1, diet.number_of_compounds + 1):
        st.calories.compound[compoundnumber] = 0
        for j in range(1, st.comp + 1):
            st.calories.compound[compoundnumber] += (
                st.calories.value[j] * diet.compound_composition[compoundnumber, j]
            )


def read_compounds() -> None:
    """Add the current meal item to its dish: amount, pH, viscosity and eating speed (compound-table lookup)."""
    st.swallow_speed = 1000
    compoundnumber = diet.meal[diet.new_meal_number]
    dishnumber = diet.dish[diet.new_meal_number]
    # NOTE(port): meal.txt uses compound numbers > 50; Pascal reads past the end of these 1..50 arrays and gets 0.
    in_range = compoundnumber < len(diet.compound_ph)
    diet.ph_dish[dishnumber] = diet.compound_ph[compoundnumber] if in_range else 0.0
    diet.liquid_viscosity_dish[dishnumber] = diet.compound_viscosity[compoundnumber] if in_range else 0.0
    st.eating_speed_dish[dishnumber] = diet.eating_speed[diet.new_meal_number]
    st.compound_dish[dishnumber, compoundnumber] += diet.meal_amount[diet.new_meal_number]
    print(
        "mealnr = ",
        diet.new_meal_number,
        "  minutes =",
        f"{st.clock.minutes:8.0f}",
        "  comp. nr. = ",
        compoundnumber,
        sep="",
    )




