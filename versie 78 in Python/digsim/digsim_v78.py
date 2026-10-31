"""Main simulation driver for DigSim version 78."""
import argparse

import numpy as np

from . import state as st
from .calc_protein_hydrogen import calc_all_ph, calc_hydrogen
from .compartment_properties import volumina
from .diet_from_file import read_compounds, read_diet
from .enums import Component as C
from .enums import Segment
from .gastric_acid_secretion import gastric_secretion
from .gastric_motility import solids_in_stomach
from .glucose_homeostasis2 import insulin_dynamics
from .hb_glycosylation import calc_hba1c
from .liver_clearance import decay_in_blood, fructose_to_glucose
from .main_part1 import (
    calc_outflows,
    caloric_densities_segments,
    calories,
    initialise,
    intestinal_secretion,
    osmotic_values_brush,
    osmotic_values_lumen,
    pancreatobiliary_secretion,
    saliva_secretion,
)
from .main_part2 import decomposition, digestion, pepsin_units
from .main_part3 import (
    absorption_in_small_intestine,
    colon,
    gut_hormones,
    intravenous,
    intubation_in_antrum,
    intubation_in_duodenum,
    mucus_layer,
    sensations,
)
from .mmc import update_mmc_status
from .model_parameters import set_model_parameters
from .state import bile, diet, glucose, hormones, stomach


def make_step() -> None:
    """Activate all processes during a time step."""
    emptying = 0.0
    decelerator = 0
    compoundnumber = 0

    update_mmc_status()
    volumina()
    calories()
    decomposition()
    digestion()
    osmotic_values_lumen()
    osmotic_values_brush()
    absorption_in_small_intestine()
    mucus_layer()
    intubation_in_antrum()
    intubation_in_duodenum()
    intravenous()
    fructose_to_glucose()
    insulin_dynamics()
    calc_hba1c()
    gut_hormones()
    decay_in_blood()
    calc_all_ph()
    saliva_secretion()
    gastric_secretion()
    pancreatobiliary_secretion()
    intestinal_secretion()
    caloric_densities_segments()
    solids_in_stomach()
    colon()

    volumina()
    calc_outflows()
    sensations()
    emptying = 0
    decelerator = st.clock.substeps
    for j in range(1, st.comp + 1):
        emptying = emptying + st.outflow[Segment.ANTRUM] * st.concentration_lumen[Segment.ANTRUM, j] * st.calories.value[j]
    for compoundnumber in range(1, diet.number_of_compounds + 1):
        emptying = emptying + st.outflow[Segment.ANTRUM] * stomach.gastric_sieve[compoundnumber] * st.concentration_compound[Segment.ANTRUM, compoundnumber] * st.calories.compound[compoundnumber]
    stomach.gastric_emptying = (decelerator * stomach.gastric_emptying + emptying) / (decelerator + 1)

    transport()


def _transport_block(
    cur: np.ndarray, nxt: np.ndarray, conc: np.ndarray, dish: np.ndarray, dishnext: np.ndarray,
    concdish: np.ndarray, cols: slice, swallow: float | None = None,
) -> None:
    """Move material between compartments 2..16 and out of the dishes (all columns at once)."""
    dt = st.clock.dt
    up = st.outflow
    mouth, corpus, antrum, duo1, colon_row = (
        Segment.MOUTH, Segment.CORPUS, Segment.ANTRUM, Segment.DUODENUM_1, Segment.COLON)
    below = slice(Segment.DUODENUM_2, colon_row)  # duodenum 2 .. ileum 6
    above = slice(duo1, colon_row - 1)            # the segment each of those receives from
    dishes = slice(1, 6)
    out_of_dishes = st.outflow_dish[dishes, None] * concdish[dishes, cols] * dt
    dishnext[dishes, cols] = dish[dishes, cols] - out_of_dishes
    inflow = out_of_dishes.sum(axis=0)
    flow_swallowed = up[mouth] * conc[mouth, cols] * dt if swallow is None else up[mouth] * swallow * conc[mouth, cols] * dt
    nxt[mouth, cols] = cur[mouth, cols] + inflow - flow_swallowed
    nxt[corpus, cols] = cur[corpus, cols] + flow_swallowed
    nxt[antrum, cols] = cur[antrum, cols]
    nxt[duo1, cols] = cur[duo1, cols] - up[duo1] * conc[duo1, cols] * dt
    if swallow is None:
        nxt[below, cols] = cur[below, cols] + up[above, None] * conc[above, cols] * dt
        nxt[below, cols] = nxt[below, cols] - up[below, None] * conc[below, cols] * dt
    else:
        nxt[below, cols] = (
            cur[below, cols]
            + (up[above, None] * conc[above, cols] - up[below, None] * conc[below, cols]) * dt
        )
    nxt[colon_row, cols] = cur[colon_row, cols] + up[colon_row - 1] * conc[colon_row - 1, cols] * dt


def transport() -> None:
    """Move liquid between neighbouring segments and dishes (outflow of each segment into the next)."""
    liq = slice(1, st.comp + 1)
    _transport_block(st.liquid, st.liquid_next, st.concentration_lumen,
                     st.liquid_dish, st.liquid_dish_next, st.concentration_dish, liq)

    cmp_ = slice(1, diet.number_of_compounds + 1)
    swallow = np.where(st.solid[cmp_] > 10, 0, 1)
    _transport_block(st.compound, st.compound_next, st.concentration_compound,
                     st.compound_dish, st.compound_dish_next, st.concentration_compound_dish,
                     cmp_, swallow)

    st.compound_dish[1:6, cmp_] = st.compound_dish_next[1:6, cmp_]
    st.liquid_dish[1:6, liq] = st.liquid_dish_next[1:6, liq]
    st.compound[2:17, cmp_] = st.compound_next[2:17, cmp_]
    st.liquid[2:17, liq] = st.liquid_next[2:17, liq]


def _write_header() -> None:
    """Write the column headings of the results file."""
    st.out.write(' minutes,')
    st.out.write(' plate, mouth, corpus, antrum, duodenum 1,duodenum 2, jejunum 1, jejunum 2, jejunum 3, ileum 1, ileum 2, ileum 3, ileum 4, ileum 5, ileum 6, colon,')
    # st.out.write('FFA upper gastric, FFA antrum, FFA duodenum 1, FFA dusodenum 2, FFA jejunum 1, FFA jejunum 2, FFA jejunum 3, FFA ileum 1, FFA ileum 2, FFA ileum 3, FFA ileum 4, FFA ileum 5, FFA ileum 6, FFA colon,')
    st.out.write(' enterocytes/100,')
    st.out.write(' total volume SI, CCK, PYY, gastric emptying, calories proximal, calories distal, calories colon, calories absorbed, gastric tone,')
    st.out.write(' fullness, hunger,')
    st.out.write(' peptides upper gastric, peptides antrum, peptides duodenum 1, peptides duodenum 2, peptides jejunum 1, peptides jejunum 2, peptides jejunum 3, peptides ileum 1, peptides ileum 2, peptides ileum 3, peptides ileum 4, peptides ileum 5, peptides ileum 6, peptides colon,')
    st.out.write(' GIP, GLP1,')
    st.out.write(' glucose upper gastric, glucose antrum, glucose duodenum 1, glucose duodenum 2, glucose jejunum 1, glucose jejunum 2, glucose jejunum 3, glucose ileum 1, glucose ileum 2, glucose ileum 3, glucose ileum 4, glucose ileum 5, glucose ileum 6, glucose colon,')
    st.out.write(' fat delivered to colon, casein delivered to colon, glucose delivered to colon,')
    st.out.write(' glucose absorbed, fatty acids absorbed, conjugated bile absorbed, fatty acid bile micelles absorbed, casein small peptides absorbed, small peptides absorbed,')
    st.out.write(' pH plate, pH mouth, pH upper gastric, pH antrum, pH duodenum 1, pH duodenum 2, pH jejunum 1, pH jejunum 2, pH jejunum 3, pH ileum 1, pH ileum 2, pH ileum 3, pH ileum 4, pH ileum 5, pH ileum 6, pH colon,')
    st.out.write(' Units pepsin in corpus, Units pepsin in antrum, solids_in_corpus, solids_in_fundus, solids_in_antrum, total solids in stomach,')
    st.out.write(' serum glucose (mMol/l), serum insulin (pMol/l), pharmaceutical absorbed (g), bloodserumconc pharmaceutical (mg/l),')
    st.out.write(' CCK emptying stress, Osmotic emptying stress, pH emptying stress,')
    st.out.write(' pH fundus, fundus,')
    st.out.write(' starch fundus, starch colon, sucrase_inhibitie[7], absorption_rate[7-3], absorption_rate[7-1], glucose in brush[7], arabinose in brush[7], sucrose in brush[7],')
    st.out.write(' sucrase_inhibitie[15], absorption_rate[15-3], absorption_rate[15-1], osmotic_value_brush_proximal, osmotic_value_brush_distal, osmotic_value_lumen_duodenum, osmotic_value_lumen_distal,')
    st.out.write(' serum_glucose_mg_dl, serum_insulin_mU_l, GLP1gastricemptying_stress, Viscousgastricemptying_stress,')
    st.out.write(' 13C serum_glucose (mMol/l), 2H serum_glucose (mMol/l), glucose RaE (Mol/min), up[4], intercellulaire glucose muscle adipose (mMol/l), intercellulaire glucose skin (mMol/l), total_glucose_absorption_rate (Mol/min), Gastric Volume (mL), HbA1C (%), Leucine absorbed (g), total SI water absorption rate (g/min),')
    st.out.write(' MMC status, Ghrelin, Motilin, bile_secretion_rate (g/min), Pancreatic fluid secretion (g/min), gastric fluid secretion (g/min), saliva secretion (g/min), ')
    st.out.write(' secreted gastric fluid (g), secreted conjugated bile salt (g), secreted panceratic fluid (g), secreted_saliva (g), gall bladder volume (g), conjugated bile salt concentration in gall bladder (g/g), duodenal conjugated bile salt concentration (g/g), gal bladder emptying rate (g/min), secteted bile (g), blood bile salt pool (g), unconjugated bile salt in colon (g), conjugated bile salt in secreted_bile (g/g), secretin (a.u.), glucagon (a.u.), glycogen (g),')
    st.out.write(' RC13/C12, total lingual lipase lipolysis rate (%/min), total gastric lipase lipolysis rate (%/min), total pancreatic lipase lipolysis rate (%/min),')
    st.out.write(' gastric emptying water, gastric emptying pharmaceutical, ')
    st.out.write(' Active units salivary amylase in mouth, Active units pancreatic amylase in duodenum, Active units pancreatic amylase at start of ileum, Active units/ml salivary amylase in mouth, Active units/ml pancreatic amylase in duodenum, Active units/ml pancreatic amylase at start of ileum, Percentual_loss_of_triglycerides_lingual, Percentual_loss_of_triglycerides_gastric, Percentual_loss_of_triglycerides_pancreatic')
    st.out.write("\n")


def _write_row() -> None:
    """Write one line of results (the state at the current minute)."""
    def w(value: float, width: int, dec: int) -> None:
        st.out.write(f"{value:{width}.{dec}f}")

    w(st.clock.minutes, 7, 5)
    st.out.write(', ')
    for i in range(1, 17):
        w(st.volume[i], 12, 8)
        st.out.write(', ')
    w(st.in_enterocytes, 7, 5)
    st.out.write(', ')
    for value in (
        st.volume_total, hormones.cck, hormones.pyy, stomach.gastric_emptying, st.calories.proximal,
        st.calories.distal, st.calories.colon, st.calories.absorbed, stomach.tone,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (st.fullness, st.hunger.level):
        w(value, 7, 5)
        st.out.write(', ')
    for i in range(Segment.CORPUS, Segment.COLON + 1):
        w(st.liquid[i, C.SMALL_PEPTIDES], 7, 5)
        st.out.write(', ')
    for value in (hormones.real_gip, hormones.real_glp1):
        w(value, 7, 5)
        st.out.write(', ')
    for i in range(Segment.CORPUS, Segment.COLON + 1):
        w(st.liquid[i, C.GLUCOSE], 7, 5)
        st.out.write(', ')
    for value in (st.fat_in_colon, st.casein_in_colon, glucose.in_colon):
        w(value, 7, 5)
        st.out.write(', ')
    for index in (C.GLUCOSE, C.FATTY_ACIDS, C.CONJUGATED_BILE_SALT, C.BILE_SALT_LCFA_COMPLEX, C.CASEIN_PEPTIDES, C.SMALL_PEPTIDES):
        w(st.absorbed[index], 7, 5)
        st.out.write(', ')
    for i in range(1, 17):
        w(st.ph[i], 12, 8)
        st.out.write(', ')
    for value in (
        pepsin_units(st.ph[Segment.CORPUS], st.liquid[Segment.CORPUS, C.PEPSIN]) / st.volume[Segment.CORPUS],
        pepsin_units(st.ph[Segment.ANTRUM], st.liquid[Segment.ANTRUM, C.PEPSIN]) / st.volume[Segment.ANTRUM],
        stomach.solids_in_corpus, stomach.solids_in_fundus, stomach.solids_in_antrum,
        stomach.total_solids_in_stomach,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (glucose.serum, glucose.serum_insulin, st.absorbed[C.PHARMACEUTICAL], st.blood_serum_conc[C.PHARMACEUTICAL]):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (hormones.cck_gastric_emptying_stress, stomach.osmotic_emptying_stress, stomach.ph_emptying_stress):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (stomach.ph_fundus, stomach.v_fundus):
        w(value, 7, 5)
        st.out.write(', ')
    w(st.liquid[Segment.CORPUS, C.HD_STARCH], 7, 5)
    st.out.write(',')
    w(st.liquid[Segment.COLON, C.HD_STARCH], 7, 5)
    st.out.write(', ')
    for value in (
        st.brush_enzymes.sucrase_inhibition[Segment.JEJUNUM_1], st.absorption_rate[Segment.JEJUNUM_1, C.GLUCOSE], st.absorption_rate[Segment.JEJUNUM_1, C.WATER],
        st.concentration_brush[7, 3], st.concentration_brush[7, 25], st.concentration_brush[7, 21],
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (
        st.brush_enzymes.sucrase_inhibition[Segment.ILEUM_6], st.absorption_rate[Segment.ILEUM_6, C.GLUCOSE], st.absorption_rate[Segment.ILEUM_6, C.WATER],
        st.osmotic.brush_duodenum, st.osmotic.brush_distal,
        st.osmotic.lumen_duodenum, st.osmotic.lumen_distal,
        glucose.serum_mg_dl, glucose.serum_insulin_mu_l, hormones.glp1_gastric_emptying_stress,
        stomach.viscous_emptying_stress, glucose.serum1, glucose.serum2,
        (glucose.g_in / (1e-2 + glucose.g1_in)) * glucose.g2_in_mol_per_min * glucose.serum1 / glucose.serum2,
        st.outflow[Segment.ANTRUM], glucose.intercellular_muscle_adipose, glucose.intercellular_skin,
        glucose.total_absorption_rate, stomach.volume,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    w(glucose.hba1c, 12, 10)
    st.out.write(', ')
    for value in (st.absorbed[70], st.total_water_absorption_rate):
        w(value, 7, 5)
        st.out.write(', ')
    st.out.write(f"{st.mmc_status:5d}, ")
    for value in (
        hormones.ghrelin, hormones.motilin, bile.total_secretion_rate, st.secretion.pancreatic_fluid,
        stomach.fluid_secretion, st.secretion.saliva_addition,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (
        stomach.total_fluid_secretion, bile.total_conjugated_salt_secretion,
        st.secretion.total_pancreatic_fluid, st.secretion.total_saliva, bile.gall_bladder_volume,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    w(bile.conjugated_salt_concentration_in_gall_bladder, 7, 5)
    st.out.write(', ')
    for value in (
        bile.duodenal_conjugated_salt_concentration, bile.gall_bladder_emptying_rate,
        bile.total_secreted,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (
        bile.blood_salt_pool, st.liquid[Segment.COLON, C.UNCONJUGATED_BILE_SALT], bile.conjugated_salt_in_secreted_bile,
        hormones.secretin, glucose.serum_glucagon,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (
        glucose.glycogen_store, st.r_c13_c12, st.lipolysis.total_lingual_rate,
        st.lipolysis.total_gastric_rate, st.lipolysis.total_pancreatic_rate,
        stomach.gastric_emptying_water, stomach.gastric_emptying_pharmaceutical,
    ):
        w(value, 7, 5)
        st.out.write(', ')
    for value in (
        st.salivary_amylase_activity[Segment.MOUTH], st.pancreatic_amylase_activity[Segment.DUODENUM_2],
        st.pancreatic_amylase_activity[Segment.ILEUM_1],
    ):
        st.out.write(str(value))
        st.out.write(', ')
    for value in (
        st.salivary_amylase_activity[Segment.MOUTH] / st.volume[Segment.MOUTH],
        st.pancreatic_amylase_activity[Segment.DUODENUM_2] / st.volume[Segment.DUODENUM_2],
        st.pancreatic_amylase_activity[Segment.ILEUM_1] / st.volume[Segment.ILEUM_1],
        st.lipolysis.percentual_loss_lingual,
        st.lipolysis.percentual_loss_gastric,
        st.lipolysis.percentual_loss_pancreatic,
    ):
        st.out.write(str(value))
        st.out.write(', ')
    st.out.write("\n")


def main(
    steps_per_minute: int | None = None,
    results_file: str = 'Results.txt',
    minutes_max: float | None = None,
    meal_file: str = 'meal.txt',
) -> None:
    """Run a simulation: read the meal file, step through time and write the results file."""
    # Pascal flushes underflow to zero silently; keep the other floating-point errors fatal.
    np.seterr(all='raise', under='ignore')

    print('DIGESTION and ABSORPTION simulation program (DigSim)')
    print('version 77, March 17, 2023')
    print('George van Aken, www.insightfoodinside.com, info@insightfoodinside.com')
    print()
    print('DigSim provides a literature-based simulation of the alimentary system, including food intake, mouth and gastrointestinal physiology and processing, and the response by the body.')
    print('Copyright (C) 2023 George van Aken')
    print()
    print('This program is free software:')
    print('you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, version 3 or any later version.')
    print()
    print('This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  ')
    print('See theGNU Affero General Public License for more details.')
    print()
    print('You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.')

    st.stop = 0
    if st.stop == 0:
        print()
        print()
        st.out = open(results_file, 'w')
        st.comp = 110
        set_model_parameters()
        read_diet(meal_file)
        if minutes_max is not None:
            st.clock.minutes_max = minutes_max
        st.counter = 0
        initialise()

        diet.new_meal_number = 0
        if steps_per_minute is None:
            print('Enter the time resolution for calculation in steps per minute (6 - 600)) = ', end='', flush=True)
            st.clock.substeps = int(input())
        else:
            st.clock.substeps = int(steps_per_minute)
        print()
        print()
        st.clock.dt = 1 / st.clock.substeps
        st.prefactor = 1
        st.prefactor = 1 / st.prefactor
        volumina()
        print('INITIALISING')

        # Reaching the fasted steady state
        st.clock.minutes = -500
        st.stars_per_line = 60
        dish = 1
        st.swallow_speed = 50
        st.ph[dish] = 7
        diet.liquid_viscosity_dish[dish] = 1
        st.liquid_dish[dish, 1] = 250
        st.liquid_dish[dish, 2] = 1
        st.liquid_dish[dish, 3] = 50
        st.liquid_dish[dish, 4] = 10
        st.liquid_dish[dish, 80] = st.liquid[1, 80] + calc_hydrogen(
            st.ph[1], st.liquid[1, 61], st.liquid[1, 66], st.liquid[1, 65],
            st.liquid[1, 67], st.liquid[1, 4], st.liquid[1, 68],
        )
        st.eating_speed_dish[1] = 50
        # Integer step counting avoids float drift that adds an extra step.
        for step in range(1, 500 * st.clock.substeps + 1):
            st.clock.minutes = -500 + step / st.clock.substeps
            st.counter = st.counter + 1
            make_step()
            if st.counter == st.stars_per_line * st.clock.substeps:
                print('*', end='', flush=True)
                st.counter = 0
        glucose.serum2_eq = glucose.serum2
        bile.blood_salt_pool = 0.008
        print('   DONE')

        st.calories.colon = 0
        st.calories.absorbed = 0
        stomach.total_fluid_secretion = 0
        st.secretion.total_pancreatic_fluid = 0
        bile.total_conjugated_salt_secretion = 0
        st.secretion.total_saliva = 0
        _write_header()

        st.counter = 0
        st.star_counter = 0
        st.clock.minutes = 0
        step = 0
        while st.clock.minutes < st.clock.minutes_max + 1:
            # Pascal's accumulated float time is slightly above the meal time at the first step; '>=' reproduces that.
            if st.clock.minutes >= diet.meal_time[diet.new_meal_number + 1] and diet.new_meal_number < diet.number_of_meals:
                print()
                diet.new_meal_number = diet.new_meal_number + 1
                read_compounds()

            if st.counter == st.prefactor * st.clock.substeps:
                st.star_counter = st.star_counter + 1
                if st.star_counter == st.stars_per_line:
                    print('*', end='', flush=True)
                    st.star_counter = 0

                colon_liq = st.liquid[Segment.COLON]
                st.fat_in_colon = colon_liq[C.FAT] + colon_liq[C.FATTY_ACIDS] + colon_liq[C.BILE_SALT_LCFA_COMPLEX]
                st.casein_in_colon = colon_liq[C.CASEIN] + colon_liq[C.CASEINPEP] + colon_liq[C.CASEINPAN] + colon_liq[C.CASEIN_PEPTIDES]
                glucose.in_colon = colon_liq[C.GLUCOSE] + colon_liq[C.TRACER_GLUCOSE_1] + colon_liq[C.TRACER_GLUCOSE_2]
                _write_row()
                st.counter = 0
            step += 1
            st.clock.minutes = step / st.clock.substeps
            st.counter = st.counter + 1
            make_step()

        st.out.close()


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the command-line arguments."""
    parser = argparse.ArgumentParser(
        prog='digsim', description='DigSim digestion and absorption simulation (Python port of DigSim v78)')
    parser.add_argument('-n', '--steps-per-minute', type=int, default=None,
                        help='time resolution; use at least 20 (6 is unstable)')
    parser.add_argument('--meal', default='meal.txt', help='meal definition file (default: meal.txt)')
    parser.add_argument('--results', default='Results.txt', help='output file (default: Results.txt)')
    parser.add_argument('--minutes-max', type=float, default=None, help='stop after this many simulated minutes')
    return parser.parse_args(argv)


def cli(argv: list[str] | None = None) -> None:
    """Command-line entry point."""
    args = _parse_args(argv)
    main(args.steps_per_minute, args.results, args.minutes_max, args.meal)


if __name__ == '__main__':
    cli()






