unit globalvariables;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils;

type
   locationcomposition = array[1..16,1..110] of real;
   solidslocationcomposition = array[0..99,1..16,1..1100] of real;
   solidsparameter = array[0..99,1..16] of real;
   MMC_timer = array[0..4] of real; {describes the duration of the different MMC statussen}
   segment = array [1..16] of real;
   component = array[1..110] of real;

var     liquid, liquidnext, liquiddish, liquiddishnext, compound, compounddish, compoundnext, compounddishnext, absorption_rate, enterocyte, concentration_lumen,
                            concentration_dish, concentration_compound, concentration_compound_dish, concentration_brush, VmaxAbsorption, KmAbsorption ,entering_serum,
                            sucrase_transformation, maltase_transformation, T, Mesentric, concentration_cell : locationcomposition;
        solids,solidsnext: solidslocationcomposition;
        solidspH, solidsviscosity, solidsGprime, solidsStrain: solidsparameter;
        hulp, minutes, minutesmax, calories_proximal, calories_distal, calories_colon, calories_absorbed, calories_brush, volumetotal, in_mucus, Transition, emptying_cycle,
                          Vfunduscorpus, pHfundus, Vfundus, viscosity, fraction, isotonic, CCKgastricEmptying_stress, OsmoticgastricEmptying_stress, pHgastricEmptying_stress,
                          Viscousgastricemptying_stress, Vcorpus, Gastric_Volume, emptying_stress, solids_in_corpus, solids_in_fundus, solids_in_antrum, total_solids_in_stomach,
                          cyclus, maxantralvolume, fraction_retropelled, fraction_emptying, rateCCKproductionref, totalrateCCKrelease, totalratePYYrelease, totalrateGLP1release,
                          totalrateGIPrelease, totalrateGhrelinrelease, totalrateMotilinrelease, totalrateSecretinrelease, gall_bladder_volume, gall_bladder_emptying_rate,
                          conjugated_bile_salt_concentration_in_gall_bladder,
                          duodenal_conjugated_bile_salt_concentration, conjugated_bile_salt_absorption_rate, blood_bile_salt_pool, conjugated_Bile_salt_in_secreted_bile,
                          solid_volume, large_solids, surface_solids, swallow_solid: real;
        i, j, teller, substeps, s, smax, sterretjes, starcounter, Stop, tel : integer;
        up, updish, volume, volumedish, CCKenterocytes, I_cell, L_cell, K_cell, G_cell, Mo_Cell,S_cell, PYYenterocytes, GIPenterocytes, GLP1enterocytes, Gastrinenterocytes,
                          Motilinenterocytes, Secretinenterocytes, Maxtotalabsorption_rate, caloric_density_local, caloric_density_mucus, volume_mucus, osmotic_value_brush,
                          osmotic_value_lumen, delta, pH,
                          liquidviscosityplate, hydrogen, surface_factor, sucrase_inhibitie, Maltase_inhibitie, Ia, Ic, glucose_absorption_rate,Mesentric_CCK, Mesentric_PYY,
                          Mesentric_GIP, Mesentric_GLP1, Mesentric_Ghrelin, Mesentric_Motilin, Mesentric_Secretin, eating_speed_dish , hypertrophy, salivary_amylase_activity, pancreatic_amylase_activity, pepsin_activity, pancreatin_activity : segment;
        dt, BasisGastricEmptyingRate, BasisSITransitRate, CCKEffectOnGastricEmptying, CCKgastricEmptying, OsmoticgastricEmptying, pHgastricEmptying, GLP1gastricemptying_stress,
                          ViscousgastricEmptying,fundus_viscosity, corpus_viscosity, antral_viscosity, Gastric_tone, segmentpH, volume_mucus_layer, bicarbonate_secretion, HbA1c,
                          HbA1c_halflife, Glycosylation_rate_constant : real;
        CCKdegradationconstant, PYYdegradationconstant, GIPdegradationconstant,GLP1degradationconstant, Ghrelindegradationconstant, Motilindegradationconstant, Secretindegradationconstant,
                          CCKattenuationSI, PYYattenuationSI, GLP1attenuationSI, GIPattenuationSI, GhrelinattenuationSI, MotilinattenuationSI, SecretinattenuationSI, LinguallipolysisCST, GastriclipolysisCST,
                          PancreaticlipolysisCST, starttime,  osmotic_value_brush_duodenum, osmotic_value_lumen_duodenum, osmotic_value_brush_distal, osmotic_value_lumen_distal,
                          osmcont, swallow_speed, total_linguallipaselipolysis_rate, total_gastriclipaselipolysis_rate, total_pancreaticlipaselipolysis_rate : real;
        CCKtrans, wateruitwisselsnelheid, hunger_SI_calories, hunger_SI_absorbable_calories, hunger_detectable_nutrients, hunger_intestinal_volume, totalcaloricabsorption_rate,
                          hunger_totalcaloricabsorptionrate, maxcaloricabsorptioncapacity: real;
        gastric_emptying, in_enterocytes, PYY, CCK, GIP, GLP1, Ghrelin, Motilin, Secretin, fullness, hunger, G2in_molpermin   : real;
        G0in,G1in,G2in, Gluc0, Gluc1,Gluc2,glycogene_store0,glycogene_store1,glycogene_store2, G0_increment, G1_increment, G2_increment,serum_glucose0_mg_dl, serum_glucose1_mg_dl,
                          serum_glucose2_mg_dl, serum_glucose0, serum_glucose1, serum_glucose2,serum_glucose2_eq, R_C13_C12:real;
        fat_in_colon,casein_in_colon,glucose_in_colon,protein_in_colon,total_solids_corpus, total_solids_antrum, viscosity_liquids_in_corpus, incretin, time_of_day, referencetime,
                          intercellulair_glucose_muscle_adipose_mg_dl, intercellulair_glucose_muscle_adipose, intercellulair_glucose_skin_mg_dl, intercellulair_glucose_skin,
                          total_glucose_absorption_rate, total_bile_secretion_rate, conjugated_bile_salt_secretion_rate, pancreatic_fluid_secretion, gastric_fluid_secretion,
                          total_conjugated_bile_salt_secretion,
                          total_bile_secreted, total_pancreatic_fluid_secretion, total_gastric_fluid_secretion, saliva_addition, total_saliva_secretion  :real;
        D, caloric_value, density, CCKresponsefactor, PYYresponsefactor, GIPresponsefactor, GLP1responsefactor, Ghrelinresponsefactor, Motilinresponsefactor, Secretinresponsefactor,
                          proteolytic_susceptibility, absorbed, vloksamenstelling, solids_corpus, solids_antrum, solids_fundus, samenstelling_syneresevloeistof,
                          totalabsorption_rate, bloodserumconc, inverse_serumhalftime, fundus, gall_bladder, Blood, calories_compound, gastric_sieve, solid, concentration_fundus : component;
        realGLP1, RealGIP,prefactor, totalwaterabsorptionrate, gastric_emptying_water, gastric_emptying_pharmaceutical, percentual_loss_of_triglycerides_lingual, percentual_loss_of_triglycerides_gastric, percentual_loss_of_triglycerides_pancreatic  : real;
        , pancreatic_fluid_secretion_rate :real;
        MMC_status : integer;
        MMC_time : MMC_timer;
        YY,MM,DD, YYE, MME, DDE: word;
        out : text;
        G, Gluc, Fruc, f1:real;   {Gluc, Fruc = amount of glucose, Fructose in mg in plasma and intercellular space}{G = concentration of glucose in mg/dl}
            M, MAX:real;
        govern: string[30];governnumber, comp: integer;
        concentration_compound_fundus, compound_fundus : component;

implementation

end.

