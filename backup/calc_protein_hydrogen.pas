unit calc_Protein_hydrogen;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables, calc_BSA_hydrogen, calc_Casein_hydrogen, calc_Mucin_hydrogen, calc_BLG_hydrogen;

function Protein_hydrogen(pH, Protein:real):real;
function Ber_hydrogen(var pH, casein, BSA, Mucin, BLG, Protein: real):real;
Function  ber_pH(var Volume, hydrogen, casein, BSA, Mucin, BLG, Protein: real):double;
Procedure ber_alle_pHs;


implementation
Function Protein_hydrogen(pH, Protein:real):real;  {als BSA; gebaseerd op Curvale, 2009 (Buffer capacity of bovine serum albumin)}  {het gaat hier kennelijk om het aantal molen H= die je aan heet aangegeven aantal gram eiwit moet toevoegen}
var correction_to_IEP:real;
begin
  correction_to_IEP := 0.0267;   {IEP is 4.7 }
  case trunc(pH) of
    1: Protein_hydrogen := 0.14 - (pH-trunc(pH))*0.01;
    2: Protein_hydrogen := 0.13 - (pH-trunc(pH))*0.05;
    3: Protein_hydrogen := 0.08 - (pH-trunc(pH))*0.037;
    4: Protein_hydrogen := 0.043 - (pH-trunc(pH))*0.018;
    5: Protein_hydrogen := 0.025 - (pH-trunc(pH))*0.02;
    6: Protein_hydrogen := 0.005 - (pH-trunc(pH))*0.015;
    7: Protein_hydrogen := -0.01 - (pH-trunc(pH))*0.01;
    8: Protein_hydrogen := -0.02 - (pH-trunc(pH))*0.005;
    9: Protein_hydrogen := -0.025 - (pH-trunc(pH))*0.01;
    10: Protein_hydrogen := -0.035 - (pH-trunc(pH))*0.06;
    11: Protein_hydrogen := -0.095 - (pH-trunc(pH))*0.105;
    12: Protein_hydrogen := -0.20 -  (pH-trunc(pH))*1;
  end;
  Protein_hydrogen := Protein*(Protein_hydrogen- correction_to_IEP)/(15E-4*66463)  {-0.001*exp(-2.303*pH)+0.001*exp(-2.303*4.7)  {deze term zou afhankelijk moeten zijn van de concentratie van eiwit concentratie in de titratiecurve} ;}
  {dit is het aantal mol waterstof ionen dat nodig is om de pH te laten afwijken van de IEP van het eiwit, dus het aantal molen dat gebonden is aan het eiwit}
end;

function Ber_hydrogen(var pH, casein, BSA, Mucin, BLG, Protein: real):real;   {berekend het aantal molen H+ dat bij de gegeven pH gebonden is aan het ewit; de hoeveelheid waterstof ionen ten opzichte van de isoelectrische punten van de aanwezige eiwitten}
begin
  Ber_hydrogen := casein_hydrogen(pH, casein)+ BSA_hydrogen(pH, BSA)+ Mucin_hydrogen(pH, Mucin)+BLG_hydrogen(pH, BLG)+Protein_hydrogen(pH, Protein);
end;

Function  ber_pH(var Volume, hydrogen, casein, BSA, Mucin, BLG, Protein: real):double;  {berekend de pH voor een gegeven hoeveelheid molen watertof ionen gebonden aan eiwit in een eiwit mengsel}
var pH,Min, Max, Guess_hydrogen :real;
const Tol  = 0.001;
begin
    MIN := 0;
    Max := 10+Tol;
      while Max-Min > Tol do
      begin
       pH:= (MAX+MIN)/2;
       Guess_hydrogen := Ber_hydrogen (pH, casein, BSA, Mucin, BLG, Protein);
       If Guess_hydrogen> hydrogen - Volume*0.001*(exp(-2.302585093*pH)-0.0000001) then MIN := pH else MAX:= pH;

      end;
    ber_pH:= (MAX+MIN)/2;
end;

Procedure ber_alle_pHs;
var hulp1, hulp2, hulp3, hulp4,hulp5,hulp6: real;
begin
  for i := 2 to 15 do
   begin
    hulp1 := liquid[i,61]+liquid[i,62]{+liquid[i,63]+liquid[i,64]+liquid[i,65]}+liquid[i,63]+liquid[i,64]+liquid[i,65];
    pH[i]:=ber_pH(Volume[i],liquid[i,80], hulp1 ,liquid[i,66],liquid[i,55], liquid[i,67], liquid[i,4]);
   end;
   hulp1 := fundus[61]+fundus[62]{+fundus[63]+fundus[64]+fundus[65]}{+solids_fundus[61]+solids_fundus[62]}{+solids_fundus[63]+solids_fundus[64]+solids_fundus[65]}+fundus[63]+fundus[64]+fundus[65];
   hulp2 := fundus[80]+solids_fundus[80];
   hulp3 := fundus[66]+solids_fundus[66];
   hulp4 := fundus[55]+solids_fundus[55];
   hulp5 := fundus[67]+solids_fundus[67];
   hulp6 := fundus[4]+solids_fundus[4];
   pHfundus:=ber_pH(Vfundus, hulp2, hulp1 ,hulp3,hulp4, hulp5, hulp6);
   {writeln(pH[1]:10, pH[2]:10);}
end;
end.

