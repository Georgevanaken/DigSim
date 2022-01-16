unit calc_BSA_hydrogen;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables;

function BSA_hydrogen(pH, BSA:real):real;

implementation
Function BSA_hydrogen(pH, BSA:real):real;  {gebaseerd op Curvale, 2009 (Buffer capacity of bovine serum albumin)}
var correction_to_IEP:real;
begin
  correction_to_IEP := 0.0267;   {IEP is 4.7 }
  case trunc(pH) of
    1: BSA_hydrogen := 0.14 - (pH-trunc(pH))*0.01;
    2: BSA_hydrogen := 0.13 - (pH-trunc(pH))*0.05;
    3: BSA_hydrogen := 0.08 - (pH-trunc(pH))*0.037;
    4: BSA_hydrogen := 0.043 - (pH-trunc(pH))*0.018;
    5: BSA_hydrogen := 0.025 - (pH-trunc(pH))*0.02;
    6: BSA_hydrogen := 0.005 - (pH-trunc(pH))*0.015;
    7: BSA_hydrogen := -0.01 - (pH-trunc(pH))*0.01;
    8: BSA_hydrogen := -0.02 - (pH-trunc(pH))*0.005;
    9: BSA_hydrogen := -0.025 - (pH-trunc(pH))*0.01;
    10: BSA_hydrogen := -0.035 - (pH-trunc(pH))*0.06;
    11: BSA_hydrogen := -0.095 - (pH-trunc(pH))*0.105;
    12: BSA_hydrogen := -0.20 -  (pH-trunc(pH))*1;
  end;
  BSA_hydrogen := BSA*(BSA_hydrogen- correction_to_IEP)/(15E-4*66463);
  {dit is het aantal mol watersof ionen dat nodig is om de pH te laten afwijken van de IEP}
end;
end.

