unit calc_citric_acid_hydrogen;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables;

function citric_acid_hydrogen(pH, citric_acid:real):real;

implementation
Function citric_acid_hydrogen(pH, citric_acid:real):real;  {gebaseerd op Beltry 2019, considerations on the derivative of the titration curve of weak diacids}
var correction_to_IEP:real;
begin
  correction_to_IEP := 0;   {IEP is 000 }
  case trunc(pH) of                                                       {first on a mole-mole basis}
    1: citric_acid_hydrogen := 3 - (pH-trunc(pH))*0;
    2: citric_acid_hydrogen := 3 - (pH-trunc(pH))*0.5;
    3: citric_acid_hydrogen := 2.5 - (pH-trunc(pH))*0.5;
    4: citric_acid_hydrogen := 2 - (pH-trunc(pH))*0.6;
    5: citric_acid_hydrogen := 1.4 - (pH-trunc(pH))*0.6;
    6: citric_acid_hydrogen := 0.8 - (pH-trunc(pH))*0.5;
    7: citric_acid_hydrogen := 0.3 - (pH-trunc(pH))*0.3;
    8: citric_acid_hydrogen := 0 - (pH-trunc(pH))*0;
    9: citric_acid_hydrogen := 0 - (pH-trunc(pH))*0;
    10: citric_acid_hydrogen := 0 - (pH-trunc(pH))*0;
    11: citric_acid_hydrogen := 0 - (pH-trunc(pH))*0;
    12: citric_acid_hydrogen := 0 - (pH-trunc(pH))*0;
  end;
  citric_acid_hydrogen:= citric_acid*(-citric_acid_hydrogen- correction_to_IEP)/192.124 {= molecular mass};
  {dit is het aantal mol waterstof ionen dat nodig is om de pH te laten afwijken van de IEP}
end;
end.

