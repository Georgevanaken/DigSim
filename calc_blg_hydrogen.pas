unit calc_BLG_hydrogen;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables;

function BLG_hydrogen(pH, BLG:real):real;

implementation
Function BLG_hydrogen(pH, BLG:real):real;  {gebaseerd op Tanford and Nozaki, Physico-chemical comparison of b-lactoglobulin s A and B, J. Biological Chem., 1959, 234(11), 2874-2877)}
var correction_to_IEP:real;
begin
  correction_to_IEP := 0;   {IEP is 000 }
  case trunc(pH) of
    1: BLG_hydrogen := 0.002173913 - (pH-trunc(pH))*9.23913E-05;
    2: BLG_hydrogen := 0.002081522 - (pH-trunc(pH))*0.000383152;
    3: BLG_hydrogen := 0.00169837 - (pH-trunc(pH))*0.001263587;
    4: BLG_hydrogen := 0.000434783 - (pH-trunc(pH))*0.000543478;
    5: BLG_hydrogen := -0.000108696 - (pH-trunc(pH))*0.000271739;
    6: BLG_hydrogen := -0.000380435 - (pH-trunc(pH))*0.00026087;
    7: BLG_hydrogen := -0.000641304 - (pH-trunc(pH))*0.000228261;
    8: BLG_hydrogen := -0.000869565 - (pH-trunc(pH))*0.000108696;
    9: BLG_hydrogen := -0.000978261 - (pH-trunc(pH))*5.43478E-05;
    10: BLG_hydrogen := -0.1032609 - (pH-trunc(pH))*0;
    11: BLG_hydrogen := -0.1032609 - (pH-trunc(pH))*0;
    12: BLG_hydrogen := -0.1032609 - (pH-trunc(pH))*0;
  end;
  BLG_hydrogen:= BLG*(BLG_hydrogen- correction_to_IEP);
  {dit is het aantal mol waterstof ionen dat nodig is om de pH te laten afwijken van de IEP}
end;
end.

