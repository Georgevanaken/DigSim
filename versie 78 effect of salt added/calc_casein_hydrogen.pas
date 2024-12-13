unit calc_Casein_hydrogen;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables;

function casein_hydrogen(pH, casein:real):real;

implementation

Function casein_hydrogen(pH, casein:real):real;  {gebaseerd op PERTZOFF and CARPENTER, 1932 (zuur en base addities) en Erdos "the effect of salts on the isoelectric point of casein"}
var correction_to_IEP, bound_hydrogen:real;
begin
  correction_to_IEP := -0.00002;   {IEP op 4.5 gesteld; dit getal corrigeert de titratiecurve van PERTZOFF and CARPENTER, 1932}
  if pH < 2 then bound_hydrogen := 0.0016 + (2-pH)*(0.002-0.0016) else
  if pH < 2.5 then bound_hydrogen := 0.0009 + (2.5-pH)*(0.0016-0.0009)*2 else
  if pH < 3 then bound_hydrogen := 0.0004 + (3-pH)*(0.0009-0.0004)*2 else
  if pH < 3.5 then bound_hydrogen := 0.0001 + (3.5-pH)*(0.0004-0.0001)*2 else
  if pH < 4 then bound_hydrogen := 0.00005 + (4-pH)*(0.0001-0.00005)*2 else
  if pH < 4.5 then bound_hydrogen := 0.00002 + (4.5-pH)*(0.00005-0.00002)*2 else
  if pH < 5 then bound_hydrogen := 0 + (5-pH)*(0.00002-0)*2 else
  if pH < 6 then bound_hydrogen := -0.00002 + (6-pH)*(0+ 0.00002) else
  if pH < 7 then bound_hydrogen := -0.0003 + (7-pH)*(-0.00002+ 0.0003) else
  if pH < 8 then bound_hydrogen := -0.0007 + (8-pH)*(-0.0003+ 0.0007) else
  if pH < 9 then bound_hydrogen := -0.0008 + (9-pH)*(-0.0007+ 0.0008) else
                 bound_hydrogen := -0.0011 + (10-pH)*(-0.0008+ 0.0011) ;
  casein_hydrogen := casein*(bound_hydrogen + correction_to_IEP);
  {dit is het aantal mol watersof ionen dat nodig is om de pH te laten afwijken van de IEP}
end;
end.

