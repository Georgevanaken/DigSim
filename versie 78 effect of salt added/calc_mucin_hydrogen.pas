unit calc_Mucin_hydrogen;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables;

function Mucin_hydrogen(pH, Mucin:real):real;

implementation

Function Mucin_hydrogen(pH, Mucin:real):real;  {gebaseerd op Leon Li, Oliver Lieleg, Katharaina Ribbeck, and Jongyoon Han, Microfluidic in vitro model for quantitative study of stomachmucinacid barrier function, 15th International Conference on Miniaturized Systems for Chemistry and Life Sciences, October 2-6, 2011, Seattle, Washington, USA, page 771, Figure C)}
var correction_to_IEP:real;
begin
  correction_to_IEP := 0.0005; {IEP is between 2 and 3 }
  case trunc(pH) of
    1: Mucin_hydrogen := 0.001333333 - (pH-trunc(pH))*0.000566667;
    2: Mucin_hydrogen := 0.000766667 - (pH-trunc(pH))*0.000466667;
    3: Mucin_hydrogen := 0.0003 - (pH-trunc(pH))*0.00016;
    4: Mucin_hydrogen := 0.00014 - (pH-trunc(pH))*0.00012667;
    5: Mucin_hydrogen := 0.00001333 - (pH-trunc(pH))*4.663E-05;
    6: Mucin_hydrogen := -0.0000333 - (pH-trunc(pH))*0.00003;
    7: Mucin_hydrogen := -0.0000633 - (pH-trunc(pH))*0.00001;
    8: Mucin_hydrogen := -0.0000733 - (pH-trunc(pH))*0.000030033;
    9: Mucin_hydrogen := -0.000103333 - (pH-trunc(pH))*4.3334-05;
    10: Mucin_hydrogen := -0.000146667 - (pH-trunc(pH))*0.000186666;
    11: Mucin_hydrogen := -0.000333333 - (pH-trunc(pH))*0.000333333;
    12: Mucin_hydrogen := -0.001 - (pH-trunc(pH))*0.000666667;
  end;
  Mucin_hydrogen := Mucin*(Mucin_hydrogen- correction_to_IEP);
  {dit is het aantal mol waterstof ionen dat nodig is om de pH te laten afwijken van de IEP}
end;
end.

