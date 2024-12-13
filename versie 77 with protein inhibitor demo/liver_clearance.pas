unit Liver_clearance;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables;

procedure decay_in_blood;
procedure Initialize_blood_serum_conc;
procedure fructose_to_glucose;

implementation
procedure Initialize_blood_serum_conc;
 begin
   for j:= 1 to 80 do
     begin
       bloodserumconc[j] := 0;
       inverse_serumhalftime[j] := 0;
     end;
   inverse_serumhalftime[5]:= 1/60{30}; {inverse of the decaytime in minutes}
 end;

procedure decay_in_blood;
  begin
    {general}
     for j:= 1 to 80 do bloodserumconc[j]:= bloodserumconc[j] + (1000*totalabsorption_rate[j]/5-bloodserumconc[j]*inverse_serumhalftime[j])*dt; {serum conc in mg/l}
  end;

Procedure fructose_to_glucose;
var  decay_fructose:real;
  begin
    decay_fructose := 0.5*Fruc*dt;
    Fruc := Fruc- decay_fructose; {assumption is that fructose is removed from the blood by the liver in 2 minutes, fruc in mg}
    Gluc:=Gluc + 1*decay_fructose; {assumption is that a part of the fructose is transformed by the liver into glucose; gluc in mg}
  end;

end.

