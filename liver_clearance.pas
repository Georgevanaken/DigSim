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
   inverse_serumhalftime[5]:= 1/60{30}; {inverse van decaytijd in minuten}
 end;

procedure decay_in_blood;
  begin
    {general}
     for j:= 1 to 80 do bloodserumconc[j]:= bloodserumconc[j] + (1000*totalabsorption_rate[j]/5-bloodserumconc[j]*inverse_serumhalftime[j])*dt; {serum conc in mg/l}
  end;

Procedure fructose_to_glucose;
var  decay_fructose:real;
  begin
    decay_fructose := 0.1*Fruc*dt;
    Fruc := Fruc- decay_fructose; {aanname is dat fructose in 2 minuut door de lever uit het bloed verwijderd wordt, fruc in mg}
    Gluc:=Gluc + 1*decay_fructose; {aanname is dat een deel van de fructose door de lever opgezet wordt in glucose; gluc in mg}
  end;

end.

