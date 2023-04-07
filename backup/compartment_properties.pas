unit Compartment_properties;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables, diet_from_file;

procedure volumina;

implementation
var dishnumber:integer;               {dishnumer relates to the dishes or glasses (currently 5)from which the meal can be eaten}
procedure volumina;         {this procedure calculates the volumes of each compartment}
begin
 for dishnumber := 1 to 5 do
 begin
       volumedish[dishnumber] := 1E-10;  {to avoid division by 0 in the calculation of concentrations}
       for j:= 1 to comp do
        begin
        volumedish[dishnumber] := volumedish[dishnumber] + liquiddish[dishnumber,j]/density[j];
        end;
       for compoundnumber:= 1 to number_of_compounds do
        begin
        volumedish[dishnumber] := volumedish[dishnumber] + compounddish[dishnumber,compoundnumber]/1;  {density compounds nog yet introduced (currently set to 1 g/ml)}
        end;

      if Volumedish[dishnumber]>1E-20 then
          begin
            for j:= 1 to comp do concentration_dish[dishnumber,j]:= liquiddish[dishnumber,j]/volumedish[dishnumber];    {only dissolved parts}
            for compoundnumber := 1 to number_of_compounds do concentration_compound_dish[dishnumber,compoundnumber]:= compounddish[dishnumber,compoundnumber]/volumedish[dishnumber];
          end
        else
          begin
             Volumedish[dishnumber]:= 1E-20;
             for j:= 1 to comp do concentration_dish[dishnumber,j]:= 0;
             for compoundnumber := 1 to number_of_compounds do concentration_compound_dish[dishnumber,compoundnumber]:= 0;
          end;

 end;



 for i:= 2 to 16 do
  begin
       volume[i] := 1E-5{10};  {to avoid division by 0 in the calculation of concentrations}
       for j:= 1 to comp do
        begin
        volume[i] := volume[i] + liquid[i,j]/density[j];
        {for s:= 1 to smax do Volume[i] := volume[i] + solids[s,i,j]/density[j];}
        end;
       for compoundnumber:= 1 to number_of_compounds do
        begin
        volume[i] := volume[i] + compound[i,compoundnumber]/compounddensity[compoundnumber];
        end;

      if Volume[i]>1E-5{10} then
          begin
            for j:= 1 to comp do concentration_lumen[i,j]:= liquid[i,j]/volume[i];    {only dissolved parts} { if concentration_lumen[i,j]>1 then writeln('!'); }
            for compoundnumber := 1 to number_of_compounds do concentration_compound[i,compoundnumber]:= compound[i,compoundnumber]/volume[i];
          end
        else
          begin
             Volume[j] := 1E-5{10};
             for j:= 1 to comp do concentration_lumen[i,j]:= 1E-10;
             for compoundnumber := 1 to number_of_compounds do concentration_compound[i,compoundnumber]:= 1E-10;
          end;
  end;
  in_enterocytes := 0;
  volumetotal:= 0;
  Vfundus:= 0;
  for j:= 1 to comp do Vfundus := Vfundus + fundus[j]/density[j];
  for compoundnumber:= 1 to number_of_compounds do Vfundus := Vfundus + compoundfundus[compoundnumber]/1;  {density compounds nog yet introduced (currently set to 1 g/ml)}
  if Vfundus>1E-5{10} then
          begin
            for j:= 1 to comp do concentration_fundus[j]:= fundus[j]/Vfundus;    {only dissolved parts}
            for compoundnumber := 1 to number_of_compounds do concentration_compound_fundus[compoundnumber]:= compound_fundus[compoundnumber]/Vfundus;
          end
        else
          begin
             {Volume[j]}Vfundus:= 1E-5{20};
             for j:= 1 to comp do concentration_fundus[j]:= 0;
             for compoundnumber := 1 to number_of_compounds do concentration_compound_fundus[compoundnumber]:= 0;
          end;

  FOR i:=5 to 15 do volumetotal := volumetotal + volume[i];
end;
end.

