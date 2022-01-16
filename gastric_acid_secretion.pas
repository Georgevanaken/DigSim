unit gastric_acid_secretion ;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables, calc_mucin_hydrogen;

procedure gastric_secretion;



implementation


Procedure gastric_secretion;
var  mole_HCl, richt_pH, basis, mucin_secretion, pH_mucus, gastricacidproduction_rate   : real;
begin
      pH_mucus := 2;
      basis := 0.1{0.1}{2}{5} ; {basis secretie snelheid}
      if MMC_status = 3 then basis := 1; {het idee hierachter is dat fase 3 de schoonmaak fase is en er dan extra maagsap afgescheiden wordt}
      richt_pH := 1;  {het idee is dat onder standaard condities (Konturec 1974) de zuuradditie proportioneel is met het verschil tussen de antrum pH en 1, in werkelijkheid stimuleren ook peptiden (via gastrin) en antral distension de release}
      mole_HCl := 0.1*(10/15000)*(pH[4]-richt_pH)/(5.5-richt_pH);   {volgens Konturec 1974, zal fysiologisch door gastrine gestuurd zijn}
      if mole_HCl <0 then mole_HCl := 0;
      If pH[3] < 4 then mole_HCl := mole_HCl*(exp(pH[4]-1)-1{2})/exp(pH[4]-1{2}); {voorlopige rem reagerend op een te lage pH, beschrijft het effect van secretine}
      {if mole_HCl > 10/15000 then mole_HCL:= 0.1/15000;}
      {secretie vindt plaats in antrum en corpus, dus compartiment 3}
      {******}

      gastric_fluid_secretion :=  (mole_HCl*10000{20000} + basis);    {0.1 of 0.05 M HCl +basis}
      gastric_fluid_secretion := gastric_fluid_secretion + Volume[4]/(Volume[4]+5); {effect distensie van antrum}
      gastric_fluid_secretion := gastric_fluid_secretion + liquid[4,62]/(liquid[4,62]+1E-5) + liquid[4,64]/(liquid[4,64]+1E-5) + liquid[4,65]/(liquid[4,65]+1E-5);  {peptiden gevormd foor pepsine}

      liquid[3,1] := liquid[3,1] +  gastric_fluid_secretion*dt; {uitgaande van 0.1 M HCl in maagsap + 10ml/min zuurloos basis} {pH is echter 1.5-3.5, stel 0.01M} { dagelijks ongeveer 1.5 l}
      liquid[3,80] := liquid[3,80] + mole_HCl*dt; {additie maagsap aan de corpus}
      liquid[3,49] := liquid[3,49] + 0.001*gastric_fluid_secretion*dt; {gokwaarde gastric lipase}
            mucin_secretion := 0.01{0.02}*gastric_fluid_secretion*dt; {gok: 2% mucus eiwit in gastric secretion}
            mucin_secretion := mucin_secretion + 0.01*dt; {gokwaarde voor continue afscheiding mucine door maagwand}
      liquid[3,55] := liquid[3,55] + mucin_secretion;
      liquid[3,80] := liquid[3,80] + Mucin_hydrogen(pH_mucus,mucin_secretion); {correctie naar IEP van mucus}
      liquid[3,40] := liquid[3,40] + 0.005*gastric_fluid_secretion*dt; {gok: 0.5% pepsine in gastric secretion}



      total_gastric_fluid_secretion := total_gastric_fluid_secretion + gastriC_fluid_secretion*dt;
end;

end.

