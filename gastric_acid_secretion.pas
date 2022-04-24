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
      basis := 0.1{0.1}{2}{5} ; {basis secretion rate}
      if MMC_status = 3 then basis := 1; {the idea behind this is that MMC3 is the cleaning phase and extra gastric juice is secreted. Exactly how much extra remains to be investigated}
      richt_pH := 1;  {the idea is that under standard conditions (Konturec 1974) acid addition is proportional to the difference between the anthral pH and 1, in reality also peptides (through gastrin) and antral distension stimulate acid release, but this should still be modelled}
      mole_HCl := 0.1*(10/15000)*(pH[4]-richt_pH)/(5.5-richt_pH);   {according to Konturec 1974, this wil be controlled by gastrin. This gastrin control still has to be modelled}
      if mole_HCl <0 then mole_HCl := 0;
      If pH[3] < 4 then mole_HCl := mole_HCl*(exp(pH[4]-1)-1{2})/exp(pH[4]-1{2}); {temporary brake on acid production, reacting to a too low pH in the stomach or duodenum, as actually controlled through secretin release by the duodenum}
      {if mole_HCl > 10/15000 then mole_HCL:= 0.1/15000;}
      {secretion takes place in the anthrum and corpus, hence compartiment 3}
      {******}

      gastric_fluid_secretion :=  (mole_HCl*10000{20000} + basis);    {0.1 of 0.05 M HCl +basis}
      gastric_fluid_secretion := gastric_fluid_secretion + Volume[4]/(Volume[4]+5); {effect of antrum distension}
      gastric_fluid_secretion := gastric_fluid_secretion + liquid[4,62]/(liquid[4,62]+1E-5) + liquid[4,64]/(liquid[4,64]+1E-5) + liquid[4,65]/(liquid[4,65]+1E-5);  {peptiden formed by pepsin}

      liquid[3,1] := liquid[3,1] +  gastric_fluid_secretion*dt; {assuming 0.1 M HCl in gastric secretion + 10ml/min neutral basis secretion} {pH is however 1.5-3.5, whcih would correspond to 0.01M HCl} { daily about 1.5 l}
      liquid[3,80] := liquid[3,80] + mole_HCl*dt; {addition gastric secrete to corpus}
      liquid[3,49] := liquid[3,49] + 0.001*gastric_fluid_secretion*dt; {estimate for gastric lipase}
            mucin_secretion := 0.01{0.02}*gastric_fluid_secretion*dt; {guess: 2% mucus protein in gastric secretion}
            mucin_secretion := mucin_secretion + 0.01*dt; {guess value for continuous secretion of mucin by gastic wall}
      liquid[3,55] := liquid[3,55] + mucin_secretion;
      liquid[3,80] := liquid[3,80] + Mucin_hydrogen(pH_mucus,mucin_secretion); {correction to IEP of mucus}
      liquid[3,40] := liquid[3,40] + 0.005*gastric_fluid_secretion*dt; {guess: 0.5% pepsine in gastric secretion}



      total_gastric_fluid_secretion := total_gastric_fluid_secretion + gastriC_fluid_secretion*dt;
end;

end.

