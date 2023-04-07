unit gastric_motility;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables, diet_from_file,calc_Protein_hydrogen, compartment_properties;

Procedure pressures;
procedure solids_in_stomach;
procedure gastric_compartment_viscosities;
procedure gastric_redistribution;

implementation
var targetcorpuscomp: component;

Procedure pressures;
begin
     {Gastric_tone :=  0.023*(1-Transition)*(Volume[3] + Volume[4])+0.5; }{maagdruk, mbar} {Tack}
     Gastric_Volume := Vfunduscorpus + Volume[4];
     Gastric_Tone := 3+ 0.0177*(Gastric_Volume + total_solids_corpus)*(1-0.5*Transition); {gastric tone in mBar} {Janssen}
     if Gastric_tone < 3 then Gastric_tone:=3;
     {Gastric compliance = 58 +/- 6 mL/mm Hg (Tack 1998); first detection and discomfort starts at 4.1 +/- 0.5 and 12.1 +/- 0.4 mm Hg, respectively}
     {Janssen 2011: 3mBar + saline: 0.0177 mBar/ml; nutrient solution: 0.00887 mBar/ml}
end;

procedure solids_in_stomach;
begin
     solids_in_antrum := 0;
     solids_in_corpus:= 0;
     solids_in_fundus:= 0;
     for compoundnumber:= 1 to number_of_compounds do
       begin
         if (compoundviscosity[compoundnumber]>100*viscosity_liquids_in_corpus) and (solid[compoundnumber] = 0) then solid[compoundnumber]:= 10;       {highly viscous liquids a emulsified to large droplets in chyme}
           if solid[compoundnumber] <> 0 then
           begin
           solids_in_antrum:= solids_in_antrum+ compound[4,compoundnumber] ;
           solids_in_corpus:= solids_in_corpus+ compound[3,compoundnumber] ;
           solids_in_fundus:= solids_in_fundus+ compoundfundus[compoundnumber];
           end;
       end;
     total_solids_in_stomach := solids_in_antrum+solids_in_corpus+solids_in_fundus;
end;

procedure gastric_compartment_viscosities;
var  phi, solids_in_corpus:real;
const phimax =0.55;
      eta = 2.5;
begin
  phi := solids_in_corpus/Volume[3];
  viscosity_liquids_in_corpus := 5;      {correct for emulsions and polymers in solution}
  for compoundnumber := 1 to number_of_compounds do if compound[3,compoundnumber]>1E-40 then viscosity_liquids_in_corpus := viscosity_liquids_in_corpus + compoundviscosity[compoundnumber]*(1/(1+Volume[3]*compounddensity[compoundnumber]/compound[3,compoundnumber]))  ;
  if phi>0.5 then
        begin
          corpus_viscosity := viscosity_liquids_in_corpus;
          for compoundnumber:= 1 to number_of_compounds do
             corpus_viscosity := corpus_viscosity + solid[compoundnumber]*compound[3,compoundnumber]*compoundviscosity[compoundnumber];
          corpus_viscosity := corpus_viscosity/solids_in_corpus;
        end

     else corpus_viscosity := viscosity_liquids_in_corpus*exp(-eta*phimax*LN(1-phi/phimax));   {Krieger-Dougherty equation}
end;

procedure gastric_redistribution;
var liquidscorpusfundus, targetcorpus,solidscorpusfundus: component;
   x, fraction_retropelled, retropulsionvolume, vlokkingssnelheid_corpus, vlokkingssnelheid_antrum, oplossnelheid_vlokken_corpus,oplossnelheid_vlokken_antrum,
        syneresesnelheid, verdeler, mengfactor, mengsnelheid, pHvlokken, hulp, bolus_oppervlak, kineticfactor, change, discost, target_Vcorpus, Corpus,
        up_corpus_to_fundus, down_fundus_to_corpus, target_antralvolume, antral, transfer_efficiency_fundus_corpus: real;
{In the model the stomach is described as a lign-up of the corpus and anthrum compartments between the oesophagus and duodenum. In addition, coupled to the corpus is the fundus
compartment, which operates as an expandable storage vessel for extra liquid and solid volumes. These three compartments can exchange their contents, by which the the volumes of these compartments
adjust to the operating pressure within these compartments. The operating pressure equals the gastric pressure or tone.
Gastric pressure is determined by the volume of product of the fundic volume and a proportionality constant that varies with the CCK bloodlevel. In the fed state, CCK is increased, whcih reduces the
poportionality constant and hence the gastric tone. Due to this the gastric volume (at constant gastric tone) adjusts to the volume of the food and the detected nutrient content in the proximal small intestin,
as messaged by blood CCK.
A cyclic peristaltic movement (set at 3 times per minute) from fundus to corpus to antrum pushes the gastric content toward or into the antrum. The speed of this transport is set to be dependent of the viscosity of the material in the emptying compartment.
After each cycle, the antrum completely contracts and repels the antral content partially back into the corpus and partially into the duodenum. The ratio of retropulsion and emptying into the duedenum corresponds to the regulation of gastric emptying.
This ration is determined by a number of factors, including the antral pressure, viscosity of the material in the antrum, and feedback regulational signals from the duodenum, which monitor duodenal pH, osmotic pressure, the presence of CCK- and GIP-releasing food materials,
each of which contribute to reduce the the passage through the pyloric valve in an additive way}
begin
  volumina;
  antral := 0.1{0};
cyclus := cyclus+ dt;
if cyclus > 0.3333 then begin cyclus := 0; tel := 0; up[4]:=0 end;  {Gastric-peristaltic cycle set to 3 times per minute}
if cyclus< 0.1111 then
    begin

  {1. persistalsis from fundus to corpus}
      Vfunduscorpus := Vfundus+Volume[3];  {corresponds to the combined volume of liquid nd solids in fundus en corpus}
      transfer_efficiency_fundus_corpus := 1;
      down_fundus_to_corpus:= transfer_efficiency_fundus_corpus*Vfundus/0.1111;

      begin
        for j:= 1 to comp do
                  begin
                       liquid[3,j]:= liquid[3,j]+down_fundus_to_corpus*dt*fundus[j]/Vfundus;
                       fundus[j] := fundus[j] - down_fundus_to_corpus*dt*fundus[j]/Vfundus;
                  end;
        for compoundnumber:= 1 to number_of_compounds do
                  begin
                  compound[3,compoundnumber]:= compound[3,compoundnumber]+down_fundus_to_corpus*dt*compoundfundus[compoundnumber]/Vfundus;
                  compoundfundus[compoundnumber] := compoundfundus[compoundnumber] - down_fundus_to_corpus*dt*compoundfundus[compoundnumber]/Vfundus;
                  end;
        volumina;
      end;
     end
else
 if cyclus < 0.2222 then {retropulsion from corpus to fundus and flow from corpus to antrum}
     begin
     Corpus := Volume[3];
     up_corpus_to_fundus :={10*(Volume[3]/target_Vcorpus - 1)}1*Volume[3]/0.1111;

         begin
              for j:= 1 to comp do
                  begin
                  change :=  liquid[3,j]*dt*up_corpus_to_fundus/Volume[3];
                  liquid[3,j]:= liquid[3,j]-change;
                  fundus[j] := fundus[j]+change;
                  end;

              for compoundnumber:= 1 to number_of_compounds do
                  begin
                  change :=  compound[3,compoundnumber]*dt*up_corpus_to_fundus/Volume[3];
                  compound[3,compoundnumber]:= compound[3,compoundnumber]-change;
                  compoundfundus[compoundnumber] := compoundfundus[compoundnumber]+ change;
                  end;
              volumina;
         end ;

    {3. Transport of corpus to antrum by peristalsis}
      pressures;
      gastric_compartment_viscosities;
      maxantralvolume := 60{240}{240}{120};
      target_antralvolume :=  5+maxantralvolume*gastric_tone/(1{5}+gastric_tone); {target antral volume}
      up[3] := Corpus*((1-Volume[4]/target_antralvolume)/0.1111)*0.1;   {it has been demonstrated that this transport is slowed by the viscosity: Benini 1995, but take care: this relates to the viscosity in the gastric lumen, after dilution and interaction with gastric juice!}
         if MMC_status = 3 then Up[3] := 2*up[3];    {MMC = Migrating Motor Complex}

        begin
          for j:= 1 to comp do
              begin
              change := up[3]*liquid[3,j]*dt/volume[3];
              liquid[3,j] := liquid[3,j] - change;
              liquid[4,j] := liquid[4,j] + change;
              end;
          for compoundnumber := 1 to number_of_compounds do
              begin
              change := up[3]*dt*compound[3,compoundnumber]/volume[3];
              compound[3,compoundnumber] := compound[3,compoundnumber] - change;
              compound[4,compoundnumber] := compound[4,compoundnumber] + change;
              end;
          Volumina;
        end
    end
else
{if cyclus < 0.3333 then}  {emptying and retropulsion from antrum to corpus}
    begin

      {emptying}
    isotonic := 0.05; {0.05 glucose/ml}
    osmcont:=2{2};

      CCKgastricEmptying_stress :=  2 *totalrateCCKrelease + 1*CCK; {CCKEffectOnGastricEmptying*CCK}
                                if MMC_Status = 2 then CCKgastricEmptying_stress := 0.5*CCKgastricEmptying_stress; {situation during MMC phase 2}
                                if MMC_Status = 3 then CCKgastricEmptying_stress := 0*CCKgastricEmptying_stress; {situation during MMC phase 3}
      OsmoticgastricEmptying_stress := osmcont*osmotic_value_brush_duodenum/isotonic; {osmotic effect}
                                if MMC_Status = 2 then OsmoticgastricEmptying_stress := 0.5*OsmoticgastricEmptying_stress; {situation during MMC phase 2}
                                if MMC_Status = 3 then OsmoticgastricEmptying_stress := 0*OsmoticgastricEmptying_stress; {situation during MMC phase 3}
      pHgastricEmptying_stress := 5*(abs(6-pH[5])/6);     {pH effect}  {prefactor 20}
                                if MMC_Status = 2 then pHgastricEmptying_stress := 0.5*pHgastricEmptying_stress; {situation during MMC phase 2}
                                if MMC_Status = 3 then pHgastricEmptying_stress := 0*pHgastricEmptying_stress; {situation during MMC phase 3}
      GLP1gastricemptying_stress := 10*GLP1;    {400 gives the same geeft inhibition at constant glucose dosage as described by Brener and Mc Hugh, Regulation of gastric emptying of glucose, gastroenterology 1983}
                                if MMC_Status = 2 then GLP1gastricEmptying_stress := 0.5*GLP1gastricEmptying_stress; {situation during MMC phase 2}
                                if MMC_Status = 3 then GLP1gastricEmptying_stress := 0*GLP1gastricEmptying_stress; {situation during MMC phase 3}
      Viscousgastricemptying_stress := 0.21*corpus_viscosity/gastric_tone;        {Poisueille flow:   flow rate = deltaP.pi.r^4/(8.eta.L)}
                                if MMC_Status = 2 then ViscousgastricEmptying_stress := 0.5*ViscousgastricEmptying_stress; {situation during MMC phase 2}
                                if MMC_Status = 3 then ViscousgastricEmptying_stress := 0*ViscousgastricEmptying_stress; {situation during MMC phase 3}
      if pHgastricEmptying_stress <0 then pHgastricEmptying_stress:=0;
      emptying_stress := CCKgastricemptying_stress + Osmoticgastricemptying_stress+ pHgastricemptying_stress+ GLP1gastricemptying_stress +Viscousgastricemptying_stress;
           emptying_stress:= 0.2*emptying_stress; {correction 0.2 based on data by Guido Camps, thesis}

           fraction_emptying :=1/(1+emptying_stress);{the mechanism is that the antrum volume empties, modulated by feedback signals}

     if Volume[4]>0.1 then
          begin
          up[4]:= 1*Volume[4]*fraction_emptying/0.1111;                      {maxima in up[4] en retropulsion volume 3*per minute emptying on the basis of nutrient feedback}
          retropulsionvolume:=  1*Volume[4]*(1-fraction_emptying)/0.1111;
          {if MMC_status = 3 then up[4]:= up[4]*2;}
          end
      else
          begin
            up[4]:=0;
            retropulsionvolume:= 0;
          end;

    {4. description of the retropulsion of antrum to corpus and Pyloric passage}
    {emptying via variaties gastric sieve: 0: liquid;  1: particle diameter < 2mm; 2: particle diameter > 2 mm}

          for j:= 1 to comp do
              begin
              Change := up[4]*liquid[4,j]*dt/volume[4];
              liquid[3,j] := liquid[3,j] + retropulsionvolume*dt*liquid[4,j]/volume[4];
              liquid[4,j] := liquid[4,j] - change- retropulsionvolume*dt*liquid[4,j]/volume[4];
              liquid[5,j] := liquid[5,j] + change;
              end;
                   j:= 1; gastric_emptying_water := up[4]*liquid[4,j]*dt/volume[4]; {gastric emptying of water}
                   j:= 5; gastric_emptying_pharmaceutical := up[4]*liquid[4,j]*dt/volume[4]; {gastric emptying of pharmaceutical}
          for compoundnumber := 1 to number_of_compounds do
              begin

               gastric_sieve[compoundnumber] := exp(-solid[Compoundnumber]); {this expression gives a number to the particle size dependence of gastric sieving of solid particles by the Pylorus}
               if MMC_status = 3 then  gastric_sieve[compoundnumber] := 1;    {no sieving for MMC = 3}
              Change := up[4]*gastric_sieve[compoundnumber]*compound[4,compoundnumber]*dt/volume[4];
              compound[3,compoundnumber] := compound[3,compoundnumber] + retropulsionvolume*dt*compound[4,compoundnumber]/volume[4];
              compound[4,compoundnumber] := compound[4,compoundnumber] - change - retropulsionvolume*dt*compound[4,compoundnumber]/volume[4];
              compound[5,compoundnumber] := compound[5,compoundnumber] + Change;

              end;

          volumina;
    end;







{2a. description of the redistribution between fundus and corpus through slow sedimentation or creaming in the liquid medium, this part needs to be updated}
if 1=1 then    {on-off switch for this redistribution; this is the most basic description}
  begin
     for j:= 1 to comp do liquidscorpusfundus[j]:= liquid[3,j]+fundus[j] ;               {liquidscorpusfundus[j]is all component j in corpus and fundus}
     for j:= 1 to comp do targetcorpus[j] := liquid[3,j];                                {targetcorpus[j] is the ultimate amount of component j (which is all components)in the corpus}
      targetcorpus[2] := 0{1}*liquidscorpusfundus[2];    {creaming of fat}                  {the ultimate "target" amount after creaming of the fat is all that was previously present in corpus and fundus}
                                                                                             {"0" becasue no fat will remain in corpus after all fat has creamed to fundus}
     for j:= 1 to comp do
       begin
         x := (targetcorpus[j]-liquid[3,j])*10*dt/corpus_viscosity;                {this describes the speed of creaming of fat; droplet size not taken into account here as they are supposed to be large lumps that cream very fast; the factor 10 corrsponds to a high limiting value that the calculational cycle can still deal with; similar to below for compound creaming/sedimentation}
         fundus[j]:= fundus[j]-x;                                                   {if x<0 then the fat will cream out of the corpus into the fundus}
         liquid[3,j]:=liquid[3,j]+x;                                                {an improvement would be to describe emulsion droplets as fat-containing compound with size solid[compoundnuber]}
         Vfundus:= Vfundus -x/density[j];
         volume[3]:= volume[3]+x/density[j];
       end;
  end;

 if 1=1 then    {on-off switch for compound creaming; this is an improved description: for example emulsion droplets can be modelled as compounds that cream}
  begin
     x := (compounddensity[compoundnumber]-1); {x is ∆ρ; "1" temporarily corresponds to the density of the gastric liquid phase}{x<1 implies creaming}
     x:= x*0.6*10{temporary area interface between fundus and corpus in cm2}*solid[compoundnumber]*9.81*SQR(solid[compoundnumber])*dt/corpus_viscosity;  {Stokes equation for creaming/sedimentation}


     for compoundnumber:= 1 to number_of_compounds do
       begin

        if x>0 then         {the compound sediments from fundus to antrum}
          begin
            compound[3,compoundnumber]:= compound[3,compoundnumber]+x*compoundfundus[compoundnumber]/Vfundus;
            compoundfundus[compoundnumber]:= compoundfundus[compoundnumber] - x*compoundfundus[compoundnumber]/Vfundus;
          end;
        if x<0 then          {the compound creams up from antrum to fundus}
          begin
            compoundfundus[compoundnumber]:= compoundfundus[compoundnumber] - x*compound[3,compoundnumber]/Volume[3];
            compound[3,compoundnumber]:= compound[3,compoundnumber]+x*compound[3,compoundnumber]/Volume[3];
           end;
       end;
     volumina;
  end;

  {end 2a.}

{2b. description of the redistribution of compounds between corpus en antrum by slow sedimetation or creaming in the liquid medium; this part needs to be updated}
 if 1=1 then    {on-off switch compound redistribution}
  begin
     x := (compounddensity[compoundnumber]-1); {x is ∆ρ; "1" temporarily corresponds to the density of the gastric liquid phase}{for x<0 the compound creams up from 4 to 3}
     x:= x*0.6*10{temporary area of the interface between corpus and antrum in cm2}*solid[compoundnumber]*9.81*SQR(solid[compoundnumber])*dt/(corpus_viscosity);  {Stokes equation for creaming/sedimentation}
     for compoundnumber:= 1 to number_of_compounds do
       begin

        if x>0 then              {the compound sediments from 3 to 4}
          begin
            compound[4,compoundnumber]:= compound[4,compoundnumber]+x*compound[3,compoundnumber]/Volume[3];
            compound[3,compoundnumber]:= compound[3,compoundnumber] - x*compound[3,compoundnumber]/Volume[3];
          end;
        if x<0 then                        {the compound creams up from 4 to 3}
          begin
            compound[3,compoundnumber]:= compound[3,compoundnumber] - x*compound[4,compoundnumber]/Volume[4];
            compound[4,compoundnumber]:= compound[4,compoundnumber]+x*compound[4,compoundnumber]/Volume[4];
           end;
       end;
     volumina;
  end;

  {end 2b.}


  {3. introduction of viscosity variations in corpus en antrum}
   if 1 = 1{0}{1} then     {schakelaar***************************************************}
    begin
     fundus_viscosity := 2 ;
     corpus_viscosity := 2 {10 mPa.s};    {viscosity of the corpus volume}
     antral_viscosity := 2 {10 mPa.s};    {viscosity of the antral volume}
    end;
  {end 3.}


  {4. solid parts in lumen of corpus and lumen of antrum follow different kinetics; therefore the model describes solids_corpus and solids_antrum}
    {solids_corpus and solids_antrum have to be calculated in a separate routine, for examples on the basis of control by acidification}
    {temporary description transfer to solids}

    if 1 = 1{1} then     {switch}
    begin
     ber_alle_pHs;
    for j:= 1 to comp do  vloksamenstelling[j]:= 0;
      pHvlokken:= 4.5;
      vloksamenstelling[61] := 1;          {vloksamenstelling bepaald welke componenten bij de uitvlokking betrokken zijn; nu caseine}
      vloksamenstelling[62] := 1;
      hulp := vloksamenstelling[61] +vloksamenstelling[62];
      vloksamenstelling[80] := 1*Ber_hydrogen(pHvlokken,hulp ,vloksamenstelling[55],vloksamenstelling[66],vloksamenstelling[67], vloksamenstelling[4], vloksamenstelling[68]);
    {4a. corpus solids}
      vlokkingssnelheid_corpus := 0;
      Kineticfactor := 3000{3E3};
      if {((pH[3] < 6) and (pH[3] > 3{3}) and} (Volume[3] <>0){)} then vlokkingssnelheid_corpus := kineticfactor*EXP(-8*ABS(pH[3]-pHvlokken));

      for j:= 1 to comp do
        begin
        change := vlokkingssnelheid_corpus*vloksamenstelling[j]*liquid[3,j]*dt;
        if change > liquid[3,j] then change := liquid[3,j] else
           begin
           solids_corpus[j] := solids_corpus[j] + change ;
           liquid[3,j] := liquid[3,j]- change;
           end;
        end;

      oplossnelheid_vlokken_corpus := 0; {flocks slowly redissolve, for example below the isolectric pH }
      discost := 1;
      {if pH[3] > 5 then oplossnelheid_vlokken_corpus := discost*(EXP(pH[3]-5)-1);
      If pH[3] < 4 then oplossnelheid_vlokken_corpus := discost*(EXP(4-pH[3])-1);}
      hulp:= {EXP(}10*ABS(pH[3]-pHvlokken){-1};
      oplossnelheid_vlokken_corpus := discost*hulp/(1+hulp);

      total_solids_corpus := 0;
      for j:= 1 to comp do total_solids_corpus := total_solids_corpus+ solids_corpus[j];
      if total_solids_corpus> 0.00001 then bolus_oppervlak := exp(0.666666666666*LN(total_solids_corpus)) else bolus_oppervlak := 0{1*total_solids_corpus};
      if total_solids_corpus >0 then
        begin
          for j:= 1 to comp do
          begin
          hulp := oplossnelheid_vlokken_corpus*bolus_oppervlak*vloksamenstelling[j]*dt/total_solids_corpus;
          if hulp > solids_corpus[j] then hulp:= solids_corpus[j];
          solids_corpus[j] := solids_corpus[j] - hulp;
          liquid[3,j] := liquid[3,j]+ hulp;
          end;
        end;


      if 1 = 0{1} then     {switch}
      begin
      syneresesnelheid := 0{0.01};    {syneresis of the flocks}
      for j:= 1 to comp do samenstelling_syneresevloeistof[j] := 0;
        samenstelling_syneresevloeistof[1]:= 1;
      for j:= 1 to comp do
        begin
        solids_corpus[j] := solids_corpus[j] - samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_corpus[j]*dt;
        liquid[3,j] := liquid[3,j]+ samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_corpus[j]*dt;
        end;
      end;

      solids_in_corpus:= 0;
      for j:= 1 to comp do solids_in_corpus := solids_in_corpus + solids_corpus[j];
      {end 4a.}

  {4b. antrum solids caseine}
       vlokkingssnelheid_antrum := 0;
      Kineticfactor := 1E3{1};
      if ((pH[4] < 5.0) and (pH[4] > 4) and (Volume[4] <>0)) then vlokkingssnelheid_antrum := kineticfactor*EXP(-ABS(pH[4]-pHvlokken));
      change:= 0;
      for j:= 1 to comp do
        begin
        change := vlokkingssnelheid_antrum*vloksamenstelling[j]*(liquid[4,j]/Volume[4])*dt;
        solids_antrum[j] := solids_antrum[j] + change ;
        liquid[4,j] := liquid[4,j]- change;
        end;


      total_solids_antrum:= 0;
      for j:= 1 to comp do total_solids_antrum := total_solids_antrum + solids_antrum[j];
      hulp:= solids_antrum[80];
      solids_antrum[80] := Vloksamenstelling[80]*total_solids_antrum;
      liquid[4,80] := liquid[4,80]- solids_antrum[80]+hulp;

      syneresesnelheid := 0{0.01};    {syneresis of the flocks}
      for j:= 1 to comp do samenstelling_syneresevloeistof[j] := 0;
        samenstelling_syneresevloeistof[1]:= 1;
      for j:= 1 to comp do
        begin
        solids_antrum[j] := solids_antrum[j] - samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_antrum[j]*dt;
        liquid[4,j] := liquid[4,j]+ samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_antrum[j]*dt;
        end;
      {end 4b.}
  {4c: all flocks move to the corpus}
     For j:=  1 to comp do
       begin
       solids_corpus[j]:= solids_corpus[j] + solids_antrum[j];
       solids_antrum[j] := 0;
       end;
     total_solids_corpus:= 0;
     for j:= 1 to comp do total_solids_corpus := total_solids_corpus + solids_corpus[j];

     total_solids_antrum := 0;
  {end 4.}
  end;

      total_solids_corpus:= 0;
      for j:= 1 to comp do total_solids_corpus := total_solids_corpus + solids_corpus[j]+solids_fundus[j];
      hulp:= solids_corpus[80]+solids_fundus[80];
      solids_corpus[80] := Vloksamenstelling[80]*total_solids_corpus;

      liquid[3,80] := liquid[3,80]{+fundus[80]}- solids_corpus[80]-solids_fundus[80]+hulp;

end;

end.

