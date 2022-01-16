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
     Gastric_Tone := 3+ 0.0177*(Gastric_Volume {+ total_solids_corpus})*(1-0.5*Transition); {maagdruk, mBar} {Janssen}
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
         if (compoundviscosity[compoundnumber]>100*viscosity_liquids_in_corpus) and (solid[compoundnumber] = 0) then solid[compoundnumber]:= 10;       {dikke vloiestoffen worden geemulgeerd tot grote druppels in de chyme}
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
  viscosity_liquids_in_corpus := 5;      {corrigeren voor emulsies en polymer in oplossing}
  for compoundnumber := 1 to number_of_compounds do if compound[3,compoundnumber]>1E-40 then viscosity_liquids_in_corpus := viscosity_liquids_in_corpus + compoundviscosity[compoundnumber]*(1/(1+Volume[3]*compounddensity[compoundnumber]/compound[3,compoundnumber]))  ;
  if phi>0.5 then
        begin
          corpus_viscosity := viscosity_liquids_in_corpus;
          for compoundnumber:= 1 to number_of_compounds do
             corpus_viscosity := corpus_viscosity + solid[compoundnumber]*compound[3,compoundnumber]*compoundviscosity[compoundnumber];
          corpus_viscosity := corpus_viscosity/solids_in_corpus;
        end

     else corpus_viscosity := viscosity_liquids_in_corpus*exp(-eta*phimax*LN(1-phi/phimax));   {Krieger-Dougherty}
end;

procedure gastric_redistribution;
var liquidscorpusfundus, targetcorpus,solidscorpusfundus: component;
   x, fraction_retropelled, retropulsionvolume, vlokkingssnelheid_corpus, vlokkingssnelheid_antrum, oplossnelheid_vlokken_corpus,oplossnelheid_vlokken_antrum,
        syneresesnelheid, verdeler, mengfactor, mengsnelheid, pHvlokken, hulp, bolus_oppervlak, kineticfactor, change, discost, target_Vcorpus, Corpus,
        up_corpus_to_fundus, down_fundus_to_corpus, target_antralvolume, antral, transfer_efficiency_fundus_corpus: real;
{de beschrijving van de maag is gebaseerd op een maagmodel dat bestaat uit een corpus en een antrum die in de lijn staan tussen de slokdarm en het duodenum.
Gekoppeld aan het corpus is het fundus, dat meer een voorraadvat vormt voor extra volume en waar ook vast materiaal tijdelijk opgeslagen wordt}
{De drie compartimenten kunnen met elkaar uitwisselen, waarbij de gemiddelde volumina zich aanpassen naar de heersende gemiddelde druk in de compartimenten, de maagdruk.}
{De maagdruk wordt bepaald door het volume van de fundus X een constante, waarbij de constante can veranderen met het CCK niveau in het bloed: het daalt in de fed state}
{een peristaltische beweging van corpus naar antrum stuwt het voedsel naar het antrum, waarna het antrum zich geheel samentrekt en daarbij het voedsel geheel in het
duodenum loost, dan wel slechts gedeeltelijk vanwege de regulering van maaglediging, Het deel dat eventueel niet in het duodenum leegt komt terug in het antrum}
{Het transport van fundus naar antrum is viscositeits afhankelijk}
begin
  volumina;
  antral := 0.1{0};
cyclus := cyclus+ dt;
if cyclus > 0.3333 then begin cyclus := 0; tel := 0; up[4]:=0 end;  {maag-peristaltische cyclus gesteld op 3 maal per minuut}
if cyclus< 0.1111 then
    begin

  {1. persistaltiek van fundus naar corpus}
      Vfunduscorpus := Vfundus+Volume[3];  {beschrijft de gezamelijke inhoud vast en vloeibaar van fundus en corpus, beginwaarde}
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
 if cyclus < 0.2222 then {terugstroming corpus naar fundus en stroming van corpus naar antrum}
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

    {3. Transport van corpus naar antrum door peristaltiek}
      pressures;
      gastric_compartment_viscosities;
      maxantralvolume := 60{240}{240}{120};
      target_antralvolume :=  5+maxantralvolume*gastric_tone/(1{5}+gastric_tone); {target antral volume}
      up[3] := Corpus*((1-Volume[4]/target_antralvolume)/0.1111)*0.1;   {het is aangetoond dat dit transport vertraagd wordt door de viscositeit: Benini 1995, maar pas op: het gaat heir om de viscositeit in de maag, na verdunnign met maagsap!}
         if MMC_status = 3 then Up[3] := 2*up[3];

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
{if cyclus < 0.3333 then}  {emptying en retropulsie van antrum naar corpus}
    begin

      {emptying}
    isotonic := 0.05; {0.05 glucose/ml}
    osmcont:=2{2};

      CCKgastricEmptying_stress :=  2 *totalrateCCKrelease + 1*CCK; {CCKEffectOnGastricEmptying*CCK}
                                if MMC_Status = 2 then CCKgastricEmptying_stress := 0.5*CCKgastricEmptying_stress; {situatie tijdens MMC fase 2}
                                if MMC_Status = 3 then CCKgastricEmptying_stress := 0*CCKgastricEmptying_stress; {situatie tijdens MMC fase 3}
      OsmoticgastricEmptying_stress := osmcont*osmotic_value_brush_duodenum/isotonic; {osmotic effect}
                                if MMC_Status = 2 then OsmoticgastricEmptying_stress := 0.5*OsmoticgastricEmptying_stress; {situatie tijdens MMC fase 2}
                                if MMC_Status = 3 then OsmoticgastricEmptying_stress := 0*OsmoticgastricEmptying_stress; {situatie tijdens MMC fase 3}
      pHgastricEmptying_stress := 1*(abs(6-pH[5])/6);     {pH effect}  {prefactor 20}
                                if MMC_Status = 2 then pHgastricEmptying_stress := 0.5*pHgastricEmptying_stress; {situatie tijdens MMC fase 2}
                                if MMC_Status = 3 then pHgastricEmptying_stress := 0*pHgastricEmptying_stress; {situatie tijdens MMC fase 3}
      GLP1gastricemptying_stress := 10*GLP1;    {400 geeft zelfde remming by constante glucose oplossing dosering als Brener and Mc Hugh regulation of gastric emptying of glucose, gastroenterology 1983}
                                if MMC_Status = 2 then GLP1gastricEmptying_stress := 0.5*GLP1gastricEmptying_stress; {situatie tijdens MMC fase 2}
                                if MMC_Status = 3 then GLP1gastricEmptying_stress := 0*GLP1gastricEmptying_stress; {situatie tijdens MMC fase 3}
      Viscousgastricemptying_stress := 0.21*corpus_viscosity/gastric_tone;        {Poisueille flow:   debiet = deltaP.pi.r^4/(8.eta.L)}
                                if MMC_Status = 2 then ViscousgastricEmptying_stress := 0.5*ViscousgastricEmptying_stress; {situatie tijdens MMC fase 2}
                                if MMC_Status = 3 then ViscousgastricEmptying_stress := 0*ViscousgastricEmptying_stress; {situatie tijdens MMC fase 3}
      if pHgastricEmptying_stress <0 then pHgastricEmptying_stress:=0;
      emptying_stress := CCKgastricemptying_stress + Osmoticgastricemptying_stress+ pHgastricemptying_stress+ GLP1gastricemptying_stress +Viscousgastricemptying_stress;
           emptying_stress:= 0.2*emptying_stress; {correctie 0.2 op basis van data van Guido Camps}

           fraction_emptying :=1/(1+emptying_stress);{het mechanisme is dat het antrum volume leegt, gemodereerd door feedback}

     if Volume[4]>0.1 then
          begin
          up[4]:= 1*Volume[4]*fraction_emptying/0.1111;                      {maxima in up[4] en retropulsie volume 3*per minuut emptyen op basis van nutrient feedback}
          retropulsionvolume:=  1*Volume[4]*(1-fraction_emptying)/0.1111;
          {if MMC_status = 3 then up[4]:= up[4]*2;}
          end
      else
          begin
            up[4]:=0;
            retropulsionvolume:= 0;
          end;

    {4. beschrijving van de retropulsie van het antrum naar het corpus en de Pylorus passage}
    {emptying via variaties gastric sieve: 0: vloeibaar;  1: deeltjes < 2mm; 2: deeltjes > 2 mm}

          for j:= 1 to comp do
              begin
              Change := up[4]*liquid[4,j]*dt/volume[4];
              liquid[3,j] := liquid[3,j] + retropulsionvolume*dt*liquid[4,j]/volume[4];
              liquid[4,j] := liquid[4,j] - change- retropulsionvolume*dt*liquid[4,j]/volume[4];
              liquid[5,j] := liquid[5,j] + change;
              end;
                   {maaglediging water} j:= 1; gastric_emptying_water := up[4]*liquid[4,j]*dt/volume[4]; {maaglediging water}
                   {maaglediging pharmaceutical} j:= 5; gastric_emptying_pharmaceutical := up[4]*liquid[4,j]*dt/volume[4]; {maaglediging pharmaceutical}
          for compoundnumber := 1 to number_of_compounds do
              begin

               gastric_sieve[compoundnumber] := exp(-solid[Compoundnumber]);
               if MMC_status = 3 then  gastric_sieve[compoundnumber] := 1;
              Change := up[4]*gastric_sieve[compoundnumber]*compound[4,compoundnumber]*dt/volume[4];
              compound[3,compoundnumber] := compound[3,compoundnumber] + retropulsionvolume*dt*compound[4,compoundnumber]/volume[4];
              compound[4,compoundnumber] := compound[4,compoundnumber] - change - retropulsionvolume*dt*compound[4,compoundnumber]/volume[4];
              compound[5,compoundnumber] := compound[5,compoundnumber] + Change;

              end;

          volumina;
    end;







{2a. beschrijving van herverdeling tussen fundus en corpus door langzame sedimetatie of oproming in de vloeibare toestand}
if 1=0 then    {aan-uit schakelaar componenten oproming}
  begin
     for j:= 1 to comp do liquidscorpusfundus[j]:= liquid[3,j]+fundus[j] ;
     for j:= 1 to comp do targetcorpus[j] := liquid[3,j];
      targetcorpus[2] := 1*liquidscorpusfundus[2];    {opromen van vet in de vloeibare fase}

     for j:= 1 to comp do
       begin
         x := (targetcorpus[j]-liquid[3,j])*10*dt/corpus_viscosity;
         fundus[j]:= fundus[j]-x;
         liquid[3,j]:=liquid[3,j]+x;
         Vfundus:= Vfundus -x/density[j];
         volume[3]:= volume[3]+x/density[j];
       end;
  end;

 if 1=0 then    {aan-uit schakelaar compound oproming}
  begin
     x := (compounddensity[compoundnumber]-1); {x bepaald wat de uiteindelijke toestand wordt na oproming of sedimentatie}
     x:= (x/ABS(x+1))*10*0.001*solid[compoundnumber]*9.81*SQR(solid[compoundnumber])*dt/(18*corpus_viscosity/1000);  {Stokes verrgelijking voor oproming/sedimentatie over 10 cm}
     for compoundnumber:= 1 to number_of_compounds do
       begin

        if x>0 then
          begin
            compound[3,compoundnumber]:= compound[3,compoundnumber]+x*compoundfundus[compoundnumber];
            compoundfundus[compoundnumber]:= compoundfundus[compoundnumber] - x*compoundfundus[compoundnumber];
          end;
        if x<0 then
          begin
            compound[3,compoundnumber]:= compound[3,compoundnumber]-x*compoundfundus[compoundnumber];
            compoundfundus[compoundnumber]:= compoundfundus[compoundnumber] + x*compoundfundus[compoundnumber];
           end;
       end;
     volumina;
  end;

  {einde 2a.}

{2b. beschrijving van herverdeling tussen corpus en antrum door langzame sedimetatie of oproming in de vloeibare toestand}
 if 1=1 then    {aan-uit schakelaar compound oproming}
  begin
     x := (compounddensity[compoundnumber]-1); {x bepaald wat de uiteindelijke toestand wordt na oproming of sedimentatie}
     x:= (x/ABS(x+1))*20*0.001*solid[compoundnumber]*9.81*SQR(solid[compoundnumber])*dt/(18*corpus_viscosity/1000);  {Stokes verrgelijking voor oproming/sedimentatie over 5 cm}
     for compoundnumber:= 1 to number_of_compounds do
       begin

        if x>0 then
          begin
            compound[4,compoundnumber]:= compound[4,compoundnumber]+x*compound[3,compoundnumber];
            compound[3,compoundnumber]:= compound[3,compoundnumber] - x*compound[3,compoundnumber];
          end;
        if x<0 then
          begin
            compound[4,compoundnumber]:= compound[4,compoundnumber]-x*compound[3,compoundnumber];
            compound[3,compoundnumber]:= compound[3,compoundnumber] + x*compound[3,compoundnumber];
           end;
       end;
     volumina;
  end;

  {einde 2b.}


  {3. invoeren viscositeit variaties in corpus en antrum}
   if 1 = 0{1}{0}{1} then     {schakelaar***************************************************}
    begin
     fundus_viscosity := 2 ;
     corpus_viscosity := 2 {10 mPa.s};    {viscosity of the corpus volume}
     antral_viscosity := 2 {10 mPa.s};    {viscosity of the antral volume}
    end;
  {einde 3.}


  {4. vaste delen in het corpus en antrum lumen volgen een andere kinetiek; we beschrijven daarom solids_corpus en solids_antrum}
    {solids_corpus en solids_antrum moeten nog in een aparte routine berekend worden, bijvoorbeeld op basis van verzuring}
    {voorlopige beschrijving transfer naar solids}
    total_solids_corpus:= 0;
    if 1 = 0{1} then     {schakelaar}
    begin
    for j:= 1 to comp do  vloksamenstelling[j]:= 0;
      pHvlokken:= 4.5;
      vloksamenstelling[61] := 1;
      vloksamenstelling[62] := 1;
      hulp := vloksamenstelling[61] +vloksamenstelling[62];
      vloksamenstelling[80] := 1*Ber_hydrogen(pHvlokken,hulp ,vloksamenstelling[55],vloksamenstelling[66],vloksamenstelling[67], vloksamenstelling[4]);vlokkingssnelheid_corpus := 1{0.01}; {deze zal in werkelijkheid varieren met bijvoorbeeld de pH}
    {4a. corpus solids caseine}
      vlokkingssnelheid_corpus := 0;
      Kineticfactor := 3E3;
      if ((pH[3] < 5.0) and (pH[3] > 4) and (Volume[3] <>0)) then vlokkingssnelheid_corpus := kineticfactor*EXP(-ABS(pH[3]-pHvlokken));
      change :=0;
      for j:= 1 to comp do
        begin
        change := vlokkingssnelheid_corpus*vloksamenstelling[j]*(liquid[3,j]/Volume[3])*dt;
        solids_corpus[j] := solids_corpus[j] + change ;
        liquid[3,j] := liquid[3,j]- change;
        end;

      oplossnelheid_vlokken_corpus := 0; {de vlokken lossen weer op, bijvoorbeeld beneden de isolectrische pH }
      discost := 0.01;
      if pH[3] > 5 then oplossnelheid_vlokken_corpus := discost*(EXP(pH[3]-5)-1);
      If pH[3] < 4 then oplossnelheid_vlokken_corpus := discost*(EXP(4-pH[3])-1);
      total_solids_corpus := 0;
      for j:= 1 to comp do total_solids_corpus := total_solids_corpus+ solids_corpus[j];
      if total_solids_corpus> 0.00001 then bolus_oppervlak := exp(0.666666666666*LN(total_solids_corpus)) else bolus_oppervlak := 1*total_solids_corpus;
      if total_solids_corpus >0 then
        begin
          for j:= 1 to comp do
          begin
          hulp := oplossnelheid_vlokken_corpus*bolus_oppervlak*vloksamenstelling[j]/total_solids_corpus*dt;
          solids_corpus[j] := solids_corpus[j] - hulp;
          if solids_corpus[j] <0 then solids_corpus[j] := solids_corpus[j] + hulp else {**************************}
          liquid[3,j] := liquid[3,j]+ hulp;
          end;
        end;



      syneresesnelheid := 0{0.01};    {synerese van de vlokken}
      for j:= 1 to comp do samenstelling_syneresevloeistof[j] := 0;
        samenstelling_syneresevloeistof[1]:= 1;
      for j:= 1 to comp do
        begin
        solids_corpus[j] := solids_corpus[j] - samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_corpus[j]*dt;
        liquid[3,j] := liquid[3,j]+ samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_corpus[j]*dt;
        end;
      total_solids_corpus:= 0;
      for j:= 1 to comp do total_solids_corpus := total_solids_corpus + solids_corpus[j];
      {einde 4a.}

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

      syneresesnelheid := 0{0.01};    {synerese van de vlokken}
      for j:= 1 to comp do samenstelling_syneresevloeistof[j] := 0;
        samenstelling_syneresevloeistof[1]:= 1;
      for j:= 1 to comp do
        begin
        solids_antrum[j] := solids_antrum[j] - samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_antrum[j]*dt;
        liquid[4,j] := liquid[4,j]+ samenstelling_syneresevloeistof[j]*syneresesnelheid*solids_antrum[j]*dt;
        end;
      {einde 4b.}
  {4c: alle vlokken gaan naar corpus}
     For j:=  1 to comp do
       begin
       solids_corpus[j]:= solids_corpus[j] + solids_antrum[j];
       solids_antrum[j] := 0;
       end;
     total_solids_corpus:= 0;
     for j:= 1 to comp do total_solids_corpus := total_solids_corpus + solids_corpus[j];

     total_solids_antrum := 0;
  {einde 4.}
  end;

      total_solids_corpus:= 0;
      for j:= 1 to comp do total_solids_corpus := total_solids_corpus + solids_corpus[j]+solids_fundus[j];
      hulp:= solids_corpus[80]+solids_fundus[80];
      solids_corpus[80] := Vloksamenstelling[80]*total_solids_corpus;

      liquid[3,80] := liquid[3,80]{+fundus[80]}- solids_corpus[80]-solids_fundus[80]+hulp;

end;

end.

