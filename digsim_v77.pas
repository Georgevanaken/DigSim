Program digsim_v77;
uses crt, sysutils, gastric_acid_secretion, globalvariables, calc_Mucin_hydrogen,
  calc_Casein_hydrogen, calc_BSA_hydrogen, calc_BLG_hydrogen, calc_Protein_hydrogen,
  SubcellularInsulin_Pederson, Liver_clearance,diet_from_file,
  Brush_border_enzymes, Model_Parameters, Gastric_motility,
  Compartment_properties, Glucose_homeostasis2, Hb_glycosylation, MMC;


 {segments "i" describes all parts of the digestive system. Originally i=1 was the dish, but later this was replaced by a larger number of dishes.
 It was chosen to break up the length of duodenum in two (before and after the Sphincter of Oddi), the jejunum in three, the ileum in 6,
 and the colon in 4 to adress the differences in length and variation in function within each kind.
 fundus and gall bladder were added later

 The segments are:

 2 = mouth
 3 = corpus
 4 = antrum
 5 = duodenum 1
 6 = duodenum 2
 7 = jejunum 1
 8 = jejunum 2
 9 = jejunum 3
 10 = ileum 1
 11 = ileum 2
 12 = ileum 3
 13 = ileum 4
 14 = ileum 5
 15 = ileum 6
 16 = ascending colon
 17 = transverse colon
 18 = decending colon
 19 = rectum
 20 = stool
 fundus
 gall bladder

 The components are in principle all molecular species. The number of components are unlimited. Unless stated differently, the amounts are in grams,
 except for the amounts enzymes which are expressed in units.
 The components are numbered as follows:
 1 = water
 2 = fat
 3 = glucose
 4 = protein (a generalized protein in dissolved state)
 5 = pharmaceutical component (the activity has not been defined yet)
 6 = pharmaeutical product inactive form (defined a inert component)
 7 = diglyceride (no distinguishment has yet been made among the group of diglycerides regarding fatty acid composition or sn position)
 8 = monoglyceride (no distinguishment has yet been made among the group of monoglycerides regarding fatty acid composition or sn position)
 10 = HDstarch (rapid digestible starch)
 11 = SDstarch  (slow digestible starch)
 12 = Rstarch   (resistant starch)
 20 = fatty acids {no distinguishment has yet been regarding fatty acid composition)
 21 = sucrose
 22 = fructose
 23 = maltose
 24 = tracer maltose (13C tracer)
 25 = arabinose
 26 = xylose
 27 = acarbose
 34 = tracer glucose 1 (13C tracer)
 35 = tracer glucose 2  (2H tracer)
 36 = tracer sucrose (13C sucrose)
 37 = tracer HDstarch
 38 = tracer SDstarch
 39 = tracer Rstarch
 40 = pepsin
 41 = pancreatic proteases {no further subdivision between the pancreatic proteases has yet been made}
 48 = lingual lipase
 49 = gastric lipase
 50 = pancreatic lipase
 51 = conjugated bile salt
 52 = bile phospholipids
 53 = bile lysophospholipids
 54 = bile salt LCFA complex   (in reality this may comprise many species, this has not yet been included}
 55 = mucus protein, isoelectric in grams
 56 = pancreatische fosfolipase
 57 = salivary_amylase
 58 = pancreatic_amylase
 59 = unconjugated bile salt
 60 = CCK releasing protein
 61 = casein
 62 = caseinpep
 63 = caseinpan
 64 = caseinpeptides
 65 = small peptides
 66 = BSA
 67 = BLG

 70 = Leucine
 71 = Leucine tracer
 72 = I_c      competitive alpa amylase inhibitor
 73 = I_a      acompetitive alpha amylase inhibitor
 74 = I_nc     non competitive alpha amylase inhibitor

 78 = NaCl
 79 = HCO3-
 80 = H+ ionen (in moles)}




procedure initialiseer;        {this procedure defines the starting conditions for the model calculations}
begin
MMC_status := 0;
MMC_time[0] := 0; MMC_time[1]:= 0; MMC_time[2] := 0; MMC_time[3]:= 0; MMC_time[4]:= 0;
emptying_cycle := 0;
glycogene_store0 := 10; glycogene_store1:= 0; glycogene_store2:=0; glycogene_store := 100;
cyclus := 0;
blood_bile_salt_pool := 0.008 ;
compoundnumber := 0;     {this countee registers how many meals have been administered}
s:=0; {this counter administers the number of solid phases}
gastric_tone:= 0;        {gastric pressure initially set to 0}
osmotic_value_lumen_duodenum := 0.05;
bicarbonate_secretion:= 0;
corpus_viscosity:= 10;
antral_viscosity := 10;
for dishnumber := 1 to 5 do liquid[dishnumber,80]:= 01E-7;
calories_colon := 0;
calories_absorbed:= 0;
for i:= 2 to 16 do      {liquid, the amount of liquid and dissolved components i and compounds j unequal to 1 are set to 0. For component 1 (water) the amount is set to an arbitrary low value 1 to avoid devision by 0 in the calcualtion of concentrations}
 begin
      liquid[i,1] := 1;
      liquidnext[i,1] := 1;     {next liquid is the same as liquid, for the subsequent timestep}
      fundus[1]:= 1;
      absorption_rate[i,1]:= 0;     {absorption describes the efflux of componentents out of the lumen}
      concentration_brush[i,1] := 1;
      for compoundnumber:= 1 to number_of_compounds do    {all dissolved compounds are set to 0 for all segments}
       begin
          compound[i,compoundnumber] := 0;
          compoundfundus[compoundnumber]:= 0;
       end;
      for j:= 2 to comp do
       begin
       liquid[i,j] := 0;         {liquid describes a temporary state of amounts (as masses) of materials j in compartment i in the lumen of the GI tract}
       liquidnext[i,j] := 0;     {next liquid is the same as liquid, for the subsequent timestep}
       fundus[j]:= 0;
       for s:= 1 to 99 do
        begin
        solids[s,i,j]:= 0;
        solidsnext[s,i,j]:= 0;
        end;
       absorption_rate[i,j]:= 0;     {absorption describes the efflux of componentent i in grams/min out of the lumen of compartment j}
       enterocyte[i,j] := 0;    {amount of component i in the enterocyte layer of compartment j}
       concentration_brush[i,j] := 0;  {amount of component i in the brush border of compartment j}
       end;

 end;

  liquid[3,1]:= 25; liquid[4,1]:= 10; {the gastric corpus and antrum initially contain an arbitrary small amount of water, to avoid divisions by 0}
    for j := 1 to comp do begin solids_corpus[j] := 0; solids_antrum[j] := 0; solids_fundus[j] := 0; end; solids_in_fundus:= 0; solids_in_corpus := 0 ; solids_in_antrum:= 0; total_solids_in_stomach := 0;
  for i:= 5 to 15 do liquid[i,50] := 0.01; {basis lipase(component 50) activity in the lumen of small intestinal compartments 5 - 15}
  liquid[6,51]:= 0.1; {initial concentration of bile salts in het second part of the duodenum}
  liquid[6,50]:= 0.012/5; {initial activity of pancreatic lipase in the second part of the duodenum, which is approx. 1/5 of the activated (fed) value; reference S. Bernbäck, L. Bläckberg, O. Hernell, Biochim. Biophys. Acta. 1001, 286 (1989)}

  Initialise_glucose_dynamics;     {here the calculational procedure is started to establish the blood serum and insulin concentrations relevant for the start up conditions}


for j:=1 to comp do
 begin
 for i:= 5 to 15 do Mesentric[i,j] := 0;  {all values for the mesentric cellular tissue are set to 0}
 Blood[j]:= 0;           {all compounds are set to 0}
 end;
for i := 2 to 15 do               {all hormonal sepcies in the mesentric layer are set to 0}
begin
   Mesentric_CCK[i]:= 0;
   Mesentric_PYY[i]:= 0;
   Mesentric_GIP[i]:=0;
   Mesentric_GLP1[i]:=0;
   Mesentric_Ghrelin[i]:= 0;
   Mesentric_Motilin[i]:=0;
end;
HbA1c:= 0;                {glycation of hemoglobin is initially set to 0}
end;


Procedure caloric_densities_segments;  {this procedure calculates the caloric value of the luminal content and brush bordersof all intestinal segments}
var calories_segment : real;
begin;
for i:= 5 to 16 do
  begin
       calories_segment := 0;
       for j:= 1 to comp do calories_segment := calories_segment + liquid[i,j]*caloric_value[j];
       for compoundnumber := 1 to number_of_compounds do
           calories_segment := calories_segment + compound[i,compoundnumber]*caloric_value[compoundnumber];
       if Volume[i] > 1E-10 then caloric_density_local[i] := calories_segment/Volume[i] else caloric_density_local[i]:= 0;
       calories_brush := 0;
       for j := 1 to comp do calories_brush := calories_brush+ concentration_brush[i,j]*caloric_value[j];
  end;
end;

Procedure osmotic_values_lumen;    {this procedure calculates the osmotic values in the lumen on the basis of small sugars present and expressed as equivalent glucose concentrations; electrolytes and amaniacids etc are still to be added}
var hulp : real;

begin
 for i:= 5 to 15 do
  begin
  osmotic_value_lumen[i] := {isotonic}{background electrolyte}0 +concentration_lumen[i,3]+concentration_lumen[i,34]+concentration_lumen[i,35]+0.5*concentration_lumen[i,21]+0.5*concentration_lumen[i,36]+concentration_lumen[i,22]+0.5*concentration_lumen[i,23]+0.5*concentration_lumen[i,24]+ concentration_lumen[i,25]+ concentration_lumen[i,26]+ 0.2*concentration_lumen[i,27];  {currently the osmotic value is defined as corresponding glucose concentrations}
  end;

       osmotic_value_lumen_duodenum :=  osmotic_value_lumen[5]+1*osmotic_value_lumen[6]+ 0{.5}*osmotic_value_lumen[7] + 0{.5}*osmotic_value_lumen[8];
       osmotic_value_lumen_distal := osmotic_value_lumen[10]+1*osmotic_value_lumen[11]+ osmotic_value_lumen[12] + osmotic_value_lumen[13]+ osmotic_value_lumen[14] + osmotic_value_lumen[15] ;

end;

Procedure osmotic_values_brush; {this procedure calculates the osmotic values in the brush border on the basis of small sugars present and expressed as equivalent glucose concentrations; electrolytes and amaniacids etc are still to be added}
var hulp, vertrager : real;

begin
 osmotic_value_brush_duodenum:= 0;
 for i:= 5 to 15 do
  begin
  osmotic_value_brush[i] := concentration_brush[i,3]+concentration_brush[i,34]+concentration_brush[i,35]+ 0.5*concentration_brush[i,21]+0.5*concentration_brush[i,36]+concentration_brush[i,22]+0.5*concentration_brush[i,23]+0.5*concentration_brush[i,24]+ concentration_brush[i,25]+ concentration_brush[i,26]+ 0.2*concentration_brush[i,27];  {nu even als de glucose concentratie gedefinieerd}
  end;
   vertrager:= 0*60*dt;
   hulp :=  1*osmotic_value_brush[5]+1*osmotic_value_brush[6]+ 0.5*osmotic_value_brush[7] + 0.2*osmotic_value_brush[8];
       osmotic_value_brush_duodenum := (osmotic_value_brush_duodenum*vertrager + hulp)/(vertrager+1);
   osmotic_value_brush_distal := osmotic_value_brush[10]+1*osmotic_value_brush[11]+ osmotic_value_brush[12] + osmotic_value_brush[13]+ osmotic_value_brush[14] + osmotic_value_brush[15] ;

end;

Procedure Calories;    {this procedure calculates the calories in the lumen of different parts of the intestines: duodenum + jejunum, ileum and colon}
begin
  calories_proximal := 0;
  calories_distal:= 0;
  calories_brush := 0;
  calories_colon:= 0;
     for j := 1 to comp do
        begin
          for i:= 5 to 8 do calories_proximal := calories_proximal+ caloric_value[j]*liquid[i,j];
          for i:= 9 to 15 do calories_distal:= calories_distal + caloric_value[j]*liquid[i,j];
          calories_colon := calories_colon + caloric_value[j]*liquid[16,j];
        end;
     for compoundnumber := 1 to number_of_compounds do
        begin
          for i:= 5 to 8 do calories_proximal := calories_proximal+compound[i,compoundnumber]*calories_compound[compoundnumber];
          for i:= 9 to 15 do calories_distal:= calories_distal + compound[i,compoundnumber]*calories_compound[compoundnumber];
          calories_colon := calories_colon + compound[16,compoundnumber]*calories_compound[compoundnumber];
        end;
end;

Function time_of_day(minutes, referencetime: real):real;
 begin
  time_of_day := {minutes+referencetime} 11*60; {this function is currently used for averaged unstimulated saliva secretion only at 11 am}
 end;

Procedure saliva_secretion;
var  mucus_added, saliva_pH, hulp1, hulp2, hulp3, amylase_added, lingual_lipase_added :real;
const mucus_content_in_saliva = 0.003; {about 0.3 % proteinn in saliva, Dawes 1972}
const amylase_content_in_saliva = 2*75; {units/ml; 1 amylase unit = 1 mg maltose in 3 min, therefore 0.333 mg maltose /min}{according to Infogest protocol}
      lingual_lipase_content_in_saliva = 0.01; {a guess}

begin

  hulp1 :=  0.1 + 0.4*sin(6.283*(time_of_day(minutes,referencetime)-11*60)/(24*60)); {circadian cycle}
  hulp2:= 5*Volume[2]/(1+Volume[2]); {food volume stimulation}
  hulp3:= 4*ABS(6-pH[2])/(1+ABS(6-pH[2])); {pH stimulation}
  saliva_addition:= (hulp1+ hulp2+ hulp3)*dt; {ml/min}       {Dawes 1972, Motoc 2003}
  liquid[2,1]:= liquid[2,1] + saliva_addition;
  mucus_added := mucus_content_in_saliva*saliva_addition;
  amylase_added := amylase_content_in_saliva*saliva_addition;              {amylase content in units/ml}
   {**** ability to secrete salivary amylase ***}amylase_added := 1*amylase_added;
  lingual_lipase_added := lingual_lipase_content_in_saliva*saliva_addition;  {lingual lipase content in saliva, units lipase in saliva still to be determined accurately}
  liquid[2,55] := liquid[2,55] + mucus_added  ;
  liquid[2,57] := liquid[2,57] + amylase_added;            {amylase in units}
  liquid[2,48] := liquid[2,48]+ lingual_lipase_added; {lingual_lipase units not yet determinded}
  saliva_pH := 6;    {set value which is actually age and sex dependent, Motoc 2003}
  liquid[2,80] := liquid[2,80] + Mucin_hydrogen(saliva_pH, mucus_added ); {correction to Isolelectric point of mucus}
  total_saliva_secretion :=total_saliva_secretion + saliva_addition;
end;


Procedure pancreabiliar_secretion;
var   richt_pH, hulp, water_rate, bile_salt_rate, bile_phospholipids_rate, filling_rate_fraction, circulated_bile_salt_rate, pancreatic_bicarbonate_secretion_rate, pancreatic_enzyme_secretion_rate, pancreatic_fluid_secretion_rate : real;
const   liver_secretion_rate = 0.3 {0.3}; max_gallbladdervolume = 50; max_filling_rate = 0.015; max_emptying_rate = 0.3; bile_content = 0.2*0.007; phospholipids_content =0.007;
         target_blood_bile_salt_pool = 0.008; {gram, based on approximately 4 micromol/l (Adams 2015 Bile Acid Determination after Standardized...), 5 l bloodd en M = 408.6 g/mol (wikipedia) }
begin
                 {ref: Howard 1991, Gall bladder emptying patterns in response to a... :gives a view on the kinetics of gall bladdder emptying}
                 {Vezina 1994, Normal Fasting Volume and Postprandial Emptying of the...}
                 {Kamerling 2004 Gallbladder volume as a biomarker for the motilin effect in healthy... : shows the effect of motilin on gall bladder emptying; the direct effect of intravenous moltilin is small}
                 {Wikipedia bile: bile salts concentration in bile is about 0.7 % and dialy bile secretion is 400-800 ml }
                 {van Deest 1968 Bile Salt and Micellar Fat Concentration in Proximal Small Bowel: normal bile salt concentrations in the proximal SI are on the order of 3 mg/g, the cmc is at about 1.7 mg/g;
                      enterohepatic circulation originating from ileal absorption delivers about 80 % of the bile salts, de novo sysnthesis by the liver approximately 20% (1/6 to 1/5 parts)}
  {liver gall excretion into biliary duct}
  blood_bile_salt_pool := blood_bile_salt_pool + 0.2*(target_blood_bile_salt_pool - blood_bile_salt_pool)*dt ; {stabilisation of bile salts in blood, 10 minuten halflife time ???}
  circulated_bile_salt_rate := 0.2*blood_bile_salt_pool; {the amount that enters bile from the blood stream}
  blood_bile_salt_pool := blood_bile_salt_pool + ( totalabsorption_rate[51]-circulated_bile_salt_rate)*dt ;  {reabsorption from colon and decrease through bile secretion}
  water_rate := (1-bile_content-phospholipids_content)*liver_secretion_rate;

  bile_salt_rate := bile_content*water_rate +circulated_bile_salt_rate;  {20% de novo + 80% recirculation }
  bile_phospholipids_rate := phospholipids_content*water_rate;
  {gall bladder filling}
  filling_rate_fraction:=(max_gallbladdervolume-gall_bladder_volume)/max_gallbladdervolume;   {a part of the diluted bile originating from the liver enters the bile sac (this fraction), the rest is delivered directly to the duodenum}
  gall_bladder[1] := gall_bladder[1] + water_rate*filling_rate_fraction*dt - gall_bladder_volume*0.01*dt  ;   {water suppletion from liver and depletion through absorption}
  gall_bladder[51]:= gall_bladder[51] + bile_salt_rate*filling_rate_fraction*dt ;    {bile_salt suppletie by liver}
  gall_bladder[52]:= gall_bladder[52] + bile_phospholipids_rate*filling_rate_fraction*dt ;       {bile phospholipids from liver}
    gall_bladder_volume := gall_bladder[1] + gall_bladder[51] + gall_bladder[52];
    conjugated_bile_salt_concentration_in_gall_bladder := gall_bladder[51]/gall_bladder_volume; {concentration bile salts in bile}
  {gallbladder emptying}
  if gall_bladder_volume> max_gallbladdervolume then gall_bladder_emptying_rate :=  0.1*(gall_bladder_volume-max_gallbladdervolume)/dt; {overflow of the gall bladder}
  gall_bladder_emptying_rate := gall_bladder_emptying_rate + 0{1} + (1*totalrateCCKrelease/(5{+totalrateCCKrelease}) + 2*exp(-duodenal_conjugated_bile_salt_concentration/0.0002));{this describes the emptying of bile into the duodenum; I related it to the neural coupling with CCK release and a feedback coupling for a too low concnetraiton in the duodenum}
  if MMC_status = 3 then  gall_bladder_emptying_rate := gall_bladder_emptying_rate + 2; {the assumption is that bile bladder contraction occurs as part of the MMC 3 motility, not directly by motilin stimulation}
  gall_bladder_emptying_rate := max_emptying_rate*gall_bladder_emptying_rate/(max_emptying_rate+gall_bladder_emptying_rate); {limits the emptying rate to a maximum}

  if gall_bladder_emptying_rate*dt > gall_bladder_volume then gall_bladder_emptying_rate := 0.9*gall_bladder_volume/dt ;  {avoiding 0-values}
   hulp := dt*gall_bladder_emptying_rate/gall_bladder_volume;
  gall_bladder[1] := gall_bladder[1]-gall_bladder[1]*hulp ;  {water in bile emptying}
  gall_bladder[51]:= gall_bladder[51]- gall_bladder[51]*hulp ;   {bile_salt emptying}
  gall_bladder[52]:= gall_bladder[52]- gall_bladder[52]*hulp ;       {bile phospholipids in bile}

  {release into duodenum}
  liquid[6,1] := liquid[6,1]+ gall_bladder[1]*hulp + (1-filling_rate_fraction)*water_rate*dt;                         {this is the part that is directly deliverd to the duodenum}
  liquid[6,51] := liquid[6,51]+ gall_bladder[51]*hulp + (1-filling_rate_fraction)*bile_salt_rate*dt;                  {this is the part that is directly deliverd by the liver to the duodenum}
  liquid[6,52] := liquid[6,52]+ gall_bladder[52]*hulp + (1-filling_rate_fraction)*bile_phospholipids_rate*dt;         {this is the part that is directly deliverd by the liver to the duodenum}
  conjugated_bile_salt_secretion_rate:=  gall_bladder[51]*hulp/dt +(1-filling_rate_fraction)*bile_salt_rate;

  liquid[6,51] := liquid[6,51]+ conjugated_bile_salt_secretion_rate*dt;
  total_conjugated_bile_salt_secretion := total_conjugated_bile_salt_secretion + conjugated_bile_salt_secretion_rate*dt;   {integral over secretion}
  total_bile_secretion_rate := (gall_bladder_emptying_rate+(1-filling_rate_fraction)*liver_secretion_rate);                {integral over secretion}
  total_bile_secreted := total_bile_secreted+ total_bile_secretion_rate*dt ;                                               {integral over secretion}
  Volumina;
  conjugated_Bile_salt_in_secreted_bile := conjugated_bile_salt_secretion_rate/(gall_bladder_emptying_rate+(1-filling_rate_fraction)*liver_secretion_rate);
  duodenal_conjugated_bile_salt_concentration := (liquid[6,51]+ liquid[6,54])/Volume[6] ;


  {pancreatic_secretion}                 {this part should be improved}
  richt_pH := 6;
  pancreatic_bicarbonate_secretion_rate := 1E-8*secretin;
  pancreatic_enzyme_secretion_rate {in Units} := 10*0.5*(totalrateCCKrelease/(1+totalrateCCKrelease)+(Volume[5]+Volume[6])/(100+Volume[5]+Volume[6])); {guess: maximum at about 1.1 ml/minute}{stimulated by CCK production} {should amount to approximately 500 ml pancreatic fluid per day}

  if MMC_status = 3 then  pancreatic_enzyme_secretion_rate := pancreatic_enzyme_secretion_rate + 1;

  liquid[6,80] := liquid[6,80]- pancreatic_bicarbonate_secretion_rate*dt; {removal H+ by bicarbonate}
  liquid[6,50] := liquid[6,50] + 0.1*pancreatic_enzyme_secretion_rate*dt; {pancreatic lipase}

  liquid[6,56] := liquid[6,56] + 0.1*pancreatic_enzyme_secretion_rate*dt; {pancreatische phospholipase}
  liquid[6,58] := liquid[6,58] + 20*75*0.1*pancreatic_enzyme_secretion_rate*dt; {pancreatic amylase in units; By Infogest this was estimated to be equal to the value for saliva}

  pancreatic_fluid_secretion_rate := 0.02 +  0.1{0.01}*pancreatic_enzyme_secretion_rate + 10*pancreatic_bicarbonate_secretion_rate;
  liquid[6,41] := liquid[6,41] + 1E-1*pancreatic_fluid_secretion_rate*dt; {guess: 10% pancreatic proteases in pancreatic fluid}


  liquid[6,1] := liquid[6,1]+ 1*pancreatic_fluid_secretion_rate*dt; {water in pancreatic secrete added to duodenum}
  liquid[6,55] := liquid[6,55] + 0.01*pancreatic_fluid_secretion_rate*dt; {guess: 1% mucus protein in pancreatic fluid}
  liquid[6,60] := liquid[6,60]+ 1E-3*pancreatic_fluid_secretion_rate*dt ;{CCK signalling protein; it is guessed that is present at 0.1% in pancreatic fluid}
  total_pancreatic_fluid_secretion:= total_pancreatic_fluid_secretion +pancreatic_fluid_secretion_rate*dt;

end;

Procedure intestinal_secretion;
var intestinal_fluids_secretion: segment;
    richt_pH, bicarbonate_secretion, mucus_added,vertrager: real;

begin
 {secretion in the first part of the duodenum: Brunner glands}
           vertrager := 0;
           richt_pH := 7;
           bicarbonate_secretion := dt*(vertrager*bicarbonate_secretion+5E-6*secretin)/(vertrager+1); {stimulated by secretin}
           liquid[5,80] := liquid[5,80]-bicarbonate_secretion;
           liquid[5,1] := liquid[5,1]+1000*bicarbonate_secretion+0.01*dt;   {1M bicarbonaat in water?}
  {****}          mucus_added:= 10*bicarbonate_secretion;   {1% mucus in Brunner gland secrete?}
                  mucus_added:= 0.01*dt;                        {0.01 g/min mucus secretion as basis?}
           liquid[i,55]:= liquid[i,55] + mucus_added;
           liquid[5,80] := liquid[5,80] + Mucin_hydrogen(richt_pH, mucus_added ); {makes the secreted mucus pH-neutral}
 {general intestinal secretion}
  For i := 6 to 15 do
           begin
            richt_pH := 8;
            bicarbonate_secretion := 1E-7*exp(richt_pH-pH[i])*dt;
            liquid[i,80] := liquid[i,80]-bicarbonate_secretion;     {********************}
            liquid[i,1] := liquid[i,1] + 100*bicarbonate_secretion+ 0.01*dt;      {guess: 0.1 ml/minute/meter intestinal fluid secretion in rest containing 0.1 M Bicarbonate}
                   mucus_added := 0.1*bicarbonate_secretion; {guess: 0.1% mucus in Brunner gland secrete?}
                   mucus_added := mucus_added + 0.01*dt;         {guess: 0.01 g/min mucus secretion as basis?}
            liquid[i,55]:= liquid[i,55] + mucus_added;
            liquid[i,80] := liquid[i,80] + Mucin_hydrogen(richt_pH, mucus_added ); {makes the secreted mucus pH neutral}
            end;
end;





Procedure ups;              {this procedure describes for each compartment the transfer of materials toward the next compartment}
var  Boost, brake:real;
const pi = 3.14159; max_mouth_volume_intake = 10; max_mouth_volume_swallow = 1.5;   {max_mouth_volume_intake is the maximum volume of food in the mouth; max_mouth_volume_swallow determines the volume of saliva or free liquid that triggers a swallowing reflex}
begin
  CCKtrans:= 1; {this is roughly the tipping point in the CCK blood level above whcih gastric tone rapidly declines and gastric emptying becomes prodominantly calcorically controlled}
  Transition:= CCK/CCKtrans; {Transition:= SQR(SQR(Transition));} Transition:=Transition/(1+Transition);

{transport from dishes}
  for dishnumber := 1 to 5 do
   begin
        if (volumedish[dishnumber] > eating_speed_dish[dishnumber]*dt) then updish[dishnumber]:= eating_speed_dish[dishnumber] else updish[dishnumber] := 0.9*Volumedish[dishnumber]/dt; {ml/min}   {speed of eating}
        if  volumedish[dishnumber] <1E-9 then updish[dishnumber]:= 0;
        if volume[2]> max_mouth_volume_intake then updish[dishnumber] := 0; {avoids mouth congestion; no more food intake with a full mouth}
   end;
{calculation of mastication time and moment of swallowing}
  solid_volume := 0;large_solids:= 0; surface_solids := 0;
  for compoundnumber := 1 to number_of_compounds do
       begin
          if solid[compoundnumber] > 5 then large_solids := large_solids + compound[2, compoundnumber]; {the massa of very large food chunks in the oral food bolus}
          if solid[compoundnumber] > 0 then solid_volume :=  solid_volume +compound[2, compoundnumber];   {the total mass of solids in the oral food bolus}
          surface_solids := surface_solids + compound[2,compoundnumber]*pi*SQR(solid[compoundnumber]/20);    {the total surface area of the oral food bolus in cm2; division by "20" to convert mm diameter to cm radius}
       end;
  up[2]:=0;
  if volume[2] > solid_volume+0.3{0.6}*surface_solids + max_mouth_volume_swallow then up[2] := volume[2]/dt;   {swallowing criterion: sufficient moistening of solid particles; a wetting layer with arbitrairy thickness of 3 mm covers the surface of the particles, but in addition the particles can also absorb water; or alternatively a maximum of 1.5 ml saliva or liquid in the mouth}
  {if (volume[2] > swallow_speed*dt) then up[2]:= swallow_speed} {else up[2] :=  1*volume[2]/dt;} {ml/min} {speed of swallowing, not more is swallowed than originally present in the mouth}

{all backward and forward transport processes in the stomach}
  gastric_redistribution;

{all transport in duodenum}
      BasisSITransitRate:= 0.05; {fraction of the voume per minute}
      up[5] := 5*BasisSITransitRate*Volume[5]; up[6] := 3*BasisSITransitRate*Volume[6];{short duodenum}

{all transport in jejunum}
      BasisSITransitRate:= 0.05; {fractie of the volume per minuut}
      for i := 7 to 9 do
          begin
          Boost := 1*Volume[i]/(Max/2+Volume[i]);
          up[i]:= BasisSITransitRate*Volume[i]*(1+Boost);
          end;

{all transport in ileum}
      for i:= 10 to 15 do
          begin
          brake:=  (1-0.5*CCKattenuationSI*CCK/(1+CCK))*(1-0.5*GIPattenuationSI*GIP/(2+GIP))*(1-0.5*PYYattenuationSI*PYY/(10+PYY))*(1-0.5*GLP1attenuationSI*GLP1/(1+GLP1));
          Boost := 1*Volume[i]/(Max/2+Volume[i]) ;
          up[i] := BasisSITransitRate*Volume[i]*(1+Boost)*brake;    {speed of ileal transport toward the subsequent compartment}
          end;
      if MMC_status = 3 then for i:= 5 to 15 do up[i] := 1+up[i]*10;
 end;


{next functions describe the pH-depence of enzymes}
function pepsin_units(var pH, pepsin:real):real;
var pepsin_act :real;
begin
  If pH < 1 then pepsin_act := 0 + 2*(pH-0.5)*(0.1-0) else
  if pH < 1.5 then pepsin_act := 0.1 + 2*(pH-1)*(0.66-0.1) else
  if pH < 2 then pepsin_act := 0.66 + 2*(pH-1.5)*(1-0.66) else
  if pH < 2.5 then pepsin_act := 1 + 2*(pH-2)*(0.82-1) else
  if pH < 3 then pepsin_act := 0.82 + 2*(pH-2.5)*(0.75-0.82) else
  if pH < 3.5 then pepsin_act := 0.75 + 2*(pH-3)*(0.73 -0.75) else
  if pH < 4 then pepsin_act := 0.73 + 2*(pH-3.5)*(0.7-0.73) else
  if pH < 4.5 then pepsin_act := 0.7 + 2*(pH-4)*(0.65-0.7) else
  if pH < 5 then pepsin_act := 0.65 + 2*(pH-4.5)*(0.4-0.65) else
  if pH < 5.5 then pepsin_act := 0.4 + 2*(pH-5)*(0.06- 0.4) else
  if pH < 6 then pepsin_act := 0.06 + 2*(pH-5.5)*(0.03-0.06) else
  if pH < 6.5 then pepsin_act := 0.03 + 2*(pH-6)*(0.01-0.03) else
  if pH < 7 then pepsin_act := 0.01 + 2*(pH-6.5)*(0.005-0.01) else
  if pH < 7.5 then pepsin_act := 0.005 + 2*(pH-7)*(0-0.005) else
   pepsin_act := 0;
  pepsin_units:=pepsin_act*pepsin*4E6;
end;

function lingual_lipase_units(var pH, lingual_lipase:real):real;              {currently set equal to pepsin; it is doubted if lingual lipase actually exists}
var lingual_lipase_act :real;
begin
  If pH < 1 then lingual_lipase_act := 0 + 2*(pH-0.5)*(0.1-0) else
  if pH < 1.5 then lingual_lipase_act := 0.1 + 2*(pH-1)*(0.66-0.1) else
  if pH < 2 then lingual_lipase_act := 0.66 + 2*(pH-1.5)*(1-0.66) else
  if pH < 2.5 then lingual_lipase_act := 1 + 2*(pH-2)*(0.82-1) else
  if pH < 3 then lingual_lipase_act := 0.82 + 2*(pH-2.5)*(0.75-0.82) else
  if pH < 3.5 then lingual_lipase_act := 0.75 + 2*(pH-3)*(0.73 -0.75) else
  if pH < 4 then lingual_lipase_act := 0.73 + 2*(pH-3.5)*(0.7-0.73) else
  if pH < 4.5 then lingual_lipase_act := 0.7 + 2*(pH-4)*(0.65-0.7) else
  if pH < 5 then lingual_lipase_act := 0.65 + 2*(pH-4.5)*(0.4-0.65) else
  if pH < 5.5 then lingual_lipase_act := 0.4 + 2*(pH-5)*(0.06- 0.4) else
  if pH < 6 then lingual_lipase_act := 0.06 + 2*(pH-5.5)*(0.03-0.06) else
  if pH < 6.5 then lingual_lipase_act := 0.03 + 2*(pH-6)*(0.01-0.03) else
  if pH < 7 then lingual_lipase_act := 0.01 + 2*(pH-6.5)*(0.005-0.01) else
  if pH < 7.5 then lingual_lipase_act := 0.005 + 2*(pH-7)*(0-0.005) else
   lingual_lipase_act := 0;
  lingual_lipase_units:=lingual_lipase_act*lingual_lipase*4E6;
end;

function gastric_lipase_units(var pH, gastric_lipase:real):real;              {currently set equal to pepsin}
var gastric_lipase_act :real;
begin
  If pH < 1 then gastric_lipase_act := 0 + 2*(pH-0.5)*(0.1-0) else
  if pH < 1.5 then gastric_lipase_act := 0.1 + 2*(pH-1)*(0.66-0.1) else
  if pH < 2 then gastric_lipase_act := 0.66 + 2*(pH-1.5)*(1-0.66) else
  if pH < 2.5 then gastric_lipase_act := 1 + 2*(pH-2)*(0.82-1) else
  if pH < 3 then gastric_lipase_act := 0.82 + 2*(pH-2.5)*(0.75-0.82) else
  if pH < 3.5 then gastric_lipase_act := 0.75 + 2*(pH-3)*(0.73 -0.75) else
  if pH < 4 then gastric_lipase_act := 0.73 + 2*(pH-3.5)*(0.7-0.73) else
  if pH < 4.5 then gastric_lipase_act := 0.7 + 2*(pH-4)*(0.65-0.7) else
  if pH < 5 then gastric_lipase_act := 0.65 + 2*(pH-4.5)*(0.4-0.65) else
  if pH < 5.5 then gastric_lipase_act := 0.4 + 2*(pH-5)*(0.06- 0.4) else
  if pH < 6 then gastric_lipase_act := 0.06 + 2*(pH-5.5)*(0.03-0.06) else
  if pH < 6.5 then gastric_lipase_act := 0.03 + 2*(pH-6)*(0.01-0.03) else
  if pH < 7 then gastric_lipase_act := 0.01 + 2*(pH-6.5)*(0.005-0.01) else
  if pH < 7.5 then gastric_lipase_act := 0.005 + 2*(pH-7)*(0-0.005) else
   gastric_lipase_act := 0;
  gastric_lipase_units:=gastric_lipase_act*gastric_lipase*4E6;
end;


procedure digestion;       {this procedure describes the digestive transformation of each compound into degested species}
var change: real;
    linguallipolysis_rate, gastriclipolysis_rate, pancreaticlipolysis_rate: component;
    salivary_amylase_inhibition, pancreatic_amylase_inhibition:component;
const pepsin_susceptibility_casein =0.01;
      pepsin_susceptibility_protein = 0.001;
      pancreatic_protease_susceptibility_casein = 1;
      pancreatic_protease_susceptibility_protein = 1;
      activity_pancreatic_protease = 1;
      pepsin_susceptibility_casein_fragments = 0.001;
      pancreatic_protease_susceptibility_casein_fragments = 1;
      salivary_amylase_susceptibility_HDstarch = 1;
      pancreatic_amylase_susceptibility_HDstarch = 1;
      salivary_amylase_susceptibility_SDstarch = 0.2;
      pancreatic_amylase_susceptibility_SDstarch = 0.2;
      salivary_amylase_susceptibility_Rstarch = 0.0;
      pancreatic_amylase_susceptibility_Rstarch = 0.0;

      AU =  1000*3.333E-4;  {1 amylase unit = 1 mg maltose/3min, based on the Infogest protocol}      {still callibrating the starch digestion parameters}
       K_ms = 5E-6 {g/ml}; {Estimated Michaelis Menten constant for salivary amylase}
       K_cs = 5E-6 {g/ml}; {competitive inibitor constant for salivary amylase}                        {50E-6 g/ml is the IC50 value for the competitive inhibitor acarbose}
       K_as = 5E-6 {g/ml}; {acompetitive inibitor constant for salivary amylase}
       K_ncs = 5E-6 {g/ml};{noncompetitive inibitor constant for salivary amylase}
       K_mp = 5E-6 {g/ml}; {Michaelis Menten constant for pancreatic amylase}
       K_cp = 5E-6 {g/ml}; {competitive inibitor constant for pancreatic amylase}
       K_ap = 5E-6 {g/ml}; {acompetitive inibitor constant for pancreatic amylase}
       K_ncp = 5E-6 {g/ml};{noncompetitive inibitor constant for pancreatic amylase}

begin
  total_linguallipaselipolysis_rate :=0;
  total_gastriclipaselipolysis_rate :=0;
  total_pancreaticlipaselipolysis_rate := 0;
  for i:= 2 to 15 do         {lumen of mouth, stomach and small intestine}
  begin

   {fat digestion}
   change := LinguallipolysisCST*liquid[i,48]*liquid[i,2]*dt/SQR(Volume[i]);   {48 = lingual lipase}  {1 TG is hydrolysed into 1 diglyceride + 1 FA}
   {if minutes> 0 then begin write(i, liquid[i,2], change); readln; end; }
          liquid[i,7] := liquid[i,7]+ 0.666666*change;                {7 = diglyceride }
          liquid[i,20] := liquid[i,20] + 0.333333*change;            {20 = fatty acid}
          liquid[i,2]  := liquid[i,2]-change;                        {2 = triglyceride}
        if liquid[i,2] > 1E-10 then linguallipolysis_rate[i] := 100*change/liquid[i,2];
        total_linguallipaselipolysis_rate := total_linguallipaselipolysis_rate + linguallipolysis_rate[i] ;

   change := LinguallipolysisCST*liquid[i,48]*liquid[i,7]*dt/SQR(Volume[i]);   {50 = pancreatic lipase}   {1 Diglyceride is hydrolysed into 1 monoglyceride + 1 FA}
   change := 0;
          liquid[i,20] := liquid[i,20]+ 0.5*change;                                   {8 = monoglyceride}
          liquid[i,8] := liquid[i,8] +0.5*change;
          liquid[i,7]  := liquid[i,7]-change;

   change := LinguallipolysisCST*liquid[i,48]*liquid[i,8]*dt/SQR(Volume[i]);      {50 = pancreatic lipase}   {1 monoglyceride is hydrolysed into 1 FA; glycerol negelcted }
   change := 0;
          liquid[i,20] := liquid[i,20]+ change;                                   {8 = monoglyceride}
          liquid[i,7]  := liquid[i,8]-change;

   change := GastriclipolysisCST*liquid[i,49]*liquid[i,2]*dt/SQR(Volume[i]);   {49 = gastric lipase}  {1 TG  is hydrolysed into 1 diglyceride + 1 FA}
          liquid[i,7] := liquid[i,7]+ 0.666666*change;                {7 = diglyceride }
          liquid[i,20] := liquid[i,20] + 0.333333*change;            {20 = fatty acid}
          liquid[i,2]  := liquid[i,2]-change;                        {2 = triglyceride}
        if liquid[i,2] > 1E-10 then gastriclipolysis_rate[i] := 100*change/liquid[i,2];
        total_gastriclipaselipolysis_rate := total_gastriclipaselipolysis_rate + gastriclipolysis_rate[i] ;

   change := PancreaticlipolysisCST*liquid[i,50]*liquid[i,2]*dt/SQR(Volume[i]);      {50 = pancreatic lipase}   {1 TG is hydrolysed into 1 diglyceride  + 1 FA}
          liquid[i,7] := liquid[i,7] +0.666666*change;
          liquid[i,20] := liquid[i,20]+ 0.333333*change;
          liquid[i,2]  := liquid[i,2]-change;
        if liquid[i,2]> 1E-10 then pancreaticlipolysis_rate[i] := 100*change/liquid[i,2];
        total_pancreaticlipaselipolysis_rate := total_pancreaticlipaselipolysis_rate + pancreaticlipolysis_rate[i];

   change := PancreaticlipolysisCST*liquid[i,50]*liquid[i,7]*dt/SQR(Volume[i]);      {50 = pancreatic lipase}   {1 Diglyceride is hydrolysed into 1 monoglyceride + 1 FA}
          liquid[i,20] := liquid[i,20]+ 0.5*change;                                   {8 = monoglyceride}
          liquid[i,8] := liquid[i,8] +0.5*change;
          liquid[i,7]  := liquid[i,7]-change;

   change := 10000*liquid[i,20]*liquid[i,51]/(0.1+liquid[i,51])*dt/SQR(Volume[i]);        {formation of bile micelles with fatty acids}
          liquid[i,54] := liquid[i,54] + change;
          liquid[i,20] := liquid[i,20]- change;
          liquid[i,51] := liquid[i,51]- change;

   change := 10000*liquid[i,8]*liquid[i,51]/(0.1+liquid[i,51])*dt/SQR(Volume[i]);        {formation of bile micellen with monoglycerides}
          liquid[i,54] := liquid[i,54] + change;
          liquid[i,8] := liquid[i,8]- change;
          liquid[i,51] := liquid[i,51]- change;

   {phospholipid digestion}
   change := 0.1*liquid[i,52]*liquid[i,56]*dt/SQR(Volume[i]);
          liquid[i, 52] := liquid[i,52]- change;   {digestion gal fosfolipiden}
          liquid[i,53] := liquid[i,53] + change; {lysofosfolipiden by hydrolysis of bile fospfolipids}
          liquid[i,20] := liquid[i,20] + change; {FA formed by hydrolysis of bile fosfolipids}

   {protein digestion} {4 = protein,  61 = casein, 62 = caseinpep, 63 = caseinpan, 64 = caseinpeptides, 65 = small peptides}
   segmentpH := pH[i];
   change := liquid[i,4]*pepsin_susceptibility_protein*pepsin_units(segmentpH,liquid[i,40])*dt/(volume[i]);
          liquid[i,65] := liquid[i,65] + change; {digestion protein by pepsin into small peptides}
          liquid[i,4] := liquid[i,4] - change;

   change := liquid[i,66]*pepsin_susceptibility_protein*pepsin_units(segmentpH,liquid[i,40])*dt/(volume[i]);
          liquid[i,65] := liquid[i,65] + change; {digestion BSA by pepsin into small peptides}
          liquid[i,66] := liquid[i,66] - change;

   change := liquid[i,67]*pepsin_susceptibility_protein*pepsin_units(segmentpH,liquid[i,40])*dt/(volume[i]);
          liquid[i,65] := liquid[i,65] + change; {digestion BLG by pepsin into small protein peptides}
          liquid[i,67] := liquid[i,67] - change;

   change := liquid[i,61]*pepsin_susceptibility_casein*pepsin_units(segmentpH,liquid[i,40])*dt/(volume[i]);
          if change > liquid [i,61] then change := 0;
          liquid[i,62] := liquid[i,62] + change; {digestion Casein by pepsin into casein peptides}
          liquid[i,61] := liquid[i,61] - change;

   change:= liquid[i,4]*liquid[i,41]*pancreatic_protease_susceptibility_protein*activity_pancreatic_protease*dt/Volume[i];
            if change > liquid [i,4] then change := 0;
            liquid[i,65] := liquid[i,65] + change; {digestion protein by pancreatic protease into small peptides}
            liquid[i,4] := liquid[i,4] - change;

   change:= liquid[i,66]*liquid[i,41]*pancreatic_protease_susceptibility_protein*activity_pancreatic_protease*dt/Volume[i];
            if change > liquid [i,66] then change := 0;
            liquid[i,65] := liquid[i,65] + change; {digestion BSA by pancreatic protease into small protein peptides}
            liquid[i,66] := liquid[i,66] - change;

   change:= liquid[i,67]*liquid[i,41]*pancreatic_protease_susceptibility_protein*activity_pancreatic_protease*dt/Volume[i];
            if change > liquid [i,67] then change := 0;
            liquid[i,65] := liquid[i,65] + change; {digestion BLG by pancreatic protease into small proteien peptides}
            liquid[i,67] := liquid[i,67] - change;

   change:= liquid[i,61]*liquid[i,41]*pancreatic_protease_susceptibility_casein*activity_pancreatic_protease*dt/Volume[i];
            if change > liquid [i,61] then change := 0;
            liquid[i,63] := liquid[i,63] + change; {digestion Casein by pancreatic protease into pancreatic casein digest}
            liquid[i,61] := liquid[i,61] - change;

   change := liquid[i,62]*pepsin_susceptibility_casein_fragments*pepsin_units(segmentpH,liquid[i,40])*dt/volume[i];
          if change > liquid [i,62] then change := 0;
          liquid[i,64] := liquid[i,64] + change;  {digestion Casein pepsic digest by pepsin into casein small peptides}
          liquid[i,62] := liquid[i,62] - change;
   change := liquid[i,63]*liquid[i,41]*pancreatic_protease_susceptibility_casein_fragments*activity_pancreatic_protease*dt/Volume[i];
          if change > liquid [i,63] then change := 0;
          liquid[i,64] := liquid[i,64] + change;  {digestion Casein pancreatic digest by pancreatic protease into casein small peptides}
          liquid[i,63] := liquid[i,63] - change;

   liquid[i,65] := liquid [i,65] + liquid[i,64];  liquid[i,64] := 0; {all small casein peptides are seen as small protein peptides]};


   {pharmaceutical degradation}   {currently turned of}
          {if pH[i] <4 then change := liquid[i,5]*0.1*(exp(-pH[i])-exp(-4))*dt Else change := 0;
          liquid[i,5] := liquid[i,5]-change;
          liquid[i,6] := liquid[i,6]+change;}

   {starch digestion}
           {72 = I_c      competitive alpa amylase inhibitor
            73 = I_a      acompetitive alpha amylase inhibitor
            74 = I_nc     non competitive alpha amylase inhibitor}
   liquid[i,57]:= liquid[i,57] -dt*liquid[i,57]*exp(5*ln(5/(5+pH[i])));    {degradation of salivary amylase at low pH}
   change := AU*liquid[i,10]*salivary_amylase_susceptibility_HDstarch*liquid[i,57]*dt/(((K_ms*(1+concentration_lumen[i,72]/K_cs) +concentration_lumen[i,10]*(1+concentration_lumen[i,73]/K_as+concentration_lumen[i,3]/0.1+concentration_lumen[i,23]/0.035)*(1+concentration_lumen[i,74]/K_ncs))));
          change := change + 20*AU*liquid[i,10]*pancreatic_amylase_susceptibility_HDstarch*liquid[i,58]*dt/(((K_mp*(1+concentration_lumen[i,72]/K_cp )+concentration_lumen[i,10]*(1+(concentration_lumen[i,73]/K_ap+concentration_lumen[i,3]/0.1+concentration_lumen[i,23]/0.035 ))*(1+concentration_lumen[i,74]/K_ncp ))));
        if change> liquid[i,10] then change:= liquid[i,10];
          liquid[i,23] := liquid[i,23] + 0.99*change;   {maltose is released from starch digestion by amylases}
          liquid[i,24] := liquid[i,24] + 0.01*change;   {13C maltose}
          liquid[i,10] := liquid[i,10] - change;
   change := AU*liquid[i,11]*salivary_amylase_susceptibility_SDstarch*liquid[i,57]*dt/(((K_ms*(1+(concentration_lumen[i,72])/K_cs )+concentration_lumen[i,11]*(1+(concentration_lumen[i,73])/K_as+concentration_lumen[i,3]/0.1+concentration_lumen[i,23]/0.035 ))*(1+concentration_lumen[i,74]/K_ncs )));
          change := change + 20*AU*liquid[i,11]*pancreatic_amylase_susceptibility_SDstarch*liquid[i,58]*dt/(((K_mp*(1+(concentration_lumen[i,72])/K_cp )+concentration_lumen[i,11]*(1+(concentration_lumen[i,73])/K_ap+concentration_lumen[i,3]/0.1+concentration_lumen[i,23]/0.035 ))*(1+concentration_lumen[i,74]/K_ncp )));
          if change> liquid[i,11] then change:= liquid[i,11];
          liquid[i,23] := liquid[i,23] + 0.99*change;   {maltose is released from starch digestion by amylases}
          liquid[i,24] := liquid[i,24] + 0.01*change;   {13C maltose}
          liquid[i,11] := liquid[i,11] - change;
   change := AU*liquid[i,12]*salivary_amylase_susceptibility_Rstarch*liquid[i,57]*dt/(((K_ms*(1+(concentration_lumen[i,72])/K_cs )+concentration_lumen[i,12]*(1+(concentration_lumen[i,73])/K_as+concentration_lumen[i,3]/0.1+concentration_lumen[i,23]/0.035 ))*(1+concentration_lumen[i,74]/K_ncs )));  {0.001 is Km in g/ml}
          change := change + 20*AU*liquid[i,12]*pancreatic_amylase_susceptibility_Rstarch*liquid[i,58]*dt/(((K_mp*(1+(concentration_lumen[i,72])/K_cp )+concentration_lumen[i,12]*(1+(concentration_lumen[i,73])/K_ap+concentration_lumen[i,3]/0.1+concentration_lumen[i,23]/0.035 ))*(1+concentration_lumen[i,74]/K_ncp )));  {0.001 is Km in g/ml}
          if change> liquid[i,12] then change:= liquid[i,12];
          liquid[i,23] := liquid[i,23] + 0.99*change;   {maltose is released from starch digestion by amylases}
          liquid[i,24] := liquid[i,24] + 0.01*change;   {13C maltose}
          liquid[i,12] := liquid[i,12] - change;
   {Volumina;}
  end;
end;


procedure decomposition;        {decomposition of compounds in other compounds and subsequently in components}
var hulp, increase:real;
    degcomp : integer;
begin
 for compoundnumber:= 1 to number_of_compounds do
 begin
    {corpus, antrum}
     for i := 3 to 4 do
      begin
         increase :=  compound[i,compoundnumber]*dt/gastric_dec_time[compoundnumber];
         hulp :=  compound[i,compoundnumber]-increase;
          degcomp := gastric_dec_comp[compoundnumber,1];
             compound[i,degcomp]:= compound[i,degcomp]+ increase*gastric_dec_frac[compoundnumber,1]/100;
          degcomp := gastric_dec_comp[compoundnumber,2];
             compound[i,degcomp]:= compound[i,degcomp]+ increase*gastric_dec_frac[compoundnumber,2]/100;
      for j:= 1 to comp do liquid[i,j]:= liquid[i,j]+ 0.01*(100-gastric_dec_frac[compoundnumber,1]-gastric_dec_frac[compoundnumber,2])*compoundcomposition[compoundnumber,j]*increase;
      compound[i,compoundnumber] := hulp;
      end;

    {mouth} i := 2;
      increase :=  compound[i,compoundnumber]*dt/oral_dec_time[compoundnumber];
         hulp :=  compound[i,compoundnumber]-increase;
           degcomp := oral_dec_comp[compoundnumber,1];
              compound[i,degcomp]:= compound[i,degcomp]+ increase*oral_dec_frac[compoundnumber,1]/100;
           degcomp := oral_dec_comp[compoundnumber,2];
              compound[i,degcomp]:= compound[i,degcomp]+ increase*oral_dec_frac[compoundnumber,2]/100;
      for j:= 1 to comp do liquid[i,j]:= liquid[i,j]+ 0.01*(100-oral_dec_frac[compoundnumber,1]-oral_dec_frac[compoundnumber,2])*compoundcomposition[compoundnumber,j]*increase;
      compound[i,compoundnumber] := hulp;

    {fundus}
     increase :=  compoundfundus[compoundnumber]*dt/gastric_dec_time[compoundnumber];
      hulp :=  compoundfundus[compoundnumber] - increase;
           degcomp := gastric_dec_comp[compoundnumber,1];
              compoundfundus[degcomp]:= compoundfundus[degcomp]+ increase*gastric_dec_frac[compoundnumber,1]/100;
           degcomp := gastric_dec_comp[compoundnumber,2];
           compoundfundus[degcomp]:= compoundfundus[degcomp]+ increase*gastric_dec_frac[compoundnumber,2]/100;
      for j:= 1 to comp do fundus[j]:= fundus[j]+ 0.01*(100-gastric_dec_frac[compoundnumber,1]-gastric_dec_frac[compoundnumber,2])*compoundcomposition[compoundnumber,j]*increase;
      compoundfundus[compoundnumber] := hulp;

    {small intestine}
      for i := 5 to 15 do
      begin
         increase :=  compound[i,compoundnumber]*dt/intestinal_dec_time[compoundnumber];
         hulp :=  compound[i,compoundnumber]-increase;
           degcomp := intestinal_dec_comp[compoundnumber,1];
              compound[i,degcomp]:= compound[i,degcomp]+ increase*intestinal_dec_frac[compoundnumber,1]/100;
           degcomp := intestinal_dec_comp[compoundnumber,2];
              compound[i,degcomp]:= compound[i,degcomp]+ increase*intestinal_dec_frac[compoundnumber,2]/100;
       for j:= 1 to comp do liquid[i,j]:= liquid[i,j]+ 0.01*(100-intestinal_dec_frac[compoundnumber,1]-intestinal_dec_frac[compoundnumber,2])*compoundcomposition[compoundnumber,j]*increase;
       compound[i,compoundnumber] := hulp;
      end;
   end;
end;

Procedure absorption_in_Small_Intestine;
var totaltargetabsorption_rate, delay, change: real;
    const R=0.8; {recovery fraction of glucose after absorption by enterocytes; part is used to fuel the enterocytes}

   begin
     totalcaloricabsorption_rate:= 0;
     G0in := 0; G1in:= 0 ; G2in:= 0;fruc:= 0;
     for i:=5 to 15 do     {these are the sections of the small intestine }
     begin
       totaltargetabsorption_rate := 0;
       for j:= 2 to comp do     {these are all components except water }
       begin
               if (j <> 3) and (j<> 34) and (j<> 35)  then
                begin
                if (concentration_brush[i,j] > 1E-40) and (j <> 80) then absorption_rate[i,j] := VmaxAbsorption[i,j]*concentration_brush[i,j]/(concentration_brush[i,j]+KmAbsorption[i,j])
                    else absorption_rate[i,j] :=0 ;
                end
              else
                begin
                   hulp:= concentration_brush[i,3]+concentration_brush[i,34]+concentration_brush[i,35];
                   if(hulp > 0)then
                   begin
                    absorption_rate[i,3] :=  VmaxAbsorption[i,3]*hulp/(hulp+KmAbsorption[i,3])+0.5*hulp/(hulp+0.1*0.0225560) ; {SGLT-1 channel };
                    absorption_rate[i,34]:= 0;               {all sorts of glucoses collectively sort as glucose (component 3) for the Weber calculation}
                    absorption_rate[i,35]:= 0;
                   end
                   else
                   begin
                      absorption_rate[i,3]:= 0;
                      absorption_rate[i,34]:= 0;
                      absorption_rate[i,35]:= 0;
                   end;
                end;
      {Weber}          {Weber equation is based on experimental evaluation on competition between macronutrients in mini Troll pigs and recalibrated for human adults by van Aken}
       totaltargetabsorption_rate := totaltargetabsorption_rate + absorption_rate[i,j];
            end;
       if (Maxtotalabsorption_rate[i] = 0) then for j:=2 to comp do absorption_rate[i,j]:= 0 else
             for j:=2 to comp do
                 begin
                      if caloric_value[j] <> 0 then absorption_rate[i,j]:= absorption_rate[i,j]*Maxtotalabsorption_rate[i]/(totaltargetabsorption_rate+Maxtotalabsorption_rate[i]);
                 end;
       {end Weber}

       If hulp <>0 then
         begin
         glucose_absorption_rate[i] := absorption_rate[i,3];
         absorption_rate[i,3]:= (concentration_brush[i,3]/(hulp+1E-3))*glucose_absorption_rate[i];        {hier worden de absorptiesnelheden weer verdeeld over de glucose soorten}
         absorption_rate[i,34]:= (concentration_brush[i,34]/(hulp+1E-3))*glucose_absorption_rate[i];
         absorption_rate[i,35]:= (concentration_brush[i,35]/(hulp+1E-3))*glucose_absorption_rate[i];
         end;
         absorption_rate[i,70] := 0.15*absorption_rate[i,64];   {Leucine is a fraction 3/20 of the protein (beta lactoglobulin?)}
        absorption_rate[i,71] := absorption_rate[i,65];

        {transfer to g/min.segment for nutrients}
        for j:= 1 to comp do
            if caloric_value[j] <> 0 then absorption_rate[i,j]:= absorption_rate[i,j]/caloric_value[j];

                            {hypertropy corrects for altered absorptive surface area by elongated or shortened vili (Vili can be shorted medical treatment or elongated by fructose as stimulus}
        absorption_rate[i,1]:= 2*{hypertrophy[i]*}(Volume[i]/(Volume[i]+100))*(0.2 -osmotic_value_brush[i])+ (30*18/180.15588)*(absorption_rate[i,3]+absorption_rate[i,34]+absorption_rate[i,35]); {klopt dit al?}

     end;
     total_glucose_absorption_rate := 0;  conjugated_bile_salt_absorption_rate:= 0;
     For i := 5 to 15 do
         total_glucose_absorption_rate := total_glucose_absorption_rate+absorption_rate[i,3]+absorption_rate[i,34]+absorption_rate[i,35] ;
     total_glucose_absorption_rate:= total_glucose_absorption_rate/180.15588; {in moles/min}
         {**** absorption hydrogen ions, currently inactive }
       { hulp := liquid[i,61]+liquid[i,62]+liquid[i,63]+liquid[i,64]+liquid[i,65];
         for i:= 5 to 15 do
         absorption_rate[i,80] := -0.002*(absorption_rate[i,64]+absorption_rate[i,65]);    {peptiden wordt in neutrale vorm geabsorpbeerd}
         absorption_rate[i,80] := absorption_rate[i,1]*Ber_hydrogen(pH[i], hulp, liquid[i,66], liquid[i,55], liquid[i,67],liquid[i,4])/liquid[i,1];}
        {**** end absorption hydrogen ions }

        {calculation of total absorption rate in the small intestines}
        for j := 1 to comp do  begin totalabsorption_rate[j]:= 0; absorbed[j] := 0;  end;
        for j := 1 to comp do
          begin
          for i:= 5 to 15 do
           begin
           totalabsorption_rate[j]:= totalabsorption_rate[j]+ absorption_rate[i,j];
           totalcaloricabsorption_rate := totalcaloricabsorption_rate  + absorption_rate[i,j]*caloric_value[j];
           absorbed[j]:= absorbed[j]+absorption_rate[i,j]*dt;
           end;
           calories_absorbed := calories_absorbed + caloric_value[j]*1  *absorbed[j];
          end;
        {**********}
     {transport to blood}
     delay := 2 ;{Delay in min} delay := 1/delay;
     for i := 5 to 15 do
      begin
        For j:= 1 to comp do
         begin
         Mesentric[i,j] := Mesentric[i,j] +  absorption_rate[i,j]*dt;  {(in (g /min) transported by upper mesentric veins)}
         change := delay*Mesentric[i,j]*dt;
              if change > Mesentric[i,j] then change :=  Mesentric[i,j] ;
           Blood[j] := Blood[j] + change;
           Mesentric[i,j]:= Mesentric[i,j]-change;
           if j=3 then G0in := G0in+ delay*Mesentric[i,3]*1000; {glucose absorption rescaled to mg/min}
           if j=34 then G1in := G1in+delay*Mesentric[i,34]*1000;
           if j=35 then G2in := G2in+delay*Mesentric[i,35]*1000;
           if j=22 then fruc := fruc + delay*Mesentric[i,22]*1000;       {total blood fructose in mg}
         end;
      end;
      G0in:= R*G0in; G1in:=R*G1in; G2in:= R*G2in; {R corrects for the usage of glucose by the absorptive enterocytes}
     totalwaterabsorptionrate := totalabsorption_rate[1];
   end;




procedure Mucus_layer;
var deler,hulp:real;
const vertrager = 5;
           begin
           deler := 1/(vertrager +1);
           sucrase;
           maltase;
           volumina;
           for i:= 5 to 15 do
            begin
            hulp :=  absorption_rate[i,1]*dt;
            liquidnext[i,1] := liquid[i,1] - absorption_rate[i,1]*dt;
            for j:= 2 to comp do
             if j <> 80 then
                 begin
                 if (j <> 3) and (j<>21) and (j<> 22) and (j<>23)and (j<> 34) and (j<>35) and (j<>36) then
                     begin
                     hulp := concentration_lumen[i,j]*(1+surface_factor[i]*drag*T[i,j]*absorption_rate[i,1])-T[i,j]*surface_factor[i]*absorption_rate[i,j];
                     concentration_brush[i,j]:=  deler*(vertrager*concentration_brush[i,j] + hulp);
                     end;
                 liquidnext[i,j] := liquid[i,j] - absorption_rate[i,j]*dt;
                                         if liquidnext[i,j]< 0 then begin liquidnext[i,j]:=0; end;
                 end;
                 liquidnext[i,80] := liquid[i,80] - absorption_rate[i,80]*dt;

 {glucose}       hulp := concentration_lumen[i,3]*(1+surface_factor[i]*T[i,3]*drag*absorption_rate[i,1])+T[i,3]*surface_factor[i]*(-absorption_rate[i,3] + sucrase_transformation[i,21]*(180.15588/342.2965) + maltase_transformation[i,23]*2*(180.15588/342.2965));
                 concentration_brush[i,3]:= deler*(vertrager*concentration_brush[i,3] + hulp); {steady state calculations: this expresion dampens the fast variations}
                 liquidnext[i,3] := liquid[i,3]+dt*(-absorption_rate[i,3]+ sucrase_transformation[i,21]*(180.15588/342.2965) + maltase_transformation[i,23]*2*(180.15588/342.2965));
 {glucose1}      hulp := concentration_lumen[i,34]*(1+surface_factor[i]*T[i,34]*drag*absorption_rate[i,1])+T[i,34]*surface_factor[i]*(-absorption_rate[i,34] + sucrase_transformation[i,36]*(180.15588/342.2965) + maltase_transformation[i,24]*2*(180.15588/342.2965));
                 concentration_brush[i,34]:= deler*(vertrager*concentration_brush[i,34] + hulp); {steady state calculations: this expresion dampens the fast variations}
                 liquidnext[i,34] := liquid[i,34]+dt*(-absorption_rate[i,34]+ sucrase_transformation[i,36]*(180.15588/342.2965) + maltase_transformation[i,24]*2*(180.15588/342.2965));
                      if liquidnext[i,34]< 1E-10 then begin liquidnext[i,34]:=1E-10; end;
 {glucose2}      hulp := concentration_lumen[i,35]*(1+surface_factor[i]*T[i,35]*drag*absorption_rate[i,1])+T[i,35]*surface_factor[i]*(-absorption_rate[i,35] {+ sucrase_transformation[i,36]*(180.15588/342.2965)} {+ maltase_transformation[i,23]*2*(180.15588/342.2965)});
                 concentration_brush[i,35]:= deler*(vertrager*concentration_brush[i,35] + hulp); {steady state calculations: this expresion dampens the fast variations}
                 {liquidnext[i,35] :=  liquid[i,35]+dt*(-absorption_rate[i,35]); }
 {sucrose}       {concentration_brush has already been calculated in the routine "brush border enzymes"}
                 liquidnext[i,21] := liquid[i,21]+dt*(-absorption_rate[i,21] -sucrase_transformation[i,21]);
 {sucrose1}      liquidnext[i,36] := liquid[i,36]+dt*(-absorption_rate[i,36] -sucrase_transformation[i,36]);
 {fructose}      hulp := liquid[i,22]+dt*(-absorption_rate[i,22]+ sucrase_transformation[i,21]*(180.15588/342.2965)+ sucrase_transformation[i,36]*(180.15588/342.2965)); if hulp<0 then absorption_rate[i,22] := 0.95*liquid[i,22]/dt + sucrase_transformation[i,21]*(180.15588/342.2965)+sucrase_transformation[i,36]*(180.15588/342.2965);
                 liquidnext[i,22] := liquid[i,22]+dt*(-absorption_rate[i,22]+ sucrase_transformation[i,21]*(180.15588/342.2965)+ sucrase_transformation[i,36]*(180.15588/342.2965));
                 hulp := concentration_lumen[i,22]*(1+surface_factor[i]*T[i,22]*drag*absorption_rate[i,1])+T[i,22]*surface_factor[i]*(-absorption_rate[i,22] + sucrase_transformation[i,21]*(180.15588/342.2965)+sucrase_transformation[i,36]*(180.15588/342.2965));
                 concentration_brush[i,22]:= deler*(vertrager*concentration_brush[i,22] + hulp); {steady state calculations: this expresion dampens the fast variations}
 {maltose}       {concentration_brush has already been calculated in the routine "brush border enzymes"}
                 liquidnext[i,23] := liquid[i,23]+dt*(-absorption_rate[i,23]-maltase_transformation[i,23]);
                 liquidnext[i,24] := liquid[i,24]+dt*(-absorption_rate[i,24]-maltase_transformation[i,24]);
 {bile salts}    liquidnext[i,51] := liquid[i,51] + absorption_rate[i,54]*dt; {conjugated bile salts are not absorbed, here a 1:1 ratio between lipids and bile salts is assumed for the micelles}
                 liquidnext[i,54] := liquid[i,54] - absorption_rate[i,54]*dt;

                 for j:= 1{2} to comp do    {all changes for all components exept water}
                 begin
                  liquid[i,j]:= liquidnext[i,j] ;
                 end;
           end;
           volumina;

        end;

procedure Colon;

    begin
    liquid[16,59]:= liquid[16,59]+liquid[16,51]; liquid[16,51] := 0; {immediate and full deconjugation by microbial pool}
    liquid[16,59]:= liquid[16,59]-0.01*liquid[16,59]*dt;   {reduction of deconjugated bile salts in colon by defeacation}
    end;

Procedure Gut_hormones;
var delay, change,hulp, hulp1, delay_entry:real;
    rateCCKrelease, ratePYYrelease, rateGIPrelease, rateGLP1release, rateGhrelinrelease, rateMotilinrelease, rateSecretinrelease, rateCCKproduction, ratePYYproduction, rateGIPproduction, rateGLP1production, rateGhrelinproduction, rateMotilinproduction, rateSecretinproduction :segment;

    begin
 {***}        delay := 0.5; delay:=1/delay;
 {***}        delay_entry := 0.5; delay_entry := 1/delay_entry;
         totalrateCCKrelease:= 0;
         totalrateGIPrelease:= 0;
         totalratePYYrelease := 0;
         totalrateGLP1release := 0;
         totalrateGhrelinrelease := 0;
         totalrateMotilinrelease := 0;
         totalrateSecretinrelease := 0;

         for i:= 5 to 15 do
            begin
             rateCCKproduction[i] :=  0;
             ratePYYproduction[i] :=  0;
             rateGIPproduction[i] :=  0;
             rateGLP1production[i] := 0;
             rateGhrelinproduction[i] := 0;
             rateMotilinproduction[i] := 0;
             rateSecretinproduction[i] := 0;
            end;

         for i := 5 to 15 do for j:= 1 to comp do
            begin
               concentration_cell[i,j] := concentration_cell[i,j] + (concentration_brush[i,j]-concentration_cell[i,j])*delay_entry*dt;
            end;

         for i:= 5 to 15 do
         begin
            hulp1 := exp(4.5-pH[i]);    {stimulation of secretin secretion below pH 4.5, for secretin secretion by S-cells, pH 4.5 is the tipping point below which secretin secretion starts: Williams http://dx.doi.org/10.3998/panc.2013.5 }
            rateSecretinproduction[i] := rateSecretinproduction[i]+ 2*S_cell[i]*hulp1/(1+hulp1);
            for j:= 1 to comp do
            begin
            rateCCKproduction[i] := rateCCKproduction[i]+CCKresponsefactor[j]*I_cell[i]*concentration_cell[i,j] + 0.5*proteolytic_susceptibility[j]*I_cell[i]*(liquid[i,j]);
            ratePYYproduction[i] := ratePYYproduction[i]+PYYresponsefactor[j]*L_cell[i]*concentration_cell[i,j];
            rateGIPproduction[i] := rateGIPproduction[i]+GIPresponsefactor[j]*K_cell[i]*concentration_cell[i,j]/(1+ concentration_cell[i,j]/1E-5);
                                 hulp := 1E2*concentration_cell[i,j]{/(1+concentration_brush[i,j]/1E-4)}; {receptor sensitivity}
            rateGLP1production[i] := rateGLP1production[i]+GLP1responsefactor[j]*L_cell[i]*5*hulp/(1+hulp);
            rateGhrelinproduction[i] := rateGhrelinproduction[i] +Ghrelinresponsefactor[j]*G_cell[i]*concentration_cell[i,j] ;
            rateMotilinproduction[i] := rateMotilinproduction[i]+Motilinresponsefactor[j]*Mo_cell[i]*concentration_cell[i,j];
            end;

          {coupling s between gut hormones}
             rateCCKproduction[i] := rateCCKproduction[i]*(1-0.2*(PYY/(1+PYY))); {Gomez 2013: PYY lowers CCK production}
             rateGIPproduction[i] := rateGIPproduction[i] ;   {neural coupling with GLP1 does not exist, see Little 2006}
             ratePYYproduction[i] := ratePYYproduction[i] +1*totalrateCCKrelease;    {neural coupling with CCK}
             rateGLP1production[i] := rateGLP1production[i] + 0*rateGIPproduction[i];  {neural coupling with GIP probably does not really exist in humans, and instead occurs through a different mechanism (Lim 2006)}

          {mesentric cells are a cell layer in the small intestines that are responsible for gut hormone accumulation and release. This needs to be checked}
           Mesentric_CCK[i]:= Mesentric_CCK[i]+rateCCKproduction[i]*dt;
           rateCCKrelease[i] := delay*Mesentric_CCK[i];
           change :=  rateCCKrelease[i]*dt;
           Mesentric_CCK[i] := Mesentric_CCK[i] - change;

           Mesentric_PYY[i]:= Mesentric_PYY[i]+ratePYYproduction[i]*dt;
           ratePYYrelease[i] := delay*Mesentric_PYY[i] ;
           change:= ratePYYrelease[i]*dt;
           Mesentric_PYY[i] := Mesentric_PYY[i] - change;

           Mesentric_GIP[i]:= Mesentric_GIP[i]+rateGIPproduction[i]*dt;
           rateGIPrelease[i] := delay*Mesentric_GIP[i];
           change:= rateGIPrelease[i]*dt;
           Mesentric_GIP[i] := Mesentric_GIP[i] - change;

           Mesentric_GLP1[i]:= Mesentric_GLP1[i]+rateGLP1production[i]*dt;
           rateGLP1release[i] := delay*Mesentric_GLP1[i];
           change:= rateGLP1release[i]*dt;
           Mesentric_GLP1[i] := Mesentric_GLP1[i] - change;

           Mesentric_Ghrelin[i]:= Mesentric_Ghrelin[i]+rateGhrelinproduction[i]*dt;
           rateGhrelinrelease[i] := delay*Mesentric_Ghrelin[i];
           change:= rateGhrelinrelease[i]*dt;
           Mesentric_Ghrelin[i] := Mesentric_Ghrelin[i] - change;

           Mesentric_Motilin[i]:= Mesentric_Motilin[i]+rateMotilinproduction[i]*dt;
           rateMotilinrelease[i] := delay*Mesentric_Motilin[i];
           change:= rateMotilinrelease[i]*dt;
           Mesentric_Motilin[i] := Mesentric_Motilin[i] - change;

           Mesentric_Secretin[i]:= Mesentric_Secretin[i]+rateSecretinproduction[i]*dt;
           rateSecretinrelease[i] := delay*Mesentric_Secretin[i];
           change:= rateSecretinrelease[i]*dt;
           Mesentric_Secretin[i] := Mesentric_Secretin[i] - change;



         totalrateCCKrelease:= totalrateCCKrelease +rateCCKrelease[i];
         totalrateGIPrelease:= totalrateGIPrelease +rateGIPrelease[i];
         totalratePYYrelease := totalratePYYrelease +ratePYYrelease[i];
         totalrateGLP1release := totalrateGLP1release +rateGLP1release[i];
         totalrateGhrelinrelease := totalrateGhrelinrelease +rateGhrelinrelease[i] {+ hunger*************}; {Ghrelin is stumulated by hunger, but the pre prandial effect is not yet added}
         totalrateMotilinrelease := totalrateMotilinrelease +rateMotilinrelease[i];
         totalrateSecretinrelease := totalrateSecretinrelease +rateSecretinrelease[i];

         end;
         totalrateGhrelinrelease := totalrateGhrelinrelease + 10;    {in reality this is also dependent on the anticipation of a meal; this has not yet been added}
         totalrateMotilinrelease := totalrateMotilinrelease + 10;
             CCK := CCK+totalrateCCKrelease*dt ;
             PYY := PYY +totalratePYYrelease*dt ;
             GIP := GIP +totalrateGIPrelease*dt  ;
             GLP1 := GLP1 +totalrateGLP1release*dt ;
             Ghrelin := Ghrelin +totalrateGhrelinrelease*dt*exp(-serum_insulin_mU_l/30)*exp(-CCK/2)*exp(-GLP1/2)*exp(-osmotic_value_lumen_duodenum/0.02)  ;   {osmotic effect roughly as described by Joost Overduin}{El Khoury D 2006 suggests that insuline stimulated Gastrin production}
             Motilin := Motilin +totalrateMotilinrelease*dt*exp(-serum_insulin_mU_l/20)*exp(-CCK/2)*exp(-GLP1/0.2)*exp(-(ABS(pH[6]-6)/6))  ;   {estimated on the basis of Saito 1980: Regulation of Motilin Secretion in the..., also glucagon has this effect}
             Secretin := Secretin +totalrateSecretinrelease*dt ;



         CCK := CCK-CCK*CCKdegradationconstant*dt;
         PYY := PYY- PYY*PYYdegradationconstant*dt;
         GIP := GIP- GIP*GIPdegradationconstant*dt;
         GLP1 := GLP1- GLP1*GLP1degradationconstant*dt;
         Ghrelin := Ghrelin- Ghrelin*Ghrelindegradationconstant*dt;
         Motilin := Motilin- Motilin*Motilindegradationconstant*dt;
         Secretin := Secretin- Secretin*Secretindegradationconstant*dt;

         realGLP1 := 0.5*GLP1;      {adjusted to the steady state value of duodenal infusion of 2 kcal/min glucose described by Marathe et al.: Relationships between gastric emptying, ... Diabetes care 36 may 2013}
         realGIP := 1*GIP;      {adjusted to the steady state value of duodenal infusion of 2 kcal/min glucose described by Marathe et al.: Relationships between gastric emptying, ... Diabetes care 36 may 2013}
        end;

procedure sensations;
var FFA, peptides, sugar :real;
begin
     fullness := (10/15)*(gastric_tone-2);
     hunger_totalcaloricabsorptionrate:= 10-25*{ln}(totalcaloricabsorption_rate/maxcaloricabsorptioncapacity);

     hunger_intestinal_volume := volume[5]+volume[6]+volume[7]+volume[8]+volume[9]+volume[10]+volume[11]+volume[12]+volume[13]+volume[14];
     hunger_intestinal_volume := 10-hunger_intestinal_volume/1000;

     hunger:= 0.8*(hunger_totalcaloricabsorptionrate+hunger_intestinal_volume)/2;
     if MMC_status = 3 then hunger:= 10;
end;

procedure intubation_in_antrum;
begin
 if (minutes >250) and (minutes< 251) then
  begin
  liquid[4,1]:= liquid[4,1]+ 0{250}*dt;
  liquid[4,3]:= liquid[4,3] + 0{50}*dt;
  end;
end;

procedure intubation_in_duodenum;
begin
  if (minutes >250) and (minutes< 300) then
  begin
  liquid[5,1]:= liquid[5,1]+ 0*dt;
  liquid[5,3]:= liquid[5,3] + 0*dt;
  end;
end;

procedure intravenous;
begin
 G2in := 1 ; {2H glucose administred intraveneously in mg/min}
 G2in_molpermin:= 0.001/180.15588;
end;

Procedure make_step;         {this procedure activates all processes during a time step}
var transfer_solids_to_antrum, emptying, change :real;
   vertrager : integer;
   compoundnumber: integer;

begin
 update_MMC_status;
 Volumina;
 calories;
 decomposition;
 digestion;
 osmotic_values_lumen;
 osmotic_values_brush;
 absorption_in_Small_Intestine;
 Mucus_layer;
 intubation_in_antrum;
 intubation_in_duodenum;
 intravenous;
 fructose_to_glucose;
 insulin_dynamics;
 calc_HbA1c;
 Gut_hormones;
 {hormones;}
 decay_in_blood;
 ber_alle_pHs;
 saliva_secretion;
 gastric_secretion;
 pancreabiliar_secretion;
 intestinal_secretion;
 caloric_densities_segments;
 solids_in_stomach;
 {gastric_compartment_viscosities; }
 colon;

 volumina;
 ups;
 sensations;
 emptying:= 0;
 vertrager:= substeps; {fast variations in gastric emptying are dampened for the graph}
 {++++++++++++++++++}for j:= 1 to comp do {begin} emptying := emptying+ up[4]*concentration_lumen[4,j]*caloric_value[j];{if (up[4]*concentration_lumen[4,j]*caloric_value[j]> 1) and (minutes >200) then  begin write(j, up[4]*concentration_lumen[4,j]*caloric_value[j]); readln; end; end; }     {gastric empting rate in calories/min}
 for compoundnumber := 1 to number_of_compounds do  emptying := emptying+ up[4]*gastric_sieve[compoundnumber]*concentration_compound[4,compoundnumber]*calories_compound[compoundnumber];
 gastric_emptying := (vertrager*gastric_emptying + emptying)/(vertrager+1);



 for j:= 1 to comp do
  begin
        {dishes}
        liquidnext[2,j] := 0;
        for dishnumber := 1 to 5 do
        begin
           liquiddishnext[dishnumber,j] := liquiddish[dishnumber,j] - updish[dishnumber]*concentration_dish[dishnumber,j]*dt ;
           liquidnext[2,j] := liquidnext[2,j] + updish[dishnumber]*concentration_dish[dishnumber,j]*dt;
        end;
        liquidnext[2,j] := liquid[2,j]+liquidnext[2,j]-up[2]*concentration_lumen[2,j]*dt;

        liquidnext[3,j] := liquid[3,j]+up[2]*concentration_lumen[2,j]*dt;
        liquidnext[4,j] := liquid[4,j];
        liquidnext[5,j] := liquid[5,j]-up[5]*concentration_lumen[5,j]*dt;

        {remaining non-gastric compartments}
        for i:= 6 to 15 do
        begin
             liquidnext[i,j] :=liquid[i,j]+up[i-1]*concentration_lumen[i-1,j]*dt;
             liquidnext[i,j] :=liquidnext[i,j]- up[i]*concentration_lumen[i,j]*dt;
        end;
        {compartment 16}
        liquidnext[16,j] := liquid[16,j] + up[15]*concentration_lumen[15,j]*dt;
  end;

  for compoundnumber:= 1 to number_of_compounds do
  begin
        {dishes}
        compoundnext[2,compoundnumber] := 0;
        for dishnumber := 1 to 5 do
        begin
           compounddishnext[dishnumber,compoundnumber] := compounddish[dishnumber,compoundnumber] - updish[dishnumber]*concentration_compound_dish[dishnumber,compoundnumber]*dt;
           compoundnext[2,compoundnumber] := compoundnext[2,compoundnumber] + updish[dishnumber]*concentration_compound_dish[dishnumber,compoundnumber]*dt;
        end;
        if solid[compoundnumber] > 10 then swallow_solid := 0 else swallow_solid := 1;           {here a maximum is set for the size of solid particles that are swallowed}
        compoundnext[2,compoundnumber] := compound[2,compoundnumber] + compoundnext[2,compoundnumber] - up[2]*swallow_solid*concentration_compound[2,compoundnumber]*dt;

        compoundnext[3,compoundnumber] := compound[3,compoundnumber] + up[2]*swallow_solid*concentration_compound[2,compoundnumber]*dt;
        compoundnext[4,compoundnumber] := compound[4,compoundnumber] ;
        compoundnext[5,compoundnumber] := compound[5,compoundnumber] -up[5]*concentration_compound[5,compoundnumber]*dt;

        {all other compartiments}
        for i:= 6 to 15 do
        begin
         compoundnext[i,compoundnumber] := compound[i,compoundnumber] + (up[i-1]*concentration_compound[i-1,compoundnumber]-up[i]*concentration_compound[i,compoundnumber])*dt;
        end;
        {compartment 16}
        compoundnext[16,compoundnumber] := compound[16,compoundnumber] + up[15]*concentration_compound[15,compoundnumber]*dt;
  end;


  {  if total_solids_corpus > 10 then transfer_solids_to_antrum := 10*dt else transfer_solids_to_antrum := total_solids_corpus*dt;
        if transfer_solids_to_antrum >0 then
          for j := 1 to comp do
          begin
             solids_corpus[j] := solids_corpus[j]-transfer_solids_to_antrum*solids_corpus[j]/total_solids_corpus;
             solids_antrum[j] := solids_antrum[j] + transfer_solids_to_antrum*solids_corpus[j]/total_solids_corpus;
           end;    }

 for dishnumber := 1 to 5 do
 begin
    for compoundnumber:= 1 to number_of_compounds do
       compounddish[dishnumber,compoundnumber] := compounddishnext[dishnumber,compoundnumber];
    for j:= 1 to comp do
       liquiddish[dishnumber,j] := liquiddishnext[dishnumber,j];
 end ;



 for i:= 2 to 16 do
 begin
    for compoundnumber:= 1 to number_of_compounds do
       compound[i,compoundnumber] := compoundnext[i,compoundnumber];
    for j:= 1 to comp do
       liquid[i,j] := liquidnext[i,j];
 end;
 {***************}
 {for i:= 2 to 16 do if (minutes >200) and (liquid[i,20] <> 0) then begin write(i, liquid[i,20]); readln; end; }



 {***************}



end;



begin
 clrscr;
 writeln('DIGESTION and ABSORPTION simulation program (DigSim)');
 writeln('version 77, March 17, 2023');
 writeln('George van Aken, www.insightfoodinside.com, info@insightfoodinside.com');
 writeln;
 writeln('DigSim provides a literature-based simulation of the alimentary system, including food intake, mouth and gastrointestinal physiology and processing, and the response by the body.');
 writeln('Copyright (C) 2023 George van Aken');
 writeln;
 writeln('This program is free software:');
 writeln('you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, version 3 or any later version.');
 writeln;
 writeln('This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  ');
 writeln('See theGNU Affero General Public License for more details.');
 writeln;
 writeln('You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.');

 YYE:= 2022; MME := 12; DDE := 31;   {einddatum}
 DecodeDate(Date,YY,MM,DD);

 Stop := 0;
 if (YY>YYE) then Stop := 1 Else
 if (YY=YYE) and (MM>MME) then Stop := 1 Else
 if (MM=MME) and (DD>DDE) then Stop := 1;

 Stop := 0; {abort end date}
 {If Stop = 1 then
    begin
        writeln('expired, renew subscription');
        readln;
    end
 else
     begin
       Write (Format ('Today is %d/%d/%d',[dd,mm,yy]));
       Writeln(Format (';  Subscription expires after %d/%d/%d',[DDE,MME,YYE]));
     end; }

 if Stop = 0 then
 begin
 writeln;writeln;
 assign(out, 'Results.txt');
 rewrite(out);
 comp := 110; {number of componenten in the calculation}
 parameter_invoer;
 read_diet;
 teller := 0 {substeps};
 initialiseer;

 newmealnumber := 0;
 write('Enter the time resolution for calculation in steps per minute (6 - 600)) = '); readln(substeps); writeln; writeln;     {substeps defines the number of calculational time steps per prefactor minutes; 60 corresponds with a time resolution of 1 second}
 dt := 1/substeps; {time step in minutes}
 prefactor:= 1; {prefactor defines the number of point per minute that is being stored} ; prefactor:= 1/prefactor;
 volumina;
 writeln('INITIALISING');

 {Reaching the fasted steady state}
 minutes := -500{start-time initialisation};
 sterretjes:= 60;
         dishnumber := 1;
         swallow_speed := 50;
           pH[dishnumber] := 7;
           liquidviscositydish[dishnumber] := 1; {for the initialization, the meal viscosities are mainly determined by the liquid phase}
           liquiddish[dishnumber,1] := 250; {water}
           liquiddish[dishnumber,2] := 1; {fat}
           liquiddish[dishnumber,3] := 50; {glucose}
           liquiddish[dishnumber,4] := 10; {protein}
           liquiddish[dishnumber,80]:= liquid[1,80]+ Ber_hydrogen(pH[1], liquid[1,61], liquid[1,66], liquid[1,65], liquid[1,67], liquid[1,4]);
           eating_speed_dish[1]:= 50;
         while minutes < 0 do
         begin
              minutes:= minutes+ dt;
              teller:= teller +1;
              make_step;
              if teller =sterretjes*substeps then
                 begin
                 write('*');
                 teller:= 0;
                 end;
         end;
         serum_glucose2_eq := serum_glucose2;
         blood_bile_salt_pool := 0.008;
         writeln('   DONE');


 {collected parameters set to 0}
 calories_colon := 0;
 calories_absorbed:= 0;
 total_gastric_fluid_secretion := 0; total_pancreatic_fluid_secretion:= 0; total_conjugated_bile_salt_secretion := 0; total_saliva_secretion := 0;
 {end collected parameters set to 0}

 write(out, ' minutes,');
 write(out, ' plate, mouth, corpus, antrum, duodenum 1,duodenum 2, jejunum 1, jejunum 2, jejunum 3, ileum 1, ileum 2, ileum 3, ileum 4, ileum 5, ileum 6, colon,');
 {write(out, 'FFA upper gastric, FFA antrum, FFA duodenum 1, FFA dusodenum 2, FFA jejunum 1, FFA jejunum 2, FFA jejunum 3, FFA ileum 1, FFA ileum 2, FFA ileum 3, FFA ileum 4, FFA ileum 5, FFA ileum 6, FFA colon,');}
 write(out, ' enterocytes/100,');
 write(out, ' total volume SI, CCK, PYY, gastric emptying, calories proximal, calories distal, calories colon, calories absorbed, gastric tone,');
 write(out, ' fullness, hunger,');
 write(out, ' peptides upper gastric, peptides antrum, peptides duodenum 1, peptides duodenum 2, peptides jejunum 1, peptides jejunum 2, peptides jejunum 3, peptides ileum 1, peptides ileum 2, peptides ileum 3, peptides ileum 4, peptides ileum 5, peptides ileum 6, peptides colon,');
 write(out, ' GIP, GLP1,');
 write(out, ' glucose upper gastric, glucose antrum, glucose duodenum 1, glucose duodenum 2, glucose jejunum 1, glucose jejunum 2, glucose jejunum 3, glucose ileum 1, glucose ileum 2, glucose ileum 3, glucose ileum 4, glucose ileum 5, glucose ileum 6, glucose colon,');
 write(out, ' fat delivered to colon, casein delivered to colon, glucose delivered to colon,');
 write(out, ' glucose absorbed, fatty acids absorbed, conjugated bile absorbed, fatty acid bile micelles absorbed, casein small peptides absorbed, small peptides absorbed,');
  write(out, ' pH plate, pH mouth, pH upper gastric, pH antrum, pH duodenum 1, pH duodenum 2, pH jejunum 1, pH jejunum 2, pH jejunum 3, pH ileum 1, pH ileum 2, pH ileum 3, pH ileum 4, pH ileum 5, pH ileum 6, pH colon,');
  write(out, ' Units pepsin in corpus, Units pepsin in antrum, solids_in_corpus, solids_in_fundus, solids_in_antrum, total solids in stomach,');
  write(out, ' serum glucose (mMol/l), serum insulin (pMol/l), pharmaceutical absorbed (g), bloodserumconc pharmaceutical (mg/l),');
  write(out, ' CCK emptying stress, Osmotic emptying stress, pH emptying stress,');
  write(out, ' pH fundus, fundus,');
  write(out, ' starch fundus, starch colon, sucrase_inhibitie[7], absorption_rate[7-3], absorption_rate[7-1], glucose in brush[7], arabinose in brush[7], sucrose in brush[7],');
  write(out, ' sucrase_inhibitie[15], absorption_rate[15-3], absorption_rate[15-1], osmotic_value_brush_proximal, osmotic_value_brush_distal, osmotic_value_lumen_duodenum, osmotic_value_lumen_distal,');
  write(out, ' serum_glucose_mg_dl, serum_insulin_mU_l, GLP1gastricemptying_stress, Viscousgastricemptying_stress,');
  write(out, ' 13C serum_glucose (mMol/l), 2H serum_glucose (mMol/l), glucose RaE (Mol/min), up[4], intercellulaire glucose muscle adipose (mMol/l), intercellulaire glucose skin (mMol/l), total_glucose_absorption_rate (Mol/min), Gastric Volume (mL), HbA1C (%), Leucine absorbed (g), total SI water absorption rate (g/min),');
  write(out, ' MMC status, Ghrelin, Motilin, bile_secretion_rate (g/min), Pancreatic fluid secretion (g/min), gastric fluid secretion (g/min), saliva secretion (g/ml), ');
  write(out, ' secreted gastric fluid (g), secreted conjugated bile salt (g), secreted panceratic fluid (g), secreted_saliva (g), gall bladder volume (g), conjugated bile salt concentration in gall bladder (g/g), duodenal conjugated bile salt concentration (g/g), gal bladder emptying rate (g/min), secteted bile (g), blood bile salt pool (g), unconjugated bile salt in colon (g), conjugated bile salt in secreted_bile (g/g), secretin (a.u.), glucagon (a.u.), glycogen (g),');
  write(out, ' RC13/C12, total lingual lipase lipolysis rate (%/min), total gastric lipase lipolysis rate (%/min), total pancreatic lipase lipolysis rate (%/min),');
  write(out, ' gastric emptying water, gastric emptying pharmaceutical, Units salivary amylkase in mouth, Units pancreatic amylase in duodenum');
  writeln(out);


 teller:= 0;
 starcounter := 0;


 while minutes < minutesmax +1 do
    begin
         if (minutes > mealtime[newmealnumber+1]) and (newmealnumber < number_of_meals) then
         begin
              writeln;
              newmealnumber:= newmealnumber+1;
              compoundreader;
         end;

         if teller = prefactor*substeps then
            begin
               starcounter:=starcounter +1;
                     if starcounter= sterretjes then
                       begin
                        write('*');
                        starcounter:=0;
                       end;
  {info on screen}
                     {writeln(Volume[2]); }
  {screen info}
                fat_in_colon:= liquid[16,2]+liquid[16,20]+liquid[16,54];
                casein_in_colon:= liquid[16,61]+liquid[16,62]+liquid[16,63]+liquid[16,64];
                glucose_in_colon:= liquid[16,3]+liquid[16,34]+liquid[16,35];
                write(out, minutes:7:5,', ');
                for i := 1 to 16 do write(out, volume[i]:12:8,', ');
                write(out,in_enterocytes:7:5,', ');
                write(out,volumetotal:7:5,', ',CCK:7:5,', ',PYY:7:5,', ',gastric_emptying:7:5,', ', calories_proximal:7:5,', ', calories_distal:7:5,', ', calories_colon:7:5,', ', calories_absorbed:7:5,', ', gastric_tone:7:5,', ');
                write(out, fullness:7:5,', ',hunger:7:5,', ');
                for i := 3 to 16 do write(out, liquid[i,65]:7:5,', ');
                write(out, RealGIP:7:5,', ', realGLP1:7:5,', ');
                for i := 3 to 16 do write(out, liquid[i,3]:7:5,', ');
                write(out, fat_in_colon:7:5,', ',casein_in_colon:7:5,', ',glucose_in_colon:7:5,', ');
                write(out, absorbed[3]:7:5,', ',absorbed[20]:7:5,', ',absorbed[51]:7:5,', ',absorbed[54]:7:5,', ',absorbed[64]:7:5,', ',absorbed[65]:7:5,', ');
                for i := 1 to 16 do write(out, pH[i]:12:8,', ');
                write(out, (pepsin_units(pH[3],liquid[3,40])/Volume[3]):7:5,', ',(pepsin_units(pH[4],liquid[4,40])/Volume[4]):7:5,', ',solids_in_corpus:7:5,', ',solids_in_fundus:7:5,', ', solids_in_antrum:7:5,', ', total_solids_in_stomach:7:5,', ');
                write(out, serum_glucose:7:5,', ',serum_insulin:7:5,', ');
                write(out, absorbed[5]:7:5,', ', bloodserumconc[5]:7:5,', ');
                write(out, CCKgastricEmptying_stress:7:5,', ', osmoticgastricemptying_stress:7:5,', ',pHgastricemptying_stress:7:5,', ');
                write(out, pHfundus:7:5,', ', Vfundus:7:5,', ');
                write(out, liquid[3,10]:7:5,',',liquid[16,10]:7:5,', ', sucrase_inhibitie[7]:7:5,', ', absorption_rate[7,3]:7:5,', ', absorption_rate[7,1]:7:5,', ', concentration_brush[7,3]:7:5,', ', concentration_brush[7,25]:7:5,', ', concentration_brush[7,21]:7:5,', ');
                write(out, sucrase_inhibitie[15]:7:5,', ', absorption_rate[15,3]:7:5,', ', absorption_rate[15,1]:7:5,', ',osmotic_value_brush_duodenum:7:5,', ',osmotic_value_brush_distal:7:5,', ',osmotic_value_lumen_duodenum:7:5,', ',osmotic_value_lumen_distal:7:5,', ');
                write(out, serum_glucose_mg_dl:7:5,', ',serum_insulin_mU_l:7:5,', ',GLP1gastricemptying_stress:7:5,', ',Viscousgastricemptying_stress:7:5,', ');
                write(out, serum_glucose1:7:5,', ',serum_glucose2:7:5,', ',(Gin/(1E-2+G1in))*G2in_molpermin*serum_glucose1/serum_glucose2:7:5,', ', up[4]:7:5,', ',intercellulair_glucose_muscle_adipose:7:5,', ', intercellulair_glucose_skin:7:5,', ', total_glucose_absorption_rate:7:5,', ', Gastric_Volume:7:5,', ', HbA1c:12:10,', ',absorbed[70]:7:5,', ',totalwaterabsorptionrate:7:5,', ');
                write(out, MMC_Status:5,', ',Ghrelin:7:5,', ',Motilin:7:5,', ', total_bile_secretion_rate:7:5,', ', pancreatic_fluid_secretion:7:5,', ', gastric_fluid_secretion:7:5,', ', saliva_addition:7:5, ', ');
                write(out, total_gastric_fluid_secretion:7:5,', ',total_conjugated_bile_salt_secretion:7:5,', ',total_pancreatic_fluid_secretion:7:5,', ', total_saliva_secretion:7:5, ', ', gall_bladder_volume:7:5, ', ',
                conjugated_bile_salt_concentration_in_gall_bladder:7:5,', ', duodenal_conjugated_bile_salt_concentration:7:5,', ', gall_bladder_emptying_rate:7:5,', ', total_bile_secreted:7:5,
                ', ', blood_bile_salt_pool:7:5,', ', liquid[16,59]:7:5,', ', conjugated_Bile_salt_in_secreted_bile:7:5,', ', secretin:7:5,', ', serum_glucagon:7:5,', ');
                write(out, glycogene_store:7:5,', ', R_C13_C12:7:5,', ', total_linguallipaselipolysis_rate:7:5,', ',total_gastriclipaselipolysis_rate:7:5,', ', total_pancreaticlipaselipolysis_rate:7:5,', ');
                write(out, gastric_emptying_water:7:5,', ', gastric_emptying_pharmaceutical:7:5,', ', liquid[2,57],', ', liquid[6,58],', ');
                writeln(out);
                teller:= 0;
            end;
         minutes:= minutes+ dt;
         teller:= teller +1;
         make_step;
        end;

        close(out);
 end;
 end.
