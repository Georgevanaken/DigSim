unit Model_Parameter_Invoer;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables;

Procedure parameter_invoer;

implementation
var i,j : integer;

Procedure parameter_invoer;
var DataFile : Text;
    TextFile : String;
begin

BasisSItransitRate := 0.02{0.02}{5};{0.05, Basis Small Intestinal Transit rate in min-1}
CCKEffectOnGastricEmptying := 2;	{2, CCK Effect On Gastric Emptying}
osmcont :=200;	{200, factor die het effect van de osmotische waarde op maaglediging bepaald}
MAX := 50; {maximaal volume van een SI segment}

for i:= 1 to 16 do hypertrophy[i] := 1;
  hypertrophy[5]:= 1; hypertrophy[6]:= 1;


 For j:= 1 to 80 do
  begin
  CCKresponsefactor[j]:= 0;
  PYYresponsefactor[j] := 0;
  GLP1responsefactor[j] := 0;
  GIPresponsefactor[j] := 0;
  Ghrelinresponsefactor[j] := 0;
  Motilinresponsefactor[j]:= 0;

  end;
 CCKresponsefactor[20] := {2000}20{0};	{20 fatty acids in min-1 CCKunit-1}
 CCKresponsefactor[54] := {2000}20{0};	{54 bile salt LCFA complex}
 CCKresponsefactor[64] := {200}20{20};	{64 caseinpeptides}
 CCKresponsefactor[65] := {100}20{20};	{65 small peptides}

 PYYresponsefactor[3]  := 100;	{3 = glucose}
  PYYresponsefactor[34]  := 100;	{34 = glucose1}
  PYYresponsefactor[35]  := 100;	{35 = glucose2}
 PYYresponsefactor[20] := 200{1}{200};	{20 fatty acids in min-1 CCKunit-1}
 PYYresponsefactor[54] := 200{1}{200};	{54 bile salt LCFA complex}
 PYYresponsefactor[64] := 100;	{64 caseinpeptides}
 PYYresponsefactor[65] := 100;	{65 small peptides}

 GLP1responsefactor[3] := 1{gestandaardiseerd op glucose response1};	{3 = glucose}         {Muller 2019 Glucagon-like peptide 1 (GLP-1)}
  GLP1responsefactor[34] := GLP1responsefactor[3];	{34 = glucose1}
  GLP1responsefactor[35] := GLP1responsefactor[3];	{35 = glucose2}
 GLP1responsefactor[20]:= 0{100};	{20 fatty acids in min-1 CCKunit-1}
 GLP1responsefactor[21]:= 0{1};	{21 absorbable sugar} GLP1responsefactor[36]:= GLP1responsefactor[21];
 GLP1responsefactor[22]:= 0{1}{35};	{22 fructose}
 GLP1responsefactor[54]:= 0{100};	{54 bile salt LCFA complex}
 GLP1responsefactor[64]:= 0.1{1}{0.1}{5};	{64 caseinpeptides}
 GLP1responsefactor[65]:= 0.1{1}{0.1}{5};	{65 small peptides}

 GIPresponsefactor[3]  := 1000;	{gestandaardiseeerd op glucose}{3 = glucose}
  GIPresponsefactor[34]  := GIPresponsefactor[3];	{34 = glucose1}
  GIPresponsefactor[35]  := GIPresponsefactor[3];	{35 = glucose2}
 GIPresponsefactor[20] := 10;	{20 fatty acids in min-1 CCKunit-1}
 GIPresponsefactor[21] := 0;	{21 absorbable sugar}
 GIPresponsefactor[22] := 5;	{22 fructose}
 GIPresponsefactor[25] := 5;	{25 arabinose}
 GIPresponsefactor[26] := 5;	{26 xylose}
 GIPresponsefactor[54] := 20;	{54 bile salt LCFA complex}

 Ghrelinresponsefactor[3]  := 0;	{gestandaardiseeerd op glucose}{3 = glucose}       {ghrelin wordt gestimuleerd door de hypothalamus, gerelateerd in trek}
  Ghrelinresponsefactor[34]  := 0;	{34 = glucose1}
  Ghrelinresponsefactor[35]  := 0;	{35 = glucose2}
 Ghrelinresponsefactor[20] := 0;	{20 fatty acids in min-1 CCKunit-1}
 Ghrelinresponsefactor[64] := 0;	{64 caseinpeptides}
 Ghrelinresponsefactor[65] := 0;	{65 small peptides}
 Ghrelinresponsefactor[25] := 0;	{25 arabinose}
 Ghrelinresponsefactor[26] := 0;	{26 xylose}
 Ghrelinresponsefactor[54] := 0;	{54 bile salt LCFA complex}

 Motilinresponsefactor[3]  := 0;	{gestandaardiseeerd op glucose}{3 = glucose}
  Motilinresponsefactor[34]  := 0;	{34 = glucose1}
  Motilinresponsefactor[35]  := 0;	{35 = glucose2}
 Motilinresponsefactor[20] := 0;	{20 fatty acids in min-1 CCKunit-1}
 Motilinresponsefactor[64] := 0;	{64 caseinpeptides}
 Motilinresponsefactor[65] := 0;	{65 small peptides}
 Motilinresponsefactor[25] := 0;	{25 arabinose}
 Motilinresponsefactor[26] := 0;	{26 xylose}
 Motilinresponsefactor[51] := 0{1};	{51 bile salt}
 Motilinresponsefactor[54] := 0{1};	{54 bile salt LCFA complex}

 Secretinresponsefactor[3]  := 0;	{gestandaardiseeerd op glucose}{3 = glucose}
  Secretinresponsefactor[34]  := 0;	{34 = glucose1}
  Secretinresponsefactor[35]  := 0;	{35 = glucose2}
 Secretinresponsefactor[20] := 0;	{20 fatty acids in min-1 CCKunit-1}
 Secretinresponsefactor[64] := 0;	{64 caseinpeptides}
 Secretinresponsefactor[65] := 0;	{65 small peptides}
 Secretinresponsefactor[25] := 0;	{25 arabinose}
 Secretinresponsefactor[26] := 0;	{26 xylose}
 Secretinresponsefactor[51] := 0;	{51 bile salt}
 Secretinresponsefactor[54] := 0;	{54 bile salt LCFA complex}



 for j := 1 to comp do proteolytic_susceptibility[j] := 0;
   proteolytic_susceptibility[4] :=1 {0.1};   {4 = protein}
   proteolytic_susceptibility[61] :=1 {0.1};   {61 = casein}
   proteolytic_susceptibility[62] :=1;   {62 = caseinpep}
   proteolytic_susceptibility[63] :=1;   {63 = caseinpan}
   proteolytic_susceptibility[64] :=0;   {64 = caseinpeptides}
   proteolytic_susceptibility[65] :=0;   {65 = small peptides}
   proteolytic_susceptibility[66] :=10 {0.1};   {66 = BSA}
   proteolytic_susceptibility[67] :=10 {0.1};   {67 = BLG}


     for i:= 1 to 15 do I_cell[i]:= 0;             {* hier wordt de verdeling van I-cellen over de dunne darm beschreven}
       I_cell[5] := 1    *hypertrophy[5];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       I_cell[6] := 1    *hypertrophy[6];            {*}
       I_cell[7] := 1    *hypertrophy[7];             {*}
       I_cell[8] := 0.2  *hypertrophy[8];             {*}
       I_cell[9] := 0.1  *hypertrophy[9];             {*}
     For i:= 1 to 15 do L_cell[i] := 0;            {* hier wordt de verdeling van L-cellen over de dunne darm beschreven}
       {for i:= 9 to 14 do L_cell[i] := 1;}{*}  {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       L_cell[5] := 0    *hypertrophy[5];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       L_cell[6] := 0.1  *hypertrophy[6];            {*}
       L_cell[7] := 0.2  *hypertrophy[7];             {*}
       L_cell[8] := 0.3  *hypertrophy[8];             {*}
       L_cell[9] := 0.4  *hypertrophy[9];
       L_cell[10] := 0.5 *hypertrophy[10];             {*}
       L_cell[11] := 0.6 *hypertrophy[11];             {*}
       L_cell[12] := 0.8 *hypertrophy[12];             {*}
       L_cell[13] := 1.0 *hypertrophy[13];             {*}
       L_cell[14] := 2   *hypertrophy[14];             {*}
       L_cell[15] := 5   *hypertrophy[15];             {*}
     For i:= 1 to 15 do K_cell[i] := 0;
       K_cell[5] := 2    *hypertrophy[5];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       K_cell[6] := 1.5  *hypertrophy[6];              {*}
       K_cell[7] := 1    *hypertrophy[7];               {*}
       K_cell[8] := 0.5  *hypertrophy[8];             {*}
       K_cell[9] := 0.2  *hypertrophy[9];             {*}
       K_cell[10] := 0.1 *hypertrophy[10];             {*}
     For i:= 1 to 15 do G_cell[i] := 0;
       G_cell[3] := 0    *hypertrophy[3];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       G_cell[4] := 0  *hypertrophy[4];              {*}
       G_cell[5] := 0    *hypertrophy[5];               {*}
       G_cell[6] := 0  *hypertrophy[6];             {*}
     For i:= 1 to 15 do Mo_cell[i] := 0;
       Mo_cell[5] := 2    *hypertrophy[5];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       Mo_cell[6] := 5  *hypertrophy[6];              {*}
       Mo_cell[7] := 1    *hypertrophy[7];               {*}
       Mo_cell[8] := 1  *hypertrophy[8];             {*}
       Mo_cell[9] := 1    *hypertrophy[9];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       Mo_cell[10] := 1  *hypertrophy[10];              {*}
       Mo_cell[11] := 1    *hypertrophy[11];               {*}
       Mo_cell[12] := 1  *hypertrophy[12];             {*}
       Mo_cell[13] := 1    *hypertrophy[13];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       Mo_cell[14] := 1  *hypertrophy[14];              {*}
       Mo_cell[15] := 1    *hypertrophy[15];               {*}
     For i:= 1 to 15 do S_cell[i] := 0;
       S_cell[5] := 1    *hypertrophy[5];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       S_cell[6] := 1  *hypertrophy[6];              {*}
       S_cell[7] := 0    *hypertrophy[7];               {*}
       S_cell[8] := 0  *hypertrophy[8];             {*}
       S_cell[9] := 0    *hypertrophy[9];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       S_cell[10] := 0  *hypertrophy[10];              {*}
       S_cell[11] := 0    *hypertrophy[11];               {*}
       S_cell[12] := 0  *hypertrophy[12];             {*}
       S_cell[13] := 0    *hypertrophy[13];             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       S_cell[14] := 0  *hypertrophy[14];              {*}
       S_cell[15] := 0    *hypertrophy[15];               {*}

 CCKdegradationconstant := 0.5 {0.5 min-1};
 PYYdegradationconstant := 0.1 {0.1 min-1};
 GIPdegradationconstant := 1/7.2 {website Drucker};
 GLP1degradationconstant := 1/2.3 {publicatie website Drucker};
 Ghrelindegradationconstant := 1/2 ;
 Motilindegradationconstant := 1/15 ; {geschat op grond van Saito 1980 Regulation of Motilin Secretion in the...}
 Secretindegradationconstant := 0.4  ; {2-3 min halfwaardetijd, Williams http://dx.doi.org/10.3998/panc.2013.5}


 CCKattenuationSI :=  0{.1}; {0.1 }      {* deze factor beschrijft de mate waarin SI transit door CCK geremd wordt}
 GIPattenuationSI := 0{1};
 PYYattenuationSI := 0.9; {0.02 }      {* deze factor beschrijft de mate waarin SI transit door PYY geremd wordt}
 GLP1attenuationSI := 0.01{2};
 GhrelinattenuationSI := 0; {0}
 MotilinattenuationSI := 2{2};
 SecretinattenuationSI := 0.5;        {gokwaarde, zou wel eens kunnen remmen}

 LinguallipolysisCST := 1500{ 100}*0.012 {0.012 min-1};  {* geschatte waarde: 0ngeveer 1.2 % per minuut (fed state Carriere)}
 GastriclipolysisCST := 0 { 100}*0.012 {0.012 min-1};  {* geschatte waarde: 0ngeveer 1.2 % per minuut (fed state Carriere)}
 PancreaticlipolysisCST := 5000{100}*0.012 {0.012 min-1};  {* geschatte waarde: 0ngeveer 1.2 % per minuut (fed state Carriere)}
 {dichtheden}
 for j := 1 to comp do density[j]:=1E20;   {enzymen nemen geen volume in dus dichthed is heel hoog}
  density[1] := 1;
  density[2]:= 0.9;
  density[3]:= 1.54;  density[34]:= 1.54; density[35]:= 1.54;
  density[4] := 1.4;
  density[5] := 2;
  density[6] := 2;
  density[10]:= 1.54;
  density[11]:= 1.54;
  density[12]:= 1.54;
  density[20]:= 0.9;
  density[21] := 1.59;density[36] := 1.59;
  density[22] := 1.69;
  density[23] := 1.69;
  density[25] := 1.69;
  density[26] := 1.69;
  density[27] := 1.69;
  density[54] := 1.69;  {?}
  density[61] := 1.2;
  density[62] := 1.2;
  density[63] := 1.2;
  density[64] := 1.2;
  density[65] := 1.2;
  density[66] := 1.2;
  density[67] := 1.2;
  density[80] := 1E6;
 {Calorische waarden}
  For j:= 1 to comp do caloric_value[j] := 0;   {calorische waarden in cal/g van de diverse componenten op 0 zetten}
  Caloric_value[2]:= 9.54;   {fat}            {calorische waarden volgens Weber en Ehrlein 1998}
  Caloric_value[3] := 3.86;  {glucose} Caloric_value[34] := 3.86;  {glucose1} Caloric_value[35] := 3.86;  {glucose2}
  caloric_value[4] := 4.76;  {protein}
  caloric_value[10] := 3.86;  {HDstarch}
  caloric_value[11] := 3.86;  {SDstarch}
  caloric_value[12] := 3.86;  {Rstarch}
  caloric_value[20] := 9.54;  {fatty acid}
  Caloric_value[21] := 3.86;  {sucrose}  Caloric_value[36] := 3.86;  {tracer sucrose}
  Caloric_value[22] := 3.86;  {fructose}
  Caloric_value[23] := 3.86;  {maltose}
  Caloric_value[25] := 0;  {arabinose}
  Caloric_value[26] := 0;  {xylose}
  Caloric_value[27] := 0;  {acarbose}
  caloric_value[54] := caloric_value[20]; {fatty acid bile complex}
  caloric_value[61] := caloric_value[4]; {casein}
  caloric_value[62] := caloric_value[4]; {caseinpep}
  caloric_value[63] := caloric_value[4]; {caseinpan}
  caloric_value[64] := caloric_value[4]; {caseinpeptides}
  caloric_value[65] := caloric_value[4]; {small peptides}
  caloric_value[66] := caloric_value[4]; {BSA}
  caloric_value[67] := caloric_value[4]; {BLG}
 {***********}

 {absorptionkinetic parameters according to Weber and Ehrlein 1998 for mini Troll pigs}
  for i:= 1 to 16 do
       for j:= 1 to comp do
        begin
         KmAbsorption[i,j]:= 0.0;  {Km each component set to 0}
         VmaxAbsorption[i,j] := 0; {maximum rate of absorption for each component set to 0}
        end;
 {***********}
  VmaxAbsorption[16,1]:= 10;        {water, in werkelijkheid een functie van de osmostische en hydrodynamische druk}
  KmAbsorption[16,1] := 1;          {nu nog alleen in dikke darm}
  wateruitwisselsnelheid:= 1;      {ter correctie calorische dichtheden}
 {***********}
  For i:= 5 to 15 do   {dunne darm}
  begin
  KmAbsorption[i,1] := 1 {1 kcal.min-1}; VmaxAbsorption[i,1]:= 100 {0 kcal.m-1.min-1} ;       {water}

  {************}
  KmAbsorption[i,20] := 0.1 {1 kcal.min-1}; VmaxAbsorption[i,20]:= 5{1} {0.5 kcal.m-1.min-1} ;      {fatty acids}

  {************}
  KmAbsorption[i,21] := 0 {1 kcal.min-1}; VmaxAbsorption[i,21]:= 0{1} {0.5 kcal.m-1.min-1} ;      {sucrose}

  {***********}

  KmAbsorption[i,54] := 0.1 {1 kcal.min-1}; VmaxAbsorption[i,54]:= 5{5}{0.5 kcal.m-1.min-1} ;       {bile salt LCFA complex}

 {***********}
  VmaxAbsorption[i,3]:= 3.48/0.2{0.2}{0.4}{0.7}{2}{0.15}{0.5}{0.7}{0.7 kcal.m-1.min-1} ;   {glucose}
  KmAbsorption[i,3]:= 1*0.022556 {10*0.022556} {0.20}{0.02} {0.07 kcal/ml};  {Olsen and Ingelfinger, jejunum; factoren 10 om te coorigeren voor metingen aan 10cm stukjes van ongeveer 1 m lange gesimuleerde segmenten}
  VmaxAbsorption[i,34]:=VmaxAbsorption[i,3];VmaxAbsorption[i,35]:=VmaxAbsorption[i,3];
  KmAbsorption[i,34]:= KmAbsorption[i,3];KmAbsorption[i,35]:= KmAbsorption[i,3];


 {***********}
  VmaxAbsorption[i,22]:= 0.5*3.48/1 {0.3}{0.5}{0.2} {kcal.m-1.min-1} ;     {fructose}
  KmAbsorption[i,22]:= 1*0.022556{0.20} {0.07 kcal.min-1};

 {***********}
  KmAbsorption[i,25]:= 1 {0.07 kcal.min-1}; VmaxAbsorption[i,25]:= 12{1}*KmAbsorption[i,25] {kcal.m-1.min-1} ;   {arabinose, diffusie}

 {***********}
  KmAbsorption[i,26]:= 1 {0.07 kcal.min-1};  VmaxAbsorption[i,26]:= 0.02*KmAbsorption[i,26]{kcal.m-1.min-1} ;     {xylose, diffusie}

 {***********}
  VmaxAbsorption[i,27]:= 0.0 {kcal.m-1.min-1} ;     {acarbose}
  KmAbsorption[i,27]:= 0.007 {0.07 kcal.min-1};



  {***********}
  VmaxAbsorption[i,4]:= 0 {0.3 kcal.m-1.min-1} ;   {protein}
  KmAbsorption[i,4]:= 0.02 {0.2 kcal.min-1};
  VmaxAbsorption[i,64]:= 0.3 {0.3 kcal.m-1.min-1} ;    {kleine caseine peptiden}
  VmaxAbsorption[i,65]:= 0.3 {0.3 kcal.m-1.min-1} ;    {peptiden}
  KmAbsorption[i,64]:= 0.02 {0.2 kcal.min-1};
  KmAbsorption[i,65]:= 0.02 {0.2 kcal.min-1};


  {***********}
   end;


 For i:= 5 to 15 do   {dunne darm, niet-macronutrienten die niet de physiologie beinvloeden  }
  begin
  VmaxAbsorption[i,5] := 1{0.01}; {pharmaceutical component in g.m-1.min-1}
  KmAbsorption[i,5] := 1; {in g.min-1}
  end;

 For i:= 10 to 15 do   {ileum}
  begin
  VmaxAbsorption[i,25] := 0*VmaxAbsorption[i,25]; {arabinose wordt niet in het ileum geabsorbeerd}
  VmaxAbsorption[i,26] := 0*VmaxAbsorption[i,26]; {xylose wordt niet in het ileum geabsorbeerd}
  KmAbsorption[i,51] := 0.1 {1 kcal.min-1}; VmaxAbsorption[i,51]:= 0.02{5}{0.5 kcal.m-1.min-1} ;       {conjugated bile salt}
  end;

 for j:= 1 to comp do
     begin
                                                            {snellere absorptie in duodenum en jejunum}
        VmaxAbsorption[5,j]:=  1{0.2}{2}*VmaxAbsorption[5,j];
        VmaxAbsorption[6,j]:= 1{0.5}{2}*VmaxAbsorption[6,j];
        VmaxAbsorption[7,j]:= 1{3}*VmaxAbsorption[7,j];
        VmaxAbsorption[8,j]:= 1{3}*VmaxAbsorption[8,j];
        VmaxAbsorption[9,j]:= 1{2}*VmaxAbsorption[9,j];
        VmaxAbsorption[10,j]:= 0.5{1}*VmaxAbsorption[10,j];
        VmaxAbsorption[11,j]:= 0.5{0.9}*VmaxAbsorption[11,j];
        VmaxAbsorption[12,j]:= 0.5{0.8}*VmaxAbsorption[12,j];
        VmaxAbsorption[13,j]:= 0.5{0.7}*VmaxAbsorption[13,j];
        VmaxAbsorption[14,j]:= 0.5{0.6}*VmaxAbsorption[14,j];
        VmaxAbsorption[15,j]:= 0.5{0.5}*VmaxAbsorption[15,j];
     end;
  for i := 5 to 15 do for j:= 1 to comp do  VmaxAbsorption[i,j]:=  VmaxAbsorption[i,j]*hypertrophy[i];

  {voorfactor in de diffusie absorptie berekening die brush border concentratie koppeld aan afvoer door absorptie en aanvoer door diffusie uit lumen}


 for i:= 5 to 16 do delta[i]:= 1*5E-2 ; {in cm; 500 micron dik}  { I T JOHNSON AND JENNIFER M GEE, Effect of gel-forming gums on the intestinal unstirred layer and sugar transport in vitro, Gut, 1981, 22, 398-403}
  for j:= 1 to comp do D[j] := 1E-20; {geen diffusie voor alle componenten door de muceuse laag}{diffucie coefficienten by 37 graden C in water}
      D[1] := 2.27E-5*60{0.10};  {water gaat de mucus laag in, in cm2min-1}
      D[3] := 9.4E-6*60;  {glucose gaat de mucus laag in, in cm2min-1}
        D[34] := D[3];D[35]:=D[3];
      D[20] := 6E-6*60{0.05}; {vetzuren gaan de mucuslaag in}
      D[21] := 7.056E-6*60{0.10};D[36] := D[21];   {sucrose gaat de mucus laag in, in cm2m-1}
      D[22] := 10E-6*60{0.10};  {fructose gaat de mucus laag in}
      D[23] := 4E-6*60{0.10};  {maltose gaat de mucus laag in}
            D[24] := D[23];  {maltose gaat de mucus laag in}
      D[25] := 10E-6*60{0.10};  {arabinose gaat de mucus laag in}
      D[26] := 10E-6*60{0.10};  {xylose gaat de mucus laag in}
      D[27] := 10E-6*60{0.10};  {acarbose gaat de mucus laag in}
      D[5] := 6E-6*60; {pharmaceutical}
      {D[50] := 0.05E-5; } {lipase}
      D[51] := 6E-6*60{0.05};  {bile salts}
      D[52] := 4E-6*60{0.05};  {bile phospholipids}
      D[53] := 4E-6*60{00.05};  {bile lysophospholipids}
      D[54] := 4E-6*60{0.05};  {bile salt LCFA complex}
      D[64] := 6.7E-6*60{0.05};   {casein kleine peptiden kunnen de mucus laag in}
      D[65] := 6.7E-6*60{0.05};   {kleine peptiden kunnen de mucus laag in}

      D[80] := 20E-5*60; {diffusie waterstof ionen gekoppeld aan een tegenion zoals chloride}

      {for j := 1 to comp do D[j]:= 0.1*D[j]; }{optioneel effect van guar gom die de permeabiliteit van de muceuze laag verminderd}



  for i:= 5 to 15 do for j:= 1 to comp do
    begin
    surface_factor[i] := 1/(3*100*100*hypertrophy[i]){2E-3};{de reden hiervoor is dat het oppervlak van de brush border door villi en microvilli sterk vergroot is en er roering optreedt door peristaltiek}
    {surface_factor[i] :=1/9000; {nieuwe schatting}  }
    T[i,j]:= delta[i]/D[j];
    end;
  For i:= 1 to 4 do Maxtotalabsorption_rate[i]:= 0;
  For i:= 5 to 9 do Maxtotalabsorption_rate[i]:= 1{*hypertrophy[i]};        {1; typische waarden in mini Troll pigs volgens Weber en Ehrlein 1998}
  For i:= 10 to 15 do Maxtotalabsorption_rate[i]:= 0.5{*hypertrophy[i]};       {1; typische waarden in mini Troll pigs volgens Weber en Ehrlein 1998}

  maxcaloricabsorptioncapacity := 0;
  for i:= 1 to 15 do maxcaloricabsorptioncapacity := maxcaloricabsorptioncapacity + Maxtotalabsorption_rate[i];
 {************}
end;

end.
