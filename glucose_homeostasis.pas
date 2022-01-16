unit Glucose_homeostasis;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables, subcellularinsulin_Pederson;
var
  serum_glucose,serum_glucose_mg_dl, serum_insulin, serum_insulin_mU_l, serum_glucagon,serum_ketons, serum_FA, serum_peptides :real;
  glycogene_store, x1, x2, x3, tp, tpf:real;

  Gin: real; {incoming glucose in mg/min through mesentric/portal vein}
  Ip: real; {amount of insulin in the plasma, in mU; 1 mU insulin = 6.67pmol =  6.67*5807.57E-12 g = 38.7 ng}{1 microgram insulin = 1/0.0387 mU insulin}
  Ii: real; {amount of insulin in the intercellular space, in mU}

const beta_cells_percentage = 1;
      alfa_cells_percentage = 1;
      Vp = 3; {3 l, distribution volume for insulin in plasma}
      Vi = 11; {11 l, effective volume of intercellular space}
      Vg = 10; {10 l, effective volume in which glucose is dissolved}
      E = 0.2; {0.2 l/min, transport rate constant of insulin between plasma and intercellular space by passive diffusion}

      ti = 100; {100 min, insulin degradation constant in intecellular space}
      td1 = 1{10}; {36 min, delay associated with the time lag between the appearance of Insulin in plasma and the inhibition of hepatic glucose production}
      td2 = 1{/10}; {36 min, delay associated with the time lag between the appearance of Insulin in plasma and the inhibition of hepatic glucose production}
      td3 = 1{/10}; {36 min, delay associated with the time lag between the appearance of Insulin in plasma and the inhibition of hepatic glucose production}
      Rm = 210{40}{50}{210}; {210 mU/min, rate constant of insulin production controlled by glucose, as used in the f1 by Tolic}
      a1 = 300; {300 mg/l, fitting parameter in f1 Tolic}
      C1 = 2000; {2000 mg/l, fitting parameter in f1 Tolic}
      Ub = 72; {72 mg/min, fitting parameter in f2, related to insulin-independent glucose utilization by brain and nerve cells}
      C2 = 144; {144 mg/l, fitting parameter in f2, related to insulin-independent glucose utilization by brain and nerve cells}
      C3 = 1000; {1000 mg/l, fitting parameter in f3, related to glucose-dependent glucose utilization by muscle and fat cells}
      U0 = 40; {40 mg/min, fitting parameter in f4, related to insulin-independent glucose utilization by muscle and fat cells}
      Um = 940; {940 mg/min, fitting parameter in f4, related to insulin-independent glucose utilization by muscle and fat cells }
      beta = 1.77; {1.77 fitting parameter in f4, related to insulin-independent glucose utilization by muscle and fat cells}
      C4 = 80; {80 mU/l, fitting parameter in f4, related to insulin-independent glucose utilization by muscle and fat cells}
      Rg = 180{160}; {180 mg/min, fitting parameter in f5, related to insulin-dependent glucose hepatic production}
      alpha = 0.29{0.4}; {0.29 l/mU, fitting parameter in f5, related to insulin-dependent glucose hepatic production}
      C5 = 26{29}; {26 mU/l, fitting parameter in f5, related to insulin-dependent glucose hepatic production}
      gamma = 5.0; {5.0 fitting constant, related to the inhibition of hepatic glucose production by hyperglycemia}
      C6 = 2.0;  {2.0 fitting constant, related to the inhibition of hepatic glucose production by hyperglycemia}
      Sb = 20; { 20 mg/min, fitting constant, related to enhanced splanchic (liver and gut) glucose uptake during hyperglycemia}
      Sm = 140; {140 mg/min, fitting constant, related to enhanced splanchic (liver and gut) glucose uptake during hyperglycemia}
      ddelta = -2.4; {-2.4 fitting constant, related to enhanced splanchic (liver and gut) glucose uptake during hyperglycemia}
      C7 = 2.0; {2.0 fitting constant, related to enhanced splanchic (liver and gut) glucose uptake during hyperglycemia}

      insulin_sensitivity = 1.0; {waarde tussen 0 en 1}
      glucose_sensitivity = 1.3{1.3}; {waarde tussen 0 en 1}

      ratio =0; {waarde tussen 0 en 1}



function f2(Gluc:real):real;
function f3(Gluc:real):real;
function f4(Ii:real):real;
function f5(x3:real):real;
function f6(Gluc, Gin:real):real;
function f7(Gluc, Gin:real):real;
procedure insulin_dynamics;
procedure glucagon_dynamics;
Procedure Initialise_glucose_dynamics;
Procedure serum_glucose_dynamics;

{procedure serum_glucose_dynamics;}

implementation


              {Gluc = mg glucose}

  function f2(Gluc:real):real;                       {glucose utilisatie door hersenen en zenuwcellen in mg/min}
    begin
      f2 := Ub*(1-exp(-glucose_sensitivity*Gluc/(C2*Vg)));               {insulin independent glucose utilization by brain and nerve cells, op basis van een fit op experimentele data van Verdonk 1981}
    end;
  function f3(Gluc:real):real;                       {glucose utilisatie door spier en vetcellen, glucose afhankelijke term, in mg/min}
    begin
      f3 := glucose_sensitivity*Gluc/(C3*Vg) ;                            {op basis van een fit van experimentele data van Rizza 1981, Verdonk 1981}
    end;
  function f4(Ii:real):real;                      {glucose utilisatie door spier en vetcellen, Insuline afhankelijke term}
    begin
      f4:=U0+ (Um-U0)/(1+exp(-beta*ln((insulin_sensitivity*Ii/C4)*(1/Vi+1/(E*ti)))));   {op basis van het werk van Rizza 1981 en een scheiding in een insuline en glucose afhankleijke term}
    end;
  function f5(x3:real):real;
    begin
      f5:=Rg/(1+exp(alpha*(x3*insulin_sensitivity/Vp-C5)));          {effect van insuline op de hepatische glucose productie op basis van werk Rizza 1981}
    end;
  function f6(Gluc,Gin:real):real;
    begin
      f6:= 1/(1+exp(gamma*(glucose_sensitivity*(Gluc/Vg+ratio*(Gin/0.57-Gluc/Vg))/C3-C6)));    {suppression of hepatic glucose production by glucose}
    end;
  function f7(Gluc, Gin:real):real;
    begin
      f7:=Sb+(Sm-Sb)*(1/(1+exp(ddelta*(glucose_sensitivity*(Gluc/Vg+ratio*(Gin/0.57-Gluc/Vg))/C3-C7)))-1/(1+exp(-ddelta*C7))); {enhanced splanchnic (liver+gut) glucose uptake at high glucose levels}
    end;

  Procedure Initialise_glucose_dynamics;
  begin
    betacelverdeling;
    serum_glucagon:= 10;{mg/dl}
    glycogene_store0 := 50000; {mg}
        glycogene_store1 := 0; {mg}
        glycogene_store2 := 0; {mg}
        glycogene_store := glycogene_store0 +glycogene_store1+glycogene_store2;
    x1 := 10;
    x2 := 10;
    x3 := 10;
    {Gluc := 120*20*Vg; }   Gluc0 := 94*10*Vg; Gluc1 := 0 ; Gluc2:= 0; Gluc := Gluc0+Gluc1+Gluc2;
    {Ip := 20*Vp;} Ip := 10*Vp;
    {Ii:= 20*Vp; }Ii:= 10*Vi;
  end;

  Procedure insulin_dynamics; {Tolic et al., J. Theor. Biol. (2000)207, 361-375}
  {input Gin: glucose input in mg/min; gluc: glucose serum kwantiteit in mg; output f1 insuline serum concentratie in mU/min}
  var Ip_increment,Ii_increment,G_increment,x1_increment,x2_increment,x3_increment, f1_Tolic, incretin_factor: real;
    begin
      tpf := 4{4}; { 4 in adjusted model, 6 min, insulin degradation constant in plasma in fasted state}
      incretin:= 1*(0.1*GIP/10+0.1*GLP1/10) {1+1.5*(1*GIP/400+1*0.5*GLP1/200)};
      incretin_factor := 1 + incretin;
      tp := tpf*(1+0.5*incretin/(1+incretin));

      Gin := G0in+G1in+G2in;
      Gluc := Gluc0 + Gluc1 + Gluc2;

      {f1 volgens Pederson}
      G:= 0.1*glucose_sensitivity*(Gluc/Vg{+Gin/1}){(Vp+Vi)} ;  {omrekening!!!!???? Gluc is in mg, G is hier in mg/dl (namelijk 0.1 * concentratie mg/l}
      {calc_f1_Pederson;} {f1 berekend volgens Pedersen et al 2008 voor beta cellen; een subcellulair model voor insuline productie uitgevoerd, levert f1}
      {f1:=incretin_factor*f1; }
      {einde f1 volgens Pederson}
      {oorspronkelijk model voor insuline productie van Tolic 2000}
      f1_Tolic := incretin_factor*Rm/(1+exp((C1-glucose_sensitivity*(Gluc/Vg+ratio*(Gin/0.57{-Gluc/Vg})))/a1));   {for the fed state: Gin is in mg/min and the glucose is desolved in volume/time  = about 9.5 ml/s = 570 ml/min : At rest, blood flow through the mesenteric artery was 6.3 ± 2.6 ml/sec and 9.5 ± 2.1 ml/sec in the carotid artery. After the test meal, mesenteric artery blood flow increased significantly (p < 0.0001) and reached maximal hyperemia (20.3 ± 7.4 ml/sec) after 45 minutes. Kurt Jäger, M.D. et al, Measurement of mesenteric blood flow by duplex scanning Journal of Vascular Surgery, 1986, Pages 462-469}
      {einde f1 volgens Tolic}

      f1:= {f1_Pederson}f1_Tolic;

      Ip_increment := f1 {in milliunits/min}-E*(Ip/Vp-Ii/Vi)-Ip/tp; {increment in plasma insuline, met transfer (E)van plasma naar interstitionele volumina en afbraak (tijdschaal tp)}
      Ii_increment := E*(Ip/Vp-Ii/Vi)-Ii/ti;           {increment in interstitiele insuline}
      if Gin>0 then
         begin
           G0_increment := G0in-(Gluc0/Gluc)*(f2(Gluc)+f3(Gluc)*F4(Ii)+f7(Gluc, Gin))+f5(x3)*f6(Gluc, Gin);   {Gin is mg glucose/min}{F2,F3,F4 en F7 hebben betrekking op glucose verwijdering uit het bloed, waarin ook tracerglucose zit, F5 en F6 op de endogeen geprodiuceerde glucose, waarin nog weinig tracer materiaal zit}
           G1_increment := G1in-(Gluc1/Gluc)*(f2(Gluc)+f3(Gluc)*F4(Ii)+f7(Gluc, Gin));   {Gin is mg glucose/min}
           G2_increment := G2in-(Gluc2/Gluc)*(f2(Gluc)+f3(Gluc)*F4(Ii)+f7(Gluc, Gin));   {Gin is mg glucose/min}
           G_increment:= G0_increment+G1_increment+G2_increment;
         end;
      x1_increment:= 3*(Ip-x1)/td1;               {x, x2, x3 beschrijven samen een derde orde afhanelijkheid in de tijd met tijdschaal td, op grond van Prager 1986}
      x2_increment:= 3*(x1-x2)/td2;
      x3_increment:= 3*(x2-x3)/td3;

      Ip:= Ip +Ip_increment*dt;  {if Ip<0.0001 then Ip := 0.0001;}
      Ii:= Ii +Ii_increment*dt;  {if Ii<0.0001 then Ii := 0.0001;}

      Gluc0:= Gluc0+ G0_increment*dt;
      Gluc1:= Gluc1+ G1_increment*dt;
      Gluc2:= Gluc2+ G2_increment*dt;
      Gluc := Gluc0+Gluc1+Gluc2;

      x1:= x1+ x1_increment*dt;
      x2:= x2+ x2_increment*dt;
      x3:= x3+ x3_increment*dt;

         serum_glucose0_mg_dl := 0.1*Gluc0/Vg;  {in mg/dl}
         serum_glucose1_mg_dl := 0.1*Gluc1/Vg;  {in mg/dl}
         serum_glucose2_mg_dl := 0.1*Gluc2/Vg;  {in mg/dl}
         serum_glucose_mg_dl :=serum_glucose0_mg_dl+serum_glucose1_mg_dl+serum_glucose2_mg_dl ;  {in mg/dl}

         serum_glucose0 := serum_glucose0_mg_dl*10/180.15588; {in mmol/l}
         serum_glucose1 := serum_glucose1_mg_dl*10/180.15588; {in mmol/l}
         serum_glucose2 := serum_glucose2_mg_dl*10/180.15588; {in mmol/l}
         serum_glucose := serum_glucose0+serum_glucose1+serum_glucose2; {in mmol/l}
      serum_insulin_mU_l := Ip/Vp;        {in milliUnits/l}
      serum_insulin := serum_insulin_mU_l*6.67; {in pmol/l}
      glucagon_dynamics;
      serum_glucose_dynamics;
    end;

  procedure glucagon_dynamics;
  var glucagon_production, glucagon_decay:real;
  const gluc_criterion = 4;
  begin
    if (gluc < gluc_criterion-0.01) then glucagon_production:= 0.01*exp(4*ln((gluc_criterion-gluc)/(2*gluc_criterion-gluc)))*alfa_cells_percentage else glucagon_production := 0;
    glucagon_decay:= serum_glucagon/5; {halflife 3-6 min}
    serum_glucagon := serum_glucagon + (glucagon_production-glucagon_decay)*dt;
  end;

  procedure serum_glucose_dynamics;
 begin
    glycogene_store:= glycogene_store0+ glycogene_store1+glycogene_store2;
    glycogene_store0 := glycogene_store0 +(0.01*gluc0*serum_insulin_mU_l/50  - 2000*(glycogene_store0/glycogene_store)*gluc/(gluc+1000)*serum_glucagon)*dt;
        glycogene_store1 := glycogene_store1 +(0.01*gluc1*serum_insulin_mU_l/50  - 2000*(glycogene_store1/glycogene_store)*gluc/(gluc+1000)*serum_glucagon)*dt;
        glycogene_store2 := glycogene_store2 +(0.01*gluc2*serum_insulin_mU_l/50  - 2000*(glycogene_store2/glycogene_store)*gluc/(gluc+1000)*serum_glucagon)*dt;
    glycogene_store:= glycogene_store0+ glycogene_store1+glycogene_store2;

     Gluc0 := Gluc0 + (2000*(glycogene_store0/glycogene_store)*gluc/(gluc+1000)*serum_glucagon-0.01*Gluc0*serum_insulin_mU_l/50) *dt;
        Gluc1 := Gluc1 + (2000*(glycogene_store1/glycogene_store)*gluc/(gluc+1000)*serum_glucagon-0.01*Gluc1*serum_insulin_mU_l/50) *dt;
        Gluc2 := Gluc2 + (2000*(glycogene_store2/glycogene_store)*gluc/(gluc+1000)*serum_glucagon-0.01*Gluc2*serum_insulin_mU_l/50) *dt;
     Gluc:= Gluc0+Gluc1+Gluc2;
  end;
end.

