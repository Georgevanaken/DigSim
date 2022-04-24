unit Glucose_homeostasis2;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils,globalvariables, subcellularinsulin_Pederson;
var
  serum_glucose,serum_glucose_mg_dl, serum_insulin, serum_insulin_mU_l, serum_glucagon,serum_ketons, serum_FA, serum_peptides, glucutili, endoglucprod, gluci, gluciskin  :real;
  glycogene_store, x1, x2, x3, tp, tpf:real;

  Gin: real; {incoming glucose in mg/min through mesentric/portal vein}
  Ip: real; {amount of insulin in the plasma, in mU; 1 mU insulin = 6.67pmol =  6.67*5807.57E-12 g = 38.7 ng}{1 microgram insulin = 1/0.0387 mU insulin}
  Ii: real; {amount of insulin in the intercellular space, in mU}

const beta_cells_percentage = 1;
      alfa_cells_percentage = 1;
      Vp = 3; {3 l, distribution volume for insulin in plasma}
      Vi = 11; {11 l, effective volume of intercellular space}
      Vg = 10; {10 l, effective volume in which glucose is dissolved}
      Vs = 1; {1 l, effectieve volume of intercellular space in skin}
      E = 0.2; {0.2 l/min, transport rate constant of insulin between plasma and intercellular space by passive diffusion}

      ti = 100; {100 min, insulin degradation constant in intecellular space}
      td1 = 3{10}{25}{10}{3}; {36 min, delay associated with the time lag between the appearance of Insulin in plasma and the inhibition of hepatic glucose production}
      td2 = 3{10}{25}{12}{3}; {36 min, delay associated with the time lag between the appearance of Insulin in plasma and the inhibition of hepatic glucose production}
      td3 = 4{10}{25}{12}{4}; {36 min, delay associated with the time lag between the appearance of Insulin in plasma and the inhibition of hepatic glucose production}
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

      insulin_sensitivity = 1.40{1.40; 0.30 voor prediabetic}{1.5}; {sensitivity of the insulin receptorsn, 1 is average}
      glucose_sensitivity = 1.20{1.20; 0.60 for pre diabetic }{1.10}{1.3}; {gevoeligheid of the glucose receptors, 1 is average}

      ratio =0.1; {value between 0 en 1}{ratio of the detection of glucose in portal blood compared to in serum blood}
      k1 = 0.1; {0.1/min; equilibrationconstant between intercellular glucose and bloodglucose in muscle and adipose tissue}
      k2 = 0.1; {0.1/min; equilibratieconstant between intercellulair glucose and bloedglucose in skin}
      k3 = 0; {0.1/min; glucose ustilisation constant in skin}

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

  function f2(Gluc:real):real;                       {glucose utilisatie by brain and nerve cells in mg/min}
    begin
      f2 := Ub*(1-exp(-glucose_sensitivity*Gluc/(C2*Vg)));               {insulin independent glucose utilization by brain and nerve cells, based on a fit of experimentele data by Verdonk 1981}
    end;
  function f3(Gluc:real):real;                       {glucose utilisatie by muscle and fat cells, glucose-dependent term, in mg/min}
    begin
      f3 := glucose_sensitivity*Gluc/(C3*Vg) ;                            {basid on a fit of experimentele data by Rizza 1981, Verdonk 1981}
    end;
  function f4(Ii:real):real;                      {glucose utilisation by muscle and fat cells, Insulin dependent term}
    begin
      f4:=U0+ (Um-U0)/(1+exp(-beta*ln((insulin_sensitivity*Ii/C4)*(1/Vi+1/(E*ti)))));   {based on Rizza 1981 and separating insulin and glucose dependent terms}
    end;
  function f5(x3:real):real;
    begin
      f5:=Rg/(1+exp(alpha*(x3*insulin_sensitivity/Vp-C5)));          {effect of insulin on the hepatic glucose production based on Rizza 1981}
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

    x1 := 10;
    x2 := 10;
    x3 := 10;
    {Gluc := 120*20*Vg; }   Gluc0 := 94*10*Vg; Gluc1 := 0 ; Gluc2:= 0; Gluc := Gluc0+Gluc1+Gluc2;
    {Ip := 20*Vp;} Ip := 10*Vp;
    {Ii:= 20*Vp; }Ii:= 10*Vi;
    gluci := 10;
    gluciskin:= 1;
  end;

  Procedure insulin_dynamics; {Tolic et al., J. Theor. Biol. (2000)207, 361-375}
  {input Gin: glucose input in mg/min; gluc: glucose serum quantity in mg; output f1 insuline serum concentratie in mU/min}
  var Ip_increment,Ii_increment,G_increment,x1_increment,x2_increment,x3_increment, f1_Tolic, incretin_factor, gluci_increment, Gluciskin_increment, Gi: real;
    begin
      tpf := 4{4}; { 4 in adjusted model, 6 min, insulin degradation constant in plasma in fasted state}
      incretin:= 5*(0.02*realGIP+0.05*realGLP1) {1+1.5*(1*GIP/400+1*0.5*GLP1/200)};   {incretine effect on insulin degradation}
      incretin_factor := 1 + incretin;
      tp := tpf*(1+0.5*incretin/(1+incretin));    {change in half lifetime of insulin degradation by incretins}

      Gin := G0in+G1in+G2in;
      Gluc := Gluc0 + Gluc1 + Gluc2;

      {f1 according to Pederson}
      {G:= 0.1*glucose_sensitivity*(Gluc/Vg{+Gin/1}){(Vp+Vi)}} ;  {conversion!!!!???? Gluc is in mg, G is here in mg/dl (to be precise: 0.1 * concentration mg/l}  {only requiredfor Pederson calculation}
      {calc_f1_Pederson;} {f1 calculated according to Pedersen et al 2008 for beta cells; a subcellular model voor insuline production delivers f1}
      {f1:=incretin_factor*f1; }
      {end f1 accroding to Pederson}
      {original model for insuline production according to Tolic 2000}
      f1_Tolic := incretin_factor*Rm/(1+exp((C1-glucose_sensitivity*(Gluc/Vg+ratio*(Gin/0.57{-Gluc/Vg})))/a1));   {for the fed state: Gin is in mg/min and the glucose is desolved in volume/time  = about 9.5 ml/s = 570 ml/min : At rest, blood flow through the mesenteric artery was 6.3 ± 2.6 ml/sec and 9.5 ± 2.1 ml/sec in the carotid artery. After the test meal, mesenteric artery blood flow increased significantly (p < 0.0001) and reached maximal hyperemia (20.3 ± 7.4 ml/sec) after 45 minutes. Kurt Jäger, M.D. et al, Measurement of mesenteric blood flow by duplex scanning Journal of Vascular Surgery, 1986, Pages 462-469}
      {end f1 according to Tolic}

      f1:= {f1_Pederson}f1_Tolic;

      Ip_increment := f1 {in milliunits/min}-E*(Ip/Vp-Ii/Vi)-Ip/tp; {increment in plasma insulin, with transfer (E) of plasma to interstitional volumes en breakdown (timescale tp)}
      Ii_increment := E*(Ip/Vp-Ii/Vi)-Ii/ti;           {increment in interstitiele insuline}
      {if Gin>0 then} {doet niets}
         begin
           gluci_increment := k1*Vi*(Gluc/Vg-gluci/Vi)-glucutili;  {dynamic of  intercellulaire glucose in muscle and adipose tissue}
           gluciskin_increment := k2*Vs*(Gluc/Vg-gluciskin/Vs)- k3*gluciskin; {dynamics of intercellular glucose in skin}
           {Gin in mg glucose/min}{dynamic of intercellular glucose}
           {F2,F3,F4 en F7 relate to glucose removal from blood, also including tracerglucose, F5 and F6 on endogeneously produced glucose with assuned very low content of tracer glucose}
           Gi := Gluc ;

           endoglucprod := f5(x3)*f6(Gi, Gin);{Gin is mg glucose/min}
           glucutili := f3(Gi)*F4(Ii)+f2(Gi)+f7(Gi, Gin);

           G0_increment := G0in-(Gluc0/Gluc)*glucutili+endoglucprod;
           G1_increment := G1in-(Gluc1/Gluc)*glucutili;
           G2_increment := G2in-(Gluc2/Gluc)*glucutili;
           G_increment:= G0_increment+G1_increment+G2_increment;
         end;
      x1_increment:= 3*(Ip-x1)/td1;               {x, x2, x3 together describe a third order dependency of time with tijdschale td, based on Prager 1986}
      x2_increment:= 3*(x1-x2)/td2;
      x3_increment:= 3*(x2-x3)/td3;

      Ip:= Ip +Ip_increment*dt;
      Ii:= Ii +Ii_increment*dt;

      Gluc0:= Gluc0+ G0_increment*dt;        {units Gluc0 is mg glucose}
      Gluc1:= Gluc1+ G1_increment*dt;
      Gluc2:= Gluc2+ G2_increment*dt;
      Gluc := Gluc0+Gluc1+Gluc2;

      Gluci := Gluci+ Gluci_increment*dt;
      Gluciskin := Gluciskin+ Gluciskin_increment*dt;

      x1:= x1+ x1_increment*dt;
      x2:= x2+ x2_increment*dt;
      x3:= x3+ x3_increment*dt;

      serum_insulin_mU_l := Ip/Vp;        {in milliUnits/l}
      serum_insulin := serum_insulin_mU_l*6.67; {in pmol/l}
      glucagon_dynamics;
      serum_glucose_dynamics;

      serum_glucose0_mg_dl := 0.1*Gluc0/Vg;  {in mg/dl}
         serum_glucose1_mg_dl := 0.1*Gluc1/Vg;  {in mg/dl}
         serum_glucose2_mg_dl := 0.1*Gluc2/Vg;  {in mg/dl}
         serum_glucose_mg_dl :=serum_glucose0_mg_dl+serum_glucose1_mg_dl+serum_glucose2_mg_dl ;  {in mg/dl}
         intercellulair_glucose_muscle_adipose_mg_dl := 0.1*Gluci/Vi; {in mg/dl}
         intercellulair_glucose_skin_mg_dl := 0.1*Gluciskin/Vs; {in mg/dl}

         serum_glucose0 := serum_glucose0_mg_dl*10/180.15588; {in mmol/l}
         serum_glucose1 := serum_glucose1_mg_dl*10/180.15588; {in mmol/l}
         serum_glucose2 := serum_glucose2_mg_dl*10/180.15588; {in mmol/l}
         serum_glucose := serum_glucose0+serum_glucose1+serum_glucose2; {in mmol/l}
         intercellulair_glucose_muscle_adipose := intercellulair_glucose_muscle_adipose_mg_dl*10/180.15588; {in mmol/l}
         intercellulair_glucose_skin := intercellulair_glucose_skin_mg_dl*10/180.15588; {in mmol/l}


        R_C13_C12 :=  1000000*serum_glucose1/serum_glucose0;

    end;

  procedure glucagon_dynamics;
  var glucagon_production, glucagon_decay, difference:real;
  const serum_glucose_criterion = 4.5;
  begin
    difference :=  serum_glucose-   serum_glucose_criterion  ;
    if difference<0 then glucagon_production:= 1*exp(-difference)*alfa_cells_percentage else glucagon_production := 0;
    glucagon_decay:= serum_glucagon/5; {halflife 3-6 min}
    serum_glucagon := serum_glucagon + (glucagon_production-glucagon_decay)*dt;
  end;

  procedure serum_glucose_dynamics;
  begin
  if glycogene_store <100 then
    begin
    glycogene_store0 := glycogene_store0 +(0.01*serum_glucose0*serum_insulin_mU_l/50  - 0.01*(glycogene_store0/glycogene_store)*serum_glucagon)*dt;
        glycogene_store1 := glycogene_store1 +(0.01*serum_glucose1*serum_insulin_mU_l/50  - 0.01*(glycogene_store1/glycogene_store)*serum_glucagon)*dt;
        glycogene_store2 := glycogene_store2 +(0.01*serum_glucose2*serum_insulin_mU_l/50  - 0.01*(glycogene_store2/glycogene_store)*serum_glucagon)*dt;
    glycogene_store:= glycogene_store0+ glycogene_store1+glycogene_store2;
    end;
     serum_glucose0 := serum_glucose0 +  0.01*(glycogene_store0/glycogene_store)*serum_glucagon*dt;
        serum_glucose1 := serum_glucose1 + 0.01*(glycogene_store1/glycogene_store)*serum_glucagon*dt;
        serum_glucose2 := serum_glucose2 +  0.01*(glycogene_store2/glycogene_store)*serum_glucagon*dt;
     serum_glucose:= serum_glucose0+serum_glucose1+serum_glucose2;
     Gluc := serum_glucose*10*Vg;
  end;
end.

