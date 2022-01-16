unit SubcellularInsulin_Pederson;
{based on Pedersen et al. 2008}
{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables;

type betacellclass = array[0..1000] of real;
var  Intermediate , RRP,Fused,
     Intermediate_increment,RRP_increment,Fused_increment, released:betacellclass;
     phi,dphi : betacellclass;
     f1_Pederson :real;


procedure calc_f1_Pederson;
procedure betacelverdeling;
function sigmoid(BASIS,MAX,POWER,x,xhalf:real):real;


implementation


const  m_secretion = 0.62; {1/min}
       fplus = 6.2; {1/min}
       Pplus = 0.021; {1/min}
       Pmin = 0.1695{0.11}; {1/min}
       k =0.47; {1/min}
       r = 3*0.0023; {1/min}
       nM = 12;
       nH = 3.3;
       c = 1.25; {microgram/min}
       M0 =3*0.020; {microgram/min}
       KmM = 250{185}; {milligram/dl}
       G50 = 150{150}; {mg/dl glucose}
       tau = 0.5{30}; {min}
       imax =500; {aantal klassen beta cellen}
       f1min = 0.01; {rekentechnische minimale insuline productie, lager geeft analytische problemen}

function sigmoid(BASIS,MAX,POWER,x,xhalf:real):real;
begin
    sigmoid := BASIS+ MAX*exp(power*LN(x))/exp(power*LN(x)+exp(power*(LN(Xhalf))));
end;

Procedure Betacelverdeling;
  var i : integer;
  begin
    for i:= 1 to imax do Phi[i] := exp(nH*LN(i))/(exp(nH*LN(imax/2))+exp(nH*LN(i)));
    dPhi[1]:= 0;
    for i:= 2 to imax do dPhi[i] := Phi[i]-phi[i-1];
    for i := 1 to imax do
        begin
             intermediate[i] := 1*dphi[i];
             RRP[i] := {1.09}1.62*dphi[i];   {microgram insuline totaal in RRP}
             Fused[i] := 0;
        end;
  end;


procedure calc_f1_Pederson;      {insuline productie op basis van Pedersen 2008, G in mg/dl}
var total, mobilization,M_increment: real;
    i: integer;
begin
  total := 0;
  if (minutes = starttime) then {initialisatie van M; G moet hier op een initiele waarde gezet worden}
    begin
    G := 40; {initiele waarde van G in mg/dl}
    M := M0+c*exp(nM*LN(G))/(exp(nM*LN(KmM))+exp(nM*LN(G)))
    end
  else

    begin

  mobilization:= M0+c*exp(nM*LN(G))/(exp(nM*LN(KmM))+exp(nM*LN(G)));    {de rol van glucose (in mg/dl) hierin is dat het via de ATP/ADP balans transport reguleerd;
    let wel: het gaat daarbij om intracellulaire glucose dit door het GLUT2 kanaal de cel is binnengekomen}
    mobilization:= mobilization{*incretin}; {intracellulair effect van incretines}
  M_increment:= (mobilization-M)/tau;
  M:= M+M_increment*dt;
    end;
  for i:= 1 to imax do
      begin
      Intermediate_increment[i]:= M*dphi[i] -(r+Pplus)*Intermediate[i]+Pmin*RRP[i];
      RRP_increment[i]:= Pplus*Intermediate[i]-Pmin*RRP[i];
      if (2*i<imax*G/G50 ) then          {G in mg/dl; transitie bij i = 0.5 imax }
         begin
         Intermediate_increment[i]:= Intermediate_increment[i]+k*Fused[i];
         RRP_increment[i] := RRP_increment[i]-fplus*RRP[i];
         Fused_increment[i]:= fplus*RRP[i]-(m_secretion+k)*Fused[i];
         released[i] := m_secretion{*incretin}*Fused[i];    {in microgram/min}
         total:= total+released[i];
         end;
      end;

   f1_Pederson := total/0.0387 {f1 hier gecorrigeerd van microgram/min naar milliUnits/min};
   f1_Pederson := f1_Pederson; {callibratiefactor, te optimaliseren}
   if f1_Pederson < f1min then f1_Pederson := f1min;


   for i := 1 to imax do
    begin
      Intermediate[i] := Intermediate[i] + intermediate_increment[i]*dt;
      RRP[i] := RRP[i] + RRP_increment[i]*dt;
      Fused[i]:= Fused[i]+Fused_increment[i]*dt;
    end;
end;

end.
