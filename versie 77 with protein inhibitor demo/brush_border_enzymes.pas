unit Brush_border_enzymes;

{$mode objfpc}{$H+}


interface

uses
  Classes, SysUtils, globalvariables;

var sucrase_presence: segment;
    maltase_presence: segment;
    Rs, R1: real;
const drag = 1;


procedure sucrase;
procedure maltase;

implementation

var sucrase_activity,maltase_activity,A,B,C,  effective_concentration, sucrose_absorption_rate, sucrose_concentration_lumen,Sucrase_transformation_overall, sucrose_concentration_brush, Maltose_absorption_rate, Maltose_concentration_lumen,Maltase_transformation_overall, Maltose_concentration_brush  : real;
const Km_sucrase = 8.56E-3;    {in gram/ml}  {based on plots in Seri et al. 1996}
      Km_maltase = 0.01;    {in gram/ml}


procedure sucrase;
begin


  for i:= 1 to 15 do sucrase_presence[i]:= 0;             {* hier wordt de verdeling van sucrase over de dunne darm beschreven}
       sucrase_presence[5] := 2;             {*}             {de preciese verdeling moet nog uit de literatuur afgeleid worden}
       sucrase_presence[6] := 1.5;             {*}
       sucrase_presence[7] := 1.2;             {*}
       sucrase_presence[8] := 1;             {*}
       sucrase_presence[9] := 1;             {*}
       sucrase_presence[10] := 0.5;             {*}
       sucrase_presence[11] := 0.5;             {*}
       sucrase_presence[12] := 0.5;             {*}
       sucrase_presence[13] := 0.5;             {*}
       sucrase_presence[14] := 0.5;             {*}
       sucrase_presence[15] := 0.5;             {*}

  For i:= 5 to 15 do
    begin
          sucrase_activity := 0.15{0.18}{0.18}{0.17}*surface_factor[i];        {g/min per darmsegment}
          {sucrase_activity := 0.04{0.18}{0.17}*surface_factor[i];}    {nu als isomaltulose}
          {acompetitive inhibition}
          Ia[i] :=1+1*concentration_brush[i,25]/6.85E-4+concentration_brush[i,26]/2E-5+0*concentration_brush[i,22]/6.85E-4; {acompetitive inhibitie door arabinose en xylose}{Michaelis Menten equation for sucrose hydrolysis by sucrase} {Ki arbinose  = 6.85E-4 g/ml; Ki acarbose = 7.1E-7 g/ml volgens Seri et al 1996}
          {competitive inhibition}
          Ic[i] :=1+concentration_brush[i,27]/1E-6;   {competitive inhibitie door acarbose}
         for j:= 1 to comp do sucrase_transformation[i,j] := 0;
            sucrose_absorption_rate:= absorption_rate[i,21]+absorption_rate[i,36];
            sucrose_concentration_lumen:= concentration_lumen[i,21]+concentration_lumen[i,36];

         Rs:= surface_factor[i]*sucrose_absorption_rate;
         R1:= surface_factor[i]*absorption_rate[i,1];
         A:= -Ia[i];
         B:= -Km_sucrase*Ic[i] + sucrose_concentration_lumen*Ia[i] + (Delta[i]/D[21])*(-sucrase_presence[i]*sucrase_activity -Rs*Ia[i] +drag*sucrose_concentration_lumen*R1*Ia[i]);
         C:= Km_sucrase*Ic[i]*sucrose_concentration_lumen+ (Delta[i]/D[21])*(-Rs*Km_sucrase*Ic[i] + drag*sucrose_concentration_lumen*Km_sucrase*Ic[i]*R1);
         sucrose_concentration_brush := (-B-SQRT(B*B-4*A*C))/(2*A);
         Sucrase_transformation_overall{in g/segment.min} := (1/surface_factor[i])*sucrase_presence[i]*sucrase_activity*sucrose_concentration_brush/(Km_sucrase*Ic[i]+sucrose_concentration_brush*Ia[i]);
         if (sucrose_concentration_brush >0) and (sucrose_concentration_lumen> 0) then
           begin
            concentration_brush[i,21]:= (concentration_lumen[i,21]/sucrose_concentration_lumen)*sucrose_concentration_brush;
            concentration_brush[i,36]:= (concentration_lumen[i,36]/sucrose_concentration_lumen)*sucrose_concentration_brush;
            sucrase_transformation[i,21]:=(concentration_brush[i,21]/sucrose_concentration_brush)*Sucrase_transformation_overall;
            sucrase_transformation[i,36]:=(concentration_brush[i,36]/sucrose_concentration_brush)*Sucrase_transformation_overall;
            sucrase_inhibitie[i]:= (Km_sucrase*Ic[i]+sucrose_concentration_brush*Ia[i])/(Km_sucrase+sucrose_concentration_brush);
           end
         else
           begin
                concentration_brush[i,21] := 0;
                concentration_brush[i,36] := 0;
                sucrase_transformation[i,21] :=0;
                sucrase_transformation[i,36] :=0;
                sucrase_inhibitie[i] := 1;
           end;


    end;                                                                               { uncompetitive inhibition:10.492 - Integrated Chemical Engineering (ICE) Topics: Biocatalysis - MIT Chemical Engineering Department         {25 = arabinose, 26 = xylose, 27 acarbose}                                 {MW acarbose = 645.61 g·mol−1}                                                                                                        ,
                                                                                       Instructor: Professor Kristala Prather - Fall 2004 - Lecture # 5, 6 – Enzyme Inhibition and Toxicity }

end;


procedure maltase;
begin
     maltase_activity := 0.2{1}{0.3}{0.25}*surface_factor[i];        {g/min per darmsegment}
     {acompetitive inhibition}
     Ia[i] :=1*(1+concentration_brush[i,26]/2E-5); {acompetitive inhibitie door arabinose en xylose}{Michaelis Menten equation for sucrose hydrolysis by sucrase} {Ki arbinose  = 6.85E-4 g/ml; Ki acarbose = 7.1E-7 g/ml volgens Seri et al 1996}
     {competitive inhibition}
     Ic[i] :=1+concentration_brush[i,27]/7E-7;   {competitive inhibitie door acarbose}

  for i:= 1 to 15 do maltase_presence[i]:= sucrase_presence[i];             {* hier wordt de verdeling van maltase over de dunne darm beschreven} {de verdeling gaat samen met die van sucrase}

  For i:= 5 to 15 do
    begin
         for j:= 1 to comp do maltase_transformation[i,j] := 0;
            Maltose_absorption_rate:= absorption_rate[i,23]+absorption_rate[i,24];
            Maltose_concentration_lumen:= concentration_lumen[i,23]+concentration_lumen[i,24];

         Rs:= surface_factor[i]*Maltose_absorption_rate;
         R1:= surface_factor[i]*absorption_rate[i,1];
         A:= -Ia[i];
         B:= -Km_maltase*Ic[i]+ Maltose_concentration_lumen*Ia[i] + (Delta[i]/D[23])*(-maltase_presence[i]*maltase_activity -Rs*Ia[i] +drag*Maltose_concentration_lumen*R1*Ia[i]);
         C:= Km_maltase*Ic[i]*Maltose_concentration_lumen+ (Delta[i]/D[23])*(-Rs*Km_maltase*Ic[i] + drag*Maltose_concentration_lumen*Km_maltase*Ic[i]*R1);
         Maltose_concentration_brush:= (-B-SQRT(B*B-4*A*C))/(2*A);
         Maltase_transformation_overall {in g/segment.min} := (1/surface_factor[i])*maltase_presence[i]*maltase_activity*Maltose_concentration_brush/(Km_maltase*Ic[i]+Maltose_concentration_brush*Ia[i]);
         if (Maltose_concentration_brush >0) and (Maltose_concentration_lumen>0) then
           begin
            concentration_brush[i,23]:= (concentration_lumen[i,23]/Maltose_concentration_lumen)*Maltose_concentration_brush;
            concentration_brush[i,24]:= (concentration_lumen[i,24]/Maltose_concentration_lumen)*Maltose_concentration_brush;
            Maltase_transformation[i,23]:=(concentration_brush[i,23]/Maltose_concentration_brush)*maltase_transformation_overall;
            Maltase_transformation[i,24]:=(concentration_brush[i,24]/Maltose_concentration_brush)*maltase_transformation_overall;
            Maltase_inhibitie[i]:= (Km_maltase*Ic[i]+Maltose_concentration_brush*Ia[i])/(Km_maltase+Maltose_concentration_brush);
           end
         else
           begin
                concentration_brush[i,23] := 0;
                concentration_brush[i,24] := 0;
                Maltase_transformation[i,23] :=0;
                Maltase_transformation[i,24] :=0;
           end;


         {***}
    end;                                                                                                                            {27 is acarbose}
end;

end.

