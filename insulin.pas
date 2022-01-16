program insulin;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, crt, globalvariables,SubcellularInsulin;
var teller:integer;
    substeps:integer;


begin

  betacelverdeling;

  assign(out, 'Insulin_Results.txt');
  rewrite(out);
  write(out, 'minutes, serum glucose, f1');
  writeln(out);
    {serum_glucose := 100.0; }{mg/dl}
    {serum_insulin:= 100;}  {mg/dl}

  minutes := 0;
  dt:= 0.01;
  teller := 0;
  substeps := Round(1/dt)div 10;
  clrscr;
  repeat
    begin
      if (minutes>5) {and (minutes<5.2)} then G := 500 else G :=0.001;
      teller:= teller+1;
      minutes:= minutes+dt;
      calc_insulin_transport;

      if (teller = substeps) then
       begin
         teller := 0;
         writeln(minutes:8:0,', ', G:12:8,', ',f1:20:15);
         write(out, minutes:12:8,',',G:12:8,',',f1:20:15);
         writeln(out);
       end;
    end
  until minutes > 100;
  {readln;}
  close(out);

  assign(out, 'betacellverdeling.txt');
  rewrite(out);
  write(out, 'i, phi[i], dphi[i]');
  writeln(out);
  for i:=0 to 100 do
       begin
         write(out, i:8,',',phi[i]:12:8,',',dphi[i]:12:8);
         writeln(out);
       end;
  close(out);
end.

