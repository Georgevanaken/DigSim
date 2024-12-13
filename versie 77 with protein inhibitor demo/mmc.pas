unit MMC;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables, glucose_homeostasis2;

procedure update_MMC_status;



implementation
var start_MMC,start_FED: Boolean;


procedure update_MMC_status;
var Criterion_start_fasted, criterion_start_fed : Boolean;
begin
   if minutes > 0 then
   begin
   {Criterion_start_fasted geeft de overgang naar fasted state aan, Criterion_start_fed geeft de overgang naar de FED state aan }
   {criteria start en stop MMC en FED state}
   if {(CCK < 2) and (pH[5] > 4) and (GLP1 < 2) and (serum_insulin_mU_l < 20) \} motilin> 40{20} then Criterion_start_fasted := True else Criterion_start_fasted := False;
   if {(CCK > 4) or (pH[5] < 3) or(GLP1 >  1) or (serum_insulin_mU_l > 25)} motilin < 40{10} then Criterion_start_fed := True else Criterion_start_fed := False;
   {end crirteria start and stop MC en FED state}

   if Criterion_start_fasted and (MMC_status = 0) then begin MMC_status := 1; MMC_time[1] := 0; end;
   if Criterion_start_fed and (MMC_status <> 0) then Begin MMC_Status := 0; MMC_time[0] := 0; end;
   case MMC_status of            {geeft de MMC status van de GIT; 0 is fed state, 1,2,3,4 de opeenvolgende fasen}
    0: MMC_time[0] := MMC_time[0] + dt;
    1: if MMC_time[1] > 50 then begin MMC_status := 2; MMC_time[1] := 0 end else MMC_time[1] := MMC_time[1] + dt;
    2: if MMC_time[2] > 20 then begin MMC_status := 3; MMC_time[2] := 0 end else MMC_time[2] := MMC_time[2] + dt;
    3: if MMC_time[3] > 7 then begin MMC_status := 4; MMC_time[3] := 0 end else MMC_time[3] := MMC_time[3] + dt;
    4: if MMC_time[4] > 10 then begin MMC_status := 1; MMC_time[4] := 0 end else MMC_time[4] := MMC_time[4] + dt;
   end;
   {if start_FED then begin MMC_status := 0; MMC_time[0] := 0; end;}
   end;
end;

end.
