unit Hb_glycosylation;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, globalvariables;
procedure calc_HbA1c;

implementation

procedure calc_HbA1c;
          begin
          HbA1c_halflife := 3*30*24*60; {3 months in minutes};
          Glycosylation_rate_constant := 1.23*1E-8;
          HbA1c := HbA1c + (Gluc*Glycosylation_rate_constant - HbA1C/(2*HbA1c_halflife))*dt;       {HbA1c in % of total Hb}
          {writeln('HbA1c = ', HbA1c:7:5);  }
          end;
end.

