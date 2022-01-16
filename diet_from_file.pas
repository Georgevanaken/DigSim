unit diet_from_file;

{$mode objfpc}{$H+}
interface
uses
  Classes, SysUtils,globalvariables, calc_Protein_hydrogen;
type
  moment = array[1..50] of real;
  dishmoment = array[1..5] of real;
  integermoment = array[1..50] of integer;
  decomposer =  array[1..50,1..2] of integer;
  fraction =  array[1..50,1..2] of real;
  composition= array[1..50,1..110] of real;
  presence = array[1..110] of Boolean;

var
  compoundnumber, mealnumber: integer;
  newmealnumber: integer;
  compound,compoundcomposition: composition;
  teller:integer;
  DataFile : Text;
  mealtime, mealamount, compoundpH, compoundviscosity, compoundGprime, compoundstrain, swallow_speed_meal, oral_dec_time, gastric_dec_time, intestinal_dec_time, compoundfundus, compounddensity  : moment;
  meal, dish, solidfundus : integermoment;
  eating_speed, pHdish, liquidviscositydish: {dish}moment;
  number_of_meals, number_of_compounds, dishnumber: integer;
  oral_dec_comp, oral_dec_frac,gastric_dec_comp, gastric_dec_frac, intestinal_dec_comp, intestinal_dec_frac: decomposer;
  active: presence;
  FileName, TextFile: String;

procedure read_diet;
Procedure compoundreader;

implementation

procedure read_diet;
var hulp :integer;
    hulpreal: real;
    C :char;
 begin
 Assign(DataFile, 'meal' + '.txt');

 Reset(DataFile); {'Reset(x)' - means open the file x}
 for compoundnumber := 0 to 50 do for j := 1 to 110 do
     begin
     compoundcomposition[compoundnumber,j]:=0;
     oral_dec_frac[compoundnumber,1]:=0;
     oral_dec_frac[compoundnumber,2]:=0;
     gastric_dec_frac[compoundnumber,1]:=0;
     gastric_dec_frac[compoundnumber,2]:=0;
     intestinal_dec_frac[compoundnumber,1]:=0;
     intestinal_dec_frac[compoundnumber,2]:=0;
     end;
 readln(DataFile, TextFile);
  C:= ' '; While C<>')' do Read (Datafile,C); {routine die characters en open spaces leest tot en met ")" }
 Readln(Datafile,referencetime);  {de echte tijd in minuten vanaf middernacht op minutes = 0}
 readln(DataFile, Textfile);
 readln(DataFile, Textfile);
 readln(DataFile, Textfile);

{ C:= ' '; While C<>')' do Read (Datafile,C);} {routine die characters en open spaces leest tot en met ")" }
 mealnumber := 0;
 minutesmax := -1;
 Repeat
    begin
    Read(DataFile,hulpreal); {tijdstip maaltijd}
     if  (hulpreal > minutesmax) and (hulpreal <> 9999) then
     begin
     minutesmax := hulpreal;
     mealnumber := mealnumber + 1;
     Read(DataFile,meal[mealnumber]); { writeln(compoundnumber); readln;}
       mealtime[mealnumber] := hulpreal;
       Read(DataFile,mealamount[mealnumber]); {g}
       Read(DataFile,dish[mealnumber]); {1,2,3,4,...; 0 is intravenous}
          dishnumber:= dish[mealnumber];
       Read(DataFile,eating_speed[mealnumber]); {ml/min}
       {writeln(mealtime[mealnumber]:5:2, mealnumber:4, meal[mealnumber]:5,mealamount[mealnumber]:10:2,dish[mealnumber]:3, eating_speed[mealnumber]:7:2); readln;  }
     end;
     end;
 until hulpreal = 9999; {eindcriterium maaltijd nummers}
 number_of_meals := mealnumber;
 minutesmax := 500+mealtime[mealnumber]; { het programma rekent 500 min door na het startmoment van de laatste maaltijd}
{        writeln(minutesmax);readln;}

 readln(DataFile, TextFile);

 number_of_compounds:= 0;
 C:= ' '; While C<>')' do Read (Datafile,C); {routine die characters en open spaces leest tot en met ")" }

 Read(DataFile,hulp);

 while hulp <> 9999 do {eindcriterium maaltijd nummers}
  begin
      number_of_compounds := hulp;
      Read(DataFile,hulp);
  end;


 C:= ' '; While C<>')' do Read (Datafile,C);
 {******} for compoundnumber := 1 to number_of_compounds do begin {writeln( compoundnumber, ' ');} Read(DataFile, compoundpH[compoundnumber]); end; readln(DataFile, textfile);

 C:= ' '; While C<>')' do Read (Datafile,C);
 for compoundnumber := 1 to number_of_compounds do Read(DataFile, solid[compoundnumber]);  readln(DataFile, textfile);

 C:= ' '; While C<>')' do Read (Datafile,C);
 for compoundnumber := 1 to number_of_compounds do Read(DataFile, compoundviscosity[compoundnumber]);  readln(DataFile, textfile);

 C:= ' '; While C<>')' do Read (Datafile,C);
 for compoundnumber := 1 to number_of_compounds do Read(DataFile, oral_dec_time[compoundnumber]);  readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, oral_dec_comp[compoundnumber,1]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, oral_dec_frac[compoundnumber,1]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, oral_dec_comp[compoundnumber,2]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, oral_dec_frac[compoundnumber,2]); readln(DataFile, textfile);

 C:= ' '; While C<>')' do Read (Datafile,C);
 for compoundnumber := 1 to number_of_compounds do Read(DataFile, gastric_dec_time[compoundnumber]);  readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, gastric_dec_comp[compoundnumber,1]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, gastric_dec_frac[compoundnumber,1]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, gastric_dec_comp[compoundnumber,2]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, gastric_dec_frac[compoundnumber,2]); readln(DataFile, textfile);

 C:= ' '; While C<>')' do Read (Datafile,C);
 for compoundnumber := 1 to number_of_compounds do Read(DataFile, intestinal_dec_time[compoundnumber]);  readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, intestinal_dec_comp[compoundnumber,1]);readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, intestinal_dec_frac[compoundnumber,1]); readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, intestinal_dec_comp[compoundnumber,2]);readln(DataFile, textfile);
     C:= ' '; While C<>')' do Read (Datafile,C);
     for compoundnumber := 1 to number_of_compounds do Read(DataFile, intestinal_dec_frac[compoundnumber,2]);readln(DataFile, textfile);

 Readln(DataFile,textfile);
 Readln(DataFile,textfile);
 Readln(DataFile,textfile);

 {inlezen compound samenstellingen}
 For j := 1 to comp do for compoundnumber := 1 to number_of_compounds do
     begin;
     compoundcomposition[compoundnumber,j] := 0;
     active[j]:= False;
     end;
 repeat
    Read(DataFile, j);
     if j<>0 then
      begin
      C:= ' '; While C<>')' do Read (Datafile,C); {routine die characters en open spaces leest tot en met ")" }
     for compoundnumber := 1 to number_of_compounds do
           begin
           Read(DataFile,compoundcomposition[compoundnumber,j]); compoundcomposition[compoundnumber,j]:=0.01*compoundcomposition[compoundnumber,j] ;
           if compoundcomposition[compoundnumber,j] <> 0 then active[j]:= True;
           end;
      {repeat read(DataFile, textfile);until EOLn(DataFile);}
      readln(DataFile, textfile);
      {readln;}
      end;
 until j=0;

 C:= ' '; While C<>')' do Read (Datafile,C); {routine die characters en open spaces leest tot en met ")" }
     for compoundnumber := 1 to number_of_compounds do
           begin
           Read(DataFile,compounddensity[compoundnumber]);
           end;
 Close(DataFile);
 {***********}
 for compoundnumber:= 1 to number_of_compounds do
 begin
      calories_compound[compoundnumber] := 0;
      for j := 1 to comp do calories_compound[compoundnumber] := calories_compound[compoundnumber] + caloric_value[j]*compoundcomposition[compoundnumber,j];

 end;

 {***********}
 end;

Procedure compoundreader;

begin
         {eating_speed := eating_speed_meal[newmealnumber]; }{writeln(pH[1]);writeln(eating_speed_compound[newcompoundnumber]);readln;}
         swallow_speed := 1000;
         compoundnumber :=  meal[newmealnumber];
           dishnumber := dish[newmealnumber];
           pHdish[dishnumber] := compoundpH[compoundnumber]; {  writeln(pH[1]);readln;}{de maaltijd pH is de pH van de vloeibare fase}
           liquidviscositydish[dishnumber] := compoundviscosity[compoundnumber]; {de maaltijd viscositeit wordt voornamelijk bepaald door de vloeibare fase}
           eating_speed_dish[dishnumber] := eating_speed[newmealnumber];
           compounddish[dishnumber,compoundnumber] := compounddish[dishnumber,compoundnumber] + mealamount[newmealnumber];
             writeln('mealnr = ', newmealnumber, '  minutes =', minutes:8:0,'  comp. nr. = ', compoundnumber);



end;

end.



