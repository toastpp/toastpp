function data = loaddatafile(filename);
%
% data = loaddatafile(filename);
% 
% Load data types in TOAST .fem format. Don't use this - better to use 
% loaddatatypes.m which calls this function.
%
% File format is a string of ascii floats, space delimited with [ and ] 
% before and after.

%-------------------------------------------------------
%  Open data file

datafile = fopen(filename,'r');
while datafile == -1;
  [filename, pathname] = uigetfile({'*.dat;*.int;*.mean;*.var', 'All data files'}, 'Data file not found. Please reselect');
  datafile = fopen(filename,'r');
end;

disp ([' - loading data file ', filename]);

%----------------------------------------------------------
% Read data. Data is all one line, so read a one line.

D = fgetl(datafile);

%----------------------------------------------------------
% Now lose leading and trailing []

data = str2num(D(2:length(D)-1));

fclose(datafile);
