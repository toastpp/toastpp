function [sourcepos, measpos, linklist] = loadqmfile(filename)

%  Will load a toast qm file 
%
%       [sourcepos, measpos,linklist] = loadqmfile(filename)
%
%  and return arrays containing positions of source optodes,
%  positions of measurement optodes and the linklist taken 
%  from the qm file, zero padded.
%
%  Adam 06/06/01
%


%-------------------------------------------------------
%  Open qm file


qmfile = fopen(filename,'r');
while qmfile == -1;
  [filename, pathname] = uigetfile({'*.qm', 'All qm files (*.qm)'}, 'QM file not found. Please reselect');
  qmfile = fopen(filename,'r');
end;
fseek(qmfile,0,'bof');

disp ([' - loading qm file ', filename]);


%-------------------------------------------------------
%  Read header and check
header = fgetl (qmfile);
if header(1:7) ~= 'QM file'
   error ([filename, ' is not a valid QM file']);
end;

% dump rest of header
dump = fgetl (qmfile);
while isletter(dump(1))
  dump = fgetl (qmfile);
  if isempty(dump) 
    dump = fgetl (qmfile);
  end
end


%-------------------------------------------------------
%  Read SourceList

sourcepos=[];

nextline = str2num(dump);
while ~isempty(nextline)
  sourcepos = [sourcepos;  nextline];
  nextline = str2num(fgetl(qmfile));
end


%-------------------------------------------------------
%  Read MeasurementList

dump = fgetl (qmfile);
while isletter(dump(1))
  dump = fgetl (qmfile);
  if isempty(dump) 
    dump = fgetl (qmfile);
  end
end

measpos=[];

nextline = str2num(dump);
while ~isempty(nextline)
  measpos = [measpos;  nextline];
  nextline = str2num(fgetl(qmfile));
end

%-----------------------------------------------------------
%   Read LinkList

dump = fgetl (qmfile);
while isletter(dump(1))
  dump = fgetl (qmfile);
  if isempty(dump) 
    dump = fgetl (qmfile);
  end
end

linklist=-1 * [ones(size(sourcepos,1),size(measpos,1))];
%ctr=1;

nextline = str2num(dump(4:end));
linklist(1,(1:length(nextline))) = nextline;
%while ~feof(qmfile)
for ctr = 2 : length(sourcepos)
  N = fgetl(qmfile);
  try
    nextline = str2num(N(4:end));
  catch
    break;
  end
  linklist(ctr,(1:length(nextline))) = nextline;
  %  ctr=ctr+1;
end
