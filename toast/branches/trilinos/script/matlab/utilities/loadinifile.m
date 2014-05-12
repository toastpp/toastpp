function A=loadinifile(filename);

%
% loads a text file into matlab
% format is 
%    field: value
%
% returns a structure of fields and values
%
% 21/02/05
%

A=struct('date',datestr(now));
fid=fopen(filename,'rt');
for line = 1:100
  clear tline L
  tline = fgetl(fid);
  if ~ischar(tline)|isempty(tline), break, end
  
  ctr=1;
  while tline(ctr) ~= ':'
    L(ctr)=tline(ctr);
    ctr=ctr+1;
  end
  ctr=ctr+1;
  while isspace(tline(ctr))
    ctr=ctr+1;
  end
  
  A=setfield(A,L,tline(ctr:end));
end
fclose(fid);
A.date=datestr(now);
