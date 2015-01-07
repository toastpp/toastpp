function saveinifile(A, filename);

%
% saves a matlab structure as a text file
% format is 
%    field: value
%
% returns a structure of fields and values
%
% 21/02/05
%


names=fieldnames(A);
fid=fopen(filename,'w+t');
for ctr = 1:length(names)
  fprintf(fid,'%s:%s',names{ctr},getfield(A,names{ctr}));
  fprintf(fid,'\n');
end
fclose(fid);

