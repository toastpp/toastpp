function [cvec] = readcplxtoastsol(Filein)
fid = fopen(Filein,'r');
vec = [];
cvec = [];
str = fgets(fid)
eof = 0;
ns = 0;
while(eof ==0)
    ns = ns +1
    fscanf(fid,'%c',1);
    i=1;
    eol = 0;
    while(eol==0)
        fscanf(fid,'%c',1);
        d = fscanf(fid,'%f %f');
        %size(d)
        vec(i,:) = [d(1) d(2)];
        i = i+1;
        str =fscanf(fid,'%c',2);
        if(str(2)==']')
            eol = 1;
        end
    end
    cvec(:,ns) = vec(:,1) + sqrt(-1)*vec(:,2);
    fscanf(fid,'%c',1);
    str = fgets(fid)
    if(str == -1)
        eof = 1
    end
end
fclose(fid);
size(vec)
