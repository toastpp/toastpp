% Returns a relative path from srcpath to file, where both srcpath and
% file contain absolute paths (or at least relative paths starting from
% the same node)

function rpath = relpath(file,srcpath)

fs = filesep;
[h1,r1]=strtok(srcpath,fs);
[h2,r2]=strtok(file,fs);

% step 1: remove identical leading components
while strcmp (h1,h2)
    srcpath = r1;
    file = r2;
    if length(srcpath) == 0
        break;
    end
    if length(file) == 0
        break;
    end
    [h1,r1] = strtok(srcpath,fs);
    [h2,r2] = strtok(file,fs);
end

% step 2: step up through the remaining directories of srcpath
rpath = '';
if length(srcpath)
    while length(srcpath)
        if length(rpath)
            rpath = [rpath fs '..'];
        else
            rpath = '..';
        end
        [h1,srcpath] = strtok(srcpath,fs);
    end
else
    % srcpath has been dissolved entirely
    if file(1) == fs
        file = file(2:end); % remove leading path delimiter
    end
end
    

% step 3: append remaining file path
rpath = [rpath file];