% This script adds the TOAST mex and script directories to the
% Matlab path.
% To add the path permanently, you can append the contents of this
% file to your startup.m file, located in our Matlab startup
% directory, or in <matlabroot>/toolbox/local.

function mtoast_install(varargin)

nargin = length(varargin);
nogui = nargin > 0 && varargin{1} == true;

if ~nogui
    fprintf(1,'\nAdding search paths for TOAST Matlab scripts and MEX files.\n\n')
end

% figure out path structure
toastdir = getenv('TOASTDIR');
if length(toastdir) == 0
    toastdir = pwd;
end
toastver = getenv('TOASTVER');
if length(toastver) == 0
    arch = computer;
    switch arch
        case 'PCWIN'
            toastver = [toastdir '\win\Win32\Release'];
        case 'PCWIN64'
            toastver = [toastdir '\win\x64\Release'];
        otherwise
            disp('Warning: could not determine location of mex files')
            disp('Please edit the matlab path manually')
    end
end

% Remove current toast references from path
p = path;
p(find(p==' ')) = '^';
p(find(p==pathsep)) = ' ';
p = textscan(p,'%s');
p = p{1};
for i = 1:size(p,1)        % restore spaces
    p{i}(find(p{i}=='^')) = ' ';
end

found_toast = false;
k = strfind(p,'toast');
for i=1:length(k)
    if length(k{i}) > 0
        if ~found_toast
            if ~nogui
                fprintf(1,'Found existing toast directories in path:\n\n');
            end
            found_toast = true;
        end
        if ~nogui
            disp(p{i});
        end
    end
end
if found_toast
    if nogui
        inp = 'Y';
    else
        inp = input('\nDelete existing toast directories? Y/N [Y]: ','s');
        if isempty(inp)
            inp = 'Y';
        end
    end
    if strcmpi(inp,'Y')
        for i=1:length(k)
            if length(k{i}) > 0
                rmpath(p{i});
            end
        end
    end
    fprintf('\n');
end

if ~nogui
    fprintf ('TOAST root directory: %s\n', toastdir);
    fprintf ('TOAST arch directory: %s\n\n', toastver);
end

% Add all directories under the script/matlab node
p = genpath([toastdir '/script/matlab/']);
p(find(p==' ')) = '^';     % protect spaces in directory names
p(find(p==pathsep)) = ' '; % replace separators with spaces
p = textscan(p,'%s');      % split path into separate elements
p = p{1};
for i = 1:size(p,1)        % restore spaces
    p{i}(find(p{i}=='^')) = ' ';
end

k1=strfind(p,'CVS');        % eliminate CVS subdirs
k2=strfind(p,'.svn');       % eliminate .svn subdirs
k = [k1,k2];

pth = '';
for i=1:size(k,1)
    if length(k{i,1}) == 0 && length(k{i,2}) == 0
        if length(pth) > 0
            pth = [pth pathsep];
        end
        pth = [pth cell2mat(p(i))];
        if ~nogui
            disp(['Adding search path ' cell2mat(p(i))])
        end
    end
end
addpath (pth);

% add the mex directory
mexp = [toastver '/mex'];
addpath (mexp);
if ~nogui
    disp(['Adding search path ' mexp])
end

if nargin == 0 || nogui==false
    % open path tool GUI to allow user to modify and save the toast paths
    fprintf(1,'\nPlease check that the toast paths are set correctly,\n')
    fprintf(1,'then save to store the paths permanently.\n')
    pathtool
end
