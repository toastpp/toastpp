% This script adds the TOAST mex and script directories to the
% Matlab path.
% To add the path permanently, you can append the contents of this
% file to your startup.m file, located in our Matlab startup
% directory, or in <matlabroot>/toolbox/local.

function mtoast_install(varargin)

nargin = length(varargin);
nogui = nargin > 0 && varargin{1} == true;

if ~nogui
    fprintf(1,'\nAdding search paths for the Toast++ toolbox.\n\n')
end

% figure out path structure
toastdir = getenv('TOASTDIR');
if length(toastdir) == 0
    [toastdir,name,ext] = fileparts(which('mtoast2_install.m'));
end
toastver = getenv('TOASTVER');
if length(toastver) == 0
    arch = computer;
    switch arch
        case 'PCWIN'
            toastver = [toastdir '\win\Win32\Release'];
        case 'PCWIN64'
            toastver = [toastdir '\win\x64\Release'];
        case 'MACI64'
            toastver = [toastdir '/darwin64'];
        otherwise
            disp('Warning: could not determine location of mex files')
            disp('Please edit the matlab path manually')
    end
end

% Remove current toast references from path
fprintf('Searching for existing toast path entries ...\n');
p = path;
p(find(p==' ')) = '^';
p(find(p==pathsep)) = ' ';
p = textscan(p,'%s');
p = p{1};
for i = 1:size(p,1)        % restore spaces
    p{i}(find(p{i}=='^')) = ' ';
end

k = strfind(p,'toast');
nfound = 0;
for i=1:length(k)
    if length(k{i}) > 0
        fprintf('Removing search path %s\n', p{i});
        rmpath(p{i});
        nfound = nfound+1;
    end
end
if nfound > 0
    fprintf('Removed %d existing toast path entries\n', nfound);
else
    fprintf('No existing toast paths found.\n');
end

if ~nogui
    fprintf ('\nTOAST root directory: %s\n', toastdir);
    fprintf ('TOAST arch directory: %s\n\n', toastver);
end

% Sanity checks: assert valid file structure
assertfile([toastdir '/mtoast2_install.m']);
assertdir([toastdir '/script']);
assertdir([toastdir '/script/matlab']);
assertdir([toastdir '/script/matlab/toast2']);
assertdir([toastdir '/script/matlab/utilities']);
assertdir([toastver '/mex2']);
assertfile([toastver '/mex2/toastmex.' mexext]);

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
mexp = [toastver '/mex2'];
addpath (mexp);
if ~nogui
    disp(['Adding search path ' mexp])
end

rmpath ([toastdir '/script/matlab/toast']); % remove original toast script directory

if nargin == 0 || nogui==false
    % open path tool GUI to allow user to modify and save the toast paths
    fprintf(1,'\nPlease check that the toast paths are set correctly,\n')
    fprintf(1,'then save to store the paths permanently.\n')
    pathtool
end


function assertfile(pathstr)
if exist(pathstr,'file') ~= 2
    error('\nToast toolbox file structure mismatch. Expected file not found:\n%s\n\nDetected toast root folder:   %s\nDetected toast binary folder: %s\n\nPlease check your toast installation.', pathstr, toastdir, toastver);
end
end

function assertdir(pathstr)
if exist(pathstr,'dir') ~= 7
    error('\nToast toolbox file structure mismatch. Expected directory not found:\n%s\n\nDetected toast root folder:   %s\nDetected toast binary folder: %s\n\nPlease check your toast installation.', pathstr, toastdir, toastver);
end
end

end