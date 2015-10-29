% Create the directory structure for the Toast++ Matlab toolbox and copy
% the relevant files.

% After running this, use Add-Ons | Package toolbox with the
% $TOASTDIR/Toast++.prj toolbox project file to actually package the
% toolbox for distribution.

% After creating the toolbox, the temporary toolbox directory can be
% deleted.

function make_toolbox

rootdir = getenv('TOASTDIR');
if isempty(rootdir)
    rootdir = pwd;
end
cd (rootdir);

% Clear any existing toolbox directory tree
[stat,mess,id] = rmdir('matlab_toolbox', 's');

% remove any pre-existing toast-related paths
display('Removing existing toast paths:');
remove_paths('toast');
remove_paths('Toast');

mkdir ('matlab_toolbox');
cd ('matlab_toolbox');
tbxdir = pwd;

% Copy all scripts
mkdir ('script');
cd ('script');
copyfile('../../script/matlab/*', '.');
movefile('html', '..');
movefile('info.xml', '..');
movefile('demos.xml', '..');

% Delete the old interface directory
rmdir('toast', 's');
cd ..

% Copy the examples directory
mkdir ('examples');
cd ('examples');
copyfile('../../examples/matlab/*', '.');
cd ..

% Create a directory for the platform-dependent binaries
switch computer
    case 'GLNXA64'
        tgtdir = 'glnxa64';
        bindir = [getenv('TOASTVER') '/mex2'];
        libdir = [getenv('TOASTVER') '/lib'];
    case 'PCWIN'
        tgtdir = 'win32';
        bindir = [rootdir '/win/Win32/Release/mex2'];
        libdir = bindir;
    case 'PCWIN64'
        tgtdir = 'win64';
        bindir = [rootdir '/win/x64/Release/mex2'];
        libdir = bindir;
    case 'MACI64'
        tgtdir = 'maci64';
        bindir = [rootdir '/darwin64'];
        libdir = bindir; % not sure
    otherwise
        error ('Unsupported platform/operating system');
end
        
mkdir (tgtdir);
mkdir (tgtdir, 'bin');
copyfile([bindir '/toastmex.*'], [tgtdir '/bin']);
copyfile([libdir '/*'], [tgtdir '/bin']);

% Add the required toast paths
addpath([tbxdir '/html']);
addpath([tbxdir '/script/demos']);
addpath([tbxdir '/script/gui']);
addpath([tbxdir '/script/utilities']);
addpath([tbxdir '/script/toast2']);
addpath([tbxdir '/' tgtdir '/bin']);

cd (rootdir);
fprintf('\nToast++ toolbox directory tree created.\n\n');
fprintf('Now package the toolbox using the "Add-Ons | Package toolbox"\n');
fprintf('option under the "Home | Resources" tab (Matlab 2014 and later).\n\n');
fprintf('Load project "Toast++.prj" and click "Package".\n\n');
fprintf('You can then delete the "matlab_toolbox" directory tree.\n\n');
fprintf('Before installing the toolbox (by double-clicking the\n');
fprintf('Toast++.mltbx file in the Matlab file manager), make sure that all\n');
fprintf('Toast-related paths are removed, e.g. by quitting Matlab and\n');
fprintf('restarting).\n');

end


% -------------------------------------------------------------------
function remove_paths(token)
p = path;
p(find(p==' ')) = '^';
p(find(p==pathsep)) = ' ';
p = textscan(p,'%s');
p = p{1};
for i = 1:size(p,1)        % restore spaces
    p{i}(find(p{i}=='^')) = ' ';
end

found_toast = false;
k = strfind(p,token);
for i=1:length(k)
    if length(k{i}) > 0
        found_toast = true;
        disp(p{i});
    end
end
if found_toast
    for i=1:length(k)
        if length(k{i}) > 0
            rmpath(p{i});
        end
    end
    fprintf('\n');
end
end

