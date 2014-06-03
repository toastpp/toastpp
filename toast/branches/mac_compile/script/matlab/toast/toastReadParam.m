function prm = toastReadParam(prmfile)
%toastReadParam       - Read a TOAST parameter file.
%
% Synopsis: prm = toastReadParam(prmfile)
%    prmfile: parameter file name (string)
%    prm:     parameter structure
%
% This function reads the toast reconstruction parameters used by
% toastRecon from a file and stores them in a structure that can
% be passed directly to toastRecon.  
  
prm.meas.qmfile = scanprm(prmfile,'QMFILE');
prm.meas.src.type = scanprm(prmfile,'SOURCETYPE');
prm.meas.src.prof = scanprm(prmfile,'SOURCEPROFILE');
prm.meas.src.width = str2num(scanprm(prmfile,'SOURCEWIDTH'));
prm.meas.det.prof = scanprm(prmfile,'MEASUREMENTPROFILE');
prm.meas.det.width = str2num(scanprm(prmfile,'MEASUREMENTWIDTH'));

% data and reference data
prm.data.freq = str2num(scanprm(prmfile,'FREQ'));
prm.data.lnampfile = scanprm(prmfile,'DATA_MOD');
prm.data.phasefile = scanprm(prmfile,'DATA_ARG');
prm.data.useref = logical(str2num(lower(scanprm(prmfile,'REFDATA'))));
if length(prm.data.useref) == 0
    prm.data.useref = false;
end
if prm.data.useref
    prm.data.ref.lnampfile = scanprm(prmfile,'REFDATA_MOD');
    prm.data.ref.phasefile = scanprm(prmfile,'REFDATA_ARG');
end

prm.fwdsolver.meshfile = scanprm(prmfile,'MESHFILE');
prm.fwdsolver.method = scanprm(prmfile,'LINSOLVER');
prm.fwdsolver.tol = str2num(scanprm(prmfile,'LINSOLVER_TOL'));
if length(prm.fwdsolver.tol) == 0
    prm.fwdsolver.tol = 1e-10;
end

tmp = scanprm(prmfile,'BASIS');
if length(tmp) > 0
    prm.solver.basis.bdim = str2num(tmp);
end
tmp = scanprm(prmfile,'GRID');
if length(tmp) > 0
    prm.solver.basis.gdim = str2num(tmp);
end
tmp = scanprm(prmfile,'SOLVER');
if length(tmp) > 0
    prm.solver.method = tmp;
end
tmp = scanprm(prmfile,'NONLIN_TOL');
if length(tmp) > 0
    prm.solver.tol = str2num(tmp);
end
tmp = scanprm(prmfile,'NONLIN_ITMAX');
if length(tmp) > 0
    prm.solver.itmax = tmp;
end
tmp = scanprm(prmfile,'DATA_SCALING');
if length(tmp) > 0
    prm.solver.dscale = tmp;
end
tmp = scanprm(prmfile,'LM_LINESEARCH');
if length(tmp) > 0
    if strcmpi(tmp,'True') == 1
        prm.solver.lsearch = true;
    else
        prm.solver.lsearch = false;
    end
end

% regularisation parameters
prm.regul.method = scanprm(prmfile,'PRIOR');
if length(prm.regul.method) == 0
    prm.regul.method = 'None';
end
if strcmpi(prm.regul.method,'None') == 0
    prm.regul.tau = str2num(scanprm(prmfile,'PRIOR_TAU'));
    switch prm.regul.method
        case 'TV'
            prm.regul.tv = struct('beta', scanprm(prmfile,'TV_BETA'));
    end
    prm.regul.prior = struct('refname',scanprm(prmfile,'PRIOR_KAPREFIMG'), ...
                             'smooth',scanprm(prmfile,'PRIOR_DIFFSCALE'), ...
                             'threshold',scanprm(prmfile,'PRIOR_PMTHRESHOLD'));
    if length(prm.regul.prior.refname) == 0 || strcmpi(prm.regul.prior.refname,'none')
        prm.regul = rmfield(prm.regul,'prior');
    end
end

% initial optical parameters
prm.initprm.mua = setinitprm(scanprm(prmfile,'RESET_MUA'));
prm.initprm.mus = setinitprm(scanprm(prmfile,'RESET_MUS'));
prm.initprm.ref = setinitprm(scanprm(prmfile,'RESET_N'));


%% ------------------------------------------------------------------
% scans the parameter file for a particular tag and returns its value
function prm = scanprm (prmfile,item)
    
prm = '';
fid = fopen(prmfile);
if fid ~= -1
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        [tag,val] = scanline(tline);
        if strcmpi(tag,item) == true
            prm = val;
            break
        end
    end
    fclose(fid);
end


%% ------------------------------------------------------------------
% splits a scan line into a tag/value pair
function [tag,val] = scanline (line)

tag = '';
val = '';
idx = find(line=='=');
if length(idx) == 1
    tag = sscanf(line(1:idx-1),'%s');
    for i=idx+1:length(line)
        if line(i) ~= ' '
            break;
        end
    end
    val = line(i:length(line));
end


%% ------------------------------------------------------------------
% defines initial parameter settings from string
function prm=setinitprm(str)
[method val] = strread(str,'%s %s');
method = upper(cell2mat(method));
prm.reset = method;
switch method
    case 'HOMOG'
        prm.val = str2num(cell2mat(val));
    case 'NIM'
        prm.nim = cell2mat(val);
        prm.nimlast = false;
    case 'NIM_LAST'
        prm.nim = cell2mat(val);
        prm.nimlast = true;
end    
