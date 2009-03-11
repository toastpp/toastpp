function cwdata = toastFwdCW(prm)
%toastFwdCW           - High-level toast diffusion forward model.
%
% Synopsis: cwdata = toastFwd(prm)
%    prm:   model parameter structure
%    cwdata: boundary data (log intensity)
%
% Calculates continuous wave (CW) boundary measurement data,
% given a mesh, optical coefficient distributions, and source and
% detector locations.
%
% Runs a toast forward solver with the parameters defined in prm.
% prm contains information about measurements, meshes,
% tolerance limits for the forward solver, etc.
%
% The calculated boundary data are returned in cwdata, in the form
% of log intensity. The data file consists of a 1-D real vector,
% containing the measurements for all source and detector
% combinations defined in the QM link list. The data are arranged
% so that the detector index changes most rapidly.
%
% See also: toastFwd, toastReadParam
  
disp('---------------------------------------')
disp('Starting forward data calculation (CW)')
disp('---------------------------------------')

% Initialisations
toastCatchErrors();
dispprm(prm,'');

% Read a TOAST mesh definition from file.
if isfield(prm,'basis') && isfield(prm.basis,'hMesh') 
    hMesh = prm.basis.hMesh;
else
    hMesh = toastReadMesh (prm.basis.meshfile);
    toastReadQM (hMesh, prm.meas.qmfile);
end
n = toastMeshNodeCount (hMesh);
dmask = toastDataLinkList (hMesh);
nqm = length(dmask);

% Set up homogeneous initial parameter estimates
mua = resetprm(prm.initprm.mua,hMesh);
mus = resetprm(prm.initprm.mus,hMesh);
ref = resetprm(prm.initprm.ref,hMesh);
kap = 1./(3*(mua+mus));

% Generate source vectors
qvec = real (toastQvec (hMesh, prm.meas.src.type, prm.meas.src.prof, prm.meas.src.width));

% Generate measurement vectors
mvec = real (toastMvec (hMesh, prm.meas.det.prof, prm.meas.det.width));

% Apply forward model: generate projection data set f[x0]
smat = real (toastSysmat (hMesh, mua, mus, ref, 0));   % system matrix

% Solve linear system: source and measurement operators
switch upper(prm.linsolver.method)
    case 'DIRECT'
        cwdata = reshape (log(mvec.' * (smat\qvec)), [], 1);
    otherwise
        nq = size(qvec,2);
        nm = size(mvec,2);
        for i=1:nq
            phi = pcg(smat,qvec(:,i),prm.linsolver.tol,1000);
            prj = log(mvec.' * phi);
            cwdata((i-1)*nm+1:i*nm,1) = prj;
        end
end
cwdata = full(cwdata(dmask)); % remove unused source-detector combinations

% Write to files
if isfield(prm,'data') && isfield(prm.data,'lnampfile')
    toastWriteRealVector(prm.data.lnampfile,cwdata);
    disp (['Log intensity data written to ' prm.data.lnampfile]);
end

end

% =============================================================
% recursively display prm structure fields
function dispprm(prm,prefix)
fields = fieldnames(prm);
for i=1:size(fields,1)
    f = getfield(prm,fields{i});
    if strcmpi(fields{i},'callback'), f = '<internal structure>'; end
    if isstruct(f)
        dispprm(f,[prefix fields{i} '.']);
    else
        fprintf (1, '%30s :', [prefix fields{i}]);
        if     ischar(f),    fprintf (1, ' %s', f);
        elseif islogical(f)
            if f == true, fprintf (1, ' true');
            else, fprintf (1, ' false'); end
        elseif isreal(f)
            if length(f) > 3
                fprintf (1, ' [%dx%d real]', size(f,1), size(f,2));
            else
                fprintf (1, ' %f', f);
            end
        elseif isinteger(f), fprintf (1, ' %d', f);
        end
        fprintf (1, '\n');
    end
end
end

%% =============================================================
function prm=resetprm(cfg,hmesh)
n = toastMeshNodeCount (hmesh);
if isnumeric(cfg) && length(cfg) == n
    prm = cfg;  % parameter vector stored directly
else
    switch upper(cfg.reset)
        case 'HOMOG'
            prm = ones(n,1) * cfg.val;
        case 'NIM'
            prm = toastReadNIM(cfg.nim);
            if length(prm) ~= n
                disp('Warning: incompatible size of NIM file')
            end
        case 'NIM_LAST'
            prm = toastReadNIM(cfg.nim,0);
            if length(prm) ~= n
                disp('Warning: incompatible size of NIM file')
            end
        otherwise
            disp('Warning: Unsupported reset method')
    end
end
end
