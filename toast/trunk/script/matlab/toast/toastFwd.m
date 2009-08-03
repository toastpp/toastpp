function [mdata pdata] = toastFwd(prm)
%toastFwd             - High-level toast diffusion forward model.
%
% Synopsis: [mdata pdata] = toastFwd(prm)
%    prm:   model parameter structure
%    mdata: boundary data (log amplitude)
%    pdata: boundary data (phase)
%
% Calculates boundary measurement data, given a mesh, optical
% coefficient distributions, and source and detector locations.
%
% Runs a toast forward solver with the parameters defined in prm.
% prm contains information about measurements, meshes,
% tolerance limits for the forward solver, etc.
%
% The calculated boundary data are returned in mdata and pdata,
% as log amplitude and phase. The data files are 1-D real vectors,
% containing the measurements for all source and detector
% combinations defined in the QM link list. The data are arranged
% so that the detector index changes most rapidly.
%
% See also: toastFwdCW, toastReadParam
  
disp('---------------------------------------')
disp('Starting forward data calculation')
disp('---------------------------------------')

% Initialisations
toastCatchErrors();

refind = 1.4;                           % refractive index

% ======================================================================
% Check consistency of parameter structure

[prm lprm] = checkprm(prm);
dispprm(prm,'');

% Set up some variables
c0 = 0.3;
cm = c0/refind;

% Read a TOAST mesh definition from file.
if isfield(prm,'fwdsolver') && isfield(prm.fwdsolver,'hmesh')
    hMesh = prm.fwdsolver.hmesh;
else
    hMesh = toastReadMesh (prm.fwdsolver.meshfile);
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
qvec = toastQvec (hMesh, prm.meas.src.type, prm.meas.src.prof, prm.meas.src.width);
nq = size(qvec,2);

% Generate measurement vectors
mvec = toastMvec (hMesh, prm.meas.det.prof, prm.meas.det.width);
nm = size(mvec,2);

% Apply forward model: generate projection data set f[x0]
proj = toastProject (hMesh, mua, mus, ref, prm.data.freq, ...
    qvec, mvec, prm.fwdsolver.method, prm.fwdsolver.tol);

% Split into log amplitude and phase components
mdata = proj(1:nqm);
pdata = proj(nqm+1:nqm*2);

% Write to files
if isfield(prm,'data')
    if isfield(prm.data,'lnampfile')
        toastWriteRealVector(prm.data.lnampfile,mdata);
        disp (['Log amplitude data written to ' prm.data.lnampfile]);
    end
    if isfield(prm.data,'phasefile')
        toastWriteRealVector(prm.data.phasefile,pdata);
        disp (['Phase data written to ' prm.data.phasefile]);
    end
end

figure;

% Display data as a function of source-detector separation
qp = toastQPos(hMesh);
mp = toastMPos(hMesh);
for i=1:size(qp,1)
    for j=1:size(mp,1)
        dst(j,i) = norm(qp(i,:)-mp(j,:));
    end
end
dst = reshape(dst,[],1);
dst = dst(dmask);

subplot(2,2,1);
plot(dst,mdata,'o'); axis tight; title('log amplitude data');
xlabel('optode separation');
subplot(2,2,2);
plot(dst,pdata,'o'); axis tight; title('phase data');
xlabel('optode separation');

% Display data as sinograms
msino = zeros(nq*nm,1);
msino(dmask) = mdata;
msino = reshape(msino,nq,nm);
subplot(2,2,3); imagesc(msino); axis equal tight; colorbar;
title('log amplitude data');
xlabel('detectors');
ylabel('sources');

psino = zeros(nq*nm,1);
psino(dmask) = pdata;
psino = reshape(psino,nq,nm);
subplot(2,2,4); imagesc(psino); axis equal tight; colorbar;
title('phase data');
xlabel('detectors');
ylabel('sources');

end

%% =============================================================
function [prm,localprm] = checkprm(prm)

% fill missing parameter entries
localprm = 0;

end

%% =============================================================
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

%% =============================================================
% inline function for display of reconstruction results (if
% prm.callback is not defined)
function disp_iter(res)

figure(1);
set(gcf,'Name','toastRecon');
subplot(2,2,1), imagesc(reshape(res.bmua,res.bdim(1),res.bdim(2))), axis equal, axis tight, colorbar;
set(gca,'Units','normalized','OuterPosition',[0 0.5 0.5 0.5]);
title('initial \mu_a');
subplot(2,2,2), imagesc(reshape(res.bmus,res.bdim(1),res.bdim(2))), axis equal, axis tight, colorbar;
set(gca,'Units','normalized','OuterPosition',[0.5 0.5 0.5 0.5]);
title('initial \mu_s');
subplot(2,1,2);
semilogy(res.of);xlabel('iteration');axis tight;title(['Objective:']);
set(gca,'Units','normalized','OuterPosition',[0 0 1 0.5]);
drawnow

drawnow;
end
