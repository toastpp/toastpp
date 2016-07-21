function toastRecon(prm)
%toastRecon           - High-level toast reconstruction function.
%
% Synopsis: toastRecon(prm)
%    prm: reconstruction parameter structure
%
% A high-level convenience function which runs a toast reconstruction with
% the parameters defined in prm.
% prm is a toastParam instance containing information about measurements,
% meshes and grids, tolerance limits for the forward and inverse solvers,
% regularisation parameters, etc.
%
% See also: toastParam
  
disp('---------------------------------------')
disp('Starting reconstruction')
disp('---------------------------------------')

% ----------------------------------------------------------------------
% Initialisations
global RES;  % a structure for collecting reconstruction results
global LPRM; % a set of local parameters
global scref;
RES.of = [];

% ----------------------------------------------------------------------
% Check consistency of parameter structure
[prm LPRM] = checkprm(prm);
dispprm(prm,'');

% ----------------------------------------------------------------------
% Set up some variables
c0 = 0.3;           % speed of light in vacuum [mm/ps]

% ----------------------------------------------------------------------
% Read a TOAST mesh definition from file.
if isfield(prm.fwdsolver,'hmesh')
    LPRM.hMesh = prm.fwdsolver.hmesh;
else
    LPRM.hMesh = toastMesh (prm.fwdsolver.meshfile);
    LPRM.hMesh.ReadQM (prm.meas.qmfile);
end
n = LPRM.hMesh.NodeCount ();
dmask = LPRM.hMesh.DataLinkList ();

% ----------------------------------------------------------------------
% Set up the mapper between FEM and solution bases
if ~isfield(prm.solver.basis,'hbasis')
    prm.solver.basis.hbasis = toastBasis(LPRM.hMesh,prm.solver.basis.bdim,'Linear');
end
LPRM.hBasis = prm.solver.basis.hbasis;
RES.bdim = LPRM.hBasis.Dims();
blen = prod(RES.bdim);

% ----------------------------------------------------------------------
% Set up homogeneous initial parameter estimates
fprintf('Initial parameter estimates:\n');
mua = resetprm(prm.initprm.mua,LPRM.hMesh);
fprintf('mua: mean=%f, std=%f\n', mean(mua), std(mua));
mus = resetprm(prm.initprm.mus,LPRM.hMesh);
fprintf('mus: mean=%f, std=%f\n', mean(mus), std(mus));
ref = resetprm(prm.initprm.ref,LPRM.hMesh);
fprintf('ref: mean=%f, std=%f\n', mean(ref), std(ref));
kap = 1./(3*(mua+mus));

% ----------------------------------------------------------------------
% Read the data
if isfield(prm.data,'lnamp')
    mdata = prm.data.lnamp;
else
    fprintf('Reading log amplitude data from %s\n', prm.data.lnampfile);
    mdata = toastReadVector(prm.data.lnampfile);
end
if isfield(prm.data,'phase')
    pdata = prm.data.phase;
else
    fprintf('Reading phase data from %s\n', prm.data.phasefile);
    pdata = toastReadVector(prm.data.phasefile);
end
LPRM.data = [mdata;pdata];
m = length(LPRM.data);
fprintf('Data space dimension: %d\n', m);

% ----------------------------------------------------------------------
% Generate source vectors
LPRM.qvec = LPRM.hMesh.Qvec (prm.meas.src.type, prm.meas.src.prof, prm.meas.src.width);
fprintf ('Source vector set up: %d sources\n', size(LPRM.qvec,2));

% ----------------------------------------------------------------------
% Generate measurement vectors
LPRM.mvec = LPRM.hMesh.Mvec (prm.meas.det.prof, prm.meas.det.width, ref);
fprintf ('Detector vector setup: %d detectors\n', size(LPRM.mvec,2));

% ----------------------------------------------------------------------
% Initial data set f[x0]
fprintf('Calculating projections from initial parameters ...\n');
proj = toastProject (LPRM.hMesh, mua, mus, ref, prm.data.freq, ...
    LPRM.qvec, LPRM.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);

% ----------------------------------------------------------------------
% Difference data setup
if isfield(prm.data,'useref') && prm.data.useref == true
    fprintf ('Using difference data\n');
    if isfield(prm.data.ref,'lnampfile')
        fprintf ('Log amplitude reference: %s\n', prm.data.ref.lnampfile);
        mref = toastReadVector(prm.data.ref.lnampfile);
        if length(mref) > 0
            mproj = proj(1:length(mdata));
            mdata = mdata-mref+mproj;
        else
            fprintf ('Log amplitude reference: file not found!\n');
            return;
        end
    end
    if isfield(prm.data.ref,'phasefile')
        fprintf ('Phase reference: %s\n', prm.data.ref.phasefile);
        pref = toastReadVector(prm.data.ref.phasefile);
        if length(pref) > 0
            pproj = proj(length(mdata)+1:end);
            pdata = pdata-pref+pproj;
        else
            fprintf ('Phase reference: file not found!\n');
            return;
        end
    end
    LPRM.data = [mdata;pdata];
end

% ----------------------------------------------------------------------
% data scaling
nqm = length(proj)/2;
sd_lnmod = ones(nqm,1);
sd_phase = ones(nqm,1);
switch prm.solver.dscale
    case 'AVG_DIFFDATA'
        sd_lnmod = sd_lnmod * norm(mdata-proj(1:nqm)); 
        sd_phase = sd_phase * norm(pdata-proj(nqm+1:nqm*2));
    case 'AVG_DATA'
        sd_lnmod = sd_lnmod * norm(proj(1:nqm));
        sd_phase = sd_phase * norm(proj(nqm+1:nqm*2));
end
LPRM.sd = [sd_lnmod;sd_phase];

% ----------------------------------------------------------------------
% initial parameter estimates in solution basis
RES.bmua = LPRM.hBasis.Map ('M->B', mua);
RES.bmus = LPRM.hBasis.Map ('M->B', mus);
bmua_itr(1,:) = RES.bmua;
bmus_itr(1,:) = RES.bmus;

bcmua = LPRM.hBasis.Map ('M->B', mua .* (c0./ref));
bckap = LPRM.hBasis.Map ('M->B', kap .* (c0./ref));

scmua = LPRM.hBasis.Map ('B->S', bcmua);
sckap = LPRM.hBasis.Map ('B->S', bckap);
scref = LPRM.hBasis.Map ('M->S', c0./ref);

x = [scmua;sckap];
logx = log(x);
p = length(x);

% ----------------------------------------------------------------------
% Initialise regularisation
%LPRM.hReg = 0;
%if isfield(prm.regul,'method') && ~strcmpi(prm.regul.method,'none')
%    LPRM.hReg = toastRegul (prm.regul, logx);
%end

% ----------------------------------------------------------------------
% Inverse solver loop
switch prm.solver.method
    case 'LM'
        solver = toastSolverGN (prm);
        solver.Solve (x);
    case 'GN_implicit'
        toastSolverGNimplicit (prm, LPRM, ref, x);
    case 'PCG'
        solver = toastSolverCG (prm);
        solver.Solve (x);
    case 'LBFGS'
        solver = toastSolverLBFGS (prm);
        solver.Solve (x);
end
        

end % of toastRecon


% =============================================================
% =============================================================
function [prm,localprm] = checkprm(prm)

if ~isa(prm,'toastParam')
    error('Expected a toastParam argument.');
end

% check for mandatory fields
if isfield(prm.solver,'basis') == false
    error('Required parameter field "solver.basis" is missing.');
end

% fill missing parameter entries
if isfield(prm.solver,'method') == false || ischar(prm.solver.method) == false || length(prm.solver.method) == 0
    prm.solver.method = 'PCG';
end
if isfield(prm.solver,'tol') == false || isnumeric(prm.solver.tol) == false
    prm.solver.tol = 1e-10;
end
if isfield(prm.solver,'itmax') == false || isnumeric(prm.solver.itmax) == false
    prm.solver.itmax = 100;
end
if isfield(prm.solver,'lsearch') == false || islogical(prm.solver.itmax) == false
    prm.solver.lsearch = true;
end
if isfield(prm.solver,'step0') == false || isnumeric(prm.solver.step0) == false
    prm.solver.step0 = 1;
end
if isfield(prm.solver,'dscale') == false || ischar(prm.solver.dscale) == false || length(prm.solver.dscale) == 0
    prm.solver.dscale = 'AVG_DIFFDATA';
end
switch upper(prm.solver.method)
    case 'PCG'
        if isfield(prm.solver,'cg') == false || isfield(prm.solver.cg,'reset') == false || isnumeric(prm.solver.cg.reset) == false
            prm.solver.cg.reset = 10;
        end
    case {'LM' 'GN_IMPLICIT'}
        if isfield(prm.solver,'krylov') == false || isfield(prm.solver.krylov,'method') == false || ischar(prm.solver.krylov.method) == false || length(prm.solver.krylov.method) == 0
            prm.solver.krylov.method = 'gmres';
        end
        if isfield(prm.solver.krylov,'tol') == false || isnumeric(prm.solver.krylov.tol) == false
            prm.solver.krylov.tol = 1e-2;
        end
        if isfield(prm.solver.krylov,'maxit') == false || isnumeric(prm.solver.krylov.maxit) == false
            prm.solver.krylov.maxit = 100;
        end
end

localprm.Himplicit = true;                  % Implicit/explicit Hessian matrix

% Parameters for display callback function
localprm.callback.iter = @solve_iter;
localprm.callback.context = prm;

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


% =============================================================
function prm=resetprm(cfg,hmesh)
n = hmesh.NodeCount();
switch upper(cfg.reset)
    case 'HOMOG'
        prm = ones(n,1) * cfg.val;
    case 'NIM'
        prm = toastNim(cfg.nim);
        if length(prm) ~= n
            disp('Warning: incompatible size of NIM file')
        end
    case 'NIM_LAST'
        prm = toastNim(cfg.nim,0);
        if length(prm) ~= n
            disp('Warning: incompatible size of NIM file')
        end
    otherwise
        disp('Warning: Unsupported reset method')
end
end


% =============================================================
% called by nonlinear solvers after each iteration

function solve_iter (prm, itr, x, err)

fprintf (1, '**** Iteration %d, objective %f\n', itr, err)

global RES LPRM
global scref;

ps = length(x)/2;
scmua = x(1:ps);
sckap = x(ps+1:2*ps);
smua  = scmua./scref;
skap  = sckap./scref;
smus  = 1./(3*skap) - smua;
RES.bmua = LPRM.hBasis.Map ('S->B', smua);
RES.bmus = LPRM.hBasis.Map ('S->B', smus);
RES.of(itr+1) = err;

if isfield(prm,'callback') && isfield(prm.callback,'iter')
    % let the calling function do the display
    if isfield(prm.callback,'request') % check requests for secondary results
        if isfield(prm.callback.request,'prior') && prm.callback.request.prior == true
            RES.kapref = LPRM.hReg.Kappa (log(x));
        end
    end
    feval(prm.callback.iter, prm, RES);
else
    % inline display
    disp_iter(RES);
end

end


% =============================================================
% inline function for display of reconstruction results (if
% prm.callback is not defined)

function disp_iter(res)

figure(1);
set(gcf,'Name','toastRecon');
if length(res.bdim) == 2
    subplot(2,2,1), imagesc(reshape(res.bmua,res.bdim(1),res.bdim(2))), axis equal, axis tight, colorbar;
    set(gca,'Units','normalized','OuterPosition',[0 0.5 0.5 0.5]);
    title('initial \mu_a');
    subplot(2,2,2), imagesc(reshape(res.bmus,res.bdim(1),res.bdim(2))), axis equal, axis tight, colorbar;
    set(gca,'Units','normalized','OuterPosition',[0.5 0.5 0.5 0.5]);
    title('initial \mu_s');
    subplot(2,1,2);
    semilogy(res.of);xlabel('iteration');axis tight;title(['Objective:']);
    set(gca,'Units','normalized','OuterPosition',[0 0 1 0.5]);
    drawnow;
else
    muamin=min(res.bmua);
    muamax=max(res.bmua);
    musmin=min(res.bmus);
    musmax=max(res.bmus);
    bmua = reshape(res.bmua,res.bdim');
    bmus = reshape(res.bmus,res.bdim');
    for i=1:4
        idx = round(res.bdim(3)/4*(i-0.5));
        subplot(4,2,i); imagesc(bmua(:,:,idx),[muamin muamax]); axis equal tight
        subplot(4,2,i+4); imagesc(bmus(:,:,idx),[musmin musmax]); axis equal tight
    end
    drawnow
    save recon_gui_res.mat res
end

end
