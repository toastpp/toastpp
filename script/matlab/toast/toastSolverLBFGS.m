function toastSolverLBFGS (prm, lprm, ref, x0)
%toastSolverLBFGS     - Limited-bandwidth BFGS iterative solver
%
%   Syntax: toastSolverLBFGS(prm, lprm, ref, x0)
%     prm:   parameter structure
%     lprm:  local function parameter structure
%     ref:   global refractive index
%     x0:    initial parameter vector
%
%   This function should normally not be called directly. Instead, it is
%   called by toastRecon when the solver.method field of the prm
%   structure is set to 'LBFGS'.
%
%   See also: toastRecon

global ITR_TERM;        % termination flag

ITR_TERM = false;
c0 = 0.3;               % speed of light
history = 10;           % number of previous steps to store
nx = length(x0);
np = nx/2;

x = x0;                 % initial parameters
logx = log(x);          % parameter transformation
logx0 = logx;

% Initialise starting parameters
proj = privProject (lprm.hMesh, lprm.hBasis, logx, ref, prm.data.freq, ...
    lprm.qvec, lprm.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);
err0 = privObjective (proj, lprm.data, lprm.sd, lprm.hReg, logx); %initial error
err  = err0;
errp = inf;
step = prm.solver.step0;
itr  = 1;

% Return initial parameters
if isfield(lprm,'callback') && isfield(lprm.callback,'iter')
    feval(lprm.callback.iter, lprm.callback.context, 0, x, err);
end

% initial gradient
g = gradient(x);

% initial inverse Hessian
H0 = ones(nx,1); % use identity matrix

mm = history;
m = 0;
upd_idx=-1;
gamma = 1.0;
p = zeros(nx,1);
S = zeros(mm,nx);
Y = zeros(mm,nx);


% LBFGS loop
while  (itr <= prm.solver.itmax) && ...
       (err > prm.solver.tol*err0) && ...
       (errp-err > prm.solver.tol) && ...
       (ITR_TERM == false)

    errp = err;

    % apply approximate inverse Hessian
    if itr == 1
        p = -(H0.*g);
    else
        p = -(H0.*g)*gamma;
        f1 = zeros(2*m,1);
        f2 = zeros(2*m,1);
        for i=1:m
            idx_i = mod (upd_idx+i,m);
            f1(i) = dot(S(idx_i+1,:),g);
            f1(i+m) = dot(Y(idx_i+1,:).*H0', g)*gamma;
        end
        f2 = D * f1;
        for i=1:m
            idx_i = mod(upd_idx+i,m);
            for j=1:nx
                p(j) = p(j) - S(idx_i+1,j) * f2(i) + Y(idx_i+1,j)*H0(j)* ...
                       f2(i+m)*gamma;
            end
        end
    end
    
    [alpha,fmin] = toastLineSearch(logx,p,step,err,@objective);
    if fmin < err
        step = alpha;
        err = fmin;
    else
        % no improvement - reset to gradient direction
        p = -H0 .* g;
        [alpha,fmin] = toastLineSearch(logx,p,step,err,@objective);
        if fmin < err
            step = alpha;
            err = fmin;
        else
            fprintf (1,'No improvement found during line search. Terminating.\n');
            return;
        end
    end
    
    % update approximate solution
    logx = logx + p*step;
    x = exp(logx);
    
    % Return solution
    if isfield(lprm,'callback') && isfield(lprm.callback,'iter')
        feval(lprm.callback.iter, lprm.callback.context, itr, x, err);
    end

    % update gradient
    g1 = gradient(x);

    proj = privProject (lprm.hMesh, lprm.hBasis, logx, ref, prm.data.freq, ...
        lprm.qvec, lprm.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);
    err = privObjective (proj, lprm.data, lprm.sd, lprm.hReg, logx);

    % update S and Y
    upd_idx = mod(upd_idx+1,mm);
    S(upd_idx+1,:) = logx-logx0;
    Y(upd_idx+1,:) = g1-g;
    
    % update hessian scaling
    gamma = dot(Y(upd_idx+1),S(upd_idx+1)) / dot(Y(upd_idx+1),Y(upd_idx+1));
    if gamma < 0
        gamma = 1.0;  % problems here
    end
    fprintf (1,'BFGS scaling set to %f\n', gamma);
    
    % grow R and D
    if m < mm
        m = m+1;
        R = zeros(m);
        D = zeros(2*m);
    end
    
    % update R
    for i=1:m
        idx_i = mod(upd_idx+i,m);
        for j=1:i-1
            R(i,j) = 0;
        end
        for j=i:m
            idx_j = mod(upd_idx+j,m);
            R(i,j) = dot (S(idx_i+1,:), Y(idx_j+1,:));
        end
    end
    
    % update D
    RI = inv(R);
    RIT = RI';
    YTY = zeros(m);
    for i=1:m
        idx_i = mod(upd_idx+i,m);
        YH = Y(idx_i+1,:) .* H0';
        for j=1:i
            idx_j = mod(upd_idx+j,m);
            YTY(i,j) = dot (YH, Y(idx_j+1,:));
            YTY(j,i) = YTY(i,j);
        end
        YTY(i,i) = dot(YH, Y(idx_i+1,:));
    end
    YTY = YTY * gamma;
    for i=1:m
        idx_i = mod(upd_idx+i,m);
        YTY(i,i) = YTY(i,i) + dot(S(idx_i+1,:),Y(idx_i+1,:));
    end
    B = YTY*RI;
    C = RIT*B;
    
    for i=1:m
        for j=1:m
            D(i,j) = C(i,j);
            D(i,j+m) = -RIT(i,j);
            D(i+m,j) = -RI(i,j);
            D(i+m,j+m) = 0;
        end
    end
    
    logx0 = logx;
    g = g1;
    
    fprintf (1, '**** LBFGS ITERATION %d, ERROR %f\n\n', itr, err);
    itr = itr+1;
end


    % =====================================================================
    % Gradient calculation: adjust input parameters
    function g = gradient(x)

    mua = toastMapSolToMesh (lprm.hBasis, x(1:np)) .* (ref/c0);
    kap = toastMapSolToMesh (lprm.hBasis, x(np+1:nx)) .* (ref/c0);
    mus = 1./(3*kap) - mua;

    g = toastGradient (lprm.hMesh, lprm.hBasis, lprm.qvec, lprm.mvec, mua, mus, ref, ...
        prm.data.freq, lprm.data, lprm.sd, prm.fwdsolver.method, prm.fwdsolver.tol);
    g = g .* x; % parameter scaling
    g = g + toastRegulGradient (lprm.hReg, logx);
    end


    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)

    proj = privProject (lprm.hMesh, lprm.hBasis, x, ref, ...
        prm.data.freq, lprm.qvec, lprm.mvec, prm.fwdsolver.method, ...
        prm.fwdsolver.tol);
    [p, p_data, p_prior] = privObjective (proj, lprm.data, lprm.sd, ...
        lprm.hReg, x);
    fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
    end


end
