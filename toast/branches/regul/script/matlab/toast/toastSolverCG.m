function toastSolverCG (prm, lprm, ref, x0)
%toastSolverCG        - Nonlinear conjugate gradient iterative solver
%
%   Syntax: toastSolverCG(prm, lprm, ref, x0)
%     prm:   parameter structure
%     lprm:  local function parameter structure
%     ref:   global refractive index
%     x0:    initial parameter vector
%
%   This function should normally not be called directly. Instead, it is
%   called by toastRecon when the solver.method field of the prm
%   structure is set to 'PCG'.
%
%   See also: toastRecon

global ITR_TERM;        % termination flag

ITR_TERM = false;
try_again = false;
c0 = 0.3;               % speed of light
p = length(x0);         % parameter space dimension
ps = p/2;               % single-parameter dimension

x = x0;                 % initial parameters
logx = log(x);          % parameter transformation

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

% Conjugate gradient loop
while ((itr <= prm.solver.itmax) && ...
       (err > prm.solver.tol*err0) && ...
       (errp-err > prm.solver.tol) && ...
       (ITR_TERM == false)) || try_again == true

    errp = err;

    % Gradient of cost function
    disp('Calculating gradient ...')
    r = gradient(x);
      
    if itr > 1
        delta_old = delta_new;
        delta_mid = r' * s;
    end
      
    % Apply PCG preconditioner
    s = r; % dummy for now
      
    if itr == 1
        dx = s;
        delta_new = r' * dx;
        delta0 = delta_new;
    else
        delta_new = r' * s;
        beta = (delta_new-delta_mid) / delta_old;
        if mod(itr, prm.solver.cg.reset) == 0 || try_again == true || beta <= 0
            disp ('Resetting CG ...');
            dx = s;
        else
            dx = s + dx*beta;
        end
    end
    delta_d = dx' * dx;

    % Line search
    if prm.solver.lsearch
        disp ('Entering line search ...')
        [step, err] = toastLineSearch (logx, dx, step, err, @objective);
    else
        [step, err] = toastStepSize (logx, dx, step*1.5, err, @objective);
    end
    
    % Add update to solution
    disp('Updating solution ...')
    logx = logx + dx*step;
    x = exp(logx);
    
        save(['iter',num2str(itr)],'x')
        
    proj = privProject (lprm.hMesh, lprm.hBasis, logx, ref, prm.data.freq, ...
        lprm.qvec, lprm.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);
    err = privObjective (proj, lprm.data, lprm.sd, lprm.hReg, logx);

    if err < errp % improvement of cost function
        % Return solution
        if isfield(lprm,'callback') && isfield(lprm.callback,'iter')
            feval(lprm.callback.iter, lprm.callback.context, itr, x, err);
        end
        itr = itr+1;
        try_again = false;
    else
        if try_again == false
            err = errp;
            try_again = true; % try again in gradient direction
        else
            try_again = false;
        end
    end
    
end % GN loop


    % =====================================================================
    % Gradient calculation: adjust input parameters
    function g = gradient(x)

    mua = toastMapSolToMesh (lprm.hBasis, x(1:ps)) .* (ref/c0);
    kap = toastMapSolToMesh (lprm.hBasis, x(ps+1:p)) .* (ref/c0);
    mus = 1./(3*kap) - mua;

    g = -toastGradient (lprm.hMesh, lprm.hBasis, lprm.qvec, lprm.mvec, mua, mus, ref, ...
        prm.data.freq, lprm.data, lprm.sd, prm.fwdsolver.method, prm.fwdsolver.tol);
    g = g .* x; % parameter scaling
    g = g - toastRegulGradient (lprm.hReg, logx);
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


end % of toastSolverCG
