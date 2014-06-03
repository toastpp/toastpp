function toastSolverGN (prm, lprm, ref, x0)
%toastSolverGN        - Gauss-Newton iterative solver
%
%   Syntax: toastSolverGN(prm, lprm, ref, x0)
%     prm:   parameter structure
%     lprm:  local function parameter structure
%     ref:   global refractive index
%     x0:    initial parameter vector
%
%   This function should normally not be called directly. Instead, it is
%   called by toastRecon when the solver.method field of the prm
%   structure is set to 'LM'.
%
%   See also: toastRecon
  
global ITR_TERM;        % termination flag

ITR_TERM = false;
c0 = 0.3;               % speed of light
m = length(lprm.data);  % data space dimension
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

% Gauss-Newton loop
while (itr <= prm.solver.itmax) && ...
      (err > prm.solver.tol*err0) && ...
      (errp-err > prm.solver.tol) && ...
      (ITR_TERM == false)

    errp = err;

    % Construct the Jacobian
    disp ('Calculating Jacobian ...')
    J = jacobian(x);
    
    % data normalisation
    for i = 1:m, J(i,:) = J(i,:) / lprm.sd(i); end;

    % parameter normalisation
    for i = 1:p, J(:,i) = J(:,i) * x(i); end;

    % Normalisation of Hessian
    psiHdiag = toastRegulHDiag (lprm.hReg, logx);
    for i = 1:p
        M(i) = sum(J(:,i) .* J(:,i));
        M(i) = M(i) + psiHdiag(i);
        M(i) = 1 ./ sqrt(M(i));
    end
    for i = 1:p, J(:,i) = J(:,i) * M(i); end;

    % Gradient of cost function
    r = J' * ((lprm.data-proj)./lprm.sd);
    r = r - toastRegulGradient (lprm.hReg, logx) .* M';
    
    if lprm.Himplicit == true
        % Update with implicit Krylov solver
        fprintf (1, 'Entering Krylov solver (tol=%g)\n', ...
             prm.solver.krylov.tol);
        if lprm.hReg
            HessReg = toastRegulHess (lprm.hReg, x);
        end
        dx = krylov(r);
         
    else
        % Update with explicit Hessian
        H = J' * J;
        lambda = 0.01;
        H = H + eye(size(H)).* lambda;
        dx = gmres (H, r, 30, prm.solver.krylov.tol, 100);
        clear H;
    end

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
    
    proj = privProject (lprm.hMesh, lprm.hBasis, logx, ref, prm.data.freq, ...
        lprm.qvec, lprm.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);
    err = privObjective (proj, lprm.data, lprm.sd, lprm.hReg, logx);

    % Return solution
    if isfield(lprm,'callback') && isfield(lprm.callback,'iter')
        feval(lprm.callback.iter, lprm.callback.context, itr, x, err);
    end

    itr = itr+1;
end % GN loop


    % =====================================================================
    % Jacobian calculation: adjust input parameters
    function J = jacobian(x)
        mua = toastMapSolToMesh (lprm.hBasis, x(1:ps)) .* (ref/c0);
        kap = toastMapSolToMesh (lprm.hBasis, x(ps+1:p)) .* (ref/c0);
        mus = 1./(3*kap) - mua;

        % Construct the Jacobian
        J = toastJacobian (lprm.hMesh, lprm.hBasis, lprm.qvec, lprm.mvec, ...
            mua, mus, ref, ...
            prm.data.freq, prm.fwdsolver.method, prm.fwdsolver.tol);
    end


    % =====================================================================
    % Krylov solver subroutine
    function dx = krylov(r)
        k_t = cputime;
        switch prm.solver.krylov.method
            case 'bicgstab'
                [dx k_flag k_res k_iter] = bicgstab(@jtjx, r, prm.solver.krylov.tol, prm.solver.krylov.maxit);
            otherwise
                [dx k_flag k_res k_iter] = gmres (@jtjx, r, 30, prm.solver.krylov.tol, prm.solver.krylov.maxit);
        end
        k_dt = cputime-k_t;
        fprintf (1, '--> iter=%0.0f(%0.0f), time=%0.1fs, res=%g\n', ...
            k_iter(1), k_iter(2), k_dt, k_res);
        clear k_t k_dt k_flag k_res k_iter
    end
    
    
    % =====================================================================
    % Callback function for matrix-vector product (called by toastKrylov)
    function b = jtjx(x)
        b = J' * (J*x);
        if lprm.hReg
            b = b + M' .* (HessReg * (M' .* x));
        end
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


end % of toastSolverGN
