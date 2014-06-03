function toastSolverLBFGS_old (prm, lprm, ref, x0)
% Limited-bandwidth BFGS iterative solver
%
% Syntax: toastSolverLBFGS(prm, lprm, ref, x0)
%
% Parameters:
%         prm [structure]:
%             parameter structure
%         lprm [structure]:
%             local function parameter structure
%         ref [real]:
%             global refractive index
%         x0 [real array]:
%             initial parameter vector
%
% Notes:  This function should normally not be called directly. Instead,
%         it is called by toastRecon when the solver.method field of the
%         prm structure is set to 'LBFGS'.
%
% See also:
%         toastRecon

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
proj = project (logx);
err0 = objective (proj, logx); %initial error
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
    
    [alpha,fmin] = toastLineSearch(logx,p,step,err,@clbk_objective);
    if fmin < err
        step = alpha;
        err = fmin;
    else
        % no improvement - reset to gradient direction
        p = -H0 .* g;
        [alpha,fmin] = toastLineSearch(logx,p,step,err,@clbk_objective);
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

    proj = project (logx);
    err = objective (proj, logx);

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
    % Generate boundary data from a solution vector
    function prj = project (logx)
        x = exp(logx);
        cm = 0.3/ref(1);
        scmua = x(1:size(x)/2);
        sckap = x(size(x)/2+1:end);
        smua = scmua/cm;
        skap = sckap/cm;
        smus = 1./(3*skap) - smua;
        mua = lprm.hBasis.Map ('S->M', smua);
        mus = lprm.hBasis.Map ('S->M', smus);
        prj = toastProject (lprm.hMesh, mua, mus, ref, prm.data.freq, ...
            lprm.qvec, lprm.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);
    end

    % =====================================================================
    % Gradient calculation: adjust input parameters
    function g = gradient(x)
        mua = lprm.hBasis.Map ('S->M', x(1:np)) .* (ref/c0);
        kap = lprm.hBasis.Map ('S->M', x(np+1:nx)) .* (ref/c0);
        mus = 1./(3*kap) - mua;

        g = toastGradient (lprm.hMesh, lprm.hBasis, lprm.qvec, lprm.mvec, mua, mus, ref, ...
            prm.data.freq, lprm.data, lprm.sd, prm.fwdsolver.method, prm.fwdsolver.tol);
        g = g .* x; % parameter scaling
        if lprm.hReg ~= 0
            g = g + lprm.hReg.Gradient (logx);
        end
    end

    % =====================================================================
    % Objective function for a set of projections
    function [ob,varargout] = objective (prj,logx)
        ob_data = sum(((lprm.data-prj)./lprm.sd).^2);
        ob_prior = 0;
        
        if lprm.hReg ~= 0
            ob_prior = lprm.hReg.Value(logx);
        end
        
        ob = ob_data+ob_prior;
        if nargout > 1
            varargout{1} = ob_data;
            if nargout > 2
                varargout{2} = ob_prior;
            end
        end
    end

    % =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function ob = clbk_objective(logx)
        prj = project (logx);
        [ob,ob_data,ob_prior] = objective (prj,logx);
        fprintf (1, '    [LH: %f, PR: %f]\n', ob_data, ob_prior);
    end


end
