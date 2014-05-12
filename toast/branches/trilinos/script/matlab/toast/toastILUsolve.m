function [x, ob_relres, iter] = toastILUsolve (A, b, p_method, s_method, relres, droptol, x0);
%toastILUsolve        - solve linear system with ILU preconditioning.
%
% Synopsis: [x, ob_relres, iter]] = toastILUsolve(A, b, p_method, s_method, relres, droptol, x0);
%
% General: this wonderful function precondion and solve a sparse symmetric system
% matrix using either matlab's precondioning tools (incomplete Choleski or
% LU), matlab's solvers (pcg, bicgstab or gmres if you insist), or
% altrernatively using ILUPACK inverse based multi-level precondioners and
% solvers (pcg, sqmr, gmres, fgmres)
%
% Input:
% A - system matrix (nxn)
% b - RHS (nxm)
% p_method - a flag choosing the preconditioner
%   'c' Matlab's incomplete Choleski
%   'l' Matlab's incomplete LU
%   'a' ILUPACK inverse-based multi - level preconditioner
% s_method - a flag choosing the solver
%   'p' Matlab's PCG
%   'b' Matlab's stabilised Bi-Conjugate Gradients
%   'm' Matlab's MinRes
%   'g' Matlab's GMRes
%   's' ILUPACK solver (could be either of the ones mentioned above, in later stages the ysytme charcterisation wil be done only once)
% relres - desired relative residual error, typically between 1e-6 to 1e-12 (optional)
% droptol - set the drop tolerance for the incomplete factorisation, typically use 1e-2 or 1e-3
% x0 - inital guess for the the solver (nxm)
%
% Output:
% x - the solution (nxm)
% ob_relres - retrieved relative residual error
% iter - number of external iterations preformed
%
% Written by Lior Horesh and Martin Schweiger 16/02/06
%--------------------------------------------------------------------------
% set defaults for testing
if ~exist('relres','var')
    relres = 1e-8; % set desired relative residual error
end
if ~exist('droptol','var')
    droptol =  1e-2;
end
if ~exist('x0','var')
    x0 = zeros(size(b));
end
nrestart = 20; % set gmres restarts sequence length
maxit = 500;
%--------------------------------------------------------------------------
% preconditioning stage
switch lower(p_method)
    case 'c'
        % factor by matlab's incomplete Choleski
        tic,[u] = cholinc(A,droptol); Mchol_p = toc;
        disp(['Matlab incomplete Cholesky preconditioning took ',num2str(Mchol_p), ' seconds'])
    case 'l'
        % factor by matlab's incomplete LU
        tic,[l,u] = luinc(A,droptol); Mlu_p = toc;
        disp(['Matlab incomplete LU preconditioning took ',num2str(Mlu_p), ' seconds'])
    otherwise
        options = AMGinit(A);
        options.droptol = droptol;
        options.nrestart = nrestart;
        options.maxit = maxit;
        if droptol<1e-2
            options.elbow = 30;
        end
        if ~isreal(A)
            options.isdefinite = false;
            options.ordering = 'rcm';
            options.restol = relres*norm(b(:,1));
        else
            options.isdefinite = true;
            options.ordering = 'amd';
            options.restol = relres;
        end
        % factor using ILUPack
        tic,[prec,options]=AMGfactor(A,options); ILU_p = toc;
        disp(['ILUPack preconditioning took ',num2str(ILU_p), ' seconds'])
end
%--------------------------------------------------------------------------
% solution stage
switch lower(s_method)
    case 'p'
        if isreal(A)
            for ii=1:size(b,2)
                tic,[x(:,ii),flag,ob_relres(ii),iter(ii)] = pcg(A,b(:,ii),relres,maxit,u',u,x0(:,ii)); MChol_Mpcg_s = toc;
                disp(['Matlab pcg with Matlab incomplete Cholesky preconditioner took ', num2str(MChol_Mpcg_s), ' seconds, relres ', num2str(ob_relres(ii)), ' iter ', num2str(iter(ii))])
            end
        else
            s_method = 'b';
            disp('come on... you cannot use PCG over this complex matrix... my godness. we will try to run bicgstab, but no promises')
        end
    case 'b'
        for ii=1:size(b,2)
            tic,[x(:,ii),flag,ob_relres(ii),iter(ii)] = bicgstab(A,b(:,ii),relres,maxit,l,u,x0(:,ii)); Mlu_Mbicg_s = toc;
            disp(['Matlab bicgstab with Matlab iLU preconditioner took ', num2str(Mlu_Mbicg_s), ' seconds, relres ', num2str(ob_relres(ii)), ' iter ', num2str(iter(ii))])
        end
    case 'm'
        for ii=1:size(b,2)
            tic,[x(:,ii),flag,ob_relres(ii),iter(ii)] = minres(A,b(:,ii),relres,maxit,l,u,x0(:,ii)); Mlu_Mminres_s = toc;
            disp(['Matlab minres with Matlab iLU preconditioner took ', num2str(Mlu_Mminres_s), ' seconds, relres ', num2str(ob_relres(ii)), ' iter ', num2str(iter(ii))])
        end
    case 'g'
        for ii=1:size(b,2)
            tic,[x(:,ii),flag,ob_relres(ii),iter(ii,:)] = gmres(A,b(:,ii),nrestart,relres, maxit,l,u,x0(:,ii));  Mlu_Mgmres_s = toc;
            disp(['Matlab gmres with Matlab iLU preconditioner took ', num2str(Mlu_Mgmres_s), ' seconds, relres ', num2str(ob_relres(ii)), ' iter ', num2str(iter(ii,:))])
        end
    otherwise
        for ii=1:size(b,2)
            options.x0 = x0(:,ii);
            tic,[x(:,ii), options] = AMGsolver(A, prec, full(b(:,ii)), options); ILU_s = toc; ob_relres(ii) = norm(A*x(:,ii)-b(:,ii))/norm(b(:,ii)); iter(ii) = options.niter;
            disp(['ILUPack solver with ILUPack preconditioner took ', num2str(ILU_s), ' seconds, relres ', num2str(ob_relres(ii)), ' iter ', num2str(options.niter)])
        end
end

%   % uncomment this section in case you would like to comapre matlab's
%   % preconditioner to ILUPACK preconditioner. here ILUPACK is used for
%   % preconditioneing, and pcg or sqmr2 is used as a solver. notice that
%   % sqmr2 need to be placed in the matlab/toolbox/matlab/sparfun folder and
%   % then rehash the toolbox
%
%    % ilupack preconditioner + matlab solver (pcg / sqmr)
%     if (p_method == 'a') & (s_method=='p')
%         tic,[x,flag,ILU_Mpcg_relres(ii),ILU_Mpcg_iter(ii)] = pcg(A,b,options.restol,options.maxit,@AMGsol2,[], options.x0,prec); ILU_Mpcg_s(ii) = toc;
%         disp(['Matlab pcg with ILUPack preconditioner took ', num2str(ILU_Mpcg_s(ii)), ' seconds, relres ', num2str(ILU_Mpcg_relres(ii)), ' iter ', num2str(ILU_Mpcg_iter(ii))])
%
%         % % %   commented until Matthias will provide a solution for the segmentation fault issue
%         % % elseif  (p_method == 'a') & (s_method~='s')
%         % %     rehash toolcache
%         % %     tic,[x,flag,ILU_Msqmr_relres(ii),ILU_Msqmr_iter(ii)] = sqmr2(A,b,options.restol,options.maxit,@AMGsol2,[], options.x0,prec);ILU_Msqmr_s(ii) = toc;
%         % %     disp(['Matlab sqmr with ILUPack preconditioner took ', num2str(ILU_Msqmr_s(ii)), 'seconds, relres ', num2str(ILU_Msqmr_relres(ii)), ' iter ', ILU_Msqmr_iter(ii)])
%
%     end
