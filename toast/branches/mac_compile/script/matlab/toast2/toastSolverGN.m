classdef toastSolverGN < toastSolver
    % Gauss-Newton iterative solver

    methods
        function this = toastSolverGN(prm)
            % Construct the solver object.
            %
            % Syntax: solver = toastSolverGN(prm)
            %
            % Parameters:
            %         prm [struct]:
            %             structure of reconstruction parameters
            this = this@toastSolver('GN',prm);
        end

        function Solve (this, x0)
            % Call the GN solver.
            %
            % Syntax: solver.Solve(x0)
            %
            % Parameters:
            %        x0:    initial parameter vector
            %
            % Notes:  This function should normally not be called directly.
            %         Instead, it is called by toastRecon when the
            %         solver.method field of the prm structure is set to
            %         'LM'.
            %
            % See also:
            %         toastSolverGN, toastSolver
            
            Solve@toastSolver(this,x0);
            
            global ITR_TERM;        % termination flag
            
            ITR_TERM = false;
            m = length(this.lprm.data);  % data space dimension
            p = length(x0);         % parameter space dimension
            
            % Initialise starting parameters
            prm = this.prm;
            x = x0;                 % initial parameters
            logx = log(x);          % parameter transformation
            
            proj = this.project (logx);
            err0 = this.objective (proj, logx); %initial error
            err  = err0;
            errp = inf;
            step = prm.solver.step0;
            itr  = 1;
            
            % Return initial parameters
            this.echo_iteration_state (0, x, err);
            
            % Gauss-Newton loop
            while (itr <= prm.solver.itmax) && ...
                    (err > prm.solver.tol*err0) && ...
                    (errp-err > prm.solver.tol) && ...
                    (ITR_TERM == false)
                
                errp = err;
                
                % Construct the Jacobian
                disp ('Calculating Jacobian ...')
                J = this.jacobian(x);
                
                % data normalisation
                for i = 1:m, J(i,:) = J(i,:) / this.lprm.sd(i); end;
                
                % parameter normalisation
                for i = 1:p, J(:,i) = J(:,i) * x(i); end;
                
                % Normalisation of Hessian
                if this.lprm.hReg ~= 0
                    psiHdiag = this.lprm.hReg.HDiag (logx);
                else
                    psiHdiag = zeros(p,1);
                end
                for i = 1:p
                    M(i) = sum(J(:,i) .* J(:,i));
                    M(i) = M(i) + psiHdiag(i);
                    M(i) = 1 ./ sqrt(M(i));
                end
                for i = 1:p, J(:,i) = J(:,i) * M(i); end;
                
                % Gradient of cost function
                r = J' * ((this.lprm.data-proj)./this.lprm.sd);
                if this.lprm.hReg ~= 0
                    r = r - this.lprm.hReg.Gradient (logx) .* M';
                end
                
                if this.lprm.Himplicit == true
                    % Update with implicit Krylov solver
                    fprintf (1, 'Entering Krylov solver (tol=%g)\n', ...
                        prm.solver.krylov.tol);
                    if this.lprm.hReg ~= 0
                        HessReg = this.lprm.hReg.Hess (x);
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
                    [step, err] = this.linesearch (logx, dx, step, err);
                else
                    [step, err] = this.stepsize (logx, dx, step*1.5, err);
                end
                
                % Add update to solution
                disp('Updating solution ...')
                logx = logx + dx*step;
                x = exp(logx);
                
                proj = this.project (logx);
                err = this.objective (proj, logx);
                
                % Return solution
                this.echo_iteration_state (itr, x, err);                
                itr = itr+1;
                
            end % GN loop
            
            
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
                if this.lprm.hReg ~= 0
                    b = b + M' .* (HessReg * (M' .* x));
                end
            end
            
        end % Solve
        
    end % methods
    
end % toastSolverGN
