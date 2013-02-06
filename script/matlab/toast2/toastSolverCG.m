classdef toastSolverCG < toastSolver
    % Nonlinear conjugate gradient iterative solver

    methods
        function this = toastSolverCG(prm)
            % Construct the solver object.
            %
            % Syntax: solver = toastSolverCG(prm)
            %
            % Parameters:
            %         prm [struct]:
            %             structure of reconstruction parameters
            this = this@toastSolver('NCG',prm);
        end
        
        function Solve (this, x0)
            % Call the NCG solver. 
            %
            % Syntax: solver.Solve(lprm,x0)
            %
            % Parameters:
            %         x0:    initial parameter vector
            %
            % Notes:  This function should normally not be called directly.
            %         Instead, it is called by toastRecon when the
            %         solver.method field of the prm structure is set to
            %         'PCG'.
            %
            % See also:
            %         toastSolverCG, toastSolver
 
            Solve@toastSolver(this,x0);
            
            global ITR_TERM;        % termination flag

            ITR_TERM = false;
            try_again = false;

            % Initialise starting parameters
            prm  = this.prm;
            x    = x0;              % initial parameters
            logx = log(x);          % parameter transformation
            
            proj = this.project (logx);
            err0 = this.objective (proj, logx); %initial error
            err  = err0;
            errp = inf;
            step = prm.solver.step0;
            itr  = 1;

            % Return initial parameters
            this.echo_iteration_state (0, x, err);
            
            % Conjugate gradient loop
            while ((itr <= prm.solver.itmax) && ...
                    (err > prm.solver.tol*err0) && ...
                    (errp-err > prm.solver.tol) && ...
                    (ITR_TERM == false)) || try_again == true
                
                errp = err;
                
                % Gradient of cost function
                disp('Calculating gradient ...')
                r = this.gradient(x);
                
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
                    [step, err] = this.linesearch (logx, dx, step, err);
                else
                    [step, err] = this.stepsize (logx, dx, step*1.5, err);
                end
                
                % Add update to solution
                disp('Updating solution ...')
                logx = logx + dx*step;
                x = exp(logx);
                
                save(['iter',num2str(itr)],'x')
                
                proj = this.project (logx);
                err = this.objective (proj, logx);
                
                if err < errp % improvement of cost function
                    % Return solution
                    this.echo_iteration_state (itr, x, err);
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
            
        end % Solve
        
    end % methods
    
end % toastSolverCG
