classdef toastSolverLBFGS < toastSolver
    % Limited-bandwidth BFGS iterative solver

    methods
        function this = toastSolverLBFGS(prm)
            % Construct the solver object.
            %
            % Syntax: solver = toastSolverLBFGS(prm)
            %
            % Parameters:
            %         prm [struct]:
            %             structure of reconstruction parameters
            this = this@toastSolver('LBFGS',prm);
        end
        
        function Solve (this, x0)
            % Call the L-BFGS solver.
            %
            % Syntax: solver.Solve(x0)
            %
            % Parameters:
            %         x0 [real array]:
            %             initial parameter vector
            %
            % Notes:  This function should normally not be called directly. Instead,
            %         it is called by toastRecon when the solver.method field of the
            %         prm structure is set to 'LBFGS'.
            %
            % See also:
            %         toastRecon
            
            Solve@toastSolver(this,x0);
            
            global ITR_TERM;        % termination flag
            
            ITR_TERM = false;
            c0 = 0.3;               % speed of light
            history = 10;           % number of previous steps to store
            nx = length(x0);
            np = nx/2;
            
            prm  = this.prm;
            x    = x0;              % initial parameters
            logx = log(x);          % parameter transformation
            logx0 = logx;
            
            % Initialise starting parameters
            proj = this.project (logx);
            err0 = this.objective (proj, logx); %initial error
            err  = err0;
            errp = inf;
            step = prm.solver.step0;
            itr  = 1;
            
            % Return initial parameters
            this.echo_iteration_state (0, x, err);
            
            % initial gradient
            g = -this.gradient(x);
            
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
                
                [alpha,fmin] = this.linesearch (logx, p, step, err);
                if fmin < err
                    step = alpha;
                    err = fmin;
                else
                    % no improvement - reset to gradient direction
                    p = -H0 .* g;
                    [alpha,fmin] = this.linesearch (logx, p, step, err);
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
                this.echo_iteration_state (itr, x, err);
                
                % update gradient
                g1 = -this.gradient(x);
                
                proj = this.project (logx);
                err = this.objective (proj, logx);
                
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
                        
        end % Solve
        
    end % methods
    
end % toastSolverLBFGS
