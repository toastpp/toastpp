classdef toastSolver < handle
    
    methods
        function obj = toastSolver(method,prm)
            if nargin > 0
                obj.method = method;
                obj.prm = obj.checkprm(prm);
                obj.lprm = obj.setup(prm);
            end
        end
        
        function SetRegul (this,prm,x0)
            % Set up a regularisation object for the solver
            %
            % Syntax: solver.SetRegul(prm)
            %
            % Parameters:
            %         prm [struct]:
            %             structure defining regularisation parameters
            %         x0 [vector]:
            %             initial parameter vector
            %
            % Notes:  If the structure passed to the solver constructor
            %         contains a 'regul' substructure, it calls
            %         SetRegul(prm.regul) directly.
            this.lprm.hReg = 0;

            
%DEBUG
%rprm.method = 'TV';
%rprm.tau = 1e-3;
%rprm.basis = prm.basis;
%rprm.prior.smooth = 0.5;
%rprm.prior.threshold = 0.25;
%rprm.prior.refimg = prm.prior.refimg;
%rprm.tv.beta = 0.01;
%this.lprm.hReg = toastRegul (rprm, x0);



            if isfield(prm,'method') && ~strcmpi(prm.method,'none')
                this.lprm.hReg = toastRegul (prm, x0);
            end
            
            
        end
        
        function Solve(obj, x0)
            global scref RES;
            RES = [];
            scref = obj.c0./obj.lprm.hBasis.Map ('M->S',obj.lprm.ref);
            obj.x0 = x0;
            obj.SetRegul(obj.prm.regul,x0);
        end
        
        function DispPrm(prm,prefix)
            % recursively display prm structure fields
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
                        else fprintf (1, ' false'); end
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
    end
    
    methods (Access=protected)
        function prj = project (obj,logx)
            x = exp(logx);
            scmua = x(1:size(x)/2);
            sckap = x(size(x)/2+1:end);
            cmua = obj.lprm.hBasis.Map ('S->M', scmua);
            ckap = obj.lprm.hBasis.Map ('S->M', sckap);
            mua = cmua .* (obj.lprm.ref/obj.c0);
            kap = ckap .* (obj.lprm.ref/obj.c0);
            mus = 1./(3*kap) - mua;
            prj = toastProject (obj.lprm.hMesh, mua, mus, obj.lprm.ref, obj.prm.data.freq, ...
                obj.lprm.qvec, obj.lprm.mvec, obj.prm.fwdsolver.method, obj.prm.fwdsolver.tol);
        end
        
        % =====================================================================
        
        function g = gradient(obj,x)
            mua = obj.lprm.hBasis.Map ('S->M', x(1:size(x)/2)) .* (obj.lprm.ref/obj.c0);
            kap = obj.lprm.hBasis.Map ('S->M', x(size(x)/2+1:end)) .* (obj.lprm.ref/obj.c0);
            mus = 1./(3*kap) - mua;
                
            g = -toastGradient (obj.lprm.hMesh, obj.lprm.hBasis, obj.lprm.qvec, obj.lprm.mvec, ...
                mua, mus, obj.lprm.ref, obj.prm.data.freq, obj.lprm.data, obj.lprm.sd, ...
                'method',obj.prm.fwdsolver.method, 'tolerance',obj.prm.fwdsolver.tol);
            g = g .* x; % parameter scaling
            if isfield(obj.lprm,'hReg') &&  obj.lprm.hReg ~= 0
                g = g - obj.lprm.hReg.Gradient (log(x));
            end
        end
            
        % =====================================================================
        
        function J = jacobian(obj,x)
            mua = obj.lprm.hBasis.Map ('S->M', x(1:size(x)/2)) .* (obj.lprm.ref/obj.c0);
            kap = obj.lprm.hBasis.Map ('S->M', x(size(x)/2+1:end)) .* (obj.lprm.ref/obj.c0);
            mus = 1./(3*kap) - mua;
            
            % Construct the Jacobian
            J = toastJacobian (obj.lprm.hMesh, obj.lprm.hBasis, ...
                obj.lprm.qvec, obj.lprm.mvec, mua, mus, obj.lprm.ref, ...
                obj.prm.data.freq, obj.prm.fwdsolver.method, ...
                obj.prm.fwdsolver.tol);
        end
            
        % =====================================================================
        
        function [ob,varargout] = objective (obj,prj,logx)
            ob_data = full(sum(((obj.lprm.data-prj)./obj.lprm.sd).^2));
            ob_prior = 0;
                
            if isfield(obj.lprm,'hReg') && obj.lprm.hReg ~= 0
                ob_prior = obj.lprm.hReg.Value(logx);
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
        
        function [s,pmin] = linesearch (obj, x0, d, s0, p0)
            [s,pmin] = toastLineSearch (x0, d, s0, p0, @obj.clbk_objective);
        end
        
        % =====================================================================
        
        function [step,err] = stepsize (obj,x0, d, s0, p0)
            maxstep = 5;
            step = s0;

            for trystep = 1:maxstep
                err = obj.clbk_objective (x0 + d*step);
                fprintf (1, '--> Step: %f, objective: %f\n', step, err);
                if err < p0
                    return;
                end
                step = step/2; % no improvement
            end
        end
        
       
        % =====================================================================
        
        function echo_iteration_state (obj, itr, x, err)
            if isfield(obj.lprm,'callback') && isfield(obj.lprm.callback,'iter')
                feval(obj.lprm.callback.iter, obj.lprm.callback.context, itr, x, err);
            end
        end
    end
    
    methods (Access=private)
        function prm = checkprm(this,prm)
            global RES scref;
            
            % Fill missing parameter entries
            if isfield(prm,'solver') == false || isfield(prm.solver,'method') == false || ischar(prm.solver.method) == false || isempty(prm.solver.method)
                prm.solver.method = 'PCG';
            end
            if isfield(prm.solver,'tol') == false || isnumeric(prm.solver.tol) == false
                prm.solver.tol = 1e-10;
            end
            if isfield(prm.solver,'itmax') == false || isnumeric(prm.solver.itmax) == false
                prm.solver.itmax = 100;
            end
            if isfield(prm.solver,'lsearch') == false || islogical(prm.solver.lsearch) == false
                prm.solver.lsearch = true;
            end
            if isfield(prm.solver,'step0') == false || isnumeric(prm.solver.step0) == false
                prm.solver.step0 = 1;
            end
            if isfield(prm.solver,'dscale') == false || ischar(prm.solver.dscale) == false || isempty(prm.solver.dscale)
                prm.solver.dscale = 'AVG_DIFFDATA';
            end
            switch upper(prm.solver.method)
                case 'PCG'
                    if isfield(prm.solver,'cg') == false || isfield(prm.solver.cg,'reset') == false || isnumeric(prm.solver.cg.reset) == false
                        prm.solver.cg.reset = 10;
                    end
                case {'LM' 'GN_IMPLICIT'}
                    if isfield(prm.solver,'krylov') == false || isfield(prm.solver.krylov,'method') == false || ischar(prm.solver.krylov.method) == false || isempty(prm.solver.krylov.method)
                        prm.solver.krylov.method = 'gmres';
                    end
                    if isfield(prm.solver.krylov,'tol') == false || isnumeric(prm.solver.krylov.tol) == false
                        prm.solver.krylov.tol = 1e-2;
                    end
                    if isfield(prm.solver.krylov,'maxit') == false || isnumeric(prm.solver.krylov.maxit) == false
                        prm.solver.krylov.maxit = 100;
                    end
            end
            if isfield(prm.regul,'basis') == false
                if isfield(prm.solver,'basis') == true && isfield(prm.solver.basis,'hbasis') == true
                    prm.regul.basis = prm.solver.basis.hbasis;
                else
                    error('Missing basis handle in control parameter structure');
                end
            end
        end
        
        function localprm = setup(this,prm)
            % Create all required data structures in the localparam
            % structure
            
            % Read mesh and QM definition
            if isfield(prm,'fwdsolver') && isfield(prm.fwdsolver,'hmesh')
                localprm.hMesh = prm.fwdsolver.hmesh;
            else
                localprm.hMesh = toastMesh (prm.fwdsolver.meshfile);
                localprm.hMesh.ReadQM (prm.meas.qmfile);
            end
            
            % Set up the mapper between FEM and solution bases
            if isfield(prm,'solver') && isfield(prm.solver,'basis') && isfield(prm.solver.basis,'hbasis')
                localprm.hBasis = prm.solver.basis.hbasis;
            else
                localprm.hBasis = toastBasis (localprm.hMesh, prm.solver.basis.bdim, 'Linear');
            end
            
            % Set up homogeneous initial parameter estimates
            fprintf('Initial parameter estimates:\n');
            mua = this.resetprm(prm.initprm.mua,localprm.hMesh);
            fprintf('mua: mean=%f, std=%f\n', mean(mua), std(mua));
            mus = this.resetprm(prm.initprm.mus,localprm.hMesh);
            fprintf('mus: mean=%f, std=%f\n', mean(mus), std(mus));
            localprm.ref = this.resetprm(prm.initprm.ref,localprm.hMesh);
            fprintf('ref: mean=%f, std=%f\n', mean(localprm.ref), std(localprm.ref));
            kap = 1./(3*(mua+mus));

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
            localprm.data = [mdata;pdata];
            m = length(localprm.data);
            fprintf('Data space dimension: %d\n', m);

            % Generate source vectors
            localprm.qvec = localprm.hMesh.Qvec (prm.meas.src.type, prm.meas.src.prof, prm.meas.src.width);
            fprintf ('Source vector set up: %d sources\n', size(localprm.qvec,2));

            % Generate measurement vectors
            localprm.mvec = localprm.hMesh.Mvec (prm.meas.det.prof, prm.meas.det.width, localprm.ref);
            fprintf ('Detector vector setup: %d detectors\n', size(localprm.mvec,2));
            
            % Initial data set f[x0]
            fprintf('Calculating projections from initial parameters ...\n');
            proj = toastProject (localprm.hMesh, mua, mus, localprm.ref, prm.data.freq, ...
                localprm.qvec, localprm.mvec, prm.fwdsolver.method, prm.fwdsolver.tol);

            % Difference data setup
            if isfield(prm.data,'useref') && prm.data.useref == true
                fprintf ('Using difference data\n');
                if isfield(prm.data.ref,'lnampfile')
                    fprintf ('Log amplitude reference: %s\n', prm.data.ref.lnampfile);
                    mref = toastReadVector(prm.data.ref.lnampfile);
                    if ~isempty(mref)
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
                    if ~isempty(pref)
                        pproj = proj(length(mdata)+1:end);
                        pdata = pdata-pref+pproj;
                    else
                        fprintf ('Phase reference: file not found!\n');
                        return;
                    end
                end
                localprm.data = [mdata;pdata];
            end

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
            localprm.sd = [sd_lnmod;sd_phase];

            % initial parameter estimates in solution basis
            RES.bmua = localprm.hBasis.Map ('M->B', mua);
            RES.bmus = localprm.hBasis.Map ('M->B', mus);

            bcmua = localprm.hBasis.Map ('M->B', mua .* (this.c0./localprm.ref));
            bckap = localprm.hBasis.Map ('M->B', kap .* (this.c0./localprm.ref));

            scmua = localprm.hBasis.Map ('B->S', bcmua);
            sckap = localprm.hBasis.Map ('B->S', bckap);
            scref = localprm.hBasis.Map ('M->S', this.c0./localprm.ref);

            x = [scmua;sckap];
            logx = log(x);

            localprm.Himplicit = true;                  % Implicit/explicit Hessian matrix
            
            % Parameters for display callback function
            localprm.callback.iter = @this.solve_iter;
            localprm.callback.context = prm;
        end
        
        function p=resetprm(this,cfg,hmesh)
            n = hmesh.NodeCount();
            switch upper(cfg.reset)
                case 'HOMOG'
                    p = ones(n,1) * cfg.val;
                case 'NIM'
                    p = toastNim(cfg.nim);
                    if length(p) ~= n
                        disp('Warning: incompatible size of NIM file')
                    end
                case 'NIM_LAST'
                    p = toastNim(cfg.nim,0);
                    if length(p) ~= n
                        disp('Warning: incompatible size of NIM file')
                    end
                otherwise
                    disp('Warning: Unsupported reset method')
            end
        end
        
        % Callback function for objective evaluation (called by toastLineSearch)
        function ob = clbk_objective(this,logx)
            prj = this.project (logx);
            [ob,ob_data,ob_prior] = this.objective (prj,logx);
            fprintf (1, '    [LH: %f, PR: %f]\n', ob_data, ob_prior);
        end
        
        function disp_iter (this, res)
            fprintf (1, 'Objective: %f\n', res.of);
        end
        
        function solve_iter (this, prm, itr, x, err)
            fprintf (1, '**** Iteration %d, objective %f\n', itr, err)

            global RES scref;

            ps = length(x)/2;
            scmua = x(1:ps);
            sckap = x(ps+1:2*ps);
            smua  = scmua./scref;
            skap  = sckap./scref;
            smus  = 1./(3*skap) - smua;
            RES.bmua = this.lprm.hBasis.Map ('S->B', smua);
            RES.bmus = this.lprm.hBasis.Map ('S->B', smus);
            RES.of(itr+1) = err;

            if isfield(prm.transient,'callback') && isfield(prm.transient.callback,'iter')
                % let the calling function do the display
                if isfield(prm.transient.callback,'request') % check requests for secondary results
                    if isfield(prm.transient.callback.request,'prior') && prm.transient.callback.request.prior == true
                        RES.kapref = this.lprm.hReg.Kappa (log(x));
                    end
                end
                feval(prm.transient.callback.iter, prm, RES);
            else
                % inline display
                this.disp_iter(RES);
            end
        end
    end
    
    properties
        method = 'None';
        prm = [];
        lprm = [];
        res = [];
        x0 = [];
        c0 = 0.3;
    end
end
