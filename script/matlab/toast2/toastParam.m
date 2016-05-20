classdef toastParam
    % A hierarchical container structure for the parameters defining a Toast
    % forward or inverse problem. This structure is used by some of the
    % high-level scripts like toastRecon, and is a convenient way of storing
    % and passing parameters by custom scripts.
    %
    % See also: toastRecon
    
    properties
        % Source/detector specification
        %
        % Can contain the following fields:
        % meas.qmfile [string]:   name of source-detector (QM) definition file
        % meas.src.type [string]: source type ['Neumann','Dirichlet']
        % meas.src.prof [string]: source profile ['Gaussian','Cosine','Point']
        % meas.src.width [real]:  source width [mm]. E.g. sigma for Gaussian profile
        % meas.det.prof [string]: detector profile ['Gaussian','Cosine','Point']
        % meas.det.width [real]:  detector width [mm]. E.g. sigma for Gaussian profile
        meas = [];
        
        % Measurement data subsection
        %
        % Can contain the following fields:
        % data.freq [real]:             modulation frequency [MHz]. 0 for steady-state
        % data.lnamp [real vector]:     log amplitude target data for each source-detector pair
        % data.lnampfile [string]:      file name containing log amplitude target data
        % data.phase [real vector]:     phase target data for each source-detector pair
        % data.phasefile [string]:      file name containing phase target data
        % data.useref [bool]:           true if reference data are defined
        % data.ref.lnamp [real vector]: log amplitude reference data for each source-detector pair
        % data.ref.lnampfile [string]:  file name containing log amplitude reference data
        % data.ref.phase [real vector]: phase reference data for each source-detector pair
        % data.ref.phasefile [string]:  file name containing phase reference data
        %
        % Note: The standard behaviour of a script interpreting a toastParam
        % instance should be as follows: If the data.lnamp field is defined,
        % a script should ignore the data.lnampfile entry. If only the data.lnampfile
        % entry is defined, a script should read the data from that file and store
        % them in data.lnamp. Similar for all other data entries.
        data = [];

        % FEM forward solver parameters
        %
        % Can contain the following fields:
        % fwdsolver.hmesh [object]:    toastMesh instance
        % fwdsolver.meshfile [string]: Mesh file name (in Toast format)
        % fwdsolver.method [string]:   linear solver method ['direct','cg','bicg','bicgstab','gmres']
        % fwdsolver.tol [real]:        solver tolerance (ignored for 'direct')
        %
        % Note: The standard behaviour of a script interpreting a toastParam
        % instance should be as follows: If the fwdsolver.hmesh field is
        % defined, fwdsolver.meshfile is ignored. If ony fwdsolver.meshfile
        % is defined, the script should read the mesh from file and store
        % it in fwdsolver.hmesh.
        fwdsolver = [];
        
        % Inverse solver parameters
        %
        % Can contain the following fields:
        % solver.method [string]:                 inverse solver method ['pcg','lbfgs','lm']
        % solver.tol [real]:                      convergence criterion
        % solver.itmax [integer]:                 iteration limit
        % solver.lsearch [bool]:                  invoke line search along search direction
        % solver.basis.bdim [dx1 integer vector]: inverse solver basis dimensions
        % solver.basis.gdim [dx1 integer vector]: intermediate high-res basis dimensions
        %
        % Note: This section is only interpreted by the inverse solver, and
        % can be left empty if the toastParam instance is passed to a
        % forward solver.
        solver = [];
        
        % (Initial) optical parameter distribution
        %
        % Can contain the following fields:
        % initprm.mua absorption parameter specs
        % initprm.mus scattering parameter specs
        % initprm.ref refractive index specs
        %
        % Note: Each of the three parameter specs can contain any of the
        % following sub-fields:
        % .nprm [real vector]: nodal parameter values
        % .bprm [real vector]: parameter values in inverse basis
        % .sprm [real vector]: parameter values in reduced solution basis
        % .reset [string]: reset method ['homog','nim']
        % .val [real]: homogeneous reset value (only used if .reset=='homog')
        % .nim [string]: name of nodal image file (only used if .reset=='nim')
        %
        % Note: The standard behaviour of a script interpreting a toastParam
        % instance should be as follows: If any of the .nprm, .bprm or
        % .sprm fields are defined, the .reset, .val. and .nim fields are
        % ignored. Otherwise, the .nprm, .bprm and .sprm fields are
        % constructed by interpreting the .reset, .val and .nim fields.
        % Mapping between .nprm, .bprm and .sprm is performed via the Map
        % method of the basis mapper instance. Not all applications may
        % require all three basis expansions. For example, a forward solver
        % application may only require the .nprm expansion.
        initprm = [];
        
        % Regularisation specification
        %
        % Can contain the following fields:
        % regul.method [string]: regularisation method ['TV','TK0','TK1','Huber','PM','QPM',Tukey']
        % regul.tau [real]: regularisation hyperparameter
        % regul.x0 [real vector]: initial solution vector
        %
        % Note: depending on the regularisation method, additional fields
        % may be required. See toastRegul for a full description.
        % See also: toastRegul
        regul = [];
        
        % Script-specific parameters
        %
        % This section can contain arbitrary user-defined fields. The
        % interpretation is up to the called script.
        user = [];
        
        transient = [];
    end
    
    methods
        function obj = toastParam(fname)
            % Creates a new parameter object.
            %
            % Syntax: prm = toastParam
            %         prm = toastParam(fname)
            %
            % Parameters:
            %         fname [string]:
            %             file name to read parameters from
            %
            % Return values:
            %         prm [object]:
            %             the new parameter object
            if nargin > 0
                obj.meas = obj.ReadMeas(fname);
                obj.data = obj.ReadData(fname);
                obj.initprm = obj.ReadInitprm(fname);
                obj.fwdsolver = obj.ReadFwdsolver(fname);
                obj.solver = obj.ReadSolver(fname);
                obj.regul = obj.ReadRegul(fname);
            end
        end
        
        function obj = Read(obj,fname)
            % Reads parameters from file
            %
            % Syntax: prm = prm.Read(fname)
            %
            % Parameters:
            %         fname [string]:
            %             file name to read parameters from
            %
            % Return values:
            %         prm [object]:
            %             parameter object containing the read values
            obj = toastParam(fname);
        end
        
        function Write(obj,fname)
            % Writes parameters to a file
            %
            % Syntax: prm.Write(fname)
            %
            % Parameters:
            %         fname [string]:
            %             file name to write parameters to
            fid = fopen(fname,'w');
            obj.WriteMeas(obj.meas,fid);
            obj.WriteData(obj.data,fid);
            obj.WriteInitprm(obj.initprm,fid);
            obj.WriteFwdsolver(obj.fwdsolver,fid);
            obj.WriteSolver(obj.solver,fid);
            obj.WriteRegul(obj.regul,fid);
            fclose(fid);
        end
        
        function Echo(obj,prm,prefix)
            % Display the parameter values on screen
            %
            % Syntax: prm.Echo
            if nargin < 2
                prm = obj;
                prefix='';
            end
            fields = fieldnames(prm);
            for i=1:size(fields,1)
                f = getfield(prm,fields{i});
                if strcmpi(fields{i},'callback'), f = '<internal structure>'; end
                if isstruct(f)
                    obj.Echo(f,[prefix fields{i} '.']);
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
        
        % ------------------------------------------------------
        
        function WriteMeas(obj,prm,fid)
            if isfield(prm,'qmfile')
                obj.saveprm(fid,'QMFILE',prm.qmfile);
            end
            if isfield(prm,'src')
                obj.WriteMeasSrc(prm.src,fid);
            end
            if isfield(prm,'det')
                obj.WriteMeasDet(prm.det,fid);
            end
        end
        
        function WriteMeasSrc(obj,prm,fid)
            if isfield(prm,'type')
                obj.saveprm(fid,'SOURCETYPE',prm.type);
            end
            if isfield(prm,'prof')
                obj.saveprm(fid,'SOURCEPROFILE',prm.prof);
            end
            if isfield(prm,'width')
                obj.saveprm(fid,'SOURCEWIDTH',num2str(prm.width));
            end
        end
        
        function WriteMeasDet(obj,prm,fid)
            if isfield(prm,'prof')
                obj.saveprm(fid,'MEASUREMENTPROFILE',prm.prof);
            end
            if isfield(prm,'width')
                obj.saveprm(fid,'MEASUREMENTWIDTH',num2str(prm.width));
            end
        end
        
        function WriteData(obj,prm,fid)
            if isfield(prm,'freq')
                obj.saveprm(fid,'FREQ',num2str(prm.freq));
            end
            if isfield(prm,'lnampfile')
                obj.saveprm(fid,'DATA_MOD',prm.lnampfile);
            end
            if isfield(prm,'phasefile')
                obj.saveprm(fid,'DATA_ARG',prm.phasefile);
            end
            if isfield(prm,'useref')
                obj.savelogical(fid,'REFDATA',prm.useref);
                if isfield(prm,'ref')
                    obj.WriteDataRef(prm.ref,fid);
                end
            end
        end
        
        function WriteDataRef(obj,prm,fid)
            if isfield(prm,'lnampfile')
                obj.saveprm(fid,'REFDATA_MOD',prm.lnampfile);
            end
            if isfield(prm,'phasefile')
                obj.saveprm(fid,'REFDATA_ARG',prm.phasefile);
            end
        end
        
        function WriteInitprm(obj,prm,fid)
            if isfield(prm,'mua')
                obj.WritePrm(prm.mua,'MUA',fid);
            end
            if isfield(prm,'mus')
                obj.WritePrm(prm.mus,'MUS',fid);
            end
            if isfield(prm,'ref')
                obj.WritePrm(prm.ref,'N',fid);
            end
        end
        
        function WritePrm(obj,prm,label,fid)
            if isfield(prm,'reset')
                switch lower(prm.reset)
                    case 'homog'
                        obj.saveprm(fid,['RESET_' label],['HOMOG ' num2str(prm.val)]);
                    case 'nim'
                        obj.saveprm(fid,['RESET_' label],['NIM ' prm.nim]);
                    case 'nim_last'
                        obj.saveprm(fid,['RESET_' label],['NIM_LAST ' prm.nim]);
                end
            end
        end
        
        function WriteFwdsolver(obj,prm,fid)
            if isfield(prm,'meshfile')
                obj.saveprm(fid,'MESHFILE',prm.meshfile);
            end
            if isfield(prm,'method')
                obj.saveprm(fid,'LINSOLVER',prm.method);
            end
            if isfield(prm,'tol') && ~isempty(prm.tol)
                obj.saveprm(fid,'LINSOLVER_TOL',prm.tol);
            end
        end
        
        function WriteSolver(obj,prm,fid)
            if isfield(prm,'basis')
                obj.WriteSolverBasis(prm.basis,fid);
            end
            if isfield(prm,'method')
                obj.saveprm(fid,'SOLVER',prm.method);
            end
            if isfield(prm,'tol')
                obj.saveprm(fid,'NONLIN_TOL',prm.tol);
            end
            if isfield(prm,'itmax')
                obj.saveprm(fid,'NONLIN_ITMAX',prm.itmax);
            end
            if isfield(prm,'dscale')
                obj.saveprm(fid,'DATA_SCALING',prm.dscale);
            end
            if isfield(prm,'lsearch')
                if prm.lsearch == true
                    obj.saveprm(fid,'LM_LINESEARCH','TRUE');
                else
                    obj.saveprm(fid,'LM_LINESEARCH','FALSE');
                end
            end
        end
        
        function WriteSolverBasis(obj,prm,fid)
            if isfield(prm,'bdim')
                obj.saveprm(fid,'BASIS',num2str(prm.bdim));
            end
            if isfield(prm,'gdim')
                obj.saveprm(fid,'GRID',num2str(prm.gdim));
            end
        end
        
        function WriteRegul(obj,prm,fid)
            if isfield(prm,'method')
                saveprm(fid,'PRIOR',prm.method);
                if strcmpi(prm.method,'NONE') == 0
                    switch upper(prm.method)
                        case 'TV'
                            if isfield(prm,'tv')
                                obj.WriteRegulTV(prm.tv,fid);
                            end
                        % more to come
                    end
                    if isfield(prm,'prior')
                        obj.WriteRegulPrior(prm.prior,fid);
                    end
                    if isfield(prm,'smooth')
                        saveprm(fid,'PRIOR_DIFFSCALE',num2str(prm.smooth));
                    end
                    if isfield(prm,'threshold')
                        saveprm(fid,'PRIOR_PMTHRESHOLD',num2str(prm.threshold));
                    end
                end
            end
        end
        
        function WriteRegulTV(obj,prm,fid)
            if isfield(prm,'beta')
                saveprm(fid,'TV_BETA',prm.beta);
            end
        end
        
        function WriteRegulPrior(obj,prm,fid)
            if isfield(prm,'refname')
                saveprm(fid,'PRIOR_KAPREFIMG',prm.refname);
            end
        end
        
        % -------------------------------------------------
        
        function prm = ReadMeas(obj,fname)
            prm.qmfile = obj.scanprm(fname,'QMFILE');
            prm.src = obj.ReadMeasSrc(fname);
            prm.det = obj.ReadMeasDet(fname);
        end
        
        function prm = ReadMeasSrc(obj,fname)
            prm.type = obj.scanprm(fname,'SOURCETYPE');
            prm.prof = obj.scanprm(fname,'SOURCEPROFILE');
            prm.width = str2double(obj.scanprm(fname,'SOURCEWIDTH'));
        end
        
        function prm = ReadMeasDet(obj,fname)
            prm.prof = obj.scanprm(fname,'MEASUREMENTPROFILE');
            prm.width = str2double(obj.scanprm(fname,'MEASUREMENTWIDTH'));
        end
        
        function prm = ReadData(obj,fname)
            prm.freq = str2double(obj.scanprm(fname,'FREQ'));
            prm.lnampfile = obj.scanprm(fname,'DATA_MOD');
            prm.phasefile = obj.scanprm(fname,'DATA_ARG');
            prm.useref = lower(obj.scanprm(fname,'REFDATA'));
            if strcmp(prm.useref,'true')
                prm.useref = true;
                prm.ref.lnampfile = obj.scanprm(fname,'REFDATA_MOD');
                prm.ref.phasefile = obj.scanprm(fname,'REFDATA_ARG');
            else
                prm.useref = false;
            end
        end
        
        function prm = ReadInitprm(obj,fname)
            prm.mua = obj.setinitprm(obj.scanprm(fname,'RESET_MUA'));
            prm.mus = obj.setinitprm(obj.scanprm(fname,'RESET_MUS'));
            prm.ref = obj.setinitprm(obj.scanprm(fname,'RESET_N'));
        end
        
        function prm = ReadFwdsolver(obj,fname)
            prm.meshfile = obj.scanprm(fname,'MESHFILE');
            prm.method = obj.scanprm(fname,'LINSOLVER');
            prm.tol = obj.scanprm(fname,'LINSOLVER_TOL');
            if isempty(prm.tol)
                prm.tol = 1e-10;
            else
                prm.tol = str2double(prm.tol);
            end
        end
        
        function prm = ReadSolver(obj,fname)
            tmp = obj.scanprm(fname,'BASIS');
            if ~isempty(tmp)
                prm.basis.bdim = str2num(tmp);
            end
            tmp = obj.scanprm(fname,'GRID');
            if ~isempty(tmp)
                prm.basis.gdim = str2num(tmp);
            end
            prm.method = obj.scanprm(fname,'SOLVER');
            tmp = obj.scanprm(fname,'NONLIN_TOL');
            if ~isempty(tmp)
                prm.tol = str2num(tmp);
            end
            tmp = obj.scanprm(fname,'NONLIN_ITMAX');
            if ~isempty(tmp)
                prm.itmax = tmp;
            end
            tmp = obj.scanprm(fname,'DATA_SCALING');
            if ~isempty(tmp)
                prm.dscale = tmp;
            end
            tmp = obj.scanprm(fname,'LM_LINESEARCH');
            if ~isempty(tmp)
                if strcmpi(tmp,'true')
                    prm.lsearch = true;
                else
                    prm.lsearch = false;
                end
            end
        end
        
        function prm = ReadRegul(obj,fname)
            prm.method = obj.scanprm(fname,'PRIOR');
            if isempty(prm.method)
                prm.method = 'None';
            end
            if ~strcmpi(prm.method,'None')
                prm.tau = str2num(obj.scanprm(fname,'PRIOR_TAU'));
                switch upper(prm.method)
                    case 'TV'
                        prm.tv = obj.ReadRegulTV(obj,fname);
                end
                prm.prior = struct('refname',obj.scanprm(fname,'PRIOR_KAPREFIMG'), ...
                    'smooth',obj.scanprm(fname,'PRIOR_DIFFSCALE'), ...
                    'threshold',obj.scanprm(fname,'PRIOR_PMTHRESHOLD'));
                if isempty(prm.prior.refname) || strcmpi(prm.prior.refname,'none')
                    prm = rmfield(prm,'prior');
                end
            end
        end
        
        function prm = ReadRegulTV(obj,fname)
            prm.beta = obj.scanprm(fname,'TV_BETA');
        end
    end
    
    methods (Access=private)
        function saveprm(obj,fid,item,value)
            fprintf (fid, '%s = %s\n', item, value);
        end
        
        function savelogical(obj,fid,item,value)
            if value
                fprintf (fid, '%s = TRUE\n', item);
            else
                fprintf (fid, '%s = FALSE\n', item);
            end
        end
        
        function p = scanprm(obj,fname,item)
            % scan parameter file for a tag and return its value
            p = '';
            fid = fopen(fname);
            if fid ~= -1
                while 1
                    tline = fgetl(fid);
                    if ~ischar(tline), break, end
                    [tag,val] = obj.scanline(tline);
                    if strcmpi(tag,item) == true
                        p = val;
                        break;
                    end
                end
                fclose(fid);
            end
        end
        
        function [tag,val] = scanline(obj,line)
            % split a scan line into a tag/value pair
            tag = '';
            val = '';
            idx = find(line=='=');
            if length(idx) == 1
                tag = sscanf(line(1:idx-1),'%s');
                for i=idx+1:length(line)
                    if line(i) ~= ' '
                        break;
                    end
                end
                val = line(i:length(line));
            end
        end
        
        function prm = setinitprm(obj,str)
            C = textscan(str,'%s %s');
            method = upper(cell2mat(C{1}));
            prm.reset = method;
            switch method
                case 'HOMOG'
                    prm.val = str2double(cell2mat(C{2}));
                case 'NIM'
                    prm.nim = cell2mat(C{2});
                    prm.nimlast = false;
                case 'NIM_LAST'
                    prm.nim = cell2mat(C{2});
                    prm.nimlast = true;
            end
        end
    end
    
end