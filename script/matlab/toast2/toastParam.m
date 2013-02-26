classdef toastParam
    % A hierarchical container for the parameters defining a Toast forward
    % or inverse problem.
    
    properties
        meas = [];
        data = [];
        fwdsolver = [];
        solver = [];
        initprm = [];
        regul = [];
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