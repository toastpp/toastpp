function toastWriteParam (prmfile, prm)
%toastWriteParam      - Write parameters to a TOAST parameter file.
%
% Synopsis: toastWriteParam(prmfile,prm)
%    prmfile: parameter file name (string)
%    prm:     parameter structure
%
% This function writes the parameters in a standard TOAST parameter
% structure to a file.
  
fid = fopen(prmfile,'w');
if isfield(prm,'meas')
    if isfield(prm.meas,'qmfile')
        saveprm(fid,'QMFILE',prm.meas.qmfile);
    end
    if isfield(prm.meas,'src')
        if isfield(prm.meas.src,'type')
            saveprm(fid,'SOURCETYPE',prm.meas.src.type);
        end
        if isfield(prm.meas.src,'prof')
            saveprm(fid,'SOURCEPROFILE',prm.meas.src.prof);
        end
        if isfield(prm.meas.src,'width')
            saveprm(fid,'SOURCEWIDTH',num2str(prm.meas.src.width));
        end
    end
    if isfield(prm.meas,'det')
        if isfield(prm.meas.det,'prof')
            saveprm(fid,'MEASUREMENTPROFILE',prm.meas.det.prof);
        end
        if isfield(prm.meas.det,'width')
            saveprm(fid,'MEASUREMENTWIDTH',num2str(prm.meas.det.width));
        end
    end
end
if isfield(prm,'data')
    if isfield(prm.data,'freq')
        saveprm(fid,'FREQ',num2str(prm.data.freq));
    end
    if isfield(prm.data,'lnampfile')
        saveprm(fid,'DATA_MOD',prm.data.lnampfile);
    end
    if isfield(prm.data,'phasefile')
        saveprm(fid,'DATA_ARG',prm.data.phasefile);
    end
    if isfield(prm.data,'ref')
        if isfield(prm.data.ref,'lnampfile')
            saveprm(fid,'REFDATA_MOD',prm.data.ref.lnampfile);
        end
        if isfield(prm.data.ref,'phasefile')
            saveprm(fid,'REFDATA_ARG',prm.data.ref.phasefile);
        end
    end
end
if isfield(prm,'fwdsolver')
    if isfield(prm.fwdsolver,'meshfile')
        saveprm(fid,'MESHFILE',prm.fwdsolver.meshfile);
    end
    if isfield(prm.fwdsolver,'method')
        saveprm(fid,'LINSOLVER',prm.fwdsolver.method);
    end
    if isfield(prm.fwdsolver,'tol')
        saveprm(fid,'LINSOLVER_TOL',prm.fwdsolver.tol);
    end
end
if isfield(prm,'solver')
    if isfield(prm.solver,'basis')
        if isfield(prm.solver.basis,'bdim')
            saveprm(fid,'BASIS',num2str(prm.solver.basis.bdim));
        end
        if isfield(prm.solver.basis,'gdim')
            saveprm(fid,'GRID',num2str(prm.solver.basis.gdim));
        end
    end
    if isfield(prm.solver,'method')
        saveprm(fid,'SOLVER',prm.solver.method);
    end
    if isfield(prm.solver,'tol')
        saveprm(fid,'NONLIN_TOL',prm.solver.tol);
    end
    if isfield(prm.solver,'itmax')
        saveprm(fid,'NONLIN_ITMAX',prm.solver.itmax);
    end
    if isfield(prm.solver,'dscale')
        saveprm(fid,'DATA_SCALING',prm.solver.dscale);
    end
    if isfield(prm.solver,'lsearch')
        if prm.solver.lsearch == true
            saveprm(fid,'LM_LINESEARCH','TRUE');
        else
            saveprm(fid,'LM_LINESEARCH','FALSE');
        end
    end
end
if isfield(prm,'regul')
    if isfield(prm.regul,'method')
        saveprm(fid,'PRIOR',prm.regul.method);
        if strcmpi(prm.regul.method,'NONE') == 0
            switch prm.regul.method
                case 'TV'
                    if isfield(prm.regul,'tv')
                        if isfield(prm.regul.tv,'beta')
                            saveprm(fid,'TV_BETA',prm.regul.tv.beta);
                        end
                    end
                % more to come
            end
            if isfield(prm.regul,'prior')
                if isfield(prm.regul.prior,'refname')
                    saveprm(fid,'PRIOR_KAPREFIMG',prm.regul.prior.refname);
                end
                if isfield(prm.regul.prior,'smooth')
                    saveprm(fid,'PRIOR_DIFFSCALE',num2str(prm.regul.prior.smooth));
                end
                if isfield(prm.regul.prior,'threshold')
                    saveprm(fid,'PRIOR_PMTHRESHOLD',num2str(prm.regul.prior.threshold));
                end
            end
        end
    end
end
if isfield(prm,'initprm')
    if isfield(prm.initprm,'mua')
        if isfield(prm.initprm.mua,'reset')
            switch prm.initprm.mua.reset
                case 'HOMOG'
                    saveprm(fid,'RESET_MUA',['HOMOG ' num2str(prm.initprm.mua.val)]);
                case 'NIM'
                    saveprm(fid,'RESET_MUA',['NIM ' prm.initprm.mua.nim]);
                case 'NIM_LAST'
                    saveprm(fid,'RESET_MUA',['NIM_LAST ' prm.initprm.mua.nim]);
            end
        end
    end
    if isfield(prm.initprm,'mus')
        if isfield(prm.initprm.mus,'reset')
            switch prm.initprm.mus.reset
                case 'HOMOG'
                    saveprm(fid,'RESET_MUS',['HOMOG ' num2str(prm.initprm.mus.val)]);
                case 'NIM'
                    saveprm(fid,'RESET_MUS',['NIM ' prm.initprm.mus.nim]);
                case 'NIM_LAST'
                    saveprm(fid,'RESET_MUS',['NIM_LAST ' prm.initprm.mus.nim]);
            end
        end
    end
    if isfield(prm.initprm,'ref')
        if isfield(prm.initprm.ref,'reset')
            switch prm.initprm.ref.reset
                case 'HOMOG'
                    saveprm(fid,'RESET_N',['HOMOG ' num2str(prm.initprm.ref.val)]);
                case 'NIM'
                    saveprm(fid,'RESET_N',['NIM ' prm.initprm.ref.nim]);
                case 'NIM_LAST'
                    saveprm(fid,'RESET_N',['NIM_LAST ' prm.initprm.ref.nim]);
            end
        end
    end
end


function saveprm (fid, item, value)
fprintf (fid, '%s = %s\n', item, value);
