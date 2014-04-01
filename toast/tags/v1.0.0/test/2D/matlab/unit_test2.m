clear all
close all

generate_ref = false;

% Make sure the toast paths are available
cwd = pwd;
cd (getenv('TOASTDIR'));
mtoast2_install(true);
cd (cwd);

if ~generate_ref
    load unit2_test.mat
end

meshdir = '../meshes/';

fprintf('toast toolbox unit tests\n\n')

%% test: toastMesh.Read, toastMesh.Data

hmesh = toastMesh([meshdir 'circle25_32.msh']);
[vtx idx eltp] = hmesh.Data();

if generate_ref
    test.ReadMeshMeshData1.vtx = vtx;
    test.ReadMeshMeshData1.idx = idx;
    test.ReadMeshMeshData1.eltp = eltp;
else
    err = max (norm (test.ReadMeshMeshData1.vtx-vtx));
    if err < 1e-8
        fprintf('PASSED ReadMesh/MeshData: vertex coordinate max displacement %f\n', err)
    else
        error('ReadMesh/MeshData: vertex coordinates displaced (%f)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.idx-idx)));
    if err==0
        fprintf('PASSED ReadMesh/MeshData: element node indices match\n')
    else
        error('ReadMesh/MeshData: element node index mismatch (%d)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.eltp-eltp)));
    if err==0
        fprintf('PASSED ReadMesh/MeshData: element type indices match\n')
    else
        error('ReadMesh/MeshData: element type index mismatch (%d)', err)
    end
end


%% test: toastDotSysmat

hmesh = toastMesh([meshdir 'circle25_32.msh']);
n = hmesh.NodeCount;
mua = ones(n,1)*0.01;
mus = ones(n,1)*2.0;
ref = ones(n,1)*1.4;
freq = 100;
smat = toastDotSysmat(hmesh,mua,mus,ref,freq);
if generate_ref
    test.DotSysmat = nonzeros(smat);
else
    err = max(norm(test.DotSysmat-nonzeros(smat)));
    if err < 1e-8
        fprintf('PASSED DotSysmat: matrix difference norm=%e\n', err);
    else
        error('DotSysmat: matrix difference norm=%e\n', err);
    end
end


%% test: toastElmat

[vtx idx, eltp] = hmesh.Data();  
nel = size(idx,1);

% Map parameters into toast format
c0 = 0.3;                              % vacuum light speed [mm/ps]
c  = c0./ref;                          % medium light speed
cmua = c.*mua;                         % absorption (c*mua)
ckap = c./(3.0*(mua+mus));             % diffusion (c*kappa)
zeta = c./(2.0.*toastDotBndterm(ref,'Keijzer')); % boundary term (c/2A)
omega = freq*2.0*pi*1e-6;              % modulation frequency [cycles/ps]

smatK = sparse(n,n);
smatC = sparse(n,n);
smatA = sparse(n,n);
smatB = sparse(n,n);
smat = sparse(n,n);

for el=1:nel
    elidx = idx(el,:);
    
    % diffusion contribution
    int_kapDD = hmesh.Elmat(el,'PDD',ckap);
    
    % absorption contribution
    int_muaFF = hmesh.Elmat(el,'PFF',cmua);

    % frequency contribution
    int_omegaFF = omega*hmesh.Elmat(el,'FF');

    % boundary contribution
    bint_zetaFF = hmesh.Elmat(el,'BndPFF',zeta);
    
    % assemble into global matrix
    smatK(elidx,elidx) = smatK(elidx,elidx) + int_kapDD;
    smatC(elidx,elidx) = smatC(elidx,elidx) + int_muaFF;
    smatA(elidx,elidx) = smatA(elidx,elidx) + bint_zetaFF;
    smatB(elidx,elidx) = smatB(elidx,elidx) + int_omegaFF;
    smat(elidx,elidx) = smat(elidx,elidx) + ...
        int_kapDD + int_muaFF + bint_zetaFF + 1i*int_omegaFF;
end

if generate_ref
    test.ElmatK = nonzeros(smatK);
    test.ElmatC = nonzeros(smatC);
    test.ElmatA = nonzeros(smatA);
    test.ElmatB = nonzeros(smatB);
else
    err = max(norm(test.ElmatK-nonzeros(smatK)));
    if err < 1e-8
        fprintf('PASSED Elmat(PDD): matrix difference norm=%e\n', err);
    else
        error('Elmat(PDD): matrix difference norm=%e\n', err);
    end
    err = max(norm(test.ElmatC-nonzeros(smatC)));
    if err < 1e-8
        fprintf('PASSED Elmat(PFF): matrix difference norm=%e\n', err);
    else
        error('Elmat(PFF): matrix difference norm=%e\n', err);
    end
    err = max(norm(test.ElmatA-nonzeros(smatA)));
    if err < 1e-8
        fprintf('PASSED Elmat(BndPFF): matrix difference norm=%e\n', err);
    else
        error('Elmat(BndPFF): matrix difference norm=%e\n', err);
    end
    err = max(norm(test.ElmatB-nonzeros(smatB)));
    if err < 1e-8
        fprintf('PASSED Elmat(FF): matrix difference norm=%e\n', err);
    else
        error('Elmat(FF): matrix difference norm=%e\n', err);
    end
    err = max(norm(test.DotSysmat-nonzeros(smat)));
    if err < 1e-8
        fprintf('PASSED Elmat: matrix difference norm=%e\n', err);
    else
        error('Elmat: matrix difference norm=%e\n', err);
    end
end


%% Output reference data
if generate_ref
    save unit2_test.mat test
    fprintf('\ntest reference data saved\n')
else
    fprintf('\nall tests passed\n')
end