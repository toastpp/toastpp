% This script demonstrates the assembly of the system matrix
% for the DOT problem

% Load a mesh
hmesh = toastMesh('../../test/2D/meshes/circle25_32.msh');

% Construct the nodal coefficient vectors (homogeneous for simplicity)
nnd = hmesh.NodeCount();
mua = ones(nnd,1)*0.01;  % absorption coefficient [1/mm]
mus = ones(nnd,1)*1;     % scattering coefficient [1/mm]
ref = ones(nnd,1)*1.4;   % refractive index
freq = 100;              % modulation frequency [MHz]

%% Direct assembly with toastDotSysmat call

% Frequency-domain problem: complex case
smat = toastDotSysmat(hmesh,mua,mus,ref,freq);

%% manual matrix assembly
% Now let's assemble the system matrix manually from individual
% element contributions

% We need the element connectivity list (idx) to map from local to global DOFs
[vtx idx, eltp] = hmesh.Data();  
nel = size(idx,1);

% Map parameters into toast format
c0 = 0.3;                              % vacuum light speed [mm/ps]
c  = c0./ref;                          % medium light speed
cmua = c.*mua;                         % absorption (c*mua)
ckap = c./(3.0*(mua+mus));             % diffusion (c*kappa)
zeta = c./(2.0.*toastDotBndterm(ref,'Keijzer')); % boundary term (c/2A)
omega = freq*2.0*pi*1e-6;              % modulation frequency [cycles/ps]

smat2 = sparse(nnd,nnd);

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
    smat2(elidx,elidx) = smat2(elidx,elidx) + ...
        int_kapDD + int_muaFF + bint_zetaFF + 1i*int_omegaFF;

end

% Compare the two matrices
ds = smat-smat2;
err = norm(nonzeros(ds));
fprintf('norm of difference matrix = %e\n', err);