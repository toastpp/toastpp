% This script demonstrates the assembly of the system matrix
% for the DOT problem

% Load a mesh
hmesh=toastReadMesh('../../test/2D/meshes/circle25_32.msh');

% Construct the nodal coefficient vectors (homogeneous for simplicity)
nnd = toastMeshNodeCount(hmesh);
mua = ones(nnd,1)*0.01;  % absorption coefficient [1/mm]
mus = ones(nnd,1)*1;     % scattering coefficient [1/mm]
ref = ones(nnd,1)*1.4;   % refractive index
freq = 100;              % modulation frequency [MHz]

%% Direct assembly with toastDotSysmat call

% Frequency-domain problem: complex case
smat = toastSysmat(hmesh,mua,mus,ref,freq);

%% manual matrix assembly from individual terms

% Map parameters into toast format
c0 = 0.3;                              % vacuum light speed [mm/ps]
c  = c0./ref;                          % medium light speed
cmua = c.*mua;                         % absorption (c*mua)
ckap = c./(3.0*(mua+mus));             % diffusion (c*kappa)
zeta = c./(2.0.*toastBndReflectionTerm(ref,'Keijzer')); % boundary term (c/2A)
omega = freq*2.0*pi*1e-6;              % modulation frequency [cycles/ps]

K = toastSysmatComponent (hmesh, 'PDD', ckap);
C = toastSysmatComponent (hmesh, 'PFF', cmua);
B = omega*toastSysmatComponent (hmesh, 'FF');
A = toastSysmatComponent (hmesh, 'BndPFF', zeta);
smat2 = K + C + A + 1i*B;

%% manual matrix assembly from individual element integrals
% Now let's assemble the system matrix manually from individual
% element contributions

% We need the element connectivity list (idx) to map from local to global DOFs
[vtx idx, eltp] = toastMeshData(hmesh);  
nel = size(idx,1);

smat3 = sparse(nnd,nnd);

for el=1:nel
    elidx = idx(el,:);
    
    % diffusion contribution
    Kel = toastElmat(hmesh,el,'PDD',ckap);
    
    % absorption contribution
    Cel = toastElmat(hmesh,el,'PFF',cmua);

    % frequency contribution
    Bel = omega*toastElmat(hmesh,el,'FF');

    % boundary contribution
    Ael = toastElmat(hmesh,el,'BndPFF',zeta);
    
    % assemble into global matrix
    smat3(elidx,elidx) = smat3(elidx,elidx) + ...
        Kel + Cel + Ael + 1i*Bel;

end

%% Compare the two matrices
ds1 = smat-smat2;
ds2 = smat-smat3;
err1 = norm(nonzeros(ds1));
err2 = norm(nonzeros(ds2));
fprintf('component-wise assembly of stiffness matrix: error=%e\n', err1);
fprintf('element-wise assembly of stiffness matrix:   error=%e\n', err2);