% This script demonstrates the assembly of the system matrix
% for the DOT problem

% Load a mesh
hmesh = toastMesh('../../test/2D/meshes/circle25_32.msh');

% Construct the nodal coefficient vectors (homogeneous for simplicity)
nnd = hmesh.NodeCount();
nel = hmesh.ElementCount();
mua = ones(nnd,1)*0.01;  % absorption coefficient [1/mm]
mus = ones(nnd,1)*1;     % scattering coefficient [1/mm]
ref = ones(nnd,1)*1.4;   % refractive index
freq = 100;              % modulation frequency [MHz]

%% Direct assembly with toastDotSysmat call

% Frequency-domain problem: complex case
smat = dotSysmat(hmesh,mua,mus,ref,freq);

%% manual matrix assembly from individual terms

% Map parameters into toast format
c0 = 0.3;                              % vacuum light speed [mm/ps]
c  = c0./ref;                          % medium light speed
cmua = c.*mua;                         % absorption (c*mua)
ckap = c./(3.0*(mua+mus));             % diffusion (c*kappa)
zeta = c./(2.0.*toastDotBndterm(ref,'Keijzer')); % boundary term (c/2A)
omega = freq*2.0*pi*1e-6;              % modulation frequency [cycles/ps]

K = hmesh.SysmatComponent ('PDD', ckap);
C = hmesh.SysmatComponent ('PFF', cmua);
B = omega*hmesh.SysmatComponent ('FF');
A = hmesh.SysmatComponent ('BndPFF', zeta);
smat2 = K + C + A + 1i*B;

%% manual matrix assembly
% Now let's assemble the system matrix manually from individual
% element contributions

smat3 = sparse(nnd,nnd);

for i=1:nel
    el = hmesh.Element(i);
    
    % diffusion contribution
    int_kapDD = el.Mat('PDD',ckap);
    
    % absorption contribution
    int_muaFF = el.Mat('PFF',cmua);

    % frequency contribution
    int_omegaFF = omega*el.Mat('FF');

    % boundary contribution
    bint_zetaFF = el.Mat('BndPFF',zeta);
    
    % assemble into global matrix
    elidx = el.Dof();  % global degrees of freedom for element nodes
    smat3(elidx,elidx) = smat3(elidx,elidx) + ...
        int_kapDD + int_muaFF + bint_zetaFF + 1i*int_omegaFF;

end

%% Compare the two matrices
ds1 = smat-smat2;
ds2 = smat-smat3;
err1 = norm(nonzeros(ds1));
err2 = norm(nonzeros(ds2));
fprintf('component-wise assembly of stiffness matrix: error=%e\n', err1);
fprintf('element-wise assembly of stiffness matrix:   error=%e\n', err2);