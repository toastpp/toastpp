function fwd2

disp('MATLAB-TOAST sample script:')
disp('Forward solver test:');
disp('Check that source input = exitance for zero absorption case');
disp('-------------------------------------------------------')

% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = {'../meshes/circle25_32.msh','../meshes/circle25_64.msh'};
qmname    = '../meshes/circle25_32x32.qm';            % source-detector file

rind = [1.0,1.5];                       % refractive index
qwidth = [2 4];                         % source width
freq = 0;                               % modulation frequency [MHz]
rad = 25;
% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
toastCatchErrors();
c0 = 0.3;

for mesh=1:2
    
    % Read a TOAST mesh definition from file.
    hMesh = toastReadMesh (meshname{mesh});
    toastReadQM (hMesh, qmname);
    n = toastMeshNodeCount (hMesh);

    for r=1:2
        refind = rind(r);
        % Set up some variables
        cm = c0/refind;
        A = toastBndReflectionTerm (refind, 'Keijzer');
        c2A = cm/(2.0*A);

        % Set up homogeneous initial parameter estimates
        mua = zeros(n,1);
        mus = ones(n,1)*1;
        ref = ones(n,1) * refind;

        for qw = 1:2
            % Generate source vectors
            qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', qwidth(qw));

            % Solve for fields
            smat = toastSysmat (hMesh, mua, mus, ref, freq);
            phi = full(smat\qvec);

            [vtx idx perm] = toastSurfData(hMesh);
            nbnd = length(perm);
            bphi = phi(perm,1);
            gamma = bphi * c2A;
            gsum = sum(gamma)/nbnd * 2*pi*rad;
    
            fprintf('---------------------------\n');
            fprintf('mesh size %d\n', n);
            fprintf('refind %f\n', refind);
            fprintf('source: Gaussian, sigma=%f\n', qwidth(qw));
            fprintf('==> total exitance = %f\n', gsum);
        end
    end
end
