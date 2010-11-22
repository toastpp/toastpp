function proj = privProject (hMesh, hBasis, logx, ref, freq, qvec, mvec, solver, tol)

% Subsititute missing parameters with defaults
if nargin < 10
    tol = 1e-10;
    if nargin < 9
        solver = 'DIRECT';
    end
end

x = exp(logx);

if hBasis == 0
    c = 0.3./ref; % speed of light
    cmua = x(1:size(x)/2);
    ckap = x(size(x)/2+1:size(x));
    mua = cmua./c;
    mus = c./(3*ckap) - mua;
else
    %c = 0.3./ref; % speed of light
    cm = 0.3/ref(1);
    scmua = x(1:size(x)/2);
    sckap = x(size(x)/2+1:size(x));
    smua = scmua/cm;
    skap = sckap/cm;
    smus = 1./(3*skap) - smua;
    mua = toastMapSolToMesh (hBasis, smua);
    mus = toastMapSolToMesh (hBasis, smus);
    %cmua = toastMapSolToMesh (hBasis, scmua);
    %ckap = toastMapSolToMesh (hBasis, sckap);
    %mua = cmua./c;
    %mus = c./(3*ckap) - mua;
end

proj = toastProject (hMesh, mua, mus, ref, freq, qvec, mvec, solver, tol);
