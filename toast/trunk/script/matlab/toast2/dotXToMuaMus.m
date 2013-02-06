function [mua,mus] = dotXToMuaMus (hBasis, x, ref)

%x = exp(logx);

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
    mua = hBasis.Map ('S->M', smua);
    mus = hBasis.Map ('S->M', smus);
    %cmua = toastMapSolToMesh (hBasis, scmua);
    %ckap = toastMapSolToMesh (hBasis, sckap);
    %mua = cmua./c;
    %mus = c./(3*ckap) - mua;
end
