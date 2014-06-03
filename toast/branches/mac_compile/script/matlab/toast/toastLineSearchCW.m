function [alpha, fmin] = toastLineSearchCW (hMesh, hBasis, hReg, p, d, xb, f0, ckap, ref, cm, qvec, mvec, data, sd)
%toastLineSearchCW    - 1-D minimisation along a given search direction.
%
% Synopsis: [alpha,fmin] = toastLineSearchCW(hMesh,hBasis,hReg,p,d,xb,f0,ckap,
%                          ref,cm,qvec,mvec,data,sd)
%    hMesh:  mesh handle
%    hBasis: handle for basis mapper
%    hReg:   regularisation handle
%    p:      current parameter array
%    d:      search direction
%    xb:     initial step size
%    f0:     initial objective value
%    ckap:   c*kappa array (speed of light times diffusion coefficient)
%    ref:    refractive index array
%    cm:     speed of light array
%    qvec:   source vectors (sparse matrix of column vectors)
%    mvec:   measurement vector (sparse matrix of column vectors)
%    data:   measurement data array
%    sd:     measurement standard deviation array
%    alpha:  final step size
%    fmin:   final objective value
%
% Continuous-wave version of 1-D line minimisation algorithm.
% Performs an approximate minimisation of the objective function Q(p+alpha*d),
% where Q(x) = ||(f(x)-y)/sd||^2, f(x) is the forward model generating boundary
% projections from parameters x, y are measurements with standard deviations
% sd, p is the initial parameter set, d, is the search direction, and alpha
% is the step size being optimised.
% This method first brackets the minimum, and then performs a quadratic
% interpolation step. d is required to be a downhill direction.
% The ckap vector is required to provide the diffusion coefficient distribution
% for the forward model, since this is not provided in the data vectors, in
% contrast to the general toastLineSearch function.

lnkap = log(ckap);
nm = length(data);
x0 = 0;
p1 = p + d*xb;
proj = privProject (hMesh, hBasis, [p1;lnkap], ref, cm, 0, qvec, mvec);
[fb fbd fbp] = privObjective (proj(1:nm), data, sd, hReg, p1);
fprintf (1, 'Step: %f, Error: %f (data), %f (prior), %f (total)\n', xb, fbd, fbp, fb);

% bracket the minimum

if fb < f0
    xm = xb;
    fm = fb;
    xb = xb*2;
    p1 = p + d*xb;
    proj = privProject (hMesh, hBasis, [p1;lnkap], ref, cm, 0, qvec, mvec);
    [fb fbd fbp] = privObjective (proj(1:nm), data, sd, hReg, p1);
    fprintf (1, 'Step: %f, Error: %f (data), %f (prior), %f (total)\n', xb, fbd, fbp, fb);

    while fb < fm
        x0 = xm; f0 = fm;
        xm = xb; fm = fb;
        xb = xb*2;
        p1 = p + d*xb;
        proj = privProject (hMesh, hBasis, [p1;lnkap], ref, cm, 0, qvec, mvec);
        [fb fbd fbp] = privObjective (proj(1:nm), data, sd, hReg, p1);
        fprintf (1, 'Step: %f, Error: %f (data), %f (prior), %f (total)\n', xb, fbd, fbp, fb);
    end
else
    xm = xb/2;
    p1 = p + d*xm;
    proj = privProject (hMesh, hBasis, [p1;lnkap], ref, cm, 0, qvec, mvec);
    [fm fmd fmp] = privObjective (proj(1:nm), data, sd, hReg, p1);
    fprintf (1, 'Step: %f, Error: %f (data), %f (prior), %f (total)\n', xb, fmd, fmp, fm);

    while fm > f0
        xb = xm; fb = fm;
        xm = xb/2;
        p1 = p + d*xm;
        proj = privProject (hMesh, hBasis, [p1;lnkap], ref, cm, 0, qvec, mvec);
        [fm fmd fmp] = privObjective (proj(1:nm), data, sd, hReg, p1);
        fprintf (1, 'Step: %f, Error: %f (data), %f (prior), %f (total)\n', xb, fmd, fmp, fm);
    end
end

% quadratic interpolation

a = ((f0-fb)/(x0-xb) - (f0-fm)/(x0-xm)) / (xb-xm);
b = (f0-fb)/(x0-xb) - a*(x0+xb);
alpha = -b/(2*a);
p1 = p + d*alpha;
proj = privProject (hMesh, hBasis, [p1;lnkap], ref, cm, 0, qvec, mvec);
[fmin fmind fminp] = privObjective (proj(1:nm), data, sd, hReg, p1);
if fmin > fm    % no improvement
    alpha = xm;
    fmin = fm;
end
fprintf (1, 'Step: %f, Error: %f (data), %f (prior), %f (total) [FINAL]\n', alpha, fmind, fminp, fmin);

