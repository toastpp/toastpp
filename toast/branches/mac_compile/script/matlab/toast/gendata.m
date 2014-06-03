global SMAT
global QVEC
global MVEC
fwdrun
phi = SMAT\QVEC;
gamma = phi'*MVEC;
lgamma = log(gamma);
lnmod = real(lgamma);
phase = imag(lgamma);
figure
imagesc(lnmod)
figure
imagesc(phase)
