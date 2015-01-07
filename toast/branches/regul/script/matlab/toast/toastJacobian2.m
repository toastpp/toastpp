function j2 = toastJacobian2 (hMesh, hBasis, qvec, mvec, mua, mus, ref, freq, solver, tol)
                               
     [dphi aphi] = toastFields (hMesh, hBasis, qvec, mvec, mua, mus, ref, freq, solver, tol);
     proj = toastProject (hMesh, mua, mus, ref, freq, qvec, mvec);
     
     for i=1:size(dphi,2)
         Gdphi = toastImageGradient (hBasis,dphi(:,i));
         for j=1:size(aphi,2)
             Gaphi = toastImageGradient(hBasis,aphi(:,j));
             pmdf_mua = -dphi(i).*aphi(j);

