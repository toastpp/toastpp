function f = toastGMRES (x, J, g, M, lambda, hReg, tol)

    n = length(M);
    RHess = toastRegulHess (hReg, x);
    dM = spdiags(M',0,n,n);
    RHess = dM*RHess*dM;
    
    f = gmres (@gmres_clbk, g, 30, tol, 100, @gmres_mfun);

    function clbk = gmres_clbk(r)
        clbk = J'*(J*r);
        clbk = clbk + RHess*r;
    end

    function mfun = gmres_mfun(r)
        mfun = r./M';
    end
end
