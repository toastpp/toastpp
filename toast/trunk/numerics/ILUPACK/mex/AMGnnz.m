function nz=AMGnnz(PREC)
% nz=AMGnnz(PREC)
%
% total number of nonzeros of the multilevel ILU
% to be consistent with MATLAB, the preconditioner is
% treated as if it were unsymmetric

nz=0;
nlev=length(PREC);
for lev=1:nlev

    if lev<nlev
       if PREC(1).issymmetric | PREC(1).ishermitian
	  nnzU=nnz(PREC(lev).L);
	  nnzF=nnz(PREC(lev).E);
       else
	  nnzU=nnz(PREC(lev).U);
	  nnzF=nnz(PREC(lev).F);
       end
       nz=nz+nnz(PREC(lev).L)+nnzU-nnz(PREC(lev).D)...
	    +nnz(PREC(lev).E)+nnzF;
       nz=nz+nnz(PREC(lev).A_H);
    else
       if issparse(PREC(lev).L)
	  if PREC(1).issymmetric | PREC(1).ishermitian
	     nnzU=nnz(PREC(lev).L);
	  else
	     nnzU=nnz(PREC(lev).U);
	  end
	  nz=nz+nnz(PREC(lev).L)+nnzU-nnz(PREC(lev).D);
       else
	  nz=nz+PREC(lev).n*PREC(lev).n;
       end
    end
       
end % for lev