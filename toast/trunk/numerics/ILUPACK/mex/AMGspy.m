function AMGspy(A,PREC)
% AMGspy(PREC)
%
% display multilevel preconditioner PREC
%
% AMGspy(A,PREC)
%
% display remapped original matrix A associated with the sequence of 
% reorderings given by PREC


if nargin==1
   PREC=A;

   n=PREC(1).n;

   A=sparse(n,n);
   spy(A);
   hold on

   nlev=length(PREC);
   sumnB=0;
   for lev=1:nlev
      
       nB=PREC(lev).nB;
       [mL,nL]=size(PREC(lev).L);
       
       if lev<nlev
	  if (mL==nL)  
 	     A=[sparse(sumnB,n);...
	        sparse(nB,sumnB) PREC(lev).L sparse(nB,n-sumnB-nB);...
	        sparse(n-sumnB-nB,n)]; 
	     spy(A,'g');pause(0.001);
	  
	     A=[sparse(sumnB+nB,n);...
	        sparse(n-sumnB-nB,sumnB) PREC(lev).E sparse(n-sumnB-nB,n-sumnB-nB)]; 
	     spy(A,'r');pause(0.001);
	     if PREC(1).issymmetric | PREC(1).ishermitian
	        A=[sparse(sumnB,n);...
		   sparse(nB,sumnB) PREC(lev).L' sparse(nB,n-sumnB-nB);...
		   sparse(n-sumnB-nB,n)]; 
	        spy(A,'b');pause(0.001);
	     
	        A=[sparse(sumnB,n);...
		   sparse(nB,sumnB+nB) PREC(lev).E';...
		   sparse(n-sumnB-nB,n)]; 
	        spy(A,'r');pause(0.001);
	     else
	        A=[sparse(sumnB,n);...
		   sparse(nB,sumnB) PREC(lev).U sparse(nB,n-sumnB-nB);...
		   sparse(n-sumnB-nB,n)]; 
	        spy(A,'b');pause(0.001);
	     
	        A=[sparse(sumnB,n);...
		   sparse(nB,sumnB+nB) PREC(lev).E';...
		   sparse(n-sumnB-nB,n)]; 
	       spy(A,'r');pause(0.001);
	     end
	  
	     A=[sparse(sumnB,n);...
	        sparse(nB,sumnB) PREC(lev).D sparse(nB,n-sumnB-nB);...
	        sparse(n-sumnB-nB,n)]; 
	     spy(A,'k');pause(0.001);
	  else
 	     A=[sparse(sumnB,n);...
	        sparse(n-sumnB,sumnB) PREC(lev).L sparse(n-sumnB,n-sumnB-nB)]; 
	     spy(A,'g');pause(0.001);
	  
	     if PREC(1).issymmetric | PREC(1).ishermitian
	        A=[sparse(sumnB,n);...
		   sparse(nB,sumnB) PREC(lev).L';...
		   sparse(n-sumnB-nB,n)]; 
	        spy(A,'b');pause(0.001);
	     else
	        A=[sparse(sumnB,n);...
		   sparse(nB,sumnB) PREC(lev).U;...
		   sparse(n-sumnB-nB,n)]; 
	        spy(A,'b');pause(0.001);
	     end
	  
	     A=[sparse(sumnB,n);...
	        sparse(nB,sumnB) PREC(lev).D sparse(nB,n-sumnB-nB);...
	        sparse(n-sumnB-nB,n)]; 
	     spy(A,'k');pause(0.001);
	  end
       else
	  if issparse(PREC(lev).L)
	     A=[sparse(sumnB,n);...
		sparse(nB,sumnB) PREC(lev).L sparse(nB,n-sumnB-nB);...
		sparse(n-sumnB-nB,n)]; 
	     spy(A,'g');pause(0.001);
	     
	     if PREC(1).issymmetric | PREC(1).ishermitian
		A=[sparse(sumnB,n);...
		   sparse(nB,sumnB) PREC(lev).L' sparse(nB,n-sumnB-nB);...
		   sparse(n-sumnB-nB,n)]; 
		spy(A,'b');pause(0.001);
	     else
		A=[sparse(sumnB,n);...
		   sparse(nB,sumnB) PREC(lev).U sparse(nB,n-sumnB-nB);...
		   sparse(n-sumnB-nB,n)]; 
		spy(A,'b');pause(0.001);
	     end
	     
	     A=[sparse(sumnB,n);...
		sparse(nB,sumnB) PREC(lev).D sparse(nB,n-sumnB-nB);...
		sparse(n-sumnB-nB,n)]; 
	     spy(A,'k');pause(0.001);
	  else
	     A=[sparse(sumnB,n);...
		sparse(nB,sumnB) tril(ones(nB,nB)) sparse(nB,n-sumnB-nB);...
		sparse(n-sumnB-nB,n)]; 
	     spy(A,'g');pause(0.001);
	     
	     A=[sparse(sumnB,n);...
		sparse(nB,sumnB) triu(ones(nB,nB)) sparse(nB,n-sumnB-nB);...
		sparse(n-sumnB-nB,n)]; 
	     spy(A,'b');pause(0.001);
	     
	     A=[sparse(sumnB,n);...
		sparse(nB,sumnB) eye(nB) sparse(nB,n-sumnB-nB);...
		sparse(n-sumnB-nB,n)]; 
	     spy(A,'k');pause(0.001);
	  end
       end
       
       sumnB=sumnB+nB;
       
       clear A;
    end % for lev
    
    sumnB=0;
    for lev=1:nlev
       nB=PREC(lev).nB;
       if lev<nlev
	  plot([sumnB+1 n],         [sumnB+nB sumnB+nB], '-k');
	  plot([sumnB+nB sumnB+nB], [sumnB+1 n],         '-k');
	  pause(0.001);
       end
       
       sumnB=sumnB+nB;
    end % for lev
 
    
    title(['ILUPACK multilevel preconditioner (' num2str(nlev) ' levels)'])
    xlabel(['nz=' num2str(AMGnnz(PREC))])
    hold off;
else % two input arguments, display A
  
   nz=nnz(A);
   n=PREC(1).n;
   p=1:n;
   q=1:n;
   
   
   nlev=length(PREC);
   sumnB=0;
   for lev=1:nlev
      
       nB=PREC(lev).nB;
       p   =PREC(lev).p;
       invq=PREC(lev).invq;
       q(invq)=1:(n-sumnB);
       q=q(1:n-sumnB);
       
       spy([sparse(sumnB,n);...
	    sparse(n-sumnB,sumnB) A(p,q)],'b'); pause(0.001)
       A=A(p(nB+1:n-sumnB),q(nB+1:n-sumnB));
       hold on;
       spy([sparse(sumnB+nB,n);...
	    sparse(n-sumnB-nB,sumnB+nB) A],'w');
       pause(0.001)
       
       % p(sumnB+1:n)=p(sumnB+PREC(lev).p);
       % q(sumnB+PREC(lev).invq)=q(sumnB+1:n);
       
       sumnB=sumnB+nB;
   end % for lev
  
   
   sumnB=0;
   for lev=1:nlev
       nB=PREC(lev).nB;
       if lev<nlev
	  plot([sumnB+1 n],         [sumnB+nB sumnB+nB], '-k');
	  plot([sumnB+nB sumnB+nB], [sumnB+1 n],         '-k');
	  pause(0.001);
       end
       
       sumnB=sumnB+nB;
   end % for lev
   
   title(['ILUPACK multilevel reordering (' num2str(nlev) ' levels)']);
   xlabel(['nz=' num2str(nz)]);
   hold off;
end % if
