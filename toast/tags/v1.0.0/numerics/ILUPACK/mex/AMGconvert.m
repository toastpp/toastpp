function PREC = AMGconvert(PREC)
%
% If the preconditioner is real symmetric and indefinite or complex Hermitian
% and indefinite, then this routine turns PREC into a positive definite preconditioner


if nargin ~=1
   error('wrong number of input arguments');
end % if 

if nargout~=1
   error('wrong number of output arguments');
end % if 

  

   
if PREC(1).isreal & PREC(1).issymmetric & ~PREC(1).isdefinite
   PREC=DSYMSPDilupackconvert(PREC);
   nlev=length(PREC);
   for i=1:nlev
       if issparse(PREC(i).L)
	  flag=-1;
	  % separate diagonal part
	  PREC(i).absdiag=sparse(size(PREC(i).L,1),size(PREC(i).L,2));
       else
	  flag=0;
       end
       j=1;
       while j<=size(PREC(i).D,1)
	     if j<size(PREC(i).D,1)
		if PREC(i).D(j+1,j)==0
		   if flag
		      PREC(i).absdiag(j,j)=abs(PREC(i).D(j,j));
		   else
		      PREC(i).D(j,j)=abs(PREC(i).D(j,j));
		   end % if
		   j=j+1;
		else
		   [V,D]=eig(full(PREC(i).D(j:j+1,j:j+1)));
		   if flag
		      PREC(i).absdiag(j:j+1,j:j+1)=V*abs(D)*V';
		   else
		      PREC(i).D(j:j+1,j:j+1)=V*abs(D)*V';
		   end % if
		   j=j+2;
		end % if-else
	     else
		if flag
		   PREC(i).absdiag(j,j)=abs(PREC(i).D(j,j));
		else
		   PREC(i).D(j,j)=abs(PREC(i).D(j,j));
		end % if
		j=j+1;
	     end
       end % while
   end % for i
elseif ~PREC(1).isreal & PREC(1).ishermitian & ~PREC(1).isdefinite
   PREC=ZHERHPDilupackconvert(PREC);
   nlev=length(PREC);
   for i=1:nlev
       if issparse(PREC(i).L)
	  flag=-1;
	  % separate diagonal part
	  PREC(i).absdiag=sparse(size(PREC(i).L,1),size(PREC(i).L,2));
       else
	  flag=0;
       end
       j=1;
       while j<=size(PREC(i).D,1)
	     if j<size(PREC(i).D,1)
		if PREC(i).D(j+1,j)==0
		   if flag
		      PREC(i).absdiag(j,j)=abs(PREC(i).D(j,j));
		   else
		      PREC(i).D(j,j)=abs(PREC(i).D(j,j));
		   end % if
		   j=j+1;
		else
		   [V,D]=eig(full(PREC(i).D(j:j+1,j:j+1)));
		   if flag
		      PREC(i).absdiag(j:j+1,j:j+1)=V*abs(D)*V';
		   else
		      PREC(i).D(j:j+1,j:j+1)=V*abs(D)*V';
		   end % if
		   j=j+2;
		end % if-else
	     else
		if flag
		   PREC(i).absdiag(j,j)=abs(PREC(i).D(j,j));
		else
		   PREC(i).D(j,j)=abs(PREC(i).D(j,j));
		end % if
		j=j+1;
	     end
       end % while
   end % for i
end % if
