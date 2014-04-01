function x = AMGsol(PREC,b)
% x = AMGsol(PREC,b)
%
% Solves Ax=b using one step of the computed ILUPACK preconditioner PREC
% 

   
% preconditioner is real
if PREC(1).isreal

   % matrix and preconditioner are symmetric
   if PREC(1).issymmetric

      % matrix and preconditioner are real SPD -> CG
      if PREC(1).isdefinite
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		x(:,i)=DSPDilupacksol(PREC, full(real(b(:,i))));
		y     =DSPDilupacksol(PREC, full(imag(b(:,i))));
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		x(:,i)=DSPDilupacksol(PREC, full(b(:,i)));
	     end % if-else
	 end % for i
      
      else % symmetric indefinite preconditioner
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		x(:,i)=DSYMilupacksol(PREC, full(real(b(:,i))));
		y     =DSYMilupacksol(PREC, full(imag(b(:,i))));
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		x(:,i)=DSYMilupacksol(PREC, full(b(:,i)));
	     end % if-else
	 end % for i
      end % if
      
   else % unsymmetric PREC
      for i=1:size(b,2)
	  if ~isreal(b(:,i))
	     x(:,i)=DGNLilupacksol(PREC, full(real(b(:,i))));
	     y     =DGNLilupacksol(PREC, full(imag(b(:,i))));
	     x(:,i)=x(:,i)+sqrt(-1)*y;
	  else
	     x(:,i)=DGNLilupacksol(PREC, full(b(:,i)));
	  end % if-else
      end % for i
   end % if
      
% complex preconditioner
else
   
   % preconditioner is Hermitian
   if PREC(1).ishermitian

      % preconditioner is HPD
      if PREC(1).isdefinite
	 for i=1:size(b,2)
	     x(:,i)=ZHPDilupacksol(PREC, full(b(:,i)));
	 end % for i
	 
      % only Hermitian preconditioner
      else
	 for i=1:size(b,2)
	     x(:,i)=ZHERilupacksol(PREC, full(b(:,i)));
	 end % for i
      end % if
      
   % complex symmetric preconditioner
   elseif PREC(1).issymmetric
      for i=1:size(b,2)
	  x(:,i)=ZSYMilupacksol(PREC, full(b(:,i)));
      end % for i
      
   else % general unsymmetric preconditioner
      for i=1:size(b,2)
	  x(:,i)=ZGNLilupacksol(PREC, full(b(:,i)));
      end % for i
   end  % if
end  % if
