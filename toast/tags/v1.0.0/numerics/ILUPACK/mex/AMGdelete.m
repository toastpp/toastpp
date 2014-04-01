function PREC=AMGdelete(PREC)
% PREC=AMGdelete(PREC)
% 
% delete preconditioner, in particular release the associated memory


if nargout~=1
   error('Output argument must be the same as the input argument');
end

if PREC(1).isreal
   if PREC(1).issymmetric
      if PREC(1).isdefinite
	 DSPDilupackdelete(PREC);      
      else
	 DSYMilupackdelete(PREC);      
      end % if-else
   else
      DGNLilupackdelete(PREC);      
   end
else
   if PREC(1).ishermitian
      if PREC(1).isdefinite
	 ZHPDilupackdelete(PREC);      
      else
	 ZHERilupackdelete(PREC);      
      end % if-else
   elseif PREC(1).issymmetric
      ZSYMilupackdelete(PREC);      
   else
      ZGNLilupackdelete(PREC);      
   end
end % if
clear PREC;
PREC=[];
