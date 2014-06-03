function [A,rhs,rhstyp]=loadhbo(filename)
% [A,rhs,rhstyp]=loadhbo(filename)
% 
% load matrix A and optionally right hand side b, initial guess x0 and
% exact solution x
%
% Input
% -----
% filename   name of the file without extension
%
% Output
% ------
% A         mxn sparse matrix
% rhs       vector of right hand side(s), initial guess(es), exact solution(s)
%           depending on `rhstyp'
% rhstyp    rhstyp(1)=='F' or rhstyp=='f'
%                 => dense right hand side(s)
%           rhstyp(2)=='G' or rhstyp=='g'
%                 => dense initial guess(es)
%           rhstyp(2)=='X' or rhstyp=='g'
%                 => dense exact solution(s)
%           depending on rhstyp it will be clear how many columns of `rhs'
%           correspond to the hand side(s), the initial guess(es) or the exact
%           solution(s)

if (nargin~=1)
   error('only one input parameters required!');
end

if nargout~=3
   error('three output parameters requried');
end

fp=fopen([filename '.rsa']);
if fp~=-1
   fclose(fp);
   [A,rhs,rhstyp]=Dloadhbo([filename '.rsa']);
   [m,n]=size(A);
   A=A-spdiags(diag(A),0,m,n)+A';
else
   fp=fopen([filename '.rua']);
   if fp~=-1
      fclose(fp);
      [A,rhs,rhstyp]=Dloadhbo([filename '.rua']);
   else
      fp=fopen([filename '.cha']);
      if fp~=-1
	 fclose(fp);
	 [A,rhs,rhstyp]=Zloadhbo([filename '.cha']);
	 [m,n]=size(A);
	 A=A-spdiags(diag(A),0,m,n)+A';
      else
	 fp=fopen([filename '.csa']);
	 if fp~=-1
	    fclose(fp);
	    [A,rhs,rhstyp]=Zloadhbo([filename '.csa']);
	    [m,n]=size(A);
	    A=A-spdiags(diag(A),0,m,n)+A.';
	 else
	    fp=fopen([filename '.cua']);
	    if fp~=-1
	       fclose(fp);
	       [A,rhs,rhstyp]=Zloadhbo([filename '.cua']);
	    else
	       error('file not found');
	    end % if-else
	 end % if-else
      end % if-else
   end % if-else
end % if-else
