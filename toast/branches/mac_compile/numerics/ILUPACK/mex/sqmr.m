function [x,flag,iter,resvec]=sqmr(A,b,tol,maxit,M1,M2,x0)
% [x,flag,iter,resvec]=sqmr(A,b,tol,maxit,M1,M2,x0)
%
% sQMR   simplified Quasi-Minimal Residual Method for symmetric matrices A
%   symmetric preconditioning.
%   x = sqmr(A,b) attempts to solve the system of linear equations A*x=b for
%   x. The n-by-n coefficient matrix A must be square and real symmetric (resp.
%   complex Hermitian or complex symmetric) and the right hand
%   side column vector b must have length n.
%
%   x = sqmr(A,b,tol) specifies the tolerance of the method. If tol not
%   specified, then sqmr uses the default, 1e-6.
% 
%   x = sqmr(A,b,tol,maxit) specifies the maximum number of iterations. If
%   maxit is not specified then sqmr uses the default, min(N,500).
% 
%   x = sqmr(A,b,tol,maxit,M) and x = sqmr(A,b,tol,maxit,M1,M2) use
%   preconditioners M or M=M1*M2 and effectively solve the system
%   inv(M)*A*x = inv(M)*b for x. M may be a function handle mfun such that 
%   mfun(x) returns M\x.
% 
%   x = sqmr(A,b,tol,maxit,M1,M2,x0) and x=sqmr(A,b,tol,maxit,M,x0) specifies 
%   the initial guess. If x0 is not specified then sqmr uses the zero vector.
% 
%   [x,flag] = sqmr(A,b,...) also returns a convergence flag:
%     0 sqmr converged to the desired tolerance tol within maxit iterations.
%       As stopping criterion the backward error ||Ax-b||/(||A|| ||x||+||b||)
%       is used
%     1 sqmr iterated maxit times but did not converge.
%     2 break down
%    -m unknown error code
% 
%   [x,flag,iter] = sqmr(A,b,...) also returns the iteration number
%   at which x was computed: 0 <= iter <= maxit.
% 
%   [x,flag,iter,resvec] = sqmr(A,b,...) also returns a vector of the
%   estimated residual norms at each iteration, including norm(b-A*x0).

if isstr(A) | isa(A,'function_handle')
   error('matrix A must be given explicitly');
end

n=size(A,1);
if nargin==2
   tol=1e-6;
elseif nargin==3
   maxit=min(n,500);
end % if

if size(b,2)~=1
   error('right hand side must be a vector');
end

% find out the meaning of the sixth input parameter
m=0;
if nargin==6
   % how many columns do we have
   m=size(M2,2);
   % initial guess is provided
   if m==1
      x0=M2;
   % M2 is really passed
   elseif m==n
      x0=zeros(n,1);
   else 
      error('initial guess must be a vector');
   end
end

if nargin<6
   x0=zeros(n,1);
end

% typeres
%  1 ||Ax_k-b|| <= tol ||Ax_0-b||
%  2 ||Ax_k-b|| <= tol ||b||
%  3 ||Ax_k-b|| <= tol (||A||+||x_k||+||b||)
typeres=3;

nrmA=1;
if typeres==3
   nrmA=sqrt(norm(A,1)*norm(A,inf));
end

choice=0;
if isreal(A)
   if norm(A-A',1)==0
      choice=1;
   end
else
   if norm(A-A',1)==0
      choice=2;
   elseif norm(A-A.',1)==0
      choice=3;
   end
end
if choice==0
   error('matrix A must be symmetric');
end

% initial guess
x=x0;
% for reverse communication, initially we formally have to pass something
% when entering sqmr for the first time, this value is ignored
drain=x0;
% vector of estimated residuals
resvec=[];
% indicate that the iteration is still going on
flag=-1;
% control flag for reverse communication
control=0;
while flag

      if choice==1
	 % real symmetric case
	 [x,src,control,iter,relres]=DSYMilupacksqmr(b,x,drain,tol,maxit,control,nrmA,typeres);
      elseif choice==2
	 % complex Hermitian case
	 [x,src,control,iter,relres]=ZHERilupacksqmr(b,x,drain,tol,maxit,control,nrmA,typeres);
      else
	 % complex symmetric case
	 [x,src,control,iter,relres]=ZSYMilupacksqmr(b,x,drain,tol,maxit,control,nrmA,typeres);
      end

      % request for matrix-vector-multiplication
      if control==1
         % when sqmr stops for the first time, it computes the initial residual
         % only after that initial step we need to extract the residual 
         % estimate
	 if iter>0
	    resvec=[resvec, relres];
	 end
         drain=A*src; 
      % request for preconditioning
      elseif control==2
	 % initially M==I, no preconditioning
         drain=src;
         % check if M=M1 is passed
         % argument no. 6 is in real life x0 and M1 is M
         if nargin==5 | (nargin==6 & m==1)
            % M1 refers to 'M', single preconditioner
  	    if isstr(M1)
	       drain=feval(M1,src);
	    elseif isa(M1,'function_handle')
	       drain=M1(src);
   	    else
	       drain=M1\src;
            end
         elseif nargin>=6
	    drain=M2\(M1\src);
	 end % if
      else
         flag=0;
      end % if
end % while

flag=-control;


