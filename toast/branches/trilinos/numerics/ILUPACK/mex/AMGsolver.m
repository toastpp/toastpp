function [x, options] = AMGsolver(A, PREC, options, b, x0)
% [x, options] = AMGsolver(A, PREC, options, b, x0)
% [x, options] = AMGsolver(A, PREC, options, b)
% [x, options] = AMGsolver(A, PREC, b, x0)
% [x, options] = AMGsolver(A, PREC, b)
%
% Solves Ax=b using ILUPACK preconditioner PREC according to the given options
% 


if nargin==3
   % [x, options] = AMGsolver(A, PREC, b)
   % shift parameter
   b=options;
   options=AMGinit(A);
elseif nargin==4 || nargin==5
   if ~isstruct(options)
      % [x, options] = AMGsolver(A, PREC, b, x0)
      % shift parameters
      x0=b;
      b=options;
      options = AMGinit(A);
   else
      % [x, options] = AMGsolver(A, PREC, options, b)
      % [x, options] = AMGsolver(A, PREC, options, b, x0)
      if isfield(options, 'isdefinite')
	 myoptions.isdefinite=options.isdefinite;
	 myoptions = AMGinit(A,myoptions);
      else
	 myoptions = AMGinit(A);
      end
   
      % complete missing data
      % [x, options] = AMGsolver(A, PREC, options, b)
      if nargin==4
	 n=size(b,1);
	 m=size(b,2);
	 if issparse(b)
	    x0=sparse(n,m);
	 else
	    x0=zeros(n,m);
	 end % if
      end

      if ~isfield(options, 'restol')
	 options.restol=myoptions.restol;
      end
      if ~isfield(options, 'maxit')
	 options.maxit=myoptions.maxit;
      end
      if ~isfield(myoptions, 'nrestart')
	 options.nrestart=myoptions.nrestart;
      end
      
      if ~isfield(options,'amg')
	 options.amg=myoptions.amg;
      end
      if ~isfield(options,'npresmoothing')
	 options.npresmoothing=myoptions.npresmoothing;
      end
      if ~isfield(options,'npostsmoothing')
	 options.npostsmoothing=myoptions.npostsmoothing;
      end
      if ~isfield(options,'ncoarse')
	 options.ncoarse=myoptions.ncoarse;
      end
      if ~isfield(options,'presmoother')
	 options.presmoother=myoptions.presmoother;
      end
      if ~isfield(options,'postsmoother')
	 options.postsmoother=myoptions.postsmoother;
      end
      if ~isfield(options, 'solver')
	 options.solver=myoptions.solver;
      end
      if ~isfield(options, 'damping')
	 options.damping=myoptions.damping;
      end
      if ~isfield(options, 'mixedprecision')
	 options.damping=myoptions.mixedprecision;
      end
   end
else
   error('wrong number of input arguments');
end % if 

if nargout~=2
   error('wrong number of output arguments');
end % if 

  
myoptions=options;

if isreal(A)
   myoptions.isreal=1;
   if norm(A-A',1)==0
      myoptions.issymmetric=1;
      myoptions.ishermitian=1;
   else
      myoptions.issymmetric=0;
      myoptions.ishermitian=0;
   end
else
   myoptions.isreal=0;
   if norm(A-A',1)==0
      myoptions.ishermitian=1;
   else
      myoptions.ishermitian=0;
   end
   if ~myoptions.ishermitian
      if norm(A-A.',1)==0
	 myoptions.issymmetric=1;
      else
	 myoptions.issymmetric=0;
      end
   else
      myoptions.issymmetric=0;
   end
end


if ~isfield(options,'isdefinite')
   myoptions.isdefinite=0;
end
   
% matrix and preconditioner are real
if myoptions.isreal & PREC(1).isreal

   % matrix and preconditioner are symmetric
   if myoptions.issymmetric & PREC(1).issymmetric

      % matrix and preconditioner are real SPD -> CG
      if myoptions.isdefinite & PREC(1).isdefinite
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		myx0=full(real(x0(:,i)));
		[x(:,i), myoptions]=DSPDilupacksolver(A,PREC,myoptions,...
		                                      full(real(b(:,i))),myx0);
		options.niter(i)=myoptions.niter;
		myx0=full(imag(x0(:,i)));
		[y,      myoptions]=DSPDilupacksolver(A,PREC,myoptions,...
		                                      full(imag(b(:,i))),myx0);
		options.niter(i)=max(options.niter(i),myoptions.niter);
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		myx0=full(x0(:,i));
		[x(:,i), myoptions]=DSPDilupacksolver(A,PREC,myoptions,...
		                                      full(b(:,i)),myx0);
		options.niter(i)=myoptions.niter;
	     end % if-else
	 end % for i
      
      % definite preconditioner but not necessarily definite matrix -> MINRES
      % <=> symmetric QMR
      elseif ~myoptions.isdefinite & PREC(1).isdefinite
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		myx0=full(real(x0(:,i)));
		[x(:,i), myoptions]=DSYMSPDilupacksolver(A,PREC,myoptions,...
		                                         full(real(b(:,i))),myx0);
		options.niter(i)=myoptions.niter;
		myx0=full(imag(x0(:,i)));
		[y,      myoptions]=DSYMSPDilupacksolver(A,PREC,myoptions,...
		                                         full(imag(b(:,i))),myx0);
		options.niter(i)=max(options.niter(i),myoptions.niter);
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		myx0=full(x0(:,i));
		[x(:,i), myoptions]=DSYMSPDilupacksolver(A,PREC,myoptions,...
		                                         full(b(:,i)),myx0);
		options.niter(i)=myoptions.niter;
	     end % if-else
	 end % for i
      
      % symmetry but no definiteness -> SQMR
      else
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		myx0=full(real(x0(:,i)));
		[x(:,i), myoptions]=DSYMilupacksolver(A,PREC,myoptions,...
		                                      full(real(b(:,i))),myx0);
		options.niter(i)=myoptions.niter;
		myx0=full(imag(x0(:,i)));
		[y,      myoptions]=DSYMilupacksolver(A,PREC,myoptions,...
		                                      full(imag(b(:,i))),myx0);
		options.niter(i)=max(options.niter(i),myoptions.niter);
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		myx0=full(x0(:,i));
		[x(:,i), myoptions]=DSYMilupacksolver(A,PREC,myoptions,...
		                                      full(b(:,i)),myx0);
		options.niter(i)=myoptions.niter;
	     end % if-else
	 end % for i
      end % if
      
   % matrix is unsymmetric but preconditioner is symmetric
   elseif ~myoptions.issymmetric & PREC(1).issymmetric

      % preconditioner is SPD -> GMRES
      if PREC(1).isdefinite
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		myx0=full(real(x0(:,i)));
		[x(:,i), myoptions]=DGNLSPDilupacksolver(A,PREC,myoptions,...
		                                         full(real(b(:,i))),myx0);
		options.niter(i)=myoptions.niter;
		myx0=full(imag(x0(:,i)));
		[y,      myoptions]=DGNLSPDilupacksolver(A,PREC,myoptions,...
		                                         full(imag(b(:,i))),myx0);
		options.niter(i)=max(options.niter(i),myoptions.niter);
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		myx0=full(x0(:,i));
		[x(:,i), myoptions]=DGNLSPDilupacksolver(A,PREC,myoptions,...
		                                         full(b(:,i)),myx0);
		options.niter(i)=myoptions.niter;
	     end % if-else
	 end % for i

      % symmetric preconditioner -> GMRES
      else
	 for i=1:size(b,2)
	     if ~isreal(b(:,i))
		myx0=full(real(x0(:,i)));
		[x(:,i), myoptions]=DGNLSYMilupacksolver(A,PREC,myoptions,...
		                                         full(real(b(:,i))),myx0);
		options.niter(i)=myoptions.niter;
		myx0=full(imag(x0(:,i)));
		[y,      myoptions]=DGNLSYMilupacksolver(A,PREC,myoptions,...
		                                         full(imag(b(:,i))),myx0);
		options.niter(i)=max(options.niter(i),myoptions.niter);
		x(:,i)=x(:,i)+sqrt(-1)*y;
	     else
		myx0=full(x0(:,i));
		[x(:,i), myoptions]=DGNLSYMilupacksolver(A,PREC,myoptions,...
		                                         full(b(:,i)),myx0);
		options.niter(i)=myoptions.niter;
	     end % if-else
	 end % for i
      end % if
      
   % in any case the preconditioner is unsymmetric -> GMRES	 
   else
      for i=1:size(b,2)
	  if ~isreal(b(:,i))
	     myx0=full(real(x0(:,i)));
	     [x(:,i), myoptions]=DGNLilupacksolver(A,PREC,myoptions,...
		                                   full(real(b(:,i))),myx0);
	     options.niter(i)=myoptions.niter;
	     myx0=full(imag(x0(:,i)));
	     [y,      myoptions]=DGNLilupacksolver(A,PREC,myoptions,...
		                                   full(imag(b(:,i))),myx0);
	     options.niter(i)=max(options.niter(i),myoptions.niter);
	     x(:,i)=x(:,i)+sqrt(-1)*y;
	  else
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=DGNLilupacksolver(A,PREC,myoptions,...
		                                   full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	  end % if-else
      end %  ifor
   end % if-elseif-else myoptions.issymmetric & PREC(1).issymmetric

% complex matrix but real preconditioner
elseif ~myoptions.isreal & PREC(1).isreal

   % matrix and preconditioner are Hermitian
   if myoptions.ishermitian & (PREC(1).ishermitian|PREC(1).issymmetric)

      % matrix and preconditioner are HPD -> CG
      if myoptions.isdefinite & PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZHPDDSPDilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      
      % definite preconditioner but not necessarily definite matrix -> MINRES
      elseif ~myoptions.isdefinite & PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZHERDSPDilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      
      % symmetry but no definiteness -> SQMR
      else
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZHERDSYMilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      end % if
      
   % matrix is complex symmetric and preconditioner is real symmetric
   elseif myoptions.issymmetric & (PREC(1).issymmetric|PREC(1).ishermitian)
      % preconditioner is real SPD -> SQMR
      if PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZSYMDSPDilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i

      % real symmetric preconditioner -> SQMR
      else
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZSYMDSYMilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      end % if-else
      
   % matrix is non-Hermitian,non-symmetric but preconditioner is real symmetric
   elseif ~myoptions.ishermitian & (PREC(1).ishermitian|PREC(1).issymmetric)

      % preconditioner is real SPD -> GMRES
      if PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZGNLDSPDilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i

      % real symmetric preconditioner -> GMRES
      else
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZGNLDSYMilupacksolver(A,PREC,myoptions,...
		                                       full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      end % if

   % matrix is non-complex symmetric but preconditioner is real symmetric
   elseif ~myoptions.issymmetric & (PREC(1).issymmetric|PREC(1).ishermitian)
      % complex symmetric preconditioner -> GMRES
      for i=1:size(b,2)
	  myx0=full(x0(:,i));
	  [x(:,i), myoptions]=ZGNLDSYMilupacksolver(A,PREC,myoptions,...
	                                            full(b(:,i)),myx0);
	  options.niter(i)=myoptions.niter;
      end % for i
       
   % in any case the preconditioner is unsymmetric -> GMRES	 
   else
      for i=1:size(b,2)
	  myx0=full(x0(:,i));
	  [x(:,i), myoptions]=ZGNLDGNLilupacksolver(A,PREC,myoptions,...
	                                            full(b(:,i)),myx0);
	  options.niter(i)=myoptions.niter;
      end % for i
   end % if-elseif-else myoptions.ishermitian & (PREC(1).ishermitian|PREC(1).issymmetric)
   
% complex preconditioner, treat matrix also as complex
else
   
   % matrix and preconditioner are Hermitian
   if myoptions.ishermitian & PREC(1).ishermitian

      % matrix and preconditioner are HPD -> CG
      if myoptions.isdefinite & PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZHPDilupacksolver(A,PREC,myoptions,...
		                                   full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      
      % definite preconditioner but not necessarily definite matrix -> MINRES
      elseif ~myoptions.isdefinite & PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZHERHPDilupacksolver(A,PREC,myoptions,...
		                                      full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      
      % symmetry but no definiteness -> SQMR
      else
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZHERilupacksolver(A,PREC,myoptions,...
		                                   full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      end % if
      
   % matrix is non-Hermitian but preconditioner is Hermitian
   elseif ~myoptions.ishermitian & PREC(1).ishermitian

      % preconditioner is HPD -> GMRES
      if PREC(1).isdefinite
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZGNLHPDilupacksolver(A,PREC,myoptions,...
		                                      full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i

      % Hermitian preconditioner -> GMRES
      else
	 for i=1:size(b,2)
	     myx0=full(x0(:,i));
	     [x(:,i), myoptions]=ZGNLHERilupacksolver(A,PREC,myoptions,...
		                                      full(b(:,i)),myx0);
	     options.niter(i)=myoptions.niter;
	 end % for i
      end % if

   % matrix and preconditioner are complex symmetric
   elseif myoptions.issymmetric & PREC(1).issymmetric
      for i=1:size(b,2)
	  myx0=full(x0(:,i));
	  [x(:,i), myoptions]=ZSYMilupacksolver(A,PREC,myoptions,...
	                                        full(b(:,i)),myx0);
	  options.niter(i)=myoptions.niter;
      end % for i
      
   % matrix is non-complex symmetric but preconditioner is complex symmetric
   elseif ~myoptions.issymmetric & PREC(1).issymmetric
      % complex symmetric preconditioner -> GMRES
      for i=1:size(b,2)
	  myx0=full(x0(:,i));
	  [x(:,i), myoptions]=ZGNLSYMilupacksolver(A,PREC,myoptions,...
	                                           full(b(:,i)),myx0);
	  options.niter(i)=myoptions.niter;
      end % for i
       
   % in any case the preconditioner is unsymmetric -> GMRES	 
   else
      for i=1:size(b,2)
	  myx0=full(x0(:,i));
	  [x(:,i), myoptions]=ZGNLilupacksolver(A,PREC,myoptions,...
	                                        full(b(:,i)),myx0);
	  options.niter(i)=myoptions.niter;
      end % for i
   end % if-elseif-else myoptions.ishermitian & PREC(1).ishermitian
end % if



% update options
options.amg           =myoptions.amg;
options.restol        =myoptions.restol;
options.maxit         =myoptions.maxit;
options.nrestart      =myoptions.nrestart;
options.npresmoothing =myoptions.npresmoothing;
options.npostsmoothing=myoptions.npostsmoothing;
options.ncoarse       =myoptions.ncoarse;
options.presmoother   =myoptions.presmoother;
options.postsmoother  =myoptions.postsmoother;
options.solver        =myoptions.solver;
options.damping       =myoptions.damping;
options.mixedprecision=myoptions.mixedprecision;
