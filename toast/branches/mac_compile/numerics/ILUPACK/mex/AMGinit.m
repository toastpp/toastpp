function options = AMGinit(A,options)
% options = AMGinit(A)
% options = AMGinit(A,options)
% 
% init structure of options to their default values for a given nxn matrix A
% 
% input
% -----
% A         nxn matrix
% options   some of the options you might have set in advance, overwritten
%           on output. In particular setting 'options.isdefinite=1' you 
%           indicate that your system is symmetric (Hermitian) positive definite
%           resulting different initial values
%
% output
% ------
% options   structure with default parameters
%
% ------------------------------------------------------------------------
%
% possible options:
%  1. options.matching 
% --------------------
%    if different from zero then maximum weight matching will be used,
%    for the real symmetric, complex Hermitian or complex symmetric the
%    associated symmetric maximum weight matching will be used.
%    `maximum weight matching' is a technique to reorder and rescale the
%    matrix such that it is becomes more diagonally dominant
%    default value: 1
%
%  2. options.ordering 
% --------------------
%    several reorderings based on |A|+|A|' are offered. Unsymmetric patterns
%    treated as if they were symmetric.
%    The orderings are repeated on any coarser level.
% 
%    `amd'     default,  Approximate Minimum Degree
%    `metisn'            METIS multilevel nested dissection by NODES
%    `metise'            METIS multilevel nested dissection by EDGES
%    `rcm'               Reverse Cuthill-McKee
%    `mmd'               Minimum Degree   
%    `amf'               Approximate Minimum Fill
%    any other           no reordering
%
%  3. options.droptol
% -------------------
%    threshold for dropping small entries during the computation of the
%    incomplete LU decomposition
%    default: 1e-2
%
%  4. options.droptolS
% --------------------
%    threshold for dropping small entries from the Schur complement
%    default: 1e-2
%
%  5. options.condest
% -------------------
%    bound for the inverse triangular factors from the incomplete LU 
%    decomposition
%    default: 100
%
%  6. options.restol
% ------------------
%    bound for the accuracy of the approximate solution of Ax=b for a given 
%    right hand side after using ILUPACK-preconditioned iterative solution
%    This tolerance refers to the BACKWARD ERROR, in the symmetric(Hermitian)
%    positive definite case to the relative error in the energy norm
%    default: sqrt(eps)
%
%  7. options.maxit
% -----------------
%    maximum number of iteration steps, before the iterative ILUPACK-
%    preconditioned solvers terminates
%    default: 500
%
%  8. options.elbow
% -----------------
%    elbow space for memory of the ILUPACK multilevel preconditioner. 
%    Since the core part of the code is FORTRAN 77, no dynamic memory
%    allocation is available, one has to estimate the memory requirement
%    in advance by multiples of the fill of the initial matrix
%    default: 10     (ten times as much fill as the initial matrix)
%
%  9. options.lfil
% ---------------
%    restrict the number of nonzeros per column in L (resp. per row in U)
%    hard to at most `lfil' entries.
%    default: n+1
%
% 10. options.lfilS
% ----------------
%    restrict the number of nonzeros per row in the approximate Schur 
%    complement hard to at most `lfilS' entries.
%    default: n+1
%
% 11. options.typetv
% ------------------
%    define, whether to include a test vector into the computations
%    If used, then (a) the test vector is also included to estimate
%                      the norm of the inverse triangular factors
%                  (b) it is used to obtain a refined fine/coarse
%                      grid partitioning (if switch is set)
%                  (c) diagonal compensation and off-diagonal lumping
%                      are added to improve the factorization
%    default:      'none'
%    alternatives: 'static'         for a fixed test vector
%                  'function_name'  for a dynamically generated test 
%                                   vector the user has to provide a 
%                                   custom routine to generate the test
%                                   vector 
%                                   format:
%                                   tvd=function_name(mat,tv_old)
%                                   On every level the associate matrix
%                                   'mat' is passed to this routine,
%                                   also the coarse grid projection of
%                                   the previous test vector 'tv_old' from
%                                   the finer grid is passed. On exit, 'tvd'
%                                   is the dynamic test vector that should be
%                                   used
%
% 12. options.tv
% --------------
%     static test vector. Ignored if options.typetv=='none'
%     If options.typetv=='function_name', then on the initial finest level,
%     options.tv is passed to 'function_name' as initial guess for 'tv_old' 
%
% 13. options.amg
% ---------------
%     type of AMG preconditioner
%     default:      'ilu'        multilevel ILU
%     alternatives: 'amli'       multilevel ILU, where on each coarser grid
%                                an inner iteration is used to solve the
%                                coarse grid system. The number of inner 
%                                interation steps is prescribed in 
%                                options.ncoarse
%                   'mg'         classical multigrid with pre and post
%                                smoothing
%
% 14. options.npresmoothing
% -------------------------
%     number of pre smoothing steps (only needed if options.amg=='mg')
%     default:      1
%
% 15. options.npostsmoothing
% --------------------------
%     number of post smoothing steps (only needed if options.amg=='mg')
%     default:      1
%
% 16. options.ncoarse
% -------------------
%     number of coarse grid correction steps (only needed if options.amg=='mg'
%                                             or options.amg=='amli')
%     default:      1
%
% 17. options.presmoother
% -------------------------
%     type of pre smoother
%     default:      'gsf'           Gauss-Seidel forward
%                   'gsb'           Gauss-Seidel backward
%                   'j'             Jacobi
%                   'ilu'           partial incomplete ILU on the F-nodes
%                   'function_name' custom routine for smoothing
%                                   d=function_name(mat,r)
%                                   Given the matrix 'mat' on each level
%                                   and the associated residual 
%                                   'r(=b-mat*x_old)', the custom smoother
%                                   is asked to compute an approximate defect
%                                   'd' such that x_new=x_old+d
%
% 18. options.postsmoother
% -------------------------
%     type of post smoother
%     default:      'gsb'           Gauss-Seidel backward
%                   'gsf'           Gauss-Seidel forward
%                   'j'             Jacobi
%                   'ilu'           partial incomplete ILU on the F-nodes
%                   'function_name' custom routine for smoothing
%                                   d=function_name(mat,r)
%                                   Given the matrix 'mat' on each level
%                                   and the associated residual 
%                                   'r(=b-mat*x_old)', the custom smoother
%                                   is asked to compute an approximate defect
%                                   'd' such that x_new=x_old+d
%
% 19. options.FCpart
% ------------------
%     preselect a partitioning into fine grid and coarse grid nodes
%     default:      'none'          
%                   'yes'           if F/C partioning is desired
%
% 20. options.typecoarse
% ----------------------
%     type of coarse grid system
%     default:      'ilu'           coarse grid system is computed from the
%                                   approximate incomplete factorization by
%                                   ignoring that entries have been discarded
%                   'amg'           Use the associated approximate interpolation
%                                   operator P and the restriction operator R
%                                   from the underlying inverse triangular
%                                   factors to generate the Galerkin type
%                                   coarse grid matrix R A P
% 21. options.damping
% ----------------------
%     damping factor if Jacobi smoothing is chosen
%
% 22. options.isdefinite
% ----------------------
%     if given on input then the matrix is assumed to be symmetric (Hermitian)
%     positive definite and the parameters will be
%
% 23. options.ind
% ---------------
%     indicate by negative signs which parts of the system refer to the
%     second (typically) zero block of a saddle point system
%
% 24. options.mixedprecision
% ---------------
%     if different from zero, then single precision for preconditioning is
%     used



n=size(A,1);

options.matching=1;
options.ordering='amd';
options.droptol =1e-2;
options.droptolS=1e-2;
options.condest=1e2;
options.restol=1e-12;
options.maxit=500;
options.elbow=10;
options.lfil =0;
options.lfilS=0;
options.typetv='none';
options.tv=zeros(n,1);
options.amg='ilu';
options.npresmoothing=1;
options.npostsmoothing=1;
options.ncoarse=1;
options.presmoother='gsf';
options.postsmoother='gsb';
options.FCpart='none';
options.typecoarse='ilu';
options.solver='gmres';
options.damping=1;
options.nrestart=0;
options.ind=zeros(n,1);
options.mixedprecision=0;


% make sure that shifted system and A have same type
if nargin==2
   if isfield(options,'shiftmatrix')     
      shiftreal=1;
      shifthermitian=1;
      shiftsymmetric=1;
      
      if ~isreal(options.shiftmatrix)
	 shiftreal=0;
      end  
      if isfield(options,'shift0')
	 if ~isreal(options.shift0)
	    shiftreal=0;
	 end
      end
      if isfield(options,'shiftmax')
	 if ~isreal(options.shiftmax)
	    shiftreal=0;
	 end
      end
      if isfield(options,'shifts')
	 if ~isreal(options.shifts)
	    shiftreal=0;
	 end
      end

      if (norm(options.shiftmatrix-options.shiftmatrix',1)==0)
	 shifthermitian=1;
	 % make sure that the shifts are real
	 if isfield(options,'shift0')
	    if ~isreal(options.shift0)
	       shifthermitian=0;
	    end
	 end
	 if isfield(options,'shiftmax')
	    if ~isreal(options.shiftmax)
	       shifthermitian=0;
	    end
	 end
	 if isfield(options,'shifts')
	    if ~isreal(options.shifts)
	       shifthermitian=0;
	    end
	 end
      else
	 shifthermitian=0;
      end
      if (norm(options.shiftmatrix-options.shiftmatrix.',1)==0)
	 shiftsymmetric=1;
      else
	 shiftsymmetric=0;
      end

      % A and shift do not match
      if isreal(A) & ~shiftreal
	 % A becomes complex hermitian
	 if norm(A-A',1)==0 & shifthermitian
	    A(1,2)=A(1,2)+norm(A,1)*eps^2*sqrt(-1);
	    A(1,2)=A(2,1)';
	 % A becomes complex symmetric
	 elseif norm(A-A',1)==0 & shiftsymmetric
	    A(1,1)=A(1,1)+norm(A,1)*eps^2*sqrt(-1);
	 % A becomes complex unsymmetric
	 else
	    A(1,2)=A(1,2)+norm(A,1)*eps^2*sqrt(-1);
	 end
      elseif isreal(A) & shiftreal
	 % A becomes unsymmetric
	 if norm(A-A',1)==0 & ~shiftsymmetric
	    A(1,2)=A(1,2)+norm(A,1)*eps^2;
	 end
      elseif ~isreal(A)
	 % A becomes complex hermitian
	 if norm(A-A',1)==0 & shifthermitian
	    % perfect
	 % A becomes complex symmetric
	 elseif norm(A-A.',1)==0 & shiftsymmetric
	    % perfect
	 % A becomes complex unsymmetric
	 else
	    A(1,2)=A(1,2)+norm(A,1)*eps^2*sqrt(-1);
	 end
      end
   end
end



if isreal(A)
   if norm(A-A',1)==0
      if isfield(options,'isdefinite')
	 if options.isdefinite
	    options=DSPDilupackinit(A,options);
	 else
	    options=DSYMilupackinit(A,options);      
	 end % if-else
      else
	 options=DSYMilupackinit(A,options);      
      end % if
   else
      options=DGNLilupackinit(A,options);      
   end
else
   if norm(A-A',1)==0
      if isfield(options,'isdefinite')
	 if options.isdefinite
	    options=ZHPDilupackinit(A,options);
	 else
	    options=ZHERilupackinit(A,options);      
	 end % if-else
      else
	 options=ZHERilupackinit(A,options);      
      end % if
   elseif norm(A-A.',1)==0
      options=ZSYMilupackinit(A,options);      
   else
      options=ZGNLilupackinit(A,options);      
   end
end % if
