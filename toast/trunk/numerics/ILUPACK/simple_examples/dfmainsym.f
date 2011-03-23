      program dfmainsym
      implicit none
      integer n,nnz
      parameter (n=10,nnz=2*n-1)
c
c     sparse matrices in compressed sparse row format
      integer         ia(n+1), ja(nnz)
      doubleprecision a(nnz)
      integer         ib(n+1), jb(nnz)
      doubleprecision b(nnz)
c
c     right hand side and approximate solution
      doubleprecision rhs(nnz)
      doubleprecision sol(nnz)
c
c     auxiliary variables
      integer         i,j,mem
c
c     external ILUPACK functions
c              DSYMAMGinit         init default parameter  
c              DSYMAMGfactor       compute multilevel ILU `PREC'
c              DSYMAMGsol          solve a single linear system with `PREC'
c              DSYMAMGsolver       solve linear system Ax=b iteratively
c                                  using the ILUPACK preconditioner
c              DSYMAMGdelete       release memory
c              DSYMAMGnnz          logical number of nonzero entries of ILU
c              DSYMAMGinfo         display multilevel structure
c
      integer  DSYMAMGfactor, DSYMAMGsolver, DSYMAMGnnz
      external DSYMAMGinit,   DSYMAMGsol, DSYMAMGdelete,
     +         DSYMAMGfactor, DSYMAMGsolver,
     +         DSYMAMGnnz, DSYMAMGinfo
c
c     ILUPACK external parameters
      integer         matching, maxit, lfil, lfilS, nrestart, ierr, 
     +                mixedprecision, ind(n)
      character       ordering*20
      doubleprecision droptol, droptolS, condest, restol, elbow
c
c     variables that cover the and pass the C-pointers
      integer*8 param, PREC
c
c
c     tridiagonal sample nxn matrix in sparse row format
c     pointer
      ia(1)=1
      do i=1,n-1
c        indices
         ja(2*i-1)=i
         ja(2*i)=i+1
c        numerical values
         a(2*i-1)= 2.0
         a(2*i)=-1.0
c        next pointer
         ia(i+1)=ia(i)+2
      end do
c     indices
      ja(2*n-1)=n
c     numerical values
      a(2*n-1)=2.0
c     next pointer
      ia(n+1)=ia(n)+1
c     END tridiagonal sample nxn matrix in sparse row format
c
c
c     compute a copy of a (there are situations where this might
c                          be required since ILUPACK alters the
c                          input matrix)
      do i=1,n+1
         ib(i)=ia(i)
      end do
      do i=1,nnz
         jb(i)=ja(i)
         b(i) =a(i)
      end do
c
c
c
c     init default parameters  
      call DSYMAMGinit(n,ia,ja,a,matching,ordering,
     +                 droptol, droptolS, condest, restol,
     +                 maxit, elbow, lfil, lfilS, nrestart, 
     +                 mixedprecision, ind)
c
c
c     now the use may vary the parameters to gain optimal performance
c
c     maximum weight matching
c     default value is different from zero, matching turned on
c     matching=1
c
c     multilevel orderings
c     'amd' (default) Approximate Minimum Degree
c     'mmd'           Minimum Degree            
c     'rcm'           Reverse Cuthill-McKee
c     'metisn'        Metis multilevel nested dissection by nodes
c     'metise'        Metis multilevel nested dissection by edges
c     'pq'            ddPQ strategy by Saad
c     ordering='metisn'
c
c     threshold for ILU, default: 1e-2
      droptol=0.1
c
c     threshold for the approximate Schur complements, default: 0.1*droptol
      droptolS=0.1*droptol
c
c     norm bound for the inverse factor L^{-1}, default: 10
      condest=5
c
c     relative error for the backward error (SPD case: relative energy
c     norm) used during the iterative solver, default: sqrt(eps)
      restol=1e-12
c
c     maximum number of iterative steps, default: 500
c     maxit=1000
c
c     elbow space factor for the relative fill w.r.t. the given matrix, 
c     computed during the ILU, default: 10
c     this value only gives a recommendation for the elbow space. ILUPACK
c     will adapt this value automatically
c     elbow=15
c
c     maximum number of nonzeros per column in L/ row in U, default: n+1
c     lfil=10 * (nnz/(1.0*n))
c
c     maximum number of nonzeros per row in the approximate Schur complement,
c     default: n+1
c     lfilS=10 * (nnz/(1.0*n))
c
c     ignored for symmetric matrices
c     nrestart=20
c
c     do you want to use a single precision preconditioner
c      mixedprecision=1
c
c     underlying block structure, only partially supported
c     default: right now turn it off!
      do i=1,n
         ind(i)=0
      end do
c
c
c     compute multilevel ILU
c     cccccccccccccccccccccc
c     Note that the initial input matrix A will be rescaled by rows and
c     by columns (powers of 2.0) and that the order in the array might have
c     been altered
c     if you do need the original matrix (ia,ja,a) in for different purposes,
c     you should use a copy (ib,jb,b) instead
c
c     compute multilevel ILU `PREC'
      ierr=DSYMAMGfactor(param,PREC,
     +                   n,ib,jb,b,matching,ordering,
     +                   droptol, droptolS, condest, restol,
     +                   maxit, elbow, lfil, lfilS, nrestart, 
     +                   mixedprecision, ind)
c
c
      if (ierr.eq.-1) then
         write (6,'(A)') 'Error. input matrix may be wrong.'
      elseif (ierr.eq.-2) then
         write (6,'(A)') 'matrix L overflow, increase elbow and retry'
      elseif (ierr.eq.-3) then
         write (6,'(A)') 'matrix U overflow, increase elbow and retry'
      elseif (ierr.eq.-4) then
         write (6,'(A)') 'Illegal value for lfil'
      elseif (ierr.eq.-5) then
         write (6,'(A)') 'zero row encountered'
      elseif (ierr.eq.-6) then
         write (6,'(A)') 'zero column encountered'
      elseif (ierr.eq.-7) then
         write (6,'(A)') 'buffers are too small'
      elseif (ierr.ne.0) then
         write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
      endif
      if (ierr.ne.0) goto 999
c
c
c
c     Just for fun: display multilevel information on the screen
      write (6,'(A,F8.2)') '   final elbow space factor=',
     +      elbow+0.005
      write (6,'(A,F8.2)') '   final condest on level 1=',
     +      condest+0.005
      write (6,'(A)')      'ILUPACK,   multilevel structure'
      call DSYMAMGinfo(param,PREC,n,ib,jb,b)
c
c     Just for fun: if you want to know the logical number of nonzeros only
      mem=DSYMAMGnnz(param,PREC)
      write (6,'(A,1P,E8.1)') 'fill-in factor nnz(LU)/nnz(A)',
     +     dble(mem)/dble(ia(n+1)-1)
c
c
c
c     solve a single system with `PREC'
c     ccccccccccccccccccccccccccccccccc
c     This might be of interest if you want to apply ILUPACK inside your
c     own iterative method without referring to the convenient driver
c
c     artificial right hand side b=A*1
      rhs(1)=1.0
      do i=2,n-1
         rhs(i)=0.0
      end do
      rhs(n)=1.0
      write (6,'(A)') 'right hand side'
      do i=1,n
         write(6,'(1P,E12.4)') rhs(i)
      enddo
c
c
c     solve a single linear system with `PREC'
      call DSYMAMGsol(param,PREC,rhs,sol,n)
c
      write (6,'(A)') 'approximate solution'
      do i=1,n
         write(6,'(1P,E12.4)') sol(i)
      enddo
c
c
c
c
c
c     ok, the preconditioner is usually not exact.
c     as convenient ALTERNATIVE, ILUPACK offers an all-in-one iterative solver
c
c
c
c
c
c     solve Ax=b  until the desired accuracy is achieved
c     cccccccccccccccccccccccccccccccccccccccccccccccccc
c     provide an initial solution, e.g. 0
      do i=1,n
         sol(i)=0
      end do
c     solve Ax=b iteratively
      ierr=DSYMAMGsolver(param,PREC,rhs,sol,
     +                   n,ib,jb,b,matching,ordering,
     +                   droptol, droptolS, condest, restol,
     +                   maxit, elbow, lfil, lfilS, nrestart, 
     +                   mixedprecision, ind)
c
c
      if (ierr.eq.-1) then
         write (6,'(A)') 'too many iteration steps'
      elseif (ierr.eq.-2) then
         write (6,'(A)') 'not enough work space'
      elseif (ierr.eq.-3) then
         write (6,'(A)') 'algorithm breaks down'
      elseif (ierr.ne.0) then
         write (6,'(A,I4)') 'solver exited with error code',ierr
      end if
      if (ierr.ne.0) goto 999

      write (6,'(I4,A)') maxit,' iterations steps needed'

      write (6,'(A)') 'approximate solution'
      do i=1,n
         write(6,'(1P,E12.4)') sol(i)
      enddo
c      
c
c
c
c
c     Finally release memory
c     cccccccccccccccccccccc
c     solve Ax=b iteratively
      call DSYMAMGdelete(param,PREC);
c
c
 999  end


