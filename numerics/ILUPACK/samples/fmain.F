#include <ilupack_fortran.h>
#include <ilupackmacros.h>

      program fmain
      implicit none
      integer n,nnz
      parameter (n=10,nnz=3*n-2)
      integer  ia(n+1), ja(nnz), i,j     
      FLOAT a(nnz)
      integer  ib(n+1), jb(nnz), mem
      FLOAT b(nnz)
      FLOAT rhs(nnz)
      FLOAT sol(nnz)

c     external ILUPACK functions: `_' in {S, D, C, Z}
c              _GNLAMGinit         init default parameter  
c              _GNLAMGfactor       compute multilevel ILU `PREC'
c              _GNLAMGsol          solve a single linear system with `PREC'
c              _GNLAMGsolver       solve linear system Ax=b iteratively
c                                  using the ILUPACK preconditioner
c              _GNLAMGdelete       release memory

      integer  ILUPACKFACTOR, ILUPACKSOLVER, ILUPACKNNZ
      external ILUPACKINIT, ILUPACKSOL, ILUPACKDELETE,
     +         ILUPACKFACTOR, ILUPACKSOLVER,
     +         ILUPACKNNZ, ILUPACKINFO

c     ILUPACK external parameters
      integer   matching, maxit, lfil, lfilS, nrestart, ierr, 
     +          mixedprecision, ind(n)
      character ordering*20
      REALS     droptol, droptolS, condest, restol, elbow
c     variables that cover the and pass the C-pointers
      INTEGER*8 param,PREC


c     tridiagonal sample nxn matrix in sparse row format
c     pointer
      ia(1)=1
      ia(2)=ia(1)+2
c     indices
      ja(1)=1
      ja(2)=2
c     numerical values
      a(1)= 2.0
      a(2)= 1.0
c     next pointer
      do i=2,n-1
c        indices
         ja(3*i-3)=i-1
         ja(3*i-2)=i
         ja(3*i-1)=i+1
c        numerical values
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
         a(3*i-3)=dcmplx(0,-1.0)
#else
         a(3*i-3)=-1.0
#endif
         a(3*i-2)= 2.0
         a(3*i-1)= 1.0
c        next pointer
         ia(i+1)=ia(i)+3
      end do
c     indices
      ja(3*n-3)=n-1
      ja(3*n-2)=n
c     numerical values
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
      a(3*n-3)=dcmplx(0,-1.0)
#else
      a(3*n-3)=-1.0
#endif
      a(3*n-2)= 2.0
c     next pointer
      ia(n+1)=ia(n)+2


c     compute a copy of a (there are situation where this might
c                          be required since ILUPACK alters the
c                          input matrix)
      do i=1,n+1
         ib(i)=ia(i)
      end do
      do i=1,nnz
         jb(i)=ja(i)
         b(i) =a(i)
      end do



c     For `_' in {S, D, C, Z}: _GNLAMGinit              
c     init default parameter  
      call ILUPACKINIT(n,ia,ja,a,matching,ordering,
     +                 droptol, droptolS, condest, restol,
     +                 maxit, elbow, lfil, lfilS, nrestart, 
     +                 mixedprecision, ind)


c     now the use may vary the parameters to gain optimal performance

c     maximum weight matching
c     default value is different from zero, matching turned on
c     matching=1

c     multilevel orderings
c     'amd' (default) Approximate Minimum Degree
c     'mmd'           Minimum Degree            
c     'rcm'           Reverse Cuthill-McKee
c     'metisn'        Metis multilevel nested dissection by nodes
c     'metise'        Metis multilevel nested dissection by edges
c     'pq'            ddPQ strategy by Saad

c     threshold for ILU, default: 1e-2
      droptol=0.1

c     threshold for the approximate Schur complements, default: 0.1*droptol
      droptolS=0.1*droptol

c     norm bound for the inverse factors L^{-1}, U^{-1}, default: 1e2
      condest=5

c     relative error for the backward error (SPD case: relative energy
c     norm) used during the iterative solver, default: sqrt(eps)
      restol=1e-12

c     maximum number of iterative steps, default: 500
c     maxit=1000

c     elbow space factor for the fill computed during the ILU, default: 10
c     elbow=15

c     maximum number of nonzeros per column in L/ row in U, default: n+1
c     lfil=10

c     maximum number of nonzeros per row in the approximate Schur complement,
c     default: n+1
c     lfilS=10

c     restart length for GMRES, default: 30
c     nrestart=20

c     do you want to use a single precision preconditioner
      mixedprecision=1

c     underlying block structure, only partially supported
c     default: right now turn it off!
      do i=1,n
         ind(i)=0
      end do


c     compute multilevel ILU
c     cccccccccccccccccccccc
c     Note that the initial input matrix A will be rescaled by rows and
c     by columns (powers of 2.0) and that the order in the array might have
c     been altered
c     if you do need the original matrix (ia,ja,a) in for differen purposes,
c     you should use a copy (ib,jb,b) instead

c     For `_' in {S, D, C, Z}: _GNLAMGfactor 
c     compute multilevel ILU `PREC'
      ierr=ILUPACKFACTOR(param,PREC,
     +                   n,ib,jb,b,matching,ordering,
     +                   droptol, droptolS, condest, restol,
     +                   maxit, elbow, lfil, lfilS, nrestart, 
     +                   mixedprecision, ind)


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



c     Just for fun: display multilevel information on the screen
c     For `_' in {S, D, C, Z}: _GNLAMGinfo
      write (6,'(A,F8.2)') '   final elbow space factor=',
     +      elbow+0.005
      write (6,'(A,F8.2)') '   final condest on level 1=',
     +      condest+0.005
      write (6,'(A)')      'ILUPACK,   multilevel structure'
      call ILUPACKINFO(param,PREC,n,ib,jb,b)

c     Just for fun: if you want to know the logical number of nonzeros only
      mem=ILUPACKNNZ(param,PREC)
      write (6,'(A,1P,E8.1)') 'fill-in factor nnz(LU)/nnz(A)',
     +     dble(mem)/dble(ia(n+1)-1)



c     solve a single system with `PREC'
c     ccccccccccccccccccccccccccccccccc
c     This might be of interest if you want to apply ILUPACK inside your
c     own iterative method without referring to the convenient driver
c     artificial right hand side b=A*1
      rhs(1)=3.0
      do i=2,n-1
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
         rhs(i)=dcmplx(3,-1.0)
#else
         rhs(i)=2.0
#endif
      end do
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
      rhs(n)=dcmplx(2,-1.0)
#else
      rhs(n)=1.0
#endif
      write (6,'(A)') 'right hand side'
      do i=1,n
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
#ifdef _SINGLE_COMPLEX_
         write(6,'(1P,E12.4,1P,E12.4)') real(rhs(i)),imag(rhs(i))
#else
         write(6,'(1P,E12.4,1P,E12.4)') dreal(rhs(i)),dimag(rhs(i))
#endif
#else
         write(6,'(1P,E12.4)') rhs(i)
#endif
      enddo


c     For `_' in {S, D, C, Z}: _GNLAMGsol
c     solve a single linear system with `PREC'
      call ILUPACKSOL(param,PREC,rhs,sol,n)


      write (6,'(A)') 'approximate solution'
      do i=1,n
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
#ifdef _SINGLE_COMPLEX_
         write(6,'(1P,E12.4,1P,E12.4)') real(sol(i)),imag(sol(i))
#else
         write(6,'(1P,E12.4,1P,E12.4)') dreal(sol(i)),dimag(sol(i))
#endif
#else
         write(6,'(1P,E12.4)') sol(i)
#endif
      enddo






c     ok, the preconditioner is usually not exact.
c     as convenient ALTERNATIVE, ILUPACK offers an all-in-one iterative solver





c     solve Ax=b  until the desired accuracy is achieved
c     cccccccccccccccccccccccccccccccccccccccccccccccccc
c     provide an initial solution, e.g. 0
      do i=1,n
         sol(i)=0
      end do
c     For `_' in {S, D, C, Z}: _GNLAMGsolver
c     solve Ax=b iteratively
      ierr=ILUPACKSOLVER(param,PREC,rhs,sol,
     +                   n,ib,jb,b,matching,ordering,
     +                   droptol, droptolS, condest, restol,
     +                   maxit, elbow, lfil, lfilS, nrestart, 
     +                   mixedprecision, ind)


      if (ierr.eq.-1) then
         write (6,'(A)') 'too many iterations'
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
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
#ifdef _SINGLE_COMPLEX_
         write(6,'(1P,E12.4,1P,E12.4)') real(sol(i)),imag(sol(i))
#else
         write(6,'(1P,E12.4,1P,E12.4)') dreal(sol(i)),dimag(sol(i))
#endif
#else
         write(6,'(1P,E12.4)') sol(i)
#endif
      enddo
      




c     Finally release memory
c     cccccccccccccccccccccc
c     For `_' in {S, D, C, Z}: _GNLAMGdelete
c     solve Ax=b iteratively
      call ILUPACKDELETE(param,PREC);


 999  end


