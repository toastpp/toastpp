#include <ilupack_fortran.h>
#include <ilupackmacros.h>

#ifdef _DOUBLE_REAL_
#define MYSYMILUPACKINIT       dsymamginit
#define MYSYMILUPACKFACTOR     dsymamgfactor
#define MYSYMILUPACKSOL        dsymamgsol
#define MYSYMILUPACKSOLVER     dsymamgsolver
#define MYSYMILUPACKDELETE     dsymamgdelete
#define MYSYMILUPACKINFO       dsymamginfo
#define MYSYMILUPACKNNZ        dsymamgnnz


#elif defined _SINGLE_REAL_
#define MYSYMILUPACKINIT       ssymamginit
#define MYSYMILUPACKFACTOR     ssymamgfactor
#define MYSYMILUPACKSOL        ssymamgsol
#define MYSYMILUPACKSOLVER     ssymamgsolver
#define MYSYMILUPACKDELETE     ssymamgdelete
#define MYSYMILUPACKINFO       ssymamginfo
#define MYSYMILUPACKNNZ        ssymamgnnz


#elif defined _SINGLE_COMPLEX_

#ifdef _COMPLEX_SYMMETRIC_
#define MYSYMILUPACKINIT       csymamginit
#define MYSYMILUPACKFACTOR     csymamgfactor
#define MYSYMILUPACKSOL        csymamgsol
#define MYSYMILUPACKSOLVER     csymamgsolver
#define MYSYMILUPACKDELETE     csymamgdelete
#define MYSYMILUPACKINFO       csymamginfo
#define MYSYMILUPACKNNZ        csymamgnnz
#else
#define MYSYMILUPACKINIT       cheramginit
#define MYSYMILUPACKFACTOR     cheramgfactor
#define MYSYMILUPACKSOL        cheramgsol
#define MYSYMILUPACKSOLVER     cheramgsolver
#define MYSYMILUPACKDELETE     cheramgdelete
#define MYSYMILUPACKINFO       cheramginfo
#define MYSYMILUPACKNNZ        cheramgnnz
#endif


#else

#ifdef _COMPLEX_SYMMETRIC_
#define MYSYMILUPACKINIT       zsymamginit
#define MYSYMILUPACKFACTOR     zsymamgfactor
#define MYSYMILUPACKSOL        zsymamgsol
#define MYSYMILUPACKSOLVER     zsymamgsolver
#define MYSYMILUPACKDELETE     zsymamgdelete
#define MYSYMILUPACKINFO       zsymamginfo
#define MYSYMILUPACKNNZ        zsymamgnnz
#else
#define MYSYMILUPACKINIT       zheramginit
#define MYSYMILUPACKFACTOR     zheramgfactor
#define MYSYMILUPACKSOL        zheramgsol
#define MYSYMILUPACKSOLVER     zheramgsolver
#define MYSYMILUPACKDELETE     zheramgdelete
#define MYSYMILUPACKINFO       zheramginfo
#define MYSYMILUPACKNNZ        zheramgnnz
#endif
#endif

      program fmainsym2
      implicit none
      integer n,nnz, i,j,k,l     
      integer, dimension(:), allocatable ::  ia, ja
      FLOAT, dimension(:), allocatable ::  a, rhs, sol
      integer  mem

c     external ILUPACK functions: `_' in {S, D, C, Z} `###' in {SYM, HER}
c              _###AMGinit         init default parameter  
c              _###AMGfactor       compute multilevel ILU `PREC'
c              _###AMGsol          solve a single linear system with `PREC'
c              _###AMGsolver       solve linear system Ax=b iteratively
c                                  using the ILUPACK preconditioner
c              _###AMGdelete       release memory

      integer  MYSYMILUPACKFACTOR, MYSYMILUPACKSOLVER, MYSYMILUPACKNNZ
      external MYSYMILUPACKINIT, MYSYMILUPACKSOL, MYSYMILUPACKDELETE,
     +         MYSYMILUPACKFACTOR, MYSYMILUPACKSOLVER,
     +         MYSYMILUPACKNNZ, MYSYMILUPACKINFO, READMTC

c     ILUPACK external parameters
      integer   matching, maxit, lfil, lfilS, nrestart, ierr
      character ordering*20
      REALS     droptol, droptolS, condest, restol, elbow
c     variables that cover the and pass the C-pointers
      INTEGER*8 param,PREC


*     ***************************************************************
*     read in complex symmetric sample matrix in CSR format
      fname='mplate.csa'
*     first call to ZREADMTC returns the size of the matrix
*     and its fill
      nmax=0
      nzmax=0
      nrhs=0
      job=0
      ierr=0
      call Zreadmtc(nmax,nzmax,job,fname,Aa,Aja,Aia,rhs,nrhs,
     *              rhstyp,nrow,ncol,nnz,title,key,type,
     *              nrhsix,rhsptr,rhsind,rhsval,ierr)
      if (ierr/=0) then
         WRITE (6,'(A,I3)')'Zreadmtc terminated with error code',ierr
         goto 999
      end if
      write(6,'(A,A,A,I5,A,I5,A,I6,A,1P,E8.1,A)')'Matrix: ',fname,
     +     ': size (',nrow,',',ncol,
     +     '), nnz=',nnz,'(',nnz/dble(nrow),' average per row)'




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



c     For `_' in {S, D, C, Z}: _###AMGinit              
c     init default parameter  
      call MYSYMILUPACKINIT(n,ia,ja,a,matching,ordering,
     +                      droptol, droptolS, condest, restol,
     +                      maxit, elbow, lfil, lfilS, nrestart)


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

c     threshold for the approximate Schur complements, default: 1e-2
c     droptolS=0.1

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




c     compute multilevel ILU
c     cccccccccccccccccccccc
c     Note that the initial input matrix A will be rescaled by rows and
c     by columns (powers of 2.0) and that the order in the array might have
c     been altered
c     if you do need the original matrix (ia,ja,a) in for differen purposes,
c     you should use a copy (ib,jb,b) instead

c     For `_' in {S, D, C, Z}: _###AMGfactor 
c     compute multilevel ILU `PREC'
      ierr=MYSYMILUPACKFACTOR(param,PREC,
     +                        n,ib,jb,b,matching,ordering,
     +                        droptol, droptolS, condest, restol,
     +                        maxit, elbow, lfil, lfilS, nrestart)


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
      write (6,'(A,F8.2)') '   final elbow space factor=',
     +      elbow+0.005
      write (6,'(A,F8.2)') '   final condest on level 1=',
     +      condest+0.005
      write (6,'(A)')      'ILUPACK,   multilevel structure'
      call MYSYMILUPACKINFO(param,PREC,n,ib,jb,b)

c     Just for fun: if you want to know the logical number of nonzeros only
      mem=MYSYMILUPACKNNZ(param,PREC)
      write (6,'(A,1P,E8.1)') 'fill-in factor nnz(LU)/nnz(A)',
     +     dble(mem)/dble(ia(n+1)-1)



c     solve a single system with `PREC'
c     ccccccccccccccccccccccccccccccccc
c     This might be of interest if you want to apply ILUPACK inside your
c     own iterative method without referring to the convenient driver
c     artificial right hand side b=A*1
#if !defined _SINGLE_REAL_ && !defined _DOUBLE_REAL_
      rhs(1)=dcmplx(2,-1.0)
      do i=2,n-1
         rhs(i)=2.0
      end do
      rhs(n)=dcmplx(2,1.0)
#else
      rhs(1)=1.0
      do i=2,n-1
         rhs(i)=0.0
      end do
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


c     For `_' in {S, D, C, Z}: _###AMGsol
c     solve a single linear system with `PREC'
      call MYSYMILUPACKSOL(param,PREC,rhs,sol,n)


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
c     For `_' in {S, D, C, Z}: _###AMGsolver
c     solve Ax=b iteratively
      ierr=MYSYMILUPACKSOLVER(param,PREC,rhs,sol,
     +                        n,ib,jb,b,matching,ordering,
     +                        droptol, droptolS, condest, restol,
     +                        maxit, elbow, lfil, lfilS, nrestart)


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
c     For `_' in {S, D, C, Z}: _###AMGdelete
c     solve Ax=b iteratively
      call MYSYMILUPACKDELETE(param,PREC);

 999  end


