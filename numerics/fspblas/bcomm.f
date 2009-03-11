c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dbcomm( transa, mb, n, kb, alpha, descra,
     *           val, bindx, bjndx, bnnz, lb, 
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface
c   dbcomm -- block coordinate matrix-matrix multiply
c  
c   C <- alpha A B + beta C
c  
c   Arguments:
c  
c   int transa	Indicates how to operate with the sparse matrix
c  		0 : operate with matrix
c  		1 : operate with transpose matrix
c  
c   int mb	Number of block rows in matrix A
c  
c   int n	Number of columns in matrix c
c  
c   int kb	Number of block columns in matrix A
c  
c   double alpha Scalar parameter
c  
c   double beta  Scalar parameter
c  
c   int descra()	Descriptor argument.  Nine element integer array
c  		descra(1) matrix structure
c  			0 : general
c  			1 : symmetric
c  			2 : Hermitian
c  			3 : Triangular
c  			4 : Skew(Anti)-Symmetric
c  			5 : Diagonal
c  		descra(2) upper/lower triangular indicator
c  			1 : lower
c  			2 : upper
c  		descra(3) main diagonal type
c  			0 : non-unit
c  			1 : unit
c  		descra(4) Array base 
c  			0 : C/C++ compatible
c  			1 : Fortran compatible
c  		descra(5) repeated indices?
c  			0 : unknown
c  			1 : no repeated indices
c  
c  
c   double val()  scalar array of length nnz containing matrix entries
c                 stored column-major within each dense block.
c  
c   int bindx()   integer array of length bnnz consisting of the
c                 block row indices of the block entries of A.
c  
c   int bjndx()   integer array of length bnnz consisting of the
c                 block column indices of the block entries of A.
c  
c   int bnnz      number of block entries
c
c   int lb        dimension of dense blocks composing A.
c  
c   double b()    rectangular array with first dimension ldb.
c  
c   double c()    rectangular array with first dimension ldc.
c  
c   double work() scratch array of length lwork.  lwork should be at least
c                 max(m,n).
c  
c       ------------ end interface description --------------
c--------------------------------------------------------------------
      implicit none
c
c     interface variables:
c
      integer transa, mb, n, kb, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
      double precision b(*), c(*)
      double precision work(*)
c
c     format specific interface variables:
c
      integer lb, bnnz
      integer bindx(*), bjndx(*)
      double precision val(*)
c
c     local variables:
c
      integer m, k, info
      character transpose
c
c     externals:
c
      external xerbla

c
c     Test input parameters:
c
      info = 0
      if ( (transa .ne. 0) .and. (transa .ne. 1) ) then
         info = 1
      else if ( mb .lt. 0 ) then
         info = 2
      else if (n .lt. 0) then
         info = 3
      else if (kb .lt. 0) then
         info = 4
      else
        m=mb*lb
        k=kb*lb
        if (transa .eq. 0) then
c         Check for consistant dimensions:
          if ( ldb .lt. k ) then 
            info = 15
          else if (ldc .lt. m) then
            info = 18
          endif
        else if (transa .eq. 1) then
c         Check for consistant dimensions:
          if ( ldb .lt. m ) then 
            info = 15
          else if (ldc .lt. k) then
            info = 18
          endif
        endif
      endif

      if ( info .ne. 0 ) then
        call xerbla('BCOMM', info)
        return
      endif


      if ( (descra(1) .ge. 0 .and. descra(1) .le. 5 ) .and.
     *      alpha .eq. 0.D0                                 ) then
c       Quick return after scaling:
        call dscal(m*n, beta, c, 1)
        return
       endif
      
      transpose = 'N'
      if ( transa .eq. 1 ) transpose = 'T'
 
c
c     Call appropriate kernel subroutine:
c

      if (descra(1) .eq. 0   .or.
     *    descra(1) .eq. 3   .or.
     *    descra(1) .eq. 5        ) then
c
c        General matrix multiply:
c
        call dbcommgk( transpose,  mb, n, kb, alpha,
     *       val, bindx, bjndx, bnnz, lb, 
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dbcommsk( mb, n, kb, alpha,
     *       val, bindx, bjndx, bnnz, lb, 
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call dbcommkk( transpose,  mb, n, kb, alpha,
     *       val, bindx, bjndx, bnnz, lb, 
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('BCOMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DBCOMMGK( trans, mb, n, kb, alpha, 
     *                     val, bindx, bjndx, bnnz,
     *                     lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, bnnz, lb, ldb, ldc, base
      integer lbsq
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, k, l
      integer brow, bcol, vb, vc, cb, bb
      integer rcol, rhscols
      double precision alpha
      double precision beta
      integer bindx(*), bjndx(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Loop through blocks:
c
      rcol = lb
      if ( n .lt. rcol ) rcol = n
      rhscols = n/rcol
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      vb = 1 - lbsq
      do 5 i=1,bnnz
c       vb = starting (point) element of A block
        vb = vb + lbsq
        brow = bindx(i)
        bcol = bjndx(i)
c       cb = starting (point) row of block
        cb = (brow-1)*lb + 1
c       vc = starting (point) column of A block
        vc = (bcol-1)*lb + 1
        bb = 1 - rcol
        do 15 l=1,rhscols
c          bb = starting (point) column of b 
           bb = bb + rcol
           k = rcol
           if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
           if ( trans .eq. 'N' ) then
              call dgemm('N','N',lb,k,lb,alpha,val(vb),lb,
     *                                        b(vc,bb),ldb,
     *                                    one,c(cb,bb),ldc)
           else
              call dgemm('T','N',lb,k,lb,alpha,val(vb),lb,
     *                                        b(cb,bb),ldb,
     *                                    one,c(vc,bb),ldc)
           endif
 15     continue
 5    continue
             
      return
      end
           
      subroutine DBCOMMSK( mb, n, kb, alpha, val, bindx, bjndx, bnnz,
     *                          lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, bnnz, lb, ldb, ldc, base
      integer lbsq
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, k, l
      integer brow, bcol, vb, vc, cb, bb
      integer rcol, rhscols
      double precision alpha
      double precision beta
      integer bindx(*), bjndx(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( lb*mb*n, beta, c(1,1), 1)
c
c     Loop through blocks:
c
      rcol = lb
      if ( n .lt. rcol ) rcol = n
      rhscols = n/rcol
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      vb = 1 - lbsq
      do 5 i=1,bnnz
c       vb = starting (point) element of A block
        vb = vb + lbsq
        brow = bindx(i)
        bcol = bjndx(i)
c       cb = starting (point) row of block
        cb = (brow-1)*lb + 1
c       vc = starting (point) column of A block
        vc = (bcol-1)*lb + 1
        bb = 1 - rcol
        do 15 l=1,rhscols
c          bb = starting (point) column of b 
           bb = bb + rcol
           k = rcol
           if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
           call dgemm('N','N',lb,k,lb,alpha,val(vb),lb,
     *                                        b(vc,bb),ldb,
     *                                    one,c(cb,bb),ldc)
           if ( brow .ne. bcol ) then
             call dgemm('T','N',lb,k,lb,alpha,val(vb),lb,
     *                                          b(cb,bb),ldb,
     *                                      one,c(vc,bb),ldc)
           endif
 15     continue
 5    continue
             
      return
      end
           
      subroutine DBCOMMKK( trans, mb, n, kb, alpha, 
     *                     val, bindx, bjndx, bnnz,
     *                     lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, bnnz, lb, ldb, ldc, base
      integer lbsq
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, k, l
      integer brow, bcol, vb, vc, cb, bb
      integer rcol, rhscols
      double precision alpha
      double precision beta
      integer bindx(*), bjndx(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( lb*mb*n, beta, c(1,1), 1)
c
c     Loop through blocks:
c
      rcol = lb
      if ( n .lt. rcol ) rcol = n
      rhscols = n/rcol
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      vb = 1 - lbsq
      do 5 i=1,bnnz
c       vb = starting (point) element of A block
        vb = vb + lbsq
        brow = bindx(i)
        bcol = bjndx(i)
c       cb = starting (point) row of block
        cb = (brow-1)*lb + 1
c       vc = starting (point) column of A block
        vc = (bcol-1)*lb + 1
        bb = 1 - rcol
        do 15 l=1,rhscols
c          bb = starting (point) column of b 
           bb = bb + rcol
           k = rcol
           if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
           if ( trans .eq. 'N' ) then
              call dgemm('N','N',lb,k,lb,alpha,val(vb),lb,
     *                                        b(vc,bb),ldb,
     *                                    one,c(cb,bb),ldc)
              call dgemm('T','N',lb,k,lb,-alpha,val(vb),lb,
     *                                        b(cb,bb),ldb,
     *                                    one,c(vc,bb),ldc)
           else
              call dgemm('N','N',lb,k,lb,-alpha,val(vb),lb,
     *                                        b(vc,bb),ldb,
     *                                    one,c(cb,bb),ldc)
              call dgemm('T','N',lb,k,lb,alpha,val(vb),lb,
     *                                        b(cb,bb),ldb,
     *                                    one,c(vc,bb),ldc)
           endif
 15     continue
 5    continue
             
      return
      end
           
