c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dbscmm( transa, mb, n, kb, alpha, descra,
     *           val, bindx, bpntrb, bpntre, lb,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface
c   dbscmm -- block sparse column matrix-matrix multiply
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
c   int bpntrb()  integer array of length kb such that 
c                 bpntrb(i)-bpntrb(1) points to location in bindx
c                 of the first block entry of the j-th block column of A.
c  
c   int bpntre()  integer array of length kb such that 
c                 bpntre(i)-bpntre(1) points to location in bindx
c                 of the last block entry of the j-th block column of A.
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
      integer lb
      integer bindx(*), bpntrb(*), bpntre(*)
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
        call xerbla('BSCMM', info)
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
        call dbscmmgk( transpose,  mb, n, kb, alpha,
     *       val, bindx, bpntrb, bpntre, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dbscmmsk( mb, n, kb, alpha,
     *       val, bindx, bpntrb, bpntre, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call dbscmmkk( transpose,  mb, n, kb, alpha,
     *       val, bindx, bpntrb, bpntre, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('BSCMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DBSCMMGK( trans, mb, n, kb, alpha, val, bindx, bpntrb,
     *                          bpntre, lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, lb, ldb, ldc, base
      integer lbsq
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, k, l, jb, je
      integer vb, vc, cb, bb
      integer rcol, rhscols
      double precision alpha
      double precision beta
      integer bindx(*), bpntrb(*), bpntre(*)
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
c     Loop through block columns:
c
      rcol = lb
      if ( n .lt. rcol ) rcol = n
      rhscols = n/rcol
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      vb = 1 - lbsq
      cb = 1 - lb
      do 5 i=1,kb
c       cb = starting (point) column of block
        cb = cb + lb
        jb = bpntrb(i)
        je = bpntre(i) - 1
        do 10 j=jb,je
c         vb = starting (point) element of A block
          vb = vb + lbsq
c         vc = starting (point) row of A block
          vc = (bindx(j) - 1)*lb + 1
          bb = 1 - rcol
          do 15 l=1,rhscols
c            bb = starting (point) column of b 
             bb = bb + rcol
             k = rcol
             if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
             if ( trans .eq. 'N' ) then
                call dgemm('N','N',lb,k,lb,alpha,val(vb),lb,
     *                                          b(cb,bb),ldb,
     *                                      one,c(vc,bb),ldc)
             else
                call dgemm('T','N',lb,k,lb,alpha,val(vb),lb,
     *                                          b(vc,bb),ldb,
     *                                      one,c(cb,bb),ldc)
             endif
 15       continue
 10     continue
 5    continue
             
      return
      end
           
      subroutine DBSCMMSK( mb, n, kb, alpha, val, bindx, bpntrb,
     *                          bpntre, lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, lb, ldb, ldc, base
      integer lbsq
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, k, l, jb, je
      integer vb, vc, cb, bb
      integer rcol, rhscols, blkind
      double precision alpha
      double precision beta
      integer bindx(*), bpntrb(*), bpntre(*)
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
c     Loop through block columns:
c
      rcol = lb
      if ( n .lt. rcol ) rcol = n
      rhscols = n/rcol
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      vb = 1 - lbsq
      cb = 1 - lb
      do 5 i=1,mb
c       cb = starting (point) column of block
        cb = cb + lb
        jb = bpntrb(i)
        je = bpntre(i) - 1
        do 10 j=jb,je
          blkind = bindx(j)
c         vb = starting (point) element of A block
          vb = vb + lbsq
c         vc = starting (point) row of A block
          vc = (blkind - 1)*lb + 1
c
c           do the explicit block multiply....
c
            bb = 1 - rcol
            do 15 l=1,rhscols
c              bb = starting (point) column of b 
               bb = bb + rcol
               k = rcol
               if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
               call dgemm('N','N',lb,k,lb,alpha,val(vb),lb,
     *                                           b(cb,bb),ldb,
     *                                      one,c(vc,bb),ldc)
               if (blkind .ne. i) then
c
c           do the implicit block multiply (with transpose of A block
c
                 call dgemm('T','N',lb,k,lb,alpha,val(vb),lb,
     *                                           b(vc,bb),ldb,
     *                                      one,c(cb,bb),ldc)
               endif
 15         continue
 10     continue
 5    continue
             
      return
      end
           
      subroutine DBSCMMKK( trans, mb, n, kb, alpha, val, bindx, bpntrb,
     *                          bpntre, lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, lb, ldb, ldc, base
      integer lbsq
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, k, l, jb, je
      integer vb, vc, cb, bb
      integer rcol, rhscols, blkind
      double precision alpha
      double precision beta
      integer bindx(*), bpntrb(*), bpntre(*)
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
c     Loop through block columns:
c
      rcol = lb
      if ( n .lt. rcol ) rcol = n
      rhscols = n/rcol
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      vb = 1 - lbsq
      cb = 1 - lb
      do 5 i=1,mb
c       cb = starting (point) column of block
        cb = cb + lb
        jb = bpntrb(i)
        je = bpntre(i) - 1
        do 10 j=jb,je
          blkind = bindx(j)
c         vb = starting (point) element of A block
          vb = vb + lbsq
c         vc = starting (point) row of A block
          vc = (blkind - 1)*lb + 1
          if ( blkind .ne. i ) then
c
c           do the explicit block multiply....
c
            bb = 1 - rcol
            do 15 l=1,rhscols
c              bb = starting (point) column of b 
               bb = bb + rcol
               k = rcol
               if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
               if ( trans .eq. 'N' ) then
                  call dgemm('N','N',lb,k,lb,alpha,val(vb),lb,
     *                                            b(cb,bb),ldb,
     *                                        one,c(vc,bb),ldc)
c
c           do the implicit block multiply (with transpose of A block
c
                  call dgemm('T','N',lb,k,lb,-alpha,val(vb),lb,
     *                                             b(vc,bb),ldb,
     *                                          one,c(cb,bb),ldc)
               else
c               Reverse sign on alpha for transpose operation:
                  call dgemm('N','N',lb,k,lb,-alpha,val(vb),lb,
     *                                            b(cb,bb),ldb,
     *                                        one,c(vc,bb),ldc)
c
c           do the implicit block multiply (with transpose of A block
c
                  call dgemm('T','N',lb,k,lb,alpha,val(vb),lb,
     *                                             b(vc,bb),ldb,
     *                                          one,c(cb,bb),ldc)
               endif
 15         continue
          endif
 10     continue
 5    continue
             
      return
      end
           
