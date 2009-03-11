c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dbdimm( transa, mb, n, kb, alpha, descra,
     *           val, blda, ibdiag, nbdiag, lb,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dbdimm -- block diagonal format matrix-matrix multiply
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
c   double val()  scalar array of length lb*lb*blda*nbdiag containing 
c                 matrix entries, stored column-major within each dense 
c                 block.
c  
c   int blda      leading block dimension of val().
c  
c   int ibdiag()  integer array of length nbdiag consisting of the
c                 corresponding indices of the non-zero block
c                 diagonals of A in val().
c  
c   int nbdiag    the number of non-zero block diagonals in A.
c  
c   double b()    rectangular array with first dimension ldb.
c  
c   double c()    rectangular array with first dimension ldc.
c  
c   double work() scratch array of length lwork.  lwork should be at least
c                 max(m,n)
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
      integer blda, nbdiag, lb
      integer ibdiag(*)
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
        call xerbla('BDIMM', info)
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
        call dbdimmgk( transpose,  mb, n, kb, alpha,
     *       val, blda, ibdiag, nbdiag, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dbdimmsk( mb, n, kb, alpha,
     *       val, blda, ibdiag, nbdiag, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call dbdimmkk( transpose,  mb, n, kb, alpha,
     *       val, blda, ibdiag, nbdiag, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('BDIMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DBDIMMGK( trans, mb, n, kb, alpha, 
     *           val, blda, ibdiag, nbdiag, lb, 
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, blda, nbdiag, lb, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, vstart, cstart, len, bstart, diag, lbsq
      double precision alpha
      double precision beta
      integer ibdiag(*)
      double precision val(blda*lb*lb,*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Loop through diagonals:
c
      do 5 i=1,nbdiag
        diag = ibdiag(i)
        if ( diag .le. 0 ) then
           vstart = -diag*lbsq+1-lbsq
           cstart = -diag*lb+1-lb
           bstart = 1-lb
           len    = blda+diag
        else
           vstart = 1-lbsq
           cstart = 1-lb
           bstart = diag*lb+1-lb
           len    = kb-diag
        endif
        do 10 j=1,len
             vstart = vstart+lbsq
             cstart = cstart+lb
             bstart = bstart+lb
             if ( trans .eq. 'N' ) then
               call dgemm('N','N',lb,n,lb,alpha,val(vstart,i),lb,
     *                                            b(bstart,1),ldb,
     *                                        one,c(cstart,1),ldc)
             else
               call dgemm('T','N',lb,n,lb,alpha,val(vstart,i),lb,
     *                                            b(cstart,1),ldb,
     *                                        one,c(bstart,1),ldc)
             endif
 10     continue
 5    continue

      return
      end

      subroutine DBDIMMSK( mb, n, kb, alpha, val, blda, ibdiag,
     *                          nbdiag, lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, blda, nbdiag, lb, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, vstart, cstart, len, bstart, diag, lbsq
      double precision alpha
      double precision beta
      integer ibdiag(*)
      double precision val(blda*lb*lb,*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( n*mb*lb, beta, c(1,1), 1)
c
c     Loop through diagonals:
c
      do 5 i=1,nbdiag
        diag = ibdiag(i)
           vstart = max(-diag*lbsq+1,1)-lbsq
           cstart = max(-diag*lb+1,1)-lb
           bstart = max(diag*lb+1,1)-lb
           len    = blda-abs(diag)
           do 10 j=1,len
             vstart = vstart+lbsq
             cstart = cstart+lb
             bstart = bstart+lb
             call dgemm('N','N',lb,n,lb,alpha,val(vstart,i),lb,
     *                                           b(bstart,1),ldb,
     *                                      one,c(cstart,1),ldc)
             if ( diag .ne. 0 ) then
               call dgemm('T','N',lb,n,lb,alpha,val(vstart,i),lb,
     *                                           b(cstart,1),ldb,
     *                                      one,c(bstart,1),ldc)
             endif
 10        continue
 5    continue

      return
      end

      subroutine DBDIMMKK( trans, mb, n, kb, alpha, val, blda, ibdiag,
     *                          nbdiag, lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, blda, nbdiag, lb, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, vstart, cstart, len, bstart, diag, lbsq
      double precision alpha
      double precision beta
      integer ibdiag(*)
      double precision val(blda*lb*lb,*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( n*mb*lb, beta, c(1,1), 1)
c
c     Loop through diagonals:
c
      do 5 i=1,nbdiag
        diag = ibdiag(i)
        if ( diag .ne. 0 ) then
           vstart = max(-diag*lbsq+1,1)-lbsq
           cstart = max(-diag*lb+1,1)-lb
           bstart = max(diag*lb+1,1)-lb
           len    = blda-abs(diag)
           do 10 j=1,len
             vstart = vstart+lbsq
             cstart = cstart+lb
             bstart = bstart+lb
             if ( trans .eq. 'N') then
               call dgemm('N','N',lb,n,lb,alpha,val(vstart,i),lb,
     *                                            b(bstart,1),ldb,
     *                                        one,c(cstart,1),ldc)
               call dgemm('T','N',lb,n,lb,-alpha,val(vstart,i),lb,
     *                                             b(cstart,1),ldb,
     *                                         one,c(bstart,1),ldc)
             else
c              Reverse the sign on alpha if working with transpose
               call dgemm('N','N',lb,n,lb,-alpha,val(vstart,i),lb,
     *                                            b(bstart,1),ldb,
     *                                        one,c(cstart,1),ldc)
               call dgemm('T','N',lb,n,lb,alpha,val(vstart,i),lb,
     *                                             b(cstart,1),ldb,
     *                                         one,c(bstart,1),ldc)
             endif
 10        continue
        endif
 5    continue

      return
      end

