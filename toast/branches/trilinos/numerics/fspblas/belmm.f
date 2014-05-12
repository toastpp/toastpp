c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dbelmm( transa, mb, n, kb, alpha, descra,
     *           val, bindx, blda, maxbnz, lb,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dbelmm -- block Ellpack format matrix-matrix multiply
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
c   double val()  scalar array of length lb*lb*blda*maxbnz containing 
c                 matrix entries, stored column-major within each dense 
c                 block.
c  
c   int bindx()   two-dimensional integer blda-by-maxbnz array such 
c                 bindx(i,:) consists of the block column indices of the 
c                 nonzero blocks in block row i, padded by the integer 
c                 value i if the number of nonzero blocks is less than 
c                 maxbnz.
c  
c   int blda      leading dimension of bindx(:,:).
c  
c   int maxbnz    max number of nonzeros blocks per row.
c  
c   int lb        row and column dimension of the dense blocks composing 
c                 val.
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
      integer blda, maxbnz, lb
      integer bindx(*)
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
        call xerbla('BELMM', info)
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
        call dbelmmgk( transpose,  mb, n, kb, alpha,
     *       val, bindx, blda, maxbnz, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dbelmmsk( mb, n, kb, alpha,
     *       val, bindx, blda, maxbnz, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call dbelmmkk( transpose,  mb, n, kb, alpha,
     *       val, bindx, blda, maxbnz, lb,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('BELMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DBELMMGK( trans, mb, n, kb, alpha, 
     *                     val, bindx, blda, maxbnz, lb,
     *                     b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, blda, maxbnz, lb, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, lbsq, cstart, vstart, bstart
      double precision alpha
      double precision beta
      integer bindx(blda,*)
      double precision val(blda*lb*lb,*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( ldc*n, beta, c(1,1), 1)
c
c     Loop through block columns:
c
      do 5 j=1,maxbnz
c
c       Loop through the block rows in this column
c
        cstart = 1-lb
        vstart = 1-lbsq
        do 10 i=1,blda
          vstart = vstart+lbsq
          cstart = cstart+lb
          bstart = (bindx(i,j)-1)*lb+1
          if ( trans .eq. 'N' ) then
            call dgemm('N','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                         b(bstart,1),ldb,
     *                                      one,c(cstart,1),ldc)
          else
            call dgemm('T','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                         b(cstart,1),ldb,
     *                                      one,c(bstart,1),ldc)
          endif
 10     continue
 5    continue
               
      return
      end

      subroutine DBELMMSK( mb, n, kb, alpha, 
     *           val, bindx, blda, maxbnz, lb,   
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, blda, maxbnz, lb, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, lbsq, cstart, vstart, bstart, blkind
      double precision alpha
      double precision beta
      integer bindx(blda,*)
      double precision val(blda*lb*lb,*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( ldc*n, beta, c(1,1), 1)
c
c     Loop through block columns:
c
      do 5 j=1,maxbnz
c
c       Loop through the block rows in this column
c
        cstart = 1-lb
        vstart = 1-lbsq
        do 10 i=1,blda
          blkind = bindx(i,j)
          vstart = vstart+lbsq
          cstart = cstart+lb
          bstart = (blkind-1)*lb+1
          if ( blkind .eq. i ) then
             call dgemm('N','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                           b(bstart,1),ldb,
     *                                      one,c(cstart,1),ldc)
          else
             call dgemm('N','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                           b(bstart,1),ldb,
     *                                      one,c(cstart,1),ldc)
             call dgemm('T','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                           b(cstart,1),ldb,
     *                                      one,c(bstart,1),ldc)
          endif
 10     continue
 5    continue
               
      return
      end

      subroutine DBELMMKK( trans, mb, n, kb, alpha, 
     *                     val, bindx, blda, maxbnz,
     *                     lb, b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, blda, maxbnz, lb, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, lbsq, cstart, vstart, bstart, blkind
      double precision alpha
      double precision beta
      integer bindx(blda,*)
      double precision val(blda*lb*lb,*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.d0
      lbsq = lb*lb
c
c     Scale c by beta:
c
      call dscal( ldc*n, beta, c(1,1), 1)
c
c     Loop through block columns:
c
      do 5 j=1,maxbnz
c
c       Loop through the block rows in this column
c
        cstart = 1-lb
        vstart = 1-lbsq
        do 10 i=1,blda
          blkind = bindx(i,j)
          vstart = vstart+lbsq
          cstart = cstart+lb
          bstart = (blkind-1)*lb+1
          if ( blkind .ne. i ) then
           if (trans .eq. 'N' ) then
             call dgemm('N','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                           b(bstart,1),ldb,
     *                                      one,c(cstart,1),ldc)
             call dgemm('T','N',lb,n,lb,-alpha,val(vstart,j),lb,
     *                                           b(cstart,1),ldb,
     *                                      one,c(bstart,1),ldc)
           else
c            Reverse sign if working with the transpose
             call dgemm('N','N',lb,n,lb,-alpha,val(vstart,j),lb,
     *                                           b(bstart,1),ldb,
     *                                      one,c(cstart,1),ldc)
             call dgemm('T','N',lb,n,lb,alpha,val(vstart,j),lb,
     *                                           b(cstart,1),ldb,
     *                                      one,c(bstart,1),ldc)

           endif
          endif
 10     continue
 5    continue
               
      return
      end

