c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine ddiamm( transa, m, n, k, alpha, descra,
     *           val, lda, idiag, ndiag,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   ddiamm -- diagonal format matrix-matrix multiply
c  
c   C <- alpha A B + beta C
c  
c   Arguments:
c  
c   int transa	Indicates how to operate with the sparse matrix
c  		0 : operate with matrix
c  		1 : operate with transpose matrix
c  
c   int m	Number of rows in matrix A
c  
c   int n	Number of columns in matrix c
c  
c   int k	Number of columns in matrix A
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
c
c   double val()  two-dimensional lda-by-ndiag array such that val(:,i)
c                 consists of non-zero elements on diagonal idiag(i)
c                 of A.  Diagonals in the lower triangular part of A
c                 are padded from the top, and those in the upper
c                 triangular part are padded from the bottom.
c  
c   int lda       leading dimension of val, must be .ge. min(m,k)
c  
c   int idiag()   integer array of length ndiag consisting of the
c                 corresponding diagonal offsets of the non-zero 
c                 diagonals of A in val.  Lower triangular diagonals 
c                 have negative offsets, the main diagonal has offset
c                 0, and upper triangular diagonals have positive offset. 
c
c   int ndiag     number of non-zero diagonals in A.
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
      integer transa, m, n, k, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
      double precision b(*), c(*)
      double precision work(*)
c
c     format specific interface variables:
c
      integer lda, ndiag
      integer idiag(*)
      double precision val(*)
c
c     local variables:
c
      integer info
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
      else if ( m .lt. 0 ) then
         info = 2
      else if (n .lt. 0) then
         info = 3
      else if (k .lt. 0) then
         info = 4
      else
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
        call xerbla('DIAMM', info)
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
        call ddiammgk( transpose,  m, n, k, alpha,
     *       val, lda, idiag, ndiag,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call ddiammsk( m, n, k, alpha,
     *       val, lda, idiag, ndiag,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call ddiammkk( transpose,  m, n, k, alpha,
     *       val, lda, idiag, ndiag,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('DIAMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DDIAMMGK( trans, m, n, k, alpha, 
     *           val, lda, idiag, ndiag, 
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, lda, ndiag, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, vstart, len, bstart, diag
      logical overdet
      double precision alpha
      double precision beta
      integer idiag(*)
      double precision val(lda,*)
      double precision b(ldb,*), c(ldc,*)
      overdet=.false.
      if (m-k .gt. 0) overdet =.true.
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Loop through diagonals:
c
      do 5 i=1,ndiag
        diag = idiag(i)
        vstart = max(-diag,0)
        bstart = max(diag,0)
        len = lda
        if (diag .lt. 0) then
          if (overdet) then
            len = len-max(-diag-(m-k),0)
          else
            len = len+diag
          endif
        elseif (diag .gt. 0) then
          if (overdet) then
            len = len-diag
          else
            len = len-max(diag-(k-m),0)
          endif
        endif
        if ( trans .eq. 'N' ) then
          do 10 l=1,n
            do 15 j=1,len
              c(vstart+j,l) = c(vstart+j,l) + alpha*b(bstart+j,l)*
     *                                     val(vstart+j,i)
 15         continue
 10       continue
        else
          do 20 l=1,n
            do 25 j=1,len
              c(bstart+j,l) = c(bstart+j,l) + alpha*b(vstart+j,l)*
     *                                     val(vstart+j,i)
 25         continue
 20       continue
        endif
 5    continue

      return
      end

      subroutine DDIAMMSK( m, n, k, alpha, val, lda, idiag,
     *                          ndiag, b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, lda, ndiag, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, vstart, len, bstart, diag
      double precision alpha
      double precision beta
      integer idiag(*)
      double precision val(lda,*)
      double precision b(ldb,*), c(ldc,*)
c
c     Scale c by beta:
c
      call dscal( n*m, beta, c(1,1), 1)
c
c     Loop through diagonals:
c
      do 5 i=1,ndiag
        diag = idiag(i)
        vstart = max(-diag,0)
        bstart = max(diag,0)
        len    = lda-abs(diag)
        if ( diag .eq. 0 ) then
           do 10 l=1,n
             do 15 j=1,len
              c(vstart+j,l) = c(vstart+j,l) + alpha*b(bstart+j,l)*
     *                                     val(vstart+j,i)
 15          continue
 10        continue
        else 
           do 20 l=1,n
             do 25 j=1,len
              c(vstart+j,l) = c(vstart+j,l) + alpha*b(bstart+j,l)*
     *                                     val(vstart+j,i)
              c(bstart+j,l) = c(bstart+j,l) + alpha*b(vstart+j,l)*
     *                                     val(vstart+j,i)
 25          continue
 20        continue
        endif
 5    continue

      return
      end

      subroutine DDIAMMKK( trans, m, n, k, alpha, val, lda, idiag,
     *                          ndiag, b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, lda, ndiag, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, vstart, len, bstart, diag
      double precision alpha
      double precision beta
      integer idiag(*)
      double precision val(lda,*)
      double precision b(ldb,*), c(ldc,*)
c
c     Scale c by beta:
c
      call dscal( n*m, beta, c(1,1), 1)
c
c     Loop through diagonals:
c
      do 5 i=1,ndiag
        diag = idiag(i)
        vstart = max(-diag,0)
        bstart = max(diag,0)
        len    = lda-abs(diag)
        if ( diag .ne. 0 ) then
         if (trans .eq. 'N' ) then
           do 20 l=1,n
             do 25 j=1,len
              c(vstart+j,l) = c(vstart+j,l) + alpha*b(bstart+j,l)*
     *                                     val(vstart+j,i)
              c(bstart+j,l) = c(bstart+j,l) - alpha*b(vstart+j,l)*
     *                                     val(vstart+j,i)
 25          continue
 20        continue
         else
c          Reverse the sign on alpha if working with the transpose
           do 30 l=1,n
             do 35 j=1,len
              c(vstart+j,l) = c(vstart+j,l) - alpha*b(bstart+j,l)*
     *                                     val(vstart+j,i)
              c(bstart+j,l) = c(bstart+j,l) + alpha*b(vstart+j,l)*
     *                                     val(vstart+j,i)
 35          continue
 30        continue
         endif
        endif
 5    continue

      return
      end
