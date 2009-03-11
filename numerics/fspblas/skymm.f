c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dskymm( transa, m, n, k, alpha, descra,
     *           val, pntr, 
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dskymm -- Skyline format matrix-matrix multiply
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
c   double val()  array contain the nonzeros of A in skyline profile form.
c                 Row-oriented if descra(2) = 1 (lower triangular), 
c                 column oriented if descra(2) = 2 (upper triangular).
c  
c   int pntr()    integer array of length m+1 (lower triangular) or
c                 k+1 (upper triangular) such that pntr(i) - pntr(1) + 1
c                 points to the location in val of the first element of
c                 the skyline profile in row (column) i.
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
      integer pntr(*)
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
        call xerbla('SKYMM', info)
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

      if (descra(1) .eq. 0) then 
c
c        General matrix multiply:
c
        info = 6
        call xerbla('SKYMM', info)
      else if (descra(1) .eq. 3   .or.
     *         descra(1) .eq. 5        ) then
        if (descra(2) .eq. 1) then
          call dskymmgk( transpose, m, n, k, alpha,
     *       val, pntr, 'L',
     *       b, ldb, beta, c, ldc, descra(4))
        else if (descra(2) .eq. 2) then
          call dskymmgk( transpose, m, n, k, alpha,
     *       val, pntr, 'U',
     *       b, ldb, beta, c, ldc, descra(4))
        else
          info = 6
        endif
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        if (descra(2) .eq. 1) then
          call dskymmsk( m, n, k, alpha,
     *       val, pntr, 'L',
     *       b, ldb, beta, c, ldc, descra(4))
        else if (descra(2) .eq. 2) then
          call dskymmsk( m, n, k, alpha,
     *       val, pntr, 'U',
     *       b, ldb, beta, c, ldc, descra(4))
        else
          info = 6
        endif
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        if (descra(2) .eq. 1) then
          call dskymmkk( transpose, m, n, k, alpha,
     *       val, pntr, 'L',
     *       b, ldb, beta, c, ldc, descra(4))
        else if (descra(2) .eq. 2) then
          call dskymmkk( transpose, m, n, k, alpha,
     *       val, pntr, 'U',
     *       b, ldb, beta, c, ldc, descra(4))
        else
          info = 6
        endif
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('SKYMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DSKYMMGK( trans, m, n, k, alpha,
     *           val, pntr, uplo, 
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, l, vb, len, bindex, cindex
      character trans, uplo
      double precision alpha
      double precision beta
      integer pntr(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision t
      double precision ddot

      if ( (trans .eq. 'N' .and. uplo .eq. 'L') .or. 
     *     (trans .eq. 'T' .and. uplo .eq. 'U')      ) then
c
c         Use dot products...
c
        do 5 i=1,m
          vb = pntr(i)
          len = pntr(i+1) - vb
          bindex = i-len+1
          do 10 l=1,n
             t = ddot(len,b(bindex,l),1,val(vb),1)
             c(i,l) = beta*c(i,l) + alpha*t
 10       continue
 5      continue
 
      else
c
c         Use daxpys...
c
c
c       Scale c by beta:
c
        call dscal( n*ldc, beta, c(1,1), 1)

        do 15 i=1,k
          vb = pntr(i)
          len = pntr(i+1) - vb
          cindex = i-len+1
          do 20 l=1,n
            call daxpy(len,alpha*b(i,l),val(vb),1,c(cindex,l),1)
 20       continue
 15     continue

      endif
 
      return
      end

      subroutine DSKYMMSK( m, n, k, alpha, val, pntr, uplo, 
     *                     b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, l, vb, len, bindex, cindex
      character uplo
      double precision alpha
      double precision beta
      integer pntr(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision t
      double precision ddot
c
c       Scale c by beta:
c
        call dscal( n*ldc, beta, c(1,1), 1)

      if ( uplo .eq. 'L' ) then
        do 5 i=1,m
          vb = pntr(i)
          len = pntr(i+1) - vb
          bindex = i-len+1
          do 10 l=1,n
            t = ddot(len,b(bindex,l),1,val(vb),1)
            c(i,l) = c(i,l) + alpha*t
            call daxpy(len-1,alpha*b(i,l),val(vb),1,c(bindex,l),1)
 10       continue
 5      continue
 
      else
        do 15 i=1,k
          vb = pntr(i)
          len = pntr(i+1) - vb
          cindex = i-len+1
          do 20 l=1,n
            t = ddot(len,b(cindex,l),1,val(vb),1)
            c(i,l) = c(i,l) + alpha*t
            call daxpy(len-1,alpha*b(i,l),val(vb),1,c(cindex,l),1)
 20       continue
 15     continue

      endif
 
      return
      end

      subroutine DSKYMMKK( trans, m, n, k, alpha, val, pntr, uplo, 
     *                     b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, l, vb, len, bindex, cindex
      character trans, uplo
      double precision alpha
      double precision beta
      integer pntr(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision t
      double precision ddot
c
c       Scale c by beta:
c
        call dscal( n*ldc, beta, c(1,1), 1)

      if ( (trans .eq. 'N' .and. uplo .eq. 'L') .or. 
     *     (trans .eq. 'T' .and. uplo .eq. 'U')      ) then
        do 5 i=1,m
          vb = pntr(i)
          len = pntr(i+1) - vb
          bindex = i-len+1
          do 10 l=1,n
            t = ddot(len-1,b(bindex,l),1,val(vb),1)
            c(i,l) = c(i,l) + alpha*t
            call daxpy(len-1,-alpha*b(i,l),val(vb),1,c(bindex,l),1)
 10       continue
 5      continue
 
      else
        do 15 i=1,k
          vb = pntr(i)
          len = pntr(i+1) - vb
          cindex = i-len+1
          do 20 l=1,n
            t = ddot(len-1,b(cindex,l),1,val(vb),1)
            c(i,l) = c(i,l) - alpha*t
            call daxpy(len-1,alpha*b(i,l),val(vb),1,c(cindex,l),1)
 20       continue
 15     continue

      endif
 
      return
      end

