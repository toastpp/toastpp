c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dcoomm( transa, m, n, k, alpha, descra,
     *           val, indx, jndx, nnz,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dcoomm -- compressed sparse row format matrix-matrix multiply
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
c   double val()  scalar array of length nnz containing matrix entries.
c  
c   int indx()    integer array of length nnz containing row indices.
c
c   int jndx()    integer array of length nnz containing column indices.
c
c   int nnz       number of non-zero elements in A.
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
      integer nnz
      integer indx(*), jndx(*)
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
        call xerbla('COOMM', info)
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
         if (transpose .eq. 'N') then
           call dcoommgk( m, n, k, alpha,
     *       val, indx, jndx, nnz,
     *       b, ldb, beta, c, ldc, descra(4))
         else
           call dcoommgk( m, n, k, alpha,
     *       val, jndx, indx, nnz,
     *       b, ldb, beta, c, ldc, descra(4))
         endif
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dcoommsk( m, n, k, alpha,
     *       val, indx, jndx, nnz,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
         if (transpose .eq. 'N') then
          call dcoommkk( m, n, k, alpha,
     *       val, indx, jndx, nnz,
     *       b, ldb, beta, c, ldc, descra(4))
         else
          call dcoommkk( m, n, k, alpha,
     *       val, jndx, indx, nnz,
     *       b, ldb, beta, c, ldc, descra(4))
         endif
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('COOMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DCOOMMGK( m, n, k, alpha, 
     *           val, indx, jndx, nnz,
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, nnz, ldb, ldc, base
      integer l,j
      double precision alpha
      double precision beta
      integer indx(*), jndx(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)

      do 5 l=1,n
         do 10 j=1,nnz
           c(indx(j),l) = c(indx(j),l)+alpha*b(jndx(j),l)*val(j)
 10      continue
 5    continue

      return
      end
         
      subroutine DCOOMMSK( m, n, k, alpha, val, indx, jndx, nnz,
     *                          b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, nnz, ldb, ldc, base
      integer l,j
      double precision alpha
      double precision beta
      integer indx(*), jndx(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
c
c     Scale c by beta:
c
      call dscal( n*m, beta, c(1,1), 1)

      do 5 l=1,n
         do 10 j=1,nnz
           if ( indx(j) .ne. jndx(j) ) then
              c(indx(j),l) = c(indx(j),l)+alpha*b(jndx(j),l)*val(j)
              c(jndx(j),l) = c(jndx(j),l)+alpha*b(indx(j),l)*val(j)
           else
              c(indx(j),l) = c(indx(j),l)+alpha*b(jndx(j),l)*val(j)
           endif
 10      continue
 5    continue

      return
      end
         
      subroutine DCOOMMKK( m, n, k, alpha, val, indx, jndx, nnz,
     *                          b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, nnz, ldb, ldc, base
      integer l,j
      double precision alpha
      double precision beta
      integer indx(*), jndx(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
c
c     Scale c by beta:
c
      call dscal( n*m, beta, c(1,1), 1)

      do 5 l=1,n
         do 10 j=1,nnz
           if ( indx(j) .ne. jndx(j) ) then
             c(indx(j),l) = c(indx(j),l)+alpha*b(jndx(j),l)*val(j)
             c(jndx(j),l) = c(jndx(j),l)-alpha*b(indx(j),l)*val(j)
           endif
 10      continue
 5    continue

      return
      end
         
