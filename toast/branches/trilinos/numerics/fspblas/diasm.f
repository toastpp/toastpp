c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine ddiasm( transa, m, n, unitd, dv, alpha, descra, 
     *           val, lda, idiag, ndiag,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   ddiasm -- diagonal format triangular solve
c  
c   C <- alpha D inv(A) B + beta C    C <- alpha D inv(A') B + beta C
c   C <- alpha inv(A) D B + beta C    C <- alpha inv(A') D B + beta C
c   
c                                      ( ' indicates matrix transpose)
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
c   int unitd	Type of scaling:
c                        1 : Identity matrix (argument dv[] is ignored)
c                        2 : Scale on left (row scaling)
c                        3 : Scale on right (column scaling)
c  
c   double alpha	Scalar parameter
c  
c   double beta 	Scalar parameter
c  
c   int descra()	Descriptor argument.  Nine element integer array
c  		descra(0) matrix structure
c  			0 : general
c  			1 : symmetric
c  			2 : Hermitian
c  			3 : Triangular
c  			4 : Skew(Anti-Symmetric
c  			5 : Diagonal
c  		descra(1) upper/lower triangular indicator
c  			1 : lower
c  			2 : upper
c  		descra(2) main diagonal type
c  			0 : non-unit
c  			1 : unit
c  		descra(3) Array base 
c  			0 : C/C++ compatible
c  			1 : Fortran compatible
c  		descra(4) repeated indices?
c  			0 : unknown
c  			1 : no repeated indices
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
      integer transa, m, n, unitd, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
      double precision dv(*)
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
      character transpose, scale, uplo, diag
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
      else if ( (unitd .lt. 1) .or. (unitd .gt. 3) ) then
         info = 4
      else if ( descra(1) .ne. 3 ) then
         info = 6
      else if ( descra(2) .lt. 1 .or. descra(2) .gt. 2) then
         info = 6
      else if ( descra(3) .lt. 0 .or. descra(3) .gt. 1) then
         info = 6
      else if ( ldb .lt. m ) then
         info = 16
      else if (ldc .lt. m) then
         info = 19
      else if (lwork .lt. m ) then
         print *,'Insufficient work space for DIASM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('DIASM', info)
         return
      endif

      if (alpha .eq. 0.0) then
c        Quick return after scaling: 
         call dscal(m*n, beta, c, 1)
         return
      endif

      transpose = 'T'
      if (transa .eq. 0) transpose = 'N'

      if (unitd .eq. 1) then
        scale = 'N'
      else if (unitd .eq. 2) then
        scale = 'L'
      else if (unitd .eq. 3) then
        scale = 'R'
      endif

      uplo = 'U'
      if ( descra(2) .eq. 1 ) uplo = 'L'
      diag = 'U'
      if ( descra(3) .eq. 0 ) diag = 'N'

c
c     Call kernel subroutine:
c

      call ddiasmk( transpose, m, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, lda, idiag, ndiag,
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DDIASMK( trans, m, n, 
     *           scale, dvr, dvl, alpha, uplo, diag, 
     *           val, lda, idiag, ndiag, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer m, n, lda, ndiag, ldb, ldc, base 
      character trans, scale, uplo, diag
      integer rcol, rhscols, rhscolb, rhscole, bl, nb
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, jb, je, index
      logical left, right, unit, nonunit, lower, upper
      double precision alpha
      double precision beta
      double precision t
      integer idiag(*)
      double precision val(lda,*)
      double precision dvr(*)
      double precision dvl(*)
      double precision b(ldb,*), c(ldc,*)
      double precision work(*)
c
c     Set some parameter flags
c
      if (diag .eq. 'U' ) then
        unit =.true.
        nonunit =.false.
      else
        unit =.false.
        nonunit =.true.
      endif
 
      left =.false.
      right =.false.
      if (scale .eq. 'L' ) then
        left =.true.
      else if (scale .eq. 'R' ) then
        right =.true.
      else if (scale .eq. 'B' ) then
        left =.true.
        right =.true.
      endif
 
      lower =.false.
      upper =.false.
      if (uplo .eq. 'L') then
         lower =.true.
      else
         upper =.true.
      endif

c
c     Calculate number of column panels, based on maxcache:
c
c     ndiag*lda + rcol*lda*2 + lda < maxcache
c
c     rcol < ( maxcache - ndiag*lda - lda ) / (2*lda)
c
      rcol = ( maxcache - ndiag*lda -lda ) / (2*lda)
      rhscols = n/rcol
c     if ( mod(n,rcol) .ne. 0 ) rhscols = rhscols + 1
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      je = ndiag - 1
      if (unit .and. idiag(ndiag) .ne. 0 ) je = ndiag
      jb = 2
      if (unit .and. idiag(1) .ne. 0 ) jb = 1
c
c     Now, loop through the rhscols block columns of c & b:
c
      do 10 bl=1,rhscols
        rhscolb = (bl - 1)*rcol + 1
        rhscole = rhscolb + rcol - 1
        if (rhscole .ge. n ) rhscole = n
        nb = rhscole - rhscolb + 1
c
c       Loop through the rcol columns in this block:
c
        do 15 l=rhscolb,rhscole

         if (trans .eq. 'N') then
          if (lower) then
c-----------------------------------------------------------------
c
c         Lower triangular:
c
          call dcopy( lda, c(1,l), 1, work(1), 1)

          do 20 i=1,m
            t = 0.D0
            do 25 j=1,je
              if ( idiag(j) .lt. -i+1 ) go to 25
              t = t + val(i,j)*c(i+idiag(j), l)
 25         continue
            if (right) then
              c(i,l) =  dvr(i)*b(i,l) - t 
            else
              c(i,l) =  b(i,l) - t 
            endif
            if (nonunit) c(i,l) = c(i,l)/val(i,ndiag)
 20       continue
          if (left) then
            do 30 i=1,m
              c(i,l) = alpha*dvl(i)*c(i,l) + beta*work(i)
 30         continue
          else
            do 35 i=1,m
              c(i,l) = alpha*c(i,l) + beta*work(i)
 35         continue
          endif

        else
c-----------------------------------------------------------------
c
c       Upper triangular:
c
          if (right) then
            do 40 i=1,m
              work(i) = dvr(i)*b(i,l)
 40         continue
          else
            call dcopy(m, b(1,l), 1, work(1), 1)
          endif
          do 45 i=m,1,-1
            if (nonunit) work(i) = work(i) / val(i,1)
            do 50 j=jb,ndiag
               index = i-idiag(j)
               work(index) = work(index) - val(index,j)*work(i)
 50         continue
 45      continue
        if (left) then
          do 55 i=1,m
             c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 55       continue
        else
          do 60 i=1,m
             c(i,l) = alpha*work(i) + beta*c(i,l)
 60       continue
        endif
 
c-----------------------------------------------------------------
        endif
       else
c       Working with the transpose:
        if (lower) then
c         transpose of lower (== upper) 
          if (right) then
            do 65 i=1,m
              work(i) = dvr(i)*b(i,l)
 65         continue
          else
            call dcopy(m, b(1,l), 1, work(1), 1)
          endif
          do 70 i=m,1,-1
            if (nonunit) work(i) = work(i)/val(i,ndiag)
            do 75 j=1,je
              if (idiag(j) .lt. -i+1) go to 75
              index = i+idiag(j)
              work(index) = work(index)-val(i,j)*work(i)
 75         continue
 70       continue
          if (left) then
            do 80 i=1,m
              c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 80         continue
          else
            do 85 i=1,m
             c(i,l) = alpha*work(i) + beta*c(i,l)
 85         continue
          endif

        else 
c        transpose of upper (== lower) 
          if (right) then
            do 90 i=1,m
              work(i) = dvr(i)*b(i,l)
 90         continue
          else
            call dcopy(m, b(1,l), 1, work(1), 1)
          endif
          do 100 i=1,m
            if (nonunit) work(i) = work(i)/val(i,1)
            do 105 j=jb,ndiag
               if (idiag(j) .gt. m-i) go to 105
               index = i+idiag(j)
               work(index) = work(index) - val(i,j)*work(i)
 105        continue
 100      continue
          if (left) then
            do 110 i=1,m
              c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 110        continue
          else
            do 115 i=1,m
             c(i,l) = alpha*work(i) + beta*c(i,l)
 115        continue
          endif
        endif
       endif
 15    continue
 10   continue

      return
      end 
