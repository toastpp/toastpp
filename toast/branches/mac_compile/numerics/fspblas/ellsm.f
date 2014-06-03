c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine dellsm( transa, m, n, unitd, dv, alpha, descra, 
     *           val, indx, lda, maxnz,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   dellsm -- Ellpack format triangular solve
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
c   double val()  two-dimensional lda-by-maxnz array such that val(i,:)
c                 consists of non-zero elements in row i of A, padded by 
c                 zero values if the row contains less than maxnz.
c  
c   int indx()    two-dimensional integer blda-by-maxbnz array such 
c                 indx(i,:) consists of the column indices of the 
c                 nonzero elements in row i, padded by the integer 
c                 value i if the number of nonzeros is less than maxnz.
c  
c   int lda       leading dimension of val and indx.
c  
c   int maxnz     max number of nonzeros elements per row.
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
      integer lda, maxnz
      integer indx(*)
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
         print *,'Insufficient work space for ELLSM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('ELLSM', info)
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

      call dellsmk( transpose, m, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, indx, lda, maxnz,
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DELLSMK( trans, m, n, 
     *           scale, dvr, dvl, alpha, uplo, diag, 
     *           val, indx, lda, maxnz,
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer m, n, lda, maxnz, ldb, ldc, base
      character trans, scale, uplo, diag
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer brow, blkrows, i, j, jb, l, rcol, rhscols, je, ind
      logical left, right, unit, nonunit, lower, upper
      integer rhscolb, rhscole, bl
      double precision alpha
      double precision beta
      integer indx(lda,*)
      double precision val(lda,*)
      double precision dvr(*)
      double precision dvl(*)
      double precision work(*)
      double precision b(ldb,*), c(ldc,*)
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

      brow = cacheline
      blkrows = m/brow
      if (blkrows*brow .ne. m ) blkrows = blkrows + 1
c
c     Calculate the number of columns for partitioning
c     b and c (rcol) by taking into acount the amount of 
c     cache required for calculation:
c     For rcol columns in block:
c       cachereq  = cacheline*maxnz + rcol*(1+maxnz) < maxcache
c            from val ---^      from c   ---^   ^--- from b
c
c     So,   rcol = (maxcache - cacheline*maxnz)/(1+maxnz)
      rcol = (maxcache - cacheline*maxnz)/(1+maxnz)
      rhscols = n/rcol
      if (rhscols*rcol .ne. n ) rhscols = rhscols+1

      jb = 2
      if (unit .and. indx(1,1) .ne. 1) jb = 1
c
c     Loop through the rhscols block columns of b & c
c
      rhscolb = 1-rcol
      do 10 bl=1,rhscols
        rhscolb = rhscolb + rcol
        rhscole = rhscolb + rcol - 1
        if ( rhscole .gt. n ) rhscole = n
c
c       Loop through the rcol columns in this block
c
        do 15 l=rhscolb, rhscole
          if (right) then
c
c           store dvr*b in work... work will contain the 
c           intermediate results of the triangular solve...
c
            do 20 i=1,m
              work(i) = dvr(i)*b(i,l)
 20         continue
          else
            call dcopy(m, b(1,l), 1, work(1), 1)
          endif
 
        if (trans .eq. 'N') then
          if (lower) then
c-----------------------------------------------------------------
c
c           Lower triangular:
c
            do 25 i=1,m
             do 30 j=1,maxnz
               if ( indx(i,j) .eq. i ) go to 35
               work(i) = work(i) - val(i,j)*work(indx(i,j))
 30          continue   
 35          if (nonunit) work(i) = work(i)/val(i,j)
 25         continue

          else
c-----------------------------------------------------------------
c
c         Upper triangular:
c
            do 40 i=m,1,-1
             j = jb
             do 45 j=jb,maxnz
               if ( indx(i,j) .le. i ) go to 50
               work(i) = work(i) - val(i,j)*work(indx(i,j))
 45          continue
 50          continue
             work(i) = work(i)/val(i,1)
 40         continue

c-----------------------------------------------------------------
          endif
       else
c       Working with the tranpose:

        if (lower) then
           do 55 i=m,1,-1
c            Find element in row i with indx i
             do 60 j=1,maxnz
                if (indx(i,j) .eq. i ) then
                  je = j-1
                  if (nonunit) work(i) = work(i)/val(i,j)
                  go to 61
                endif
 60          continue
 61          do 65 j=1,je
                work(indx(i,j)) = work(indx(i,j)) - val(i,j)*work(i)
 65          continue
 55        continue
         else
           do 70 i=1,m
              if (nonunit) work(i) = work(i)/val(i,1)
              do 75 j = jb,maxnz
                ind = indx(i,j)
                if (ind .eq. i ) go to 70
                work(ind) = work(ind) - val(i,j)*work(i)
 75           continue
 70        continue
        endif
       endif

          if (left) then
            do 80 i=1,m
               c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 80         continue
          else
            do 85 i=1,m
               c(i,l) = alpha*work(i) + beta*c(i,l)
 85         continue
          endif
 15     continue
 10   continue
               
      return
      end

