c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine dskysm( transa, m, n, unitd, dv, alpha, descra, 
     *           val, pntr,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   dskysm -- Skyline format triangular solve
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
      integer pntr(*)
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

      call dskysmk( transpose, m, n, scale, dv, dv, alpha, 
     *     val, pntr, uplo, diag, 
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DSKYSMK( trans, m, n, 
     *           scale, dvl, dvr, alpha,
     *           val, pntr, uplo, diag, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer m, n, ldb, ldc, base
      character trans, scale, uplo, diag
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, l, vb, len, cindex
      logical left, right, unit, nonunit, lower, upper
      double precision alpha
      double precision beta
      integer pntr(*)
      double precision val(*)
      double precision dvl(*)
      double precision dvr(*)
      double precision work(*)
      double precision b(ldb,*), c(ldc,*)
      double precision ddot
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



      if ( ( lower .and. trans.eq.'N') .or.
     *     ( upper .and. trans.eq.'T')      ) then
c
c         Use dot products...
c
      do 5 l=1,n
        if (right) then
          do 10 i=1,m
            work(i) = dvr(i)*b(i,l)
 10       continue
        else
          call dcopy(m,b(1,l),1,work(1),1)
        endif
        do 15 i=1,m
          vb = pntr(i)
          if (nonunit) then
            len = pntr(i+1) - vb - 1
            cindex = i-len
            work(i) = work(i) - ddot(len,work(cindex),1,val(vb),1)
            work(i) = work(i) / val(pntr(i+1)-1)
          else
            len = pntr(i+1) - vb
            cindex = i-len
            work(i) = work(i) - ddot(len,work(cindex),1,val(vb),1)
          endif
 15     continue
        if (left) then
          do 20 i=1,m
             c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 20       continue
        else
          do 21 i=1,m
             c(i,l) = alpha*work(i) + beta*c(i,l)
 21       continue
        endif
  5   continue
 
      else
c
c         Use daxpys...
c
      do 25 l=1,n
        if (right) then
          do 30 i=1,m
            work(i) = dvr(i)*b(i,l)
 30       continue
        else
          call dcopy(m,b(1,l),1,work(1),1)
        endif
        do 35 i=m,1,-1
          vb = pntr(i)
          if (nonunit) then
            len = pntr(i+1) - vb - 1
            cindex = i-len
            work(i) = work(i) / val(pntr(i+1)-1)
          else
            len = pntr(i+1) - vb
            cindex = i-len
          endif
          call daxpy(len,-work(i),val(vb),1,work(cindex),1)
 35     continue
        if (left) then
          do 40 i=1,m
            c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 40       continue
        else
          do 41 i=1,m
            c(i,l) = alpha*work(i) + beta*c(i,l)
 41       continue
        endif
 25   continue

      endif
 
      return
      end

