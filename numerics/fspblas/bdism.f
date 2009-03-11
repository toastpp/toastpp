c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine dbdism( transa, mb, n, unitd, dv, alpha, descra, 
     *           val, blda, ibdiag, nbdiag, lb,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   dbdism --  block diagonal format triangular solve
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
c   int mb	Number of block rows in matrix A
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
      integer transa, mb, n, unitd, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
      double precision dv(*)
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
      integer m, info
      character transpose, scale, uplo, diag
c
c     externals:
c
      external xerbla

      m=mb*lb
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
         print *,'Insufficient work space for BDISM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('BDISM', info)
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

      call dbdismk( transpose, mb, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, blda, ibdiag, nbdiag, lb,
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DBDISMK( trans, mb, n, 
     *           scale, dvr, dvl, alpha, uplo, diag, 
     *           val, blda, ibdiag, nbdiag, lb, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer mb, n, blda, nbdiag, lb, ldb, ldc, base 
      character trans, scale, uplo, diag
      integer rcol, rhscols, rhscolb, rhscole, bl, nb, lwork
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, jb, je, lbsq, cb, bb, db, vb, colind
      logical left, right
      logical unit, nonunit, lower, upper
      double precision alpha
      double precision beta
      integer ibdiag(*)
      double precision val(*)
      double precision dvr(*)
      double precision dvl(*)
      double precision b(ldb,*), c(ldc,*)
      double precision work(ldc,*)
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


      lbsq = lb*lb
      rcol = lb
      if ( rcol .gt. n ) rcol = n
      rhscols = n/rcol
c     if ( mod(n,rcol) .ne. 0 ) rhscols = rhscols + 1
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1

      je = nbdiag - 1
c
c     Now, loop through the rhscols block columns of c & b:
c
      do 15 bl=1,rhscols
        rhscolb = (bl - 1)*rcol + 1
        rhscole = rhscolb + rcol - 1
        if (rhscole .ge. n ) rhscole = n
        nb = rhscole - rhscolb + 1
        lwork = lb*mb*nb
        if (right) then
c
c         Assign dvr*b to work:
c
c
          db = 1 - lbsq
          bb = 1 - lb
          do 20 j=1,mb
             db = db + lbsq
             bb = bb + lb
             call dgemm('N','N',lb,nb,lb,1.D0,dvr(db),lb,
     *                                       b(bb,rhscolb),ldc,
     *                                 0.D0, work(bb,1),ldc)
 20       continue
        else
c
c        Copy b into work:
c
          call dcopy(lwork, b(1,rhscolb), 1, work(1,1), 1)
        endif

        if (lower) then
C-----------------------------------------------------------------
c
c       Lower triangular:
c
        if ( trans .eq. 'N') then
        bb = 1-lb
        do 25 i=1,mb
          bb = bb+lb
          jb = max(nbdiag-i+1,1)
          do 30 j=1,je
             if ( ibdiag(j) .lt. -i+1) go to 30
             colind = (i+ibdiag(j)-1)*lb + 1
             vb = (j-1)*lbsq*mb + (i-1)*lbsq + 1
             call dgemm('N','N', lb, nb, lb,-1.D0, val(vb), lb,
     *                                    work(colind,1), ldc,
     *                                    1.D0, work(bb,1), ldc)
 30       continue
          if (nonunit) then
            vb = (nbdiag-1)*lbsq*mb + (i-1)*lbsq + 1
            call dtrsm('L','L','N','N',lb,nb,1.D0,val(vb),lb,
     *                  work(bb,1),ldc)
          endif
 25     continue

        else

c       Working with the transpose of a lower triangular (== upper)

        do 100 i=mb,1,-1
          cb = (i-1)*lb + 1
          if (nonunit) then
            vb = (nbdiag-1)*lbsq*mb + (i-1)*lbsq + 1
            call dtrsm('L','L','T','N',lb,nb,1.D0,val(vb),lb,
     *                  work(cb,1),ldc)
          endif
          jb = max(nbdiag-i+1,1)
          do 120 j=jb,nbdiag-1
            if (ibdiag(j) .ge. 0 ) go to 120
            bb = (i+ibdiag(j) - 1)*lb + 1
            vb = (j-1)*lbsq*mb + (i-1)*lbsq + 1
            call dgemm('T','N', lb, nb, lb,-1.D0, val(vb), lb,
     *                                  work(cb,1), ldc,
     *                                  1.D0, work(bb,1), ldc)
 120      continue
 100    continue
        endif
      
        else
C-----------------------------------------------------------------
c
c       Upper triangular:
c
        if ( trans .eq. 'N') then
        bb = mb*lb + 1
        do 26 i=mb,1,-1
          bb = bb-lb
          if (nonunit) then
            vb = (i-1)*lbsq + 1
            call dtrsm('L','L','N','N',lb,nb,1.D0,val(vb),lb,
     *                work(bb,1),ldc)
          endif
          do 31 j=2,nbdiag
             if ( ibdiag(j) .gt. mb-i) go to 31
             colind = (i+ibdiag(j)-1)*lb + 1
             vb = (j-1)*lbsq*mb + (i-1)*lbsq + 1
             call dgemm('N','N', lb, nb, lb,-1.D0, val(vb), lb,
     *                                    work(colind,1), ldc,
     *                                    1.D0, work(bb,1), ldc)
 31       continue
 26     continue
 
        else
c       working with the transpose of an upper triangular (== lower)
        do 130 i=1,mb
          cb = (i-1)*lb + 1
          jb = 1
          if (nonunit) then
            vb = (i-1)*lbsq+1
            call dtrsm('L','L','T','N',lb,nb,1.D0,val(vb),lb,
     *                work(cb,1),ldc)
            jb = 2
          endif
          do 140 j=je,nbdiag  
            if (ibdiag(j) .gt. mb-i) go to 140
            bb = (i+ibdiag(j) - 1)*lb + 1
            vb = (j-1)*lbsq*mb + (i-1)*lbsq + 1
            call dgemm('T','N', lb, nb, lb,-1.D0, val(vb), lb,
     *                                  work(cb,1), ldc,
     *                                  1.D0, work(bb,1), ldc)
 140      continue
 130    continue
        endif

C-----------------------------------------------------------------
        endif

        if (left) then
          db = 1 - lbsq
          cb = 1 - lb
          do 35 j=1,mb
             db = db + lbsq
             cb = cb + lb
             call dgemm('N','N',lb,nb,lb,alpha,dvl(db),lb,
     *                                       work(cb,1),ldc,
     *                                 beta, c(cb,rhscolb),ldc)
 35       continue
        else
          call dscal(lwork, beta, c(1,rhscolb), 1)
          call daxpy(lwork, alpha, work(1,1), 1, c(1, rhscolb), 1)
        endif

 15   continue

      return
      end 
