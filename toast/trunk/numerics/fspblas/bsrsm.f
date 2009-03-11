c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine dbsrsm( transa, mb, n, unitd, dv, alpha, descra, 
     *           val, bindx, bpntrb, bpntre, lb,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   dbsrsm -- block sparse row format triangular solve
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
c   double val()  scalar array of length nnz containing matrix entries
c                 stored column-major within each dense block.
c  
c   int bindx()   integer array of length bnnz consisting of the
c                 block column indices of the block entries of A.
c  
c   int bpntrb()  integer array of length mb such that 
c                 bpntrb(i)-bpntrb(1) points to location in bindx
c                 of the first block entry of the j-th block row of A.
c  
c   int bpntre()  integer array of length mb such that 
c                 bpntre(i)-bpntre(1) points to location in bindx
c                 of the last block entry of the j-th block row of A.
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
      integer lb
      integer bindx(*), bpntrb(*), bpntre(*)
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
         print *,'Insufficient work space for BSRSM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('BSRSM', info)
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

      call dbsrsmk( transpose, mb, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, bindx, bpntrb, bpntre, lb,
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DBSRSMK( trans, mb, n, 
     *           scale, dvl, dvr, alpha, uplo, diag, 
     *           val, bindx, bpntrb, bpntre, lb, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer mb, n, lb, ldb, ldc, base
      character trans, scale, uplo, diag
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, ii, j, bl, jb, je, vb, db, cb, bb, nb, lbsq, blkind
      logical left, right, unit, nonunit, lower, upper
      integer rcol, rhscols, rhscolb, rhscole, lwork
      double precision alpha
      double precision beta
      integer bindx(*), bpntrb(*), bpntre(*)
      double precision dvl(*)
      double precision dvr(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision work(ldc, *)
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
c
c     Now, loop through the rhscols block columns of c & b:
c
      do 15 bl=1,rhscols
        rhscolb = (bl - 1)*rcol + 1
        rhscole = rhscolb + rcol - 1
        if (rhscole .ge. n ) rhscole = n
        nb = rhscole - rhscolb + 1
        lwork = lb*mb*nb

        if (lower) then
          if (trans .eq. 'N') then
c-----------------------------------------------------------------
c
c         Lower triangular:
c
          vb = 1 - lbsq
          db = 1 - lbsq
          bb = 1 - lb
          do 20 j=1,mb
              bb = bb + lb
              jb = bpntrb(j)
              je = bpntre(j)
              do 25 i=1,nb
                 call dscal(lb,0.D0,work(bb,i),1)   
 25           continue
              do 30 i=jb,je-2
                 vb = vb + lbsq
                 cb = (bindx(i) - 1)*lb + 1
                 call dgemm('N','N', lb, nb, lb, 1.D0, val(vb), lb,
     *                                           work(cb,1), ldc,
     *                                      1.D0, work(bb,1), ldc)
 30           continue
              vb = vb + lbsq
              if (unit .and.  bindx(je-1) .ne. j ) then
                 cb = (bindx(je-1) - 1)*lb + 1
                 call dgemm('N','N', lb, nb, lb, 1.D0, val(vb), lb,
     *                                           work(cb,1), ldc,
     *                                      1.D0, work(bb,1), ldc)
              endif
              if (right) then
                db = db + lbsq             
                call dgemm('N','N', lb, nb, lb, 1.D0, dvr(db),lb,
     *                                          b(bb,rhscolb), ldb,
     *                                   -1.D0, work(bb,1),ldc)
              else
                do 35 ii=1,nb
                  do 35 i=1,lb
                    work(bb-1+i,ii) =  b(bb-1+i,rhscolb-1+ii) - 
     *                                              work(bb-1+i,ii)
 35             continue
              endif
              if (nonunit) then
                call dtrsm('L','L','N','N',lb,nb,1.D0,val(vb),lb,
     *                     work(bb,1),ldc)
              endif
 20        continue
         else
c         Working with the transpose of a lower triangular (== upper)

          if (right) then
c
c           Assign dvr*b to work:
c
            db = 1 - lbsq
            bb = 1 - lb
            do 40 j=1,mb
             db = db + lbsq
             bb = bb + lb
             call dgemm('N','N',lb,nb,lb,1.D0,dvr(db),lb,
     *                                       b(bb,rhscolb),ldc,
     *                                   0.D0, work(bb,1),ldc)
 40         continue
          else
c
c          Copy b into work:
           call dcopy(mb*lb*nb, b(1,rhscolb), 1, work(1,1), 1)
          endif
          bb = mb*lb+1
          do 45 j=mb,1,-1
           bb = bb-lb
           if (nonunit) then
             jb = bpntrb(j)
             je = bpntre(j)-2
             vb = lbsq*je+1
             call dtrsm('L','U','T','N',lb,nb,1.D0,val(vb),lb,
     *                  work(bb,1),ldc)
           else
             jb = bpntrb(j)
             je = bpntre(j)-1
           endif
           vb = lbsq*(jb-2)+1
           do 50 i=jb,je
              vb = vb + lbsq
              blkind = bindx(i)
              if ( blkind .ne. j ) then
                cb = (blkind - 1)*lb + 1
                call dgemm('T','N', lb, nb, lb,-1.D0, val(vb), lb,
     *                     work(bb,1), ldc, 1.D0, work(cb,1), ldc)
              endif
 50        continue
 45       continue
        endif

        else
c-----------------------------------------------------------------
c
c       Upper triangular:
c
        if (trans .eq. 'N') then
          bb = mb*lb+1
          db = mb*lbsq+1
          do 55 j=mb,1,-1
              bb = bb - lb
              jb = bpntrb(j)
              je = bpntre(j)
              do 60 i=1,nb
                 call dscal(lb,0.D0,work(bb,i),1)   
 60           continue
              vb = lbsq*(jb-1)+1
              if (unit .and. bindx(jb) .ne. j) then
                 cb = (bindx(jb) - 1)*lb + 1
                 call dgemm('N','N', lb, nb, lb, 1.D0, val(vb), lb,
     *                                           work(cb,1), ldc,
     *                                      1.D0, work(bb,1), ldc)
              endif
              do 70 i=jb+1,je-1
                 vb = vb + lbsq
                 cb = (bindx(i) - 1)*lb + 1
                 call dgemm('N','N', lb, nb, lb, 1.D0, val(vb), lb,
     *                                           work(cb,1), ldc,
     *                                      1.D0, work(bb,1), ldc)
 70           continue
              if (right) then
                db = db - lbsq             
                call dgemm('N','N', lb, nb, lb, 1.D0, dvr(db),lb,
     *                                          b(bb,rhscolb), ldb,
     *                                     -1.D0, work(bb,1),ldc)
              else
                do 75 ii=1,nb
                  do 75 i=1,lb
                    work(bb-1+i,ii) =  b(bb-1+i,rhscolb-1+ii) - 
     *                                              work(bb-1+i,ii)
 75             continue
              endif
              if (nonunit) then
                vb = lbsq*(jb-1)+1
                call dtrsm('L','U','N','N',lb,nb,1.D0,val(vb),lb,
     *                     work(bb,1),ldc)
              endif
 55        continue
         else
c          Working with the transpose of a upper triangular ( == lower)
           if (right) then
c
c           Assign dvr*b to work:
c
            db = 1 - lbsq
            bb = 1 - lb
            do 80 j=1,mb
             db = db + lbsq
             bb = bb + lb
             call dgemm('N','N',lb,nb,lb,1.D0,dvr(db),lb,
     *                                       b(bb,rhscolb),ldc,
     *                                   0.D0, work(bb,1),ldc)
 80         continue
          else
c
c           Copy b into work:
            call dcopy(mb*lb*nb, b(1,rhscolb), 1, work(1,1), 1)
          endif
      
          vb = 1 - lbsq
          bb = 1 - lb
          do 90 j=1,mb
           bb = bb + lb
           if (nonunit) then
             jb = bpntrb(j)+1
             je = bpntre(j)-1
             vb = vb + lbsq
             call dtrsm('L','L','T','N',lb,nb,1.D0,val(vb),lb,
     *                  work(bb,1),ldc)
           else
             jb = bpntrb(j)
             je = bpntre(j)-1
           endif
           do 100 i=jb,je
              vb = vb + lbsq
              blkind = bindx(i)
              if ( blkind .ne. j ) then
                cb = (blkind - 1)*lb + 1
                call dgemm('T','N', lb, nb, lb,-1.D0, val(vb), lb,
     *                     work(bb,1), ldc, 1.D0, work(cb,1), ldc)
              endif
 100       continue
  90      continue
         endif

c-----------------------------------------------------------------
        endif

        if (left) then
          db = 1 - lbsq
          cb = 1 - lb
          do 105 j=1,mb
            db = db + lbsq
            cb = cb + lb
            call dgemm('N','N',lb,nb,lb,alpha,dvl(db),lb,
     *                                        work(cb,1),ldc,
     *                                  beta, c(cb,rhscolb),ldc)
 105      continue
        else
          call dscal(lwork, beta, c(1,rhscolb), 1)
          call daxpy(lwork, alpha, work(1,1), 1, c(1, rhscolb), 1)
        endif

 15   continue
        
      return
      end

