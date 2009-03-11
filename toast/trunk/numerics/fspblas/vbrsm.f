c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine dvbrsm( transa, mb, n, unitd, dv, alpha, descra, 
     *           val, indx, bindx, rpntr, cpntr, bpntrb, bpntre,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   dvbrsm -- variable block sparse row format triangular solve
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
c  
c   double *val	scalar array of length nnz containing matrix entries
c  
c   int *indx	integer array of length bnnz+1 such that the i-th element of
c   		indx[] points to the location in val of the (1,1) element
c  		of the i-th block entry.
c  
c   int *bindx	integer array of length bnnz consisting of the block column
c   		indices of the entries of A.
c  
c   int *rpntr	integer array of length mb+1 such that rpntr(i)-rpntr(1)
c  		is the row index of the first point row in the i-th block row. 
c  		rpntr(mb+1) is set to m+rpntr(1).
c  		Thus, the number of point rows in the i-th block row is
c  		rpntr(i+1)-rpntr(i).
c  
c   int *cpntr	integer array of length mb+1 such that cpntr(j)-cpntr(1)
c  		is the column index of the first point column in the j-th
c  		block column. cpntr(mb+1) is set to k+cpntr(1).
c  		Thus, the number of point columns in the j-th block column is
c  		cpntr(j+1)-cpntr(j).
c  
c   int *bpntrb	integer array of length mb such that bpntrb(i)-bpntrb(1)
c                points to location in bindx of the first block entry of 
c  		the j-th row of A.
c  
c   int *bpntre	integer array of length mb such that bpntre(i)-bpntrb(1)
c                points to location in bindx of the last block entry of
c  		the j-th row of A.
c  
c   double *b	rectangular array with first dimension ldb
c  
c   double *c	rectangular array with first dimension ldc
c  
c   double *work	scratch array of length lwork.  
c                lwork must be at least (n*m + max_block_size)
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
      integer indx(*), bindx(*)
      integer rpntr(*), cpntr(*), bpntrb(*), bpntre(*)
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

      m=rpntr(mb+1)-rpntr(1)
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
         print *,'Insufficient work space for VBRSM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('VBRSM', info)
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

      call dvbrsmk( transpose, mb, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, indx, bindx, rpntr, cpntr, bpntrb, bpntre, 
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DVBRSMK( trans, mb, n, 
     *           scale, dvl, dvr, alpha, uplo, diag, 
     *           val, indx, bindx, rpntr, cpntr, bpntrb, bpntre, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer mb, n, ldb, ldc, base
      character trans, scale, uplo, diag
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, ii, jj, bl,jb,je,vb,db,cb,bb,nb, bi, dblast, lwork
      logical left, right, unit, nonunit, lower, upper
      integer rcol, rhscols, rhscolb, rhscole, maxcol, rows, cols
      double precision alpha
      double precision beta
      integer indx(*), bindx(*)
      integer rpntr(*), cpntr(*), bpntrb(*), bpntre(*)
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
c
c     Find max number of rows in block:
c     and last block starting index of scaling matrix:
c
      dblast = 1
      maxcol = 0
      do 5 i=1,mb
         rows = rpntr(i+1) - rpntr(i)
         dblast = dblast + rows*rows
         if ( rows .gt. maxcol )
     *       maxcol =  rows
 5    continue


      rcol = maxcol
      if ( rcol .gt. n ) rcol = n
      rhscols = n/rcol
c     if ( mod(n,rcol) .ne. 0 ) rhscols = rhscols + 1
      if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c
c     Now, loop through the rhscols block columns of c & b:
c
      do 10 bl=1,rhscols
        rhscolb = (bl - 1)*rcol + 1
        rhscole = rhscolb + rcol - 1
        if (rhscole .ge. n ) rhscole = n
        nb = rhscole - rhscolb + 1
        lwork = ldc*nb

        if (lower) then
c-----------------------------------------------------------------
c
c       Lower triangular:
c
        if (trans .eq. 'N') then
           db = 1
           do 15 j=1,mb
             rows = rpntr(j+1) - rpntr(j)
c            bb = starting (point) row of block
             bb = rpntr(j)
             jb = bpntrb(j)
             if (nonunit) then
              je = bpntre(j)-2
             else
              je = bpntre(j)-1
             endif
             do 20 i=1,nb
              call dscal(rows,0.D0,work(bb,i),1)
 20          continue
             do 25 i=jb,je
               bi = bindx(i)
               cols = cpntr(bi+1) - cpntr(bi)
               if (bi .ne. j ) then
                 cb = cpntr(bi)
                 vb = indx(i)
                 call dgemm('N','N', rows, nb, cols, 1.D0, val(vb), 
     *                     rows, work(cb,1), ldc,
     *                     1.D0, work(bb,1), ldc)
               endif
 25          continue
             if (right) then
              call dgemm('N','N', rows, nb, rows, 1.D0, dvr(db),rows,
     *                                        b(bb,rhscolb), ldb,
     *                                   -1.D0, work(bb,1),ldc)
             else 
               do 30 jj=1,nb
                 do 30 ii=1,rows
                  work(bb-1+ii,jj) = b(bb-1+ii,rhscolb-1+jj) - 
     *                                           work(bb-1+ii,jj)
 30            continue
             endif
             
             db = db + rows*rows
             if (nonunit) then
               vb = indx(je+1)
               call dtrsm('L','L','N','N',rows,nb,1.D0,val(vb),rows,
     *                     work(bb,1),ldc)
             endif
 15       continue


        else

c         Working with a transpose lower triangular ( == upper)
        if (right) then
c
c         Assign dvr*b to work:
c
          rows = rpntr(2)-rpntr(1)
          db = 1
          bb = 1
          do 40 j=1,mb
             rows = rpntr(j+1)-rpntr(j)
             call dgemm('N','N',rows,nb,rows,1.D0,dvr(db),rows,
     *                                       b(bb,rhscolb),ldc,
     *                                   0.D0, work(bb,1),ldc)
             db = db + rows*rows
             bb = bb + rows
 40       continue
        else
c
c         Copy b into work:
          call dcopy(ldb*rcol,b(rhscolb,1),1,work(1,1),1)
        endif


          do 45 j=mb,1,-1
            cols = rpntr(j+1) - rpntr(j)
c           bb = starting (point) row of block
            bb = cpntr(j)
            jb = bpntrb(j)
            if (nonunit) then
              je = bpntre(j)-2
            else
              je = bpntre(j)-1
            endif
            if (nonunit) then
              vb = indx(je+1)
              call dtrsm('L','L','T','N',cols,nb,1.D0,val(vb),
     *             cols, work(bb,1),ldc)
            endif
            do 50 i=jb,je
             bi = bindx(i)
             if ( bi .eq. j ) go to 50
             cb = cpntr(bi)
             rows = cpntr(bi+1) - cb
             vb = indx(i)
             call dgemm('T','N', rows, nb, cols, -1.D0, val(vb), 
     *                   cols, work(bb,1), ldc,
     *                   1.D0, work(cb,1), ldc)
 50         continue
 45       continue

        endif

        else
c-----------------------------------------------------------------
c
c       Upper triangular:
c
        if ( trans .eq. 'N') then
        db = dblast
        do 55 j=mb,1,-1
           rows = rpntr(j+1) - rpntr(j)
c          bb = starting (point) row of block
           bb = rpntr(j)
           if (nonunit) then
             jb = bpntrb(j) + 1
           else
             jb = bpntrb(j)
           endif
           je = bpntre(j)-1
           do 60 i=1,nb
              call dscal(rows,0.D0,work(bb,i),1)
 60        continue
           do 65 i=jb,je
             bi = bindx(i)
             cols = cpntr(bi+1) - cpntr(bi)
             if ( bi .ne. j ) then
               cb = cpntr(bi)
               vb = indx(i)
               call dgemm('N','N', rows, nb, cols, 1.D0, val(vb), 
     *                    rows, work(cb,1), ldc,
     *                    1.D0, work(bb,1), ldc)
             endif
 65        continue
           if (right) then
             db = db - rows*rows
             call dgemm('N','N', rows, nb, rows, 1.D0, dvr(db),rows,
     *                                        b(bb,rhscolb), ldb,
     *                                   -1.D0, work(bb,1),ldc)
           else
             do 70 jj=1,nb
               do 70 ii=1,rows
                 work(bb-1+ii,jj) = b(bb-1+ii,rhscolb-1+jj) - 
     *                                           work(bb-1+ii,jj)
 70          continue
           endif

           if (nonunit) then
               vb = indx(jb-1)
               call dtrsm('L','U','N','N',rows,nb,1.D0,val(vb),
     *                  rows, work(bb,1),ldc)
           endif
 55     continue

        else
c
c       Working with the tranpose of an upper triangular ( == lower)
c
        if (right) then
c
c         Assign dvr*b to work:
c
          rows = rpntr(2)-rpntr(1)
          db = 1
          bb = 1
          do 75 j=1,mb
             rows = rpntr(j+1)-rpntr(j)
             call dgemm('N','N',rows,nb,rows,1.D0,dvr(db),rows,
     *                                       b(bb,rhscolb),ldc,
     *                                   0.D0, work(bb,1),ldc)
             db = db + rows*rows
             bb = bb + rows
 75       continue
        else
c
c         Copy b into work:
          call dcopy(ldb*rcol,b(rhscolb,1),1,work(1,1),1)
        endif

          do 80 j=1,mb
            cols = rpntr(j+1) - rpntr(j)
c           bb = starting (point) row of block
            bb = cpntr(j)
            je = bpntre(j)-1
            if (nonunit) then
              jb = bpntrb(j)+1
            else
              jb = bpntrb(j)
            endif
            if (nonunit) then
              vb = indx(jb-1)
              call dtrsm('L','L','T','N',cols,nb,1.D0,val(vb),
     *             cols, work(bb,1),ldc)
            endif
            do 85 i=jb,je
             bi = bindx(i)
             if (bi.eq.j) go to 85
             cb = cpntr(bi)
             rows = cpntr(bi+1) - cb
             vb = indx(i)
             call dgemm('T','N', rows, nb, cols, -1.D0, val(vb), 
     *                   cols, work(bb,1), ldc,
     *                   1.D0, work(cb,1), ldc)
 85         continue
 80       continue
         endif
c-----------------------------------------------------------------
        endif

        if (left) then
          db = 1
          do 90 j=1,mb
             rows = rpntr(j+1) - rpntr(j)
             cb = rpntr(j)
             call dgemm('N','N',rows,nb,rows,alpha,dvl(db),rows,
     *                                       work(cb,1),ldc,
     *                                 beta, c(cb,rhscolb),ldc)
             db = db + rows*rows
 90       continue
        else
          call dscal(lwork, beta, c(1,rhscolb), 1)
          call daxpy(lwork, alpha, work(1,1), 1, c(1, rhscolb), 1)
        endif

 10   continue
        
      return
      end
