c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine djadsm( transa, m, n, unitd, dv, alpha, descra, 
     *           val, indx, pntr, maxnz, iperm,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   djadsm -- Jagged-diagonal format triangular solve
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
c   double val()  array of length nnz consisting of entries of A.
c                 val can be viewed as a column major ordering of a    
c                 row permutation of the Ellpack representation of A, 
c                 where the Ellpack representation is permuted so that
c                 the rows are non-increasing in the number of nonzero
c                 entries.  Values added for padding in Ellpack are
c                 not included in the Jagged-Diagonal format.
c  
c   int indx()    array of length nnz consisting of the column indices
c                 of the corresponding entries in val.
c  
c   int pntr()    array of length maxnz+1, where pntr(i) - pntr(1) + 1
c                 points to the location in val of the first element
c                 in the row-permuted Ellpack represenation of A.
c  
c   int maxnz     max number of nonzeros elements per row.
c  
c   int iperm()   integer array of length m such that i = iperm(i'), 
c                 where row i in the original Ellpack representation
c                 corresponds to row i' in the permuted representation. 
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
      integer maxnz
      integer indx(*), pntr(*), iperm(*)
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
         print *,'Insufficient work space for JADSM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('JADSM', info)
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

      call djadsmk( transpose, m, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, indx, pntr, maxnz, iperm, 
     *     b, ldb, beta, c, ldc, work, lwork)
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DJADSMK(trans, m, n, 
     *           scale, dvr, dvl, alpha, uplo, diag, 
     *           val, indx, pntr, maxnz, iperm, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer m, n, maxnz, ldb, ldc, base
      character trans, scale, uplo, diag
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer brow, blkrows, i, j, l, rcol, rhscols
      integer rhscolb, rhscole, bl
      integer rowind, colind, jb, je
      logical left, right, unit, nonunit, lower, upper
      double precision alpha
      double precision beta
      integer indx(*)
      integer pntr(*)
      integer iperm(*)
      double precision val(*)
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


c
c     set up order of rows to process in and store in last
c     elements of work:
c
      do 5 i=1,m
        work(m+iperm(i)) = float(i)
 5    continue
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
             rowind = int(work(m+i))
             do 30 j=1,maxnz
               jb = pntr(j)
               je = pntr(j+1)
               if ( je - jb .lt. rowind) go to 35
               colind =  indx(jb+rowind-1)
               if ( colind .ge. i ) go to 35
               work(i) = work(i) - val(jb+rowind-1)*work(colind)
 30          continue
 35          if (nonunit) work(i) = work(i)/val(jb + rowind -1)
 25         continue

          else
c-----------------------------------------------------------------
c
c           Upper triangular:
c
            do 40 i=m,1,-1
             rowind = int(work(m+i))
             do 45 j=1,maxnz
               jb = pntr(j)
               je = pntr(j+1)
               if ( je - jb .lt. rowind) go to 50
               colind =  indx(jb+rowind-1)
               if ( colind .eq. i ) go to 45
               work(i) = work(i) - val(jb+rowind-1)*work(colind)
 45          continue
 50          if (nonunit) work(i) = work(i)/val(rowind)
 40         continue

c-----------------------------------------------------------------
          endif
        else
 
c       Working with the tranpose:
c
c 
         if (lower) then
            do 55 i=m,1,-1
             rowind = int(work(m+i))
             if (nonunit) then
               do 60 j=1,maxnz
                 jb = pntr(j)
                 je = pntr(j+1)
                 if (  (je - jb .ge. rowind) .and. 
     *                (indx(jb+rowind-1) .eq. i )   ) then
                    work(i) = work(i)/val(jb+rowind-1)
                    go to 65
                endif
 60            continue
            endif
 65            do 70 j=1,maxnz
                 jb = pntr(j)
                 je = pntr(j+1)
                 if ( je - jb .lt. rowind) go to 55 
                 colind =  indx(jb+rowind-1)
                 if ( colind .eq. i ) go to 55
                 work(colind) = work(colind) - val(jb+rowind-1)*work(i)
 70          continue
 55         continue
         else
            do 75 i=1,m
             rowind = int(work(m+i))
             if (nonunit) then
               do 80 j=1,maxnz
                 jb = pntr(j)
                 je = pntr(j+1)
                 if ( (je - jb .ge. rowind) .and. 
     *              (indx(jb+rowind-1) .eq. i )  ) then
                   work(i) = work(i)/val(jb+rowind-1)
                   go to 85
                 endif
 80            continue
             endif
 85          do 90 j=1,maxnz
               jb = pntr(j)
               je = pntr(j+1)
               if ( je - jb .lt. rowind) go to 75
               colind =  indx(jb+rowind-1)
               if (colind .eq. i) go to 90 
               work(colind) = work(colind) - val(jb+rowind-1)*work(i)
 90          continue
 75         continue
          endif
        endif

          if (left) then
            do 95 i=1,m
               c(i,l) = alpha*dvl(i)*work(i) + beta*c(i,l)
 95         continue
          else
            do 100 i=1,m
               c(i,l) = alpha*work(i) + beta*c(i,l)
 100        continue
          endif
 15     continue
 10   continue
               
      return
      end
