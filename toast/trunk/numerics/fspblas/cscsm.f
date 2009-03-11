c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine: 
c--------------------------------------------------------------------
      subroutine dcscsm( transa, m, n, unitd, dv, alpha, descra, 
     *           val, indx, pntrb, pntre,
     *           b, ldb, beta, c, ldc, work, lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c
c   Toolkit interface:
c   dcscsm -- compressed sparse column format triangular solve
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
c   double val()  scalar array of length nnz containing matrix entries.
c  
c   int indx()    integer array of length nnz containing row indices.
c
c   int pntrb()   integer array of length k such that pntrb(j)-pntrb(1)
c                 points to location in val of the first nonzero element 
c                 in column j.
c
c   int pntre()   integer array of length k such that pntre(j)-pntre(1)
c                 points to location in val of the last nonzero element 
c                 in column j.
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
      integer indx(*), pntrb(*), pntre(*)
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
         print *,'Insufficient work space for CSCSM.'
         print *,'   lwork must be at least m'
         info = 21
      endif

      if ( info .ne. 0 ) then
         call xerbla('CSCSM', info)
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

      diag = 'U'
      if ( descra(3) .eq. 0 ) diag = 'N'

c
c     Call kernel subroutine:
c

      if ( transpose .eq. 'N') then
        uplo = 'U'
        if ( descra(2) .eq. 1 ) uplo = 'L'
        call dcscsmk( m, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, indx, pntrb, pntre,
     *     b, ldb, beta, c, ldc, work, lwork)
      else
        uplo = 'U'
        if ( descra(2) .eq. 2 ) uplo = 'L'
        call dcsrsmk( m, n, scale, dv, dv, alpha, uplo, diag, 
     *     val, indx, pntrb, pntre,
     *     b, ldb, beta, c, ldc, work, lwork)
      endif
   
      return
      end



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DCSCSMK( m, n, scale, dvl, dvr, alpha, uplo, diag, 
     *           val, indx, pntrb, pntre, 
     *           b, ldb, beta, c, ldc, work, base)
      implicit none
      integer m, n, ldb, ldc, base
      character scale, uplo, diag
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, bl, nb, jb, je, len, off
      logical left, right, unit, nonunit, lower, upper
      integer rcol, rhscols, rhscolb, rhscole
      double precision alpha
      double precision beta
      double precision t
      integer indx(*), pntrb(*), pntre(*)
      double precision dvl(*)
      double precision dvr(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision work(ldc,*)
      double precision z
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
c     Calculate the number of columns for partitioning
c     b and c (rcol) by taking into acount the amount of 
c     cache required for each point column:
c     For rcol columns in block:
c         cachereq  = nnz + rcol*m*3 < maxcache
c          from val ---^      from c, b, and work
c
c     So,   rcol = (maxcache-nnz)/(m*3)
c
      rcol = (maxcache-pntre(m))/(m*3)
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
c
c       Copy c into work
c
        call dcopy(ldc*rcol,c(1,bl),1,work(1,1),1)

        if (right) then
c
c         Assign dvr*b to c:
c
          do 20 i=1,m
            z = dvr(i)
            do 25 j=rhscolb,rhscole
               c(i,j) = z*b(i,j)
 25         continue
 20       continue
        else
c
c         Assign b to c:
c
          call dcopy(m*nb, b(1,rhscolb), 1, c(1,rhscolb), 1)

        endif
c
c       Loop through the rcol columns in this block:
c
        do 30 l=rhscolb,rhscole

        if (lower) then
c-----------------------------------------------------------------
c
c       Lower triangular:
c
           do 35 j=1,m
              jb = pntrb(j)
              je = pntre(j)
              z = c(j,l)   
              if (nonunit) then
                 z = z/val(jb)
                 c(j,l) = z
              endif
              len = je -jb -1
              if (unit .and. indx(jb) .ne. j) then
                 len = je - jb
                 off = 0
              else
                 off = 1
              endif
              call daxpyi(len, -z, val(jb+off), indx(jb+off), c(1,l))
 35        continue

        else
c-----------------------------------------------------------------
c
c       Upper triangular:
c
           do 36 j=m,1,-1
              jb = pntrb(j)
              je = pntre(j)
              z = c(j,l)
              if (nonunit) then 
                 z = z/val(je-1)
                 c(j,l) = z
              endif
              len = je - jb - 1
              if (unit .and. indx(je-1) .ne. j) len = je - jb
              call daxpyi(len, -z, val(jb), indx(jb), c(1,l))
 36        continue
c-----------------------------------------------------------------
        endif

 30     continue

        if (left) then
          do 40 i=1,m
            t = alpha*dvl(i)
            do 45 j=rhscolb,rhscole
              c(i,j) = t*c(i,j) + beta*work(i,j)
 45         continue
 40       continue
       else
          do 41 i=1,m
            do 46 j=rhscolb,rhscole
              c(i,j) = alpha*c(i,j) + beta*work(i,j)
 46         continue
 41       continue
       endif

 15   continue
        
      return
      end
