c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dellmm( transa, m, n, k, alpha, descra,
     *           val, indx, lda, maxnz,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dellmm -- Ellpack format matrix-matrix multiply
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
      integer transa, m, n, k, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
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
        call xerbla('ELLMM', info)
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
        call dellmmgk( transpose,  m, n, k, alpha,
     *       val, indx, lda, maxnz,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dellmmsk( m, n, k, alpha,
     *       val, indx, lda, maxnz,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call dellmmkk( transpose,  m, n, k, alpha,
     *       val, indx, lda, maxnz,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('ELLMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DELLMMGK( trans, m, n, k, alpha, 
     *           val, indx, lda, maxnz,
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, lda, maxnz, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer brow, blkrows, i, j, l, rcol, rhscols
      integer jb, je, rhscolb, rhscole, bl, vj
      double precision alpha
      double precision beta
      integer indx(lda,*)
      double precision val(lda,*)
      double precision b(ldb,*), c(ldc,*)
      double precision t
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
      if ( trans .eq. 'T' ) then
        call dscal(n*ldc, beta, c(1,1), 1)
      endif
c
c     Loop through block rows:
c
      jb = 1-brow
      do 5 i=1,blkrows
        jb = jb+brow
        je = jb+brow-1
        if ( je .gt. m ) je = m
c
c       Loop through the brow rows in this block
c
        do 10 j=jb,je
c
c         Loop through the rhscols block columns of b & c
c
          rhscolb = 1-rcol
          do 15 bl=1,rhscols
            rhscolb = rhscolb + rcol
            rhscole = rhscolb + rcol - 1
            if ( rhscole .gt. n ) rhscole = n
c
c           Loop through the rcol columns in this block
c
            do 20 l=rhscolb, rhscole
c
c             Loop through the maxnz entries of val:
c        
              if (trans .eq. 'N' ) then
                t = 0.D0
                do 25 vj=1,maxnz
                  t = t + val(j,vj)*b(indx(j,vj),l)
 25             continue
                c(j,l) = beta*c(j,l) + alpha*t
              else
                do 26 vj=1,maxnz
                  c(indx(j,vj),l) =  
     *                   c(indx(j,vj),l) + alpha*val(j,vj)*b(j,l)
 26             continue
              endif
 20         continue
 15       continue
 10     continue
 5    continue
               
      return
      end

      subroutine DELLMMSK( m, n, k, alpha, val, indx, lda, maxnz,
     *                          b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, lda, maxnz, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer brow, blkrows, i, j, l, rcol, rhscols
      integer jb, je, rhscolb, rhscole, bl, vj, index
      double precision alpha
      double precision beta
      integer indx(lda,*)
      double precision val(lda,*)
      double precision b(ldb,*), c(ldc,*)
      double precision t
      brow = cacheline
      blkrows = m/brow
      if (blkrows*brow .ne. m ) blkrows = blkrows + 1
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
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
c     Loop through block rows:
c
      jb = 1-brow
      do 5 i=1,blkrows
        jb = jb+brow
        je = jb+brow-1
        if ( je .gt. m ) je = m
c
c       Loop through the brow rows in this block
c
        do 10 j=jb,je
c
c         Loop through the rhscols block columns of b & c
c
          rhscolb = 1-rcol
          do 15 bl=1,rhscols
            rhscolb = rhscolb + rcol
            rhscole = rhscolb + rcol - 1
            if ( rhscole .gt. n ) rhscole = n
c
c           Loop through the rcol columns in this block
c
            do 20 l=rhscolb, rhscole
c
c             Loop through the maxnz entries of val:
c
              t = 0.D0
              do 25 vj=1,maxnz
                index = indx(j,vj)
                if ( index .eq. j ) then
                  t = t + val(j,vj)*b(index,l)
                else 
                  t = t + val(j,vj)*b(index,l)
                  c(index,l) = c(index,l) + alpha*val(j,vj)*b(j,l)
                endif
 25           continue
 26           c(j,l) = c(j,l) + alpha*t
 20         continue
 15       continue
 10     continue
 5    continue
               
      return
      end

      subroutine DELLMMKK( trans, m, n, k, alpha, val, indx, lda, maxnz,
     *                          b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, lda, maxnz, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer brow, blkrows, i, j, l, rcol, rhscols
      integer jb, je, rhscolb, rhscole, bl, vj, index
      double precision alpha
      double precision beta
      integer indx(lda,*)
      double precision val(lda,*)
      double precision b(ldb,*), c(ldc,*)
      double precision t
      brow = cacheline
      blkrows = m/brow
      if (blkrows*brow .ne. m ) blkrows = blkrows + 1
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
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
c     Loop through block rows:
c
      jb = 1-brow
      do 5 i=1,blkrows
        jb = jb+brow
        je = jb+brow-1
        if ( je .gt. m ) je = m
c
c       Loop through the brow rows in this block
c
        do 10 j=jb,je
c
c         Loop through the rhscols block columns of b & c
c
          rhscolb = 1-rcol
          do 15 bl=1,rhscols
            rhscolb = rhscolb + rcol
            rhscole = rhscolb + rcol - 1
            if ( rhscole .gt. n ) rhscole = n
c
c           Loop through the rcol columns in this block
c
            do 20 l=rhscolb, rhscole
c
c             Loop through the maxnz entries of val:
c
              if ( trans .eq. 'N' ) then
                t = 0.D0
                do 25 vj=1,maxnz
                  index = indx(j,vj)
                  if ( index .ne. j ) then
                    t = t + val(j,vj)*b(index,l)
                    c(index,l) = c(index,l) - alpha*val(j,vj)*b(j,l)
                  endif
 25             continue
 26             c(j,l) = c(j,l) + alpha*t
              else
c               Reverse the sign on alpha if working with transpose
                t = 0.D0
                do 35 vj=1,maxnz
                  index = indx(j,vj)
                  if ( index .ne. j ) then
                    t = t + val(j,vj)*b(index,l)
                    c(index,l) = c(index,l) + alpha*val(j,vj)*b(j,l)
                  endif
 35             continue
 36             c(j,l) = c(j,l) - alpha*t
              endif
 20         continue
 15       continue
 10     continue
 5    continue
               
      return
      end

