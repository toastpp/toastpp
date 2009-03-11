c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dcscmm( transa, m, n, k, alpha, descra,
     *           val, indx, pntrb, pntre,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dcscmm -- compressed sparse column format matrix-matrix multiply
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
      integer transa, m, n, k, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
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
        call xerbla('CSRMM', info)
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
          call dcscmmgk( m, n, k, alpha,
     *       val, indx, pntrb, pntre,
     *       b, ldb, beta, c, ldc, descra(4))
        else
          call dcsrmmgk( k, n, m, alpha,
     *       val, indx, pntrb, pntre,
     *       b, ldb, beta, c, ldc, descra(4))
        endif
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dcscmmsk( m, n, k, alpha,
     *       val, indx, pntrb, pntre,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        if (transpose .eq. 'N') then
          call dcscmmkk( m, n, k, alpha,
     *       val, indx, pntrb, pntre,
     *       b, ldb, beta, c, ldc, descra(4))
        else
          call dcsrmmkk( m, n, k, alpha,
     *       val, indx, pntrb, pntre,
     *       b, ldb, beta, c, ldc, descra(4))
        endif
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('CSRMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DCSCMMGK( m, n, k, alpha, val, indx, pntrb,
     *                          pntre, b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer blknnz, i, j, l, bl, jb, je
      integer bcol, blkcols, blkcolb, blkcole
      integer rcol, rhscols, rhscolb, rhscole
      double precision alpha
      double precision beta
      integer indx(*), pntrb(*), pntre(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Calculate number of block columns of A, based on CACHELINE:
c
      bcol = cacheline
      blkcols = k/bcol
c     if ( mod(k,bcol) .ne. 0 ) blkcols = blkcols + 1
      if ( blkcols*bcol .ne. k ) blkcols = blkcols + 1
				
c
c     Loop through blkcols block columns:
c
      do 5 i=1,blkcols
        blkcolb = (i-1)*bcol + 1
        blkcole = blkcolb + bcol - 1
        if (blkcole .ge. k ) blkcole = k
c
c       Count the number of nonzeros in this block column:
c
        blknnz = 0
        do 10 j=blkcolb,blkcole
          blknnz = blknnz + pntre(j) - pntrb(j)
 10     continue
c
c       Calculate the number of columns for partitioning
c       b and c (rcol) by taking into acount the amount of cache
c       required for each point column:
c       For rcol columns in block:
c         cachereq  = blknnz + rcol*bcol*(1+blknnz) < maxcache
c          from val  ---^      from b  ---^   ^--- from c
c
c       So,   rcol = (maxcache-blknnz)/(bcol*(1+blknnz))
c
        rcol = (maxcache-blknnz)/(bcol*(1+blknnz))
        rhscols = n/rcol
c       if ( mod(n,rcol) .ne. 0 ) rhscols = rhscols + 1
        if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c
c       Now, loop through the rhscols columns:
c
        do 15 bl=1,rhscols
          rhscolb = (bl - 1)*rcol + 1
          rhscole = rhscolb + rcol - 1
          if (rhscole .ge. n ) rhscole = n
c
c         Loop through the rcol columns in this block:
c
          do 20 l=rhscolb,rhscole
c
c            Loop through the bcol columns (of A) in this block:
c
             do 25 j=blkcolb,blkcole
                jb = pntrb(j)
                je = pntre(j)
                call daxpyi(je-jb, alpha*b(j,l),val(jb),indx(jb),c(1,l))
 25          continue
 20        continue
 15     continue
  5   continue
        
      return
      end

c Symmetric variant:
      subroutine DCSCMMSK( m, n, k, alpha, val, indx, pntrb,
     *                          pntre, b, ldb, beta, c, ldc, base)
	  implicit none
      integer m, n, k, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer blknnz, i, j, l, bl, jb, je
      integer bcol, blkcols, blkcolb, blkcole
      integer rcol, rhscols, rhscolb, rhscole
      double precision alpha
      double precision beta
      double precision t
      integer indx(*), pntrb(*), pntre(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision ddoti
c
c     Scale c by beta:
c
      call dscal( n*m, beta, c(1,1), 1)
c
c     Calculate number of block columns of A, based on CACHELINE:
c
      bcol = cacheline
      blkcols = m/bcol
c     if ( mod(m,bcol) .ne. 0 ) blkcols = blkcols + 1
      if ( blkcols*bcol .ne. m ) blkcols = blkcols + 1
				
c
c     Loop through blkcols block columns:
c
      do 5 i=1,blkcols
        blkcolb = (i-1)*bcol + 1
        blkcole = blkcolb + bcol - 1
        if (blkcole .ge. m ) blkcole = m
c
c       Count the number of nonzeros in this block column:
c
        blknnz = 0
        do 10 j=blkcolb,blkcole
          blknnz = blknnz + pntre(j) - pntrb(j)
 10     continue
c
c       Calculate the number of columns for partitioning
c       b and c (rcol) by taking into acount the amount of cache
c       required for each point column:
c       For rcol columns in block:
c         cachereq  = blknnz + rcol*bcol*(2+2*blknnz) < maxcache
c          from val  ---^      from b  ---^   ^--- from c
c
c       So,   rcol = (maxcache-blknnz)/(bcol*(1+blknnz))/2
c
        rcol = (maxcache-blknnz)/(bcol*(1+blknnz))/2
        rhscols = n/rcol
c       if ( mod(n,rcol) .ne. 0 ) rhscols = rhscols + 1
        if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c
c       Now, loop through the rhscols columns:
c
        do 15 bl=1,rhscols
          rhscolb = (bl - 1)*rcol + 1
          rhscole = rhscolb + rcol - 1
          if (rhscole .ge. n ) rhscole = n
c
c         Loop through the rcol columns in this block:
c
          do 20 l=rhscolb,rhscole
c
c            Loop through the bcol columns (of A) in this block:
c
             do 25 j=blkcolb,blkcole
                jb = pntrb(j)
                je = pntre(j)
                if ( indx(jb) .eq. j ) then
                  c(j,l) = c(j,l) + alpha * b(j,l) * val(jb)
                  jb = jb + 1 
                elseif (indx(je-1) .eq. j ) then
                  c(j,l) = c(j,l) + alpha * b(j,l) * val(je-1)
                  je = je - 1 
                endif
                t = ddoti(je-jb, val(jb), indx(jb), b(1,l))
                call daxpyi(je-jb, alpha*b(j,l),val(jb),indx(jb),c(1,l))
                c(j,l) = c(j,l) + alpha * t
 25          continue
 20        continue
 15     continue
  5   continue
        
      return
      end

c Skew-Symmetric variant:
      subroutine DCSCMMKK( m, n, k, alpha, val, indx, pntrb,
     *                          pntre, b, ldb, beta, c, ldc, base)
	  implicit none
      integer m, n, k, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer blknnz, i, j, l, bl, jb, je
      integer bcol, blkcols, blkcolb, blkcole
      integer rcol, rhscols, rhscolb, rhscole
      double precision alpha
      double precision beta
      double precision t
      integer indx(*), pntrb(*), pntre(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision ddoti
c
c     Scale c by beta:
c
      call dscal( n*m, beta, c(1,1), 1)
c
c     Calculate number of block columns of A, based on CACHELINE:
c
      bcol = cacheline
      blkcols = m/bcol
c     if ( mod(m,bcol) .ne. 0 ) blkcols = blkcols + 1
      if ( blkcols*bcol .ne. m ) blkcols = blkcols + 1
				
c
c     Loop through blkcols block columns:
c
      do 5 i=1,blkcols
        blkcolb = (i-1)*bcol + 1
        blkcole = blkcolb + bcol - 1
        if (blkcole .ge. m ) blkcole = m
c
c       Count the number of nonzeros in this block column:
c
        blknnz = 0
        do 10 j=blkcolb,blkcole
          blknnz = blknnz + pntre(j) - pntrb(j)
 10     continue
c
c       Calculate the number of columns for partitioning
c       b and c (rcol) by taking into acount the amount of cache
c       required for each point column:
c       For rcol columns in block:
c         cachereq  = blknnz + rcol*bcol*(2+2*blknnz) < maxcache
c          from val  ---^      from b  ---^   ^--- from c
c
c       So,   rcol = (maxcache-blknnz)/(bcol*(1+blknnz))/2
c
        rcol = (maxcache-blknnz)/(bcol*(1+blknnz))/2
        rhscols = n/rcol
c       if ( mod(n,rcol) .ne. 0 ) rhscols = rhscols + 1
        if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c
c       Now, loop through the rhscols columns:
c
        do 15 bl=1,rhscols
          rhscolb = (bl - 1)*rcol + 1
          rhscole = rhscolb + rcol - 1
          if (rhscole .ge. n ) rhscole = n
c
c         Loop through the rcol columns in this block:
c
          do 20 l=rhscolb,rhscole
c
c            Loop through the bcol columns (of A) in this block:
c
             do 25 j=blkcolb,blkcole
                jb = pntrb(j)
                je = pntre(j)
                if ( indx(jb) .eq. j ) jb = jb + 1 
                if ( indx(je-1) .eq. j ) je = je - 1 
                t = ddoti(je-jb, val(jb), indx(jb), b(1,l))
                call daxpyi(je-jb, alpha*b(j,l),val(jb),indx(jb),c(1,l))
                c(j,l) = c(j,l) - alpha * t
 25          continue
 20        continue
 15     continue
  5   continue
        
      return
      end

