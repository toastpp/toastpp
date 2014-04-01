c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine dvbrmm( transa, mb, n, kb, alpha, descra,
     *           val, indx, bindx, rpntr, cpntr, bpntrb, bpntre,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   dvbrmm -- variable block sparse row format matrix-matrix multiply
c  
c   C <- alpha A B + beta C
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
c   int kb	Number of block columns in matrix A
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
c   int *cpntr	integer array of length kb+1 such that cpntr(j)-cpntr(1)
c  		is the column index of the first point column in the j-th
c  		block column. cpntr(kb+1) is set to k+cpntr(1).
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
c   double *work	scratch array of length lwork.  lwork should be at least
c  		max(m,n)
c  
c       ------------ end interface description --------------
c--------------------------------------------------------------------
      implicit none
c
c     interface variables:
c
      integer transa, mb, n, kb, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
      double precision b(*), c(*)
      double precision work(*)
c
c     format specific interface variables:
c
      double precision val(*)
      integer indx(*), bindx(*)
      integer rpntr(*), cpntr(*), bpntrb(*), bpntre(*)
c
c     local variables:
c
      integer m, k, info
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
      else if ( mb .lt. 0 ) then
         info = 2
      else if (n .lt. 0) then
         info = 3
      else if (kb .lt. 0) then
         info = 4
      else
        m=rpntr(mb+1)-rpntr(1)
        k=cpntr(kb+1)-cpntr(1)
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
        call xerbla('VBRMM', info)
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
        call dvbrmmgk( transpose,  mb, n, kb, alpha,
     *       val, indx, bindx, rpntr, cpntr, bpntrb, bpntre,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        call dvbrmmsk( mb, n, kb, alpha,
     *       val, indx, bindx, rpntr, cpntr, bpntrb, bpntre,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        call dvbrmmkk( transpose,  mb, n, kb, alpha,
     *       val, indx, bindx, rpntr, cpntr, bpntrb, bpntre,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('VBRMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DVBRMMGK( trans,  mb, n, kb, alpha, 
     *           val, indx, bindx, rpntr, cpntr, bpntrb, bpntre, 
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer m, i, j, k, l, jb, je
      integer vb, vc, cb, bb, blkind
      integer maxcol, rows, cols,  rcol, rhscols
      double precision alpha
      double precision beta
      integer indx(*), bindx(*)
      integer rpntr(*), cpntr(*), bpntrb(*), bpntre(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.D0
      m = rpntr(mb+1) - rpntr(1)
      k = cpntr(kb+1) - cpntr(1)
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Find max number of columns in block:
c
      maxcol = 0
      do 5 i=1,kb
         if ( (cpntr(i+1) - cpntr(i) ) .gt. maxcol )
     *       maxcol =  cpntr(i+1) - cpntr(i) 
 5    continue

c
c     Loop through block rows:
c

      do 10 i=1,mb
        rows = rpntr(i+1) - rpntr(i)
c       determine number of columns of c to process with each
c       inner loop based on the cache requirement:   
c          rows*maxcol + rcol(maxcol + rows) < maxcache
c              ^                 ^        ^
c    from val /          from b / from c /
c
        rcol =  (maxcache - rows*maxcol)/(maxcol + rows)
        if ( n .lt. rcol ) rcol = n
        rhscols = n/rcol
        if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c       cb = starting (point) row of block
        cb = rpntr(i)
        jb = bpntrb(i)
        je = bpntre(i) - 1
        do 15 j=jb,je
          blkind = bindx(j)
c         vb = starting (point) element of A block
          vb = indx(j)  
c         vc = starting (point) column of A block
          vc = cpntr(blkind)
          cols = cpntr(blkind+1) - vc
          bb = 1 - rcol
          do 20 l=1,rhscols
c            bb = starting (point) column of b 
             bb = bb + rcol
             k = rcol
             if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
             if (trans .eq. 'N' ) then
               call dgemm('N','N',rows,k,cols,alpha,val(vb),rows,
     *                                             b(vc,bb),ldb,
     *                                         one,c(cb,bb),ldc)
             else
               call dgemm('T','N',cols,k,rows,alpha,val(vb),rows,
     *                                             b(cb,bb),ldb,
     *                                         one,c(vc,bb),ldc)
             endif
 20       continue
 15     continue
 10   continue
             
      return
      end
           
      subroutine DVBRMMSK( mb, n, kb, alpha, val, indx, bindx, 
     *                   rpntr, cpntr, bpntrb, bpntre, b, ldb, 
     *                   beta, c, ldc, base)
      implicit none
      integer mb, n, kb, ldb, ldc, base
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer m, i, j, k, l, jb, je
      integer vb, vc, cb, bb, blkind
      integer maxcol, rows, cols,  rcol, rhscols
      double precision alpha
      double precision beta
      integer indx(*), bindx(*)
      integer rpntr(*), cpntr(*), bpntrb(*), bpntre(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.D0
      m = rpntr(mb+1) - rpntr(1)
      k = cpntr(kb+1) - cpntr(1)
c
c     Scale c by beta:
c
      call dscal( m*n, beta, c(1,1), 1)
c
c     Find max number of columns in block:
c
      maxcol = 0
      do 5 i=1,kb
         if ( (cpntr(i+1) - cpntr(i) ) .gt. maxcol )
     *       maxcol =  cpntr(i+1) - cpntr(i) 
 5    continue

c
c     Loop through block rows:
c

      do 10 i=1,mb
        rows = rpntr(i+1) - rpntr(i)
c       determine number of columns of c to process with each
c       inner loop based on the cache requirement:   
c          rows*maxcol + rcol(maxcol + rows) < maxcache
c              ^                 ^        ^
c    from val /          from b / from c /
c
        rcol =  (maxcache - rows*maxcol)/(maxcol + rows)
        if ( n .lt. rcol ) rcol = n
        rhscols = n/rcol
        if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c       cb = starting (point) row of block
        cb = rpntr(i)
        jb = bpntrb(i)
        je = bpntre(i) - 1
        do 15 j=jb,je
          blkind = bindx(j)
c         vb = starting (point) element of A block
          vb = indx(j)  
c         vc = starting (point) column of A block
          vc = cpntr(blkind)
          cols = cpntr(blkind+1) - vc
c         if ( blkind .le. i ) then
            bb = 1 - rcol
            do 20 l=1,rhscols
c              bb = starting (point) column of b 
               bb = bb + rcol
               k = rcol
               if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
               call dgemm('N','N',rows,k,cols,alpha,val(vb),rows,
     *                                           b(vc,bb),ldb,
     *                                       one,c(cb,bb),ldc)
c              if ( blkind .lt. i ) then
               if ( blkind .ne. i ) then
                 call dgemm('T','N',cols,k,rows,alpha,val(vb),rows,
     *                                             b(cb,bb),ldb,
     *                                         one,c(vc,bb),ldc)
               endif
 20         continue
c         endif
 15     continue
 10   continue
             
      return
      end
           
      subroutine DVBRMMKK( trans, mb, n, kb, alpha, 
     *           val, indx, bindx, rpntr, cpntr, bpntrb, bpntre, 
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer mb, n, kb, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer m, i, j, k, l, jb, je
      integer vb, vc, cb, bb, blkind
      integer maxcol, rows, cols,  rcol, rhscols
      double precision alpha
      double precision beta
      integer indx(*), bindx(*)
      integer rpntr(*), cpntr(*), bpntrb(*), bpntre(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      double precision one
      one = 1.D0
      m = rpntr(mb+1) - rpntr(1)
      k = cpntr(kb+1) - cpntr(1)
c
c     Scale c by beta:
c
      call dscal( m*n, beta, c(1,1), 1)
c
c     Find max number of columns in block:
c
      maxcol = 0
      do 5 i=1,kb
         if ( (cpntr(i+1) - cpntr(i) ) .gt. maxcol )
     *       maxcol =  cpntr(i+1) - cpntr(i) 
 5    continue

c
c     Loop through block rows:
c

      do 10 i=1,mb
        rows = rpntr(i+1) - rpntr(i)
c       determine number of columns of c to process with each
c       inner loop based on the cache requirement:   
c          rows*maxcol + rcol(maxcol + rows) < maxcache
c              ^                 ^        ^
c    from val /          from b / from c /
c
        rcol =  (maxcache - rows*maxcol)/(maxcol + rows)
        if ( n .lt. rcol ) rcol = n
        rhscols = n/rcol
        if ( rhscols*rcol .ne. n ) rhscols = rhscols + 1
c       cb = starting (point) row of block
        cb = rpntr(i)
        jb = bpntrb(i)
        je = bpntre(i) - 1
        do 15 j=jb,je
          blkind = bindx(j)
c         vb = starting (point) element of A block
          vb = indx(j)  
c         vc = starting (point) column of A block
          vc = cpntr(blkind)
          cols = cpntr(blkind+1) - vc
c         if ( blkind .lt. i ) then
          if ( blkind .ne. i ) then
            bb = 1 - rcol
            do 20 l=1,rhscols
c              bb = starting (point) column of b 
               bb = bb + rcol
               k = rcol
               if ( n .lt. (bb - 1 + rcol) ) k = n + 1 - bb
               if ( trans .eq. 'N' ) then
                 call dgemm('N','N',rows,k,cols,alpha,val(vb),rows,
     *                                               b(vc,bb),ldb,
     *                                           one,c(cb,bb),ldc)
                 call dgemm('T','N',cols,k,rows,-alpha,val(vb),rows,
     *                                                 b(cb,bb),ldb,
     *                                             one,c(vc,bb),ldc)
               else
c                Reverse sign on alpha if working with transpose:
                 call dgemm('N','N',rows,k,cols,-alpha,val(vb),rows,
     *                                               b(vc,bb),ldb,
     *                                           one,c(cb,bb),ldc)
                 call dgemm('T','N',cols,k,rows,alpha,val(vb),rows,
     *                                                 b(cb,bb),ldb,
     *                                             one,c(vc,bb),ldc)
               endif
 20         continue
          endif
 15     continue
 10   continue
             
      return
      end
           
