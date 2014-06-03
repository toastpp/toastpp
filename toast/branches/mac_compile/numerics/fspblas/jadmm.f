c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine djadmm( transa, m, n, k, alpha, descra,
     *           val, indx, pntr, maxnz, iperm,
     *           b, ldb, beta, c, ldc, work, lwork)

c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   djadmm -- Jagged diagonal matrix-matrix multiply
c             (modified Ellpack)
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
      integer transa, m, n, k, ldb, ldc, lwork
      double precision alpha
      double precision beta
      integer descra(*)
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
        call xerbla('JADMM', info)
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
        call djadmmgk( transpose,  m, n, k, alpha,
     *       val, indx, pntr, maxnz, iperm,
     *       b, ldb, beta, c, ldc, descra(4))
        return
      else if (descra(1) .eq. 1  .or. 
     *         descra(1) .eq. 2        ) then
c
c       Symmetric/Hermitian  matrix multiply:
c
        if ( descra(2) .eq. 1) then
           call djadmmsk( m, n, k, alpha,
     *       val, indx, pntr, maxnz, iperm, 'L',
     *       b, ldb, beta, c, ldc, descra(4))
        else if (descra(2) .eq. 2) then
           call djadmmsk( m, n, k, alpha,
     *       val, indx, pntr, maxnz, iperm, 'U',
     *       b, ldb, beta, c, ldc, descra(4))
        else
           info = 6
           call xerbla('JADMM', info)
        endif
        return
      else if (descra(1) .eq. 4 ) then
c
c        Skew-Symmetric matrix multiply:
c
        if ( descra(2) .eq. 1) then
          call djadmmkk( transpose,  m, n, k, alpha,
     *       val, indx, pntr, maxnz, iperm, 'L',
     *       b, ldb, beta, c, ldc, descra(4))
        else if (descra(2) .eq. 2) then
          call djadmmkk( transpose,  m, n, k, alpha,
     *       val, indx, pntr, maxnz, iperm, 'L',
     *       b, ldb, beta, c, ldc, descra(4))
        else
           info = 6
           call xerbla('JADMM', info)
        endif
        return
      else
        info = 6
      endif
  
      if ( info .ne. 0 ) then
        call xerbla('JADMM', info)
        return
      endif
 
      return
      end 



c--------------------------------------------------------------------
c  Sparse BLAS kernel routine(s):
c--------------------------------------------------------------------
      subroutine DJADMMGK( trans, m, n, k, alpha, 
     *           val, indx, pntr, maxnz, iperm, 
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, maxnz, ldb, ldc, base
      character trans
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, rcol, rhscols
      integer jb, je, rhscolb, rhscole, bl, rowind, cp
      logical permute
      double precision alpha
      double precision beta
      integer indx(*)
      integer pntr(*)
      integer iperm(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      if (iperm(1) .eq. 0) then
         permute =.false.
      else
         permute =.true.
      endif
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Loop through "columns" of val:
c
      do 5 i=1,maxnz
        jb = pntr(i)
        je = pntr(i+1)
c
c     Calculate the number of columns for partitioning
c     b and c (rcol) by taking into acount the amount of 
c     cache required for calculation:
c     For rcol columns in block:
c       cachereq  =  (je-jb) + rcol*2*(je-jb) < maxcache
c           from val ---^    from c&b ---^ 
c
c     So,   rcol = (maxcache - (je-jb))/2/(je-jb)
      rcol = (maxcache - (je-jb))/2/(je-jb)
      rhscols = n/rcol
      if (rhscols*rcol .ne. n ) rhscols = rhscols+1
c
c     Loop through the rhscols block columns of b & c
c
      rhscolb = 1-rcol
      do 15 bl=1,rhscols
         rhscolb = rhscolb + rcol
         rhscole = rhscolb + rcol - 1
         if ( rhscole .gt. n ) rhscole = n
c
c        Loop through the rcol columns in this block
c
         do 20 l=rhscolb, rhscole
c
c           Loop through the "rows" of val:
c
            rowind = 1
            if (trans .eq. 'N' ) then
              if (permute) then
                do 25 j=jb,je-1
                  cp = iperm(rowind)
                  c(cp,l) = c(cp,l) + alpha*val(j)*b(indx(j),l)
                  rowind = rowind + 1
 25             continue
              else
                do 26 j=jb,je-1
                  c(rowind,l) = c(rowind,l) + alpha*val(j)*b(indx(j),l)
                  rowind = rowind + 1
 26             continue
              endif
            else
c           Working with transpose:
              if (permute) then
                do 35 j=jb,je-1
                  cp = iperm(rowind)
                  c(indx(j),l) = c(indx(j),l) + alpha*val(j)*b(cp,l)
                  rowind = rowind + 1
 35             continue
              else
                do 36 j=jb,je-1
                  c(indx(j),l) = c(indx(j),l) + alpha*val(j)*b(rowind,l)
                  rowind = rowind + 1
 36             continue
              endif
            endif
 20       continue
 15     continue
 5    continue
               
      return
      end

      subroutine DJADMMSK( m, n, k, alpha, 
     *           val, indx, pntr, maxnz, iperm, uplo,
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, maxnz, ldb, ldc, base
      character uplo
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, rcol, rhscols
      integer jb, je, rhscolb, rhscole, bl, rowind, cp, ip
      logical permute, lower
      double precision alpha
      double precision beta
      integer indx(*)
      integer pntr(*)
      integer iperm(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      if (uplo .eq. 'L' ) then 
         lower =.true.
      else
         lower =.false.
      endif
      if (iperm(1) .ne. 0) then
         permute =.true.
      else
         permute =.false.
      endif
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Loop through "columns" of val:
c
      do 5 i=1,maxnz
        jb = pntr(i)
        je = pntr(i+1)
c
c     Calculate the number of columns for partitioning
c     b and c (rcol) by taking into acount the amount of 
c     cache required for calculation:
c     For rcol columns in block:
c       cachereq  =  (je-jb) + rcol*2*(je-jb) < maxcache
c           from val ---^    from c&b ---^ 
c
c     So,   rcol = (maxcache - (je-jb))/2/(je-jb)
      rcol = (maxcache - (je-jb))/2/(je-jb)
      rhscols = n/rcol
      if (rhscols*rcol .ne. n ) rhscols = rhscols+1
c
c     Loop through the rhscols block columns of b & c
c
      rhscolb = 1-rcol
      do 15 bl=1,rhscols
         rhscolb = rhscolb + rcol
         rhscole = rhscolb + rcol - 1
         if ( rhscole .gt. n ) rhscole = n
c
c        Loop through the rcol columns in this block
c
         do 20 l=rhscolb, rhscole
c
c           Loop through the "rows" of val:
c
            rowind = 1
            if (permute) then
              do 25 j=jb,je-1
                cp = iperm(rowind)
                ip = indx(j)
                if ( (lower .and. cp .gt. ip) .or.
     *               ( (.not. lower) .and. cp .lt. ip) ) then
                  c(cp,l) = c(cp,l) + alpha*val(j)*b(ip,l)
                  c(ip,l) = c(ip,l) + alpha*val(j)*b(cp,l)
                else if (cp .eq. ip ) then
                  c(cp,l) = c(cp,l) + alpha*val(j)*b(ip,l)
                endif
                rowind = rowind + 1
 25           continue
            else
              do 26 j=jb,je-1
                ip = indx(j)
                if ( (lower .and. rowind .gt. ip) .or.
     *               ( (.not. lower) .and. rowind .lt. ip) ) then
                  c(rowind,l) = c(rowind,l) + alpha*val(j)*b(ip,l)
                  c(ip,l) = c(ip,l) + alpha*val(j)*b(rowind,l)
                else if (rowind .eq. ip ) then
                  c(rowind,l) = c(rowind,l) + alpha*val(j)*b(ip,l)
                endif
                rowind = rowind + 1
 26           continue
            endif
 20       continue
 15     continue
 5    continue
               
      return
      end

      subroutine DJADMMKK( trans, m, n, k, alpha, 
     *           val, indx, pntr, maxnz, iperm, uplo,
     *           b, ldb, beta, c, ldc, base)
      implicit none
      integer m, n, k, maxnz, ldb, ldc, base
      character trans, uplo
      integer maxcache, cacheline
      parameter (maxcache=32000, cacheline=4)
      integer i, j, l, rcol, rhscols
      integer jb, je, rhscolb, rhscole, bl, rowind, cp, ip
      logical permute, lower
      double precision alpha
      double precision beta
      integer indx(*)
      integer pntr(*)
      integer iperm(*)
      double precision val(*)
      double precision b(ldb,*), c(ldc,*)
      if (uplo .eq. 'L') then
         lower =.true.
      else
         lower =.false.
      endif
      if (iperm(1) .eq. 0) then
         permute =.false.
      else
         permute =.true.
      endif
c
c     Scale c by beta:
c
      call dscal( n*ldc, beta, c(1,1), 1)
c
c     Loop through "columns" of val:
c
      do 5 i=1,maxnz
        jb = pntr(i)
        je = pntr(i+1)
c
c     Calculate the number of columns for partitioning
c     b and c (rcol) by taking into acount the amount of 
c     cache required for calculation:
c     For rcol columns in block:
c       cachereq  =  (je-jb) + rcol*2*(je-jb) < maxcache
c           from val ---^    from c&b ---^ 
c
c     So,   rcol = (maxcache - (je-jb))/2/(je-jb)
      rcol = (maxcache - (je-jb))/2/(je-jb)
      rhscols = n/rcol
      if (rhscols*rcol .ne. n ) rhscols = rhscols+1
c
c     Loop through the rhscols block columns of b & c
c
      rhscolb = 1-rcol
      do 15 bl=1,rhscols
         rhscolb = rhscolb + rcol
         rhscole = rhscolb + rcol - 1
         if ( rhscole .gt. n ) rhscole = n
c
c        Loop through the rcol columns in this block
c
         do 20 l=rhscolb, rhscole
c
c           Loop through the "rows" of val:
c
            rowind = 1
            if (trans .eq. 'N') then
              if (permute) then
                do 25 j=jb,je-1
                  cp = iperm(rowind)
                  ip = indx(j)
                  if ( (lower .and. cp .gt. ip) .or. 
     *                 ( (.not. lower) .and.  cp .lt. ip)) then
                    c(cp,l) = c(cp,l) + alpha*val(j)*b(ip,l)
                    c(ip,l) = c(ip,l) - alpha*val(j)*b(cp,l)
                  endif
                  rowind = rowind + 1
 25             continue
              else
                do 26 j=jb,je-1
                  ip = indx(j)
                  if ( (lower .and. rowind .gt. ip) .or. 
     *                 ( (.not. lower) .and.  rowind .lt. ip)) then
c                 if ( rowind .gt. ip) then
                    c(rowind,l) = c(rowind,l) + alpha*val(j)*b(ip,l)
                    c(ip,l) = c(ip,l) - alpha*val(j)*b(rowind,l)
                  endif
                  rowind = rowind + 1
 26             continue
              endif
            else
c           Reverse the sign on alpha when working with transpose:
              if (permute) then
                do 35 j=jb,je-1
                  cp = iperm(rowind)
                  ip = indx(j)
                  if ( (lower .and. cp .gt. ip) .or. 
     *                 ( (.not. lower) .and.  cp .lt. ip)) then
c                 if ( cp .gt. ip) then
                    c(cp,l) = c(cp,l) - alpha*val(j)*b(ip,l)
                    c(ip,l) = c(ip,l) + alpha*val(j)*b(cp,l)
                  endif
                  rowind = rowind + 1
 35             continue
              else
                do 36 j=jb,je-1
                  ip = indx(j)
                  if ( (lower .and. rowind .gt. ip) .or. 
     *                 ( (.not. lower) .and.  rowind .lt. ip)) then
c                 if ( rowind .gt. ip) then
                    c(rowind,l) = c(rowind,l) - alpha*val(j)*b(ip,l)
                    c(ip,l) = c(ip,l) + alpha*val(j)*b(rowind,l)
                  endif
                  rowind = rowind + 1
 36             continue
              endif
            endif
 20       continue
 15     continue
 5    continue
               
      return
      end
