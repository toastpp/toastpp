c--------------------------------------------------------------------
c  Sparse BLAS Toolkit interface routine:
c--------------------------------------------------------------------
      subroutine DJADRP( transp, m, k, val, indx, pntr, maxnz,
     *                     iperm,work,lwork)
c--------------------------------------------------------------------
c         ------------ begin interface description ------------
c   Toolkit interface:
c   djadrp -- right permutation of a jagged diagonal matrix
c  
c   A <- A P 
c   A <- A P' 
c  
c   Arguments:
c  
c   int transa	Indicates how to operate with the _permutation_ matrix
c  		0 : operate with matrix
c  		1 : operate with transpose matrix
c  
c   int m	Number of rows in matrix A
c  
c   int k	Number of columns in matrix A
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
c   double work() scratch array of length lwork.  lwork should be at least
c                 max(m,n)
c
c   int lwork     length of work array
c  
c       ------------ end interface description --------------
c--------------------------------------------------------------------
      implicit none
c
c     interface variables:
c
      integer transp, m, k, lwork
      integer work(*)
c
c     format specific interface variables:
c
      integer maxnz
      integer indx(*)
      integer pntr(*)
      integer iperm(*)
      double precision val(*)
c
c     local variables:
c
      integer nnz,i,j,ii,pb,len,pbi,leni,it
      double precision t

      nnz = pntr(maxnz+1)-1
      if ( lwork .lt. m) then
        print *, 'DJADRP: Insufficient workspace.'
        print *, '        lwork must be .ge. m'
        stop 
      endif

      do 5 i=1,m
         work(iperm(i)) = i
 5    continue
      do 10 i=1,nnz
         indx(i) = work(indx(i))
 10   continue

      do 15 i=1,maxnz
        pb = pntr(i)
        len = pntr(i+1) - pb
        do 20 j=0,len-1
          do 25 ii=i+1,maxnz
            pbi = pntr(ii) 
            leni = pntr(ii+1)-pbi
            if ( leni .ge. j+1) then
               if ( indx(pbi+j) .lt. indx(pb+j)) then
                  t = val(pbi+j)
                  val(pbi+j) = val(pb+j)
                  val(pb+j) = t
                  it = indx(pbi+j)
                  indx(pbi+j) = indx(pb+j)
                  indx(pb+j) = it
               endif
            endif
 25       continue
 20     continue
 15   continue

      return
      end
