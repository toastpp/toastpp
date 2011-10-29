#ifdef __BORLANDC__
#include <strstrea.h>
#include <conio.h>
#include <process.h>

typedef unsigned pid_t;
#else
#include <sstream>
#include <unistd.h>

#endif
#include <time.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include <toast.h>

#include <sphints.h>
RDenseMatrix Chop(const RDenseMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      RDenseMatrix AA (n,m);
      for (int i = 0; i < n ; i++)        // chop small numbers
       for(int j = 0; j < m; j++)
        AA(i,j) = (fabs(A(i,j)) < 1e-12 ? 0 : A(i,j)); 
      return AA;
}
RDenseMatrix RealMat(const CDenseMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      RDenseMatrix AA (n,m);
      for (int i = 0; i < n ; i++)        // take real 
       for(int j = 0; j < m; j++)
        AA(i,j) =  A(i,j).re; 
      return AA;
}
RDenseMatrix ImagMat(const CDenseMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      RDenseMatrix AA (n,m);
      for (int i = 0; i < n ; i++)        // take imaginary
       for(int j = 0; j < m; j++)
        AA(i,j) =  A(i,j).im; 
      return AA;
}
RDenseMatrix Sparse2Dense(const RCompRowMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      RDenseMatrix AA (n,m);
      for (int j = 0; j < m ; j++) {

        RVector ej(m), vj(n); 
        ej[j] = 1;
        vj =  A*ej;
        for (int i = 0; i < n ; i++) {
	  AA(i,j) = vj[i];
	}
      }

      return AA;
}
CDenseMatrix Sparse2Dense(const CCompRowMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      CDenseMatrix AA (n,m);
      for (int j = 0; j < m ; j++) {

        CVector ej(m), vj(n); 
        ej[j] = 1;
        vj =  A*ej;
        for (int i = 0; i < n ; i++) {
	  AA(i,j) = vj[i];
	}
      }

      return AA;
}

/*--------------- sparse versions ------------------------*/
RCompRowMatrix Chop( RCompRowMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      const idxtype *rowptr, *colidx;
      int nnz;
      const double *valA = A.ValPtr();
      nnz = A.GetSparseStructure (&rowptr, &colidx);
      RCompRowMatrix AA (n,m);
      AA.Initialise(rowptr, colidx);
      double *valAA = AA.ValPtr();
      for (int i = 0; i < nnz; i++)
	valAA[i] = (fabs(valA[i]) < 1e-12 ? 0 : valA[i]); // chop small numbers
      //delete [] rowptr;
      //delete [] colidx;
      return AA;
}
RCompRowMatrix RealMat( CCompRowMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      int nnz;
      const idxtype *rowptr, *colidx;
      const toast::complex *valA = A.ValPtr();
      nnz = A.GetSparseStructure (&rowptr, &colidx);
      RCompRowMatrix AA (n,m);
      AA.Initialise(rowptr, colidx);
      double *valAA = AA.ValPtr();
      for (int i = 0; i < nnz; i++)
	valAA[i] = valA[i].re;
      //delete [] rowptr;
      //delete [] colidx;
      return AA;
}
RCompRowMatrix ImagMat( CCompRowMatrix& A) {
      int n = A.nRows();
      int m = A.nCols();
      int nnz;
      const idxtype *rowptr, *colidx;
      const toast::complex *valA = A.ValPtr();
      nnz = A.GetSparseStructure (&rowptr, &colidx);
      RCompRowMatrix AA (n,m);
      AA.Initialise(rowptr, colidx);
      double *valAA = AA.ValPtr();
      for (int i = 0; i < nnz; i++)
	valAA[i] = valA[i].im;
      //delete [] rowptr;
      //delete [] colidx;
      return AA;
}




