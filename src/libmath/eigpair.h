// Interface to BLZPACK routine blzdrd:
// eigenvalue and eigenvector calculation
// using Block Lanczos method
// Requires A symmetric

#ifndef __EIGPAIR_H
#define __EIGPAIR_H

#include "mathlib.h"

int Eigenpair (TMatrix<double> *A, int neig,
	       TVector<double> *eigval,
	       TDenseMatrix<double> *eigvec,
	       TVector<double> *residual = 0);

// A: symmetric matrix whose eigenpairs are to be determined
// neig: number of required eigenpairs
// eigval: vector of eigenvalues (dim = neig on exit)
// eigvec: matrix of eigenvectors (dim = neig x n on exit)
// residual: vector of residuals (dim = neig on exit, if defined)
// return value: error flag (0=ok)

int Eigenpair_low (TCompRowMatrix<double> *L,
		   TVector<double> *d,
		   int neig,
		   int &nconv,
		   TVector<double> *eigval,
		   TDenseMatrix<double> *eigvec,
		   TVector<double> *residual = 0);

// L,d: Cholesky decomposition of symmetric matrix
// neig: number of required eigenpairs
// nconv: number of converged eigenpairs on exit
// eigval: vector of eigenvalues (dim = nconv on exit)
// eigvec: matrix of eigenvectors (dim = nconv x n on exit)
// residual: vector of residuals (dim = nconv on exit, if defined)
// return value: error flag (0=ok)

#endif // !__EIGPAIR_H
