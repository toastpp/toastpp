#ifndef __ARPACK_H
#define __ARPACK_H

#include "mathlib.h"

int ARPACK_Eigenpair (TMatrix<double> *A, int neig,
	       TVector<double> *eigval,
	       TDenseMatrix<double> *eigvec,
	       TVector<double> *residual = 0);

#endif // !__ARPACK_H
