#ifndef __TOASTMEX_H
#define __TOASTMEX_H

#include "mex.h"
#include "matrix.h"
#include "mathlib.h"
//#include "toasttype.h"

#ifndef MX_API_VER
#define MX_API_VER 0
#endif

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#else
#define INDEX64
#endif

typedef enum {
    ROWVEC,
    COLVEC
} VectorOrientation;

// =========================================================================
// TOAST -> MATLAB interface conversions

void CopyVector (mxArray **array, const RVector &vec,
    VectorOrientation vo=COLVEC);
void CopyVector (mxArray **array, const CVector &vec,
    VectorOrientation vo=COLVEC);
void CopyVector (mxArray **array, const IVector &vec,
    VectorOrientation vo=COLVEC);

void CopyMatrix (mxArray **array, const RDenseMatrix &mat);
void CopyMatrix (mxArray **array, const CDenseMatrix &mat);
void CopyMatrix (mxArray **array, const RCompRowMatrix &mat);
void CopyMatrix (mxArray **array, const CCompRowMatrix &mat);
void CopyMatrix (mxArray **array, const CSymCompRowMatrix &mat);

void CopyTMatrix (mxArray **array, const RDenseMatrix &mat);
void CopyTMatrix (mxArray **array, const CCompRowMatrix &mat);


// =========================================================================
// MATLAB -> TOAST interface conversions

void CopyVector (RVector &vec, const mxArray *array);
void CopyVector (CVector &vec, const mxArray *array);

void CopyMatrix (RDenseMatrix &mat, const mxArray *array);
void CopyMatrix (IDenseMatrix &mat, const mxArray *array);

void CopyTMatrix (RDenseMatrix &mat, const mxArray *array);
void CopyTMatrix (CCompRowMatrix &mat, const mxArray *array);
void CopyTMatrix (RCompRowMatrix &mat, const mxArray *array);

// =========================================================================
// Assertion functions

void dAssert (bool cond, char *msg);
// assertion checked only in debug versions (compiled with -g)

void xAssert (bool cond, char *msg);
// assertion checked always

#endif // !__TOASTMEX_H
