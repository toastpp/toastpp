#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "mathdef.h"
#include "complex.h"
#include "vector.h"
#include "rtmatrix.h"
#include "matrix.h"
#include "sqmatrix.h"
#include "symatrix.h"
#include "bsmatrix.h"
#include "sbsmatrix.h"

// ==========================================================================
// member definitions

template<class MT>
TSparseBandSymMatrix<MT>::TSparseBandSymMatrix ()
    : TBandSymMatrix<MT>()
{
    rowfill = 0;
    colfill = 0;
    rowindex = 0;
    colindex = 0;
}

template<class MT>
TSparseBandSymMatrix<MT>::TSparseBandSymMatrix (int rc, int hb)
    : TBandSymMatrix<MT> (rc, hb)
{
    rowfill = 0;
    colfill = 0;
    rowindex = 0;
    colindex = 0;
}

template<class MT>
TSparseBandSymMatrix<MT>::TSparseBandSymMatrix
(const TSparseBandSymMatrix<MT> &A)
    : TBandSymMatrix<MT> (A)
{
    rowfill = 0;
    colfill = 0;
    rowindex = 0;
    colindex = 0;
}

template<class MT>
TSparseBandSymMatrix<MT>::~TSparseBandSymMatrix ()
{
    if (rowfill) delete []rowfill;
    if (rowindex) {
	for (int r = 0; r < rows; r++) if (rowindex[r]) delete []rowindex[r];
	delete []rowindex;
    }
    if (colfill) delete []colfill;
    if (colindex) {
	for (int r = 0; r < rows; r++) if (colindex[r]) delete []colindex[r];
	delete []colindex;
    }
}

template<class MT>
void TSparseBandSymMatrix<MT>::Register ()
{
    int r, c, i, rr;

    if (rowfill) delete []rowfill;
    if (rowindex) {
	for (r = 0; r < rows; r++) if (rowindex[r]) delete []rowindex[r];
	delete []rowindex;
    }
    if (colfill) delete[]colfill;
    if (colindex) {
	for (r = 0; r < rows; r++) if (colindex[r]) delete []colindex[r];
	delete []colindex;
    }

    rowfill = new int[rows];
    colfill = new int[rows];
    rowindex = new int*[rows];
    colindex = new int*[rows];

    for (r = 0; r < rows; r++) {
	rowfill[r] = 0;
	for (c = 0; c < hband; c++)
	    if (data[r][c] != 0) rowfill[r]++;
	rowindex[r] = new int[rowfill[r]];
	for (c = i = 0; c < hband; c++)
	    if (data[r][c] != 0) rowindex[r][i++] = c;

	colfill[r] = 0;
	for (c = hband-2; c >= 0; c--) {
	    if ((rr = r-c+hband-1) >= rows) break;
	    if (data[rr][c] != 0) colfill[r]++;
	}
	colindex[r] = new int[colfill[r]];
	for (i = 0, c = hband-2; c >= 0; c--) {
	    if ((rr = r-c+hband-1) >= rows) break;
	    if (data[rr][c] != 0) colindex[r][i++] = c;
	}
    }
}

template<class MT>
void TSparseBandSymMatrix<MT>::Ax (const TVector<MT> &x, TVector<MT> &b) const
{
    int r, c, i, k;

    for (r = 0; r < rows; r++) {
	MT linesum = MT(0);
	for (i = 0; i < rowfill[r]; i++) {
	    c = rowindex[r][i];
	    linesum += data[r][c] * x[c-hband+r+1];
	}
	for (i = 0; i < colfill[r]; i++) {
	    c = colindex[r][i];
	    linesum += data[r-c+hband-1][c] * x[r-c+hband-1];
	}
	b[r] = linesum;
    }
}

template<class MT>
void CHdecomp (TSparseBandSymMatrix<MT> &A)
{
    int i, j, k, l, wj;
    int n = A.Dim(ROW);
    int w = A.hband - 1;
    MT x;

    for (i = 0; i < n; i++) {
	TVector<MT> &Ai = A[i];
	for (x = (MT)0, j = 0; j < A.rowfill[i]; j++) {
	    k = A.rowindex[i][j];
	    x += Ai[k] * Ai[k];
	}
	//for (x = (MT)0, j = 0; j < w; j++) x += Ai[j] * Ai[j];
	xASSERT(Ai[w] > x, Matrix not positive definite.);
	Ai[w] = sqrt (Ai[w] - x);
	for (j = 1; j <= w && i+j < n; j++) {
	    x = (MT)0;
	    wj = w-j;
	    TVector<MT> &Aij = A[i+j];
	    for (l = 0; l < A.rowfill[i+j]; l++)
		if ((k = A.rowindex[i+j][l]) < wj) x += Aij[k] * Ai[k+j];
	    //for (k = wj-1; k >= 0; k--) x += Aij[k] * Ai[k+j];
	    Aij[wj] = (Aij[wj] - x) / Ai[w];
	}
    }
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class TSparseBandSymMatrix<double>;
template class TSparseBandSymMatrix<float>;
template class TSparseBandSymMatrix<complex>;
template class TSparseBandSymMatrix<int>;

template void CHdecomp (RSparseBandSymMatrix &A);
template void CHdecomp (FSparseBandSymMatrix &A);

#endif // NEED_EXPLICIT_INSTANTIATION

