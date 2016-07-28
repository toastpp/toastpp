#include "mex.h"
//#include "ilupack.h"
#include "mathlib.h"
#include "felib.h"
#include "mexutil.h"

using namespace std;

// ============================================================================
// ============================================================================
// PART 1: TOAST -> MATLAB conversions
// ============================================================================
// ============================================================================


// ============================================================================
// Copy a dense RVector from TOAST to MATLAB format

void CopyVector (mxArray **array, const RVector &vec, VectorOrientation vo)
{
    int i, m = vec.Dim();
    mxArray *tmp;

    // note: these seem transposed, but produce correct result
    if (vo==ROWVEC) tmp = mxCreateDoubleMatrix (1, m, mxREAL);
    else            tmp = mxCreateDoubleMatrix (m, 1, mxREAL);
    double *pr = mxGetPr (tmp);

    for (i = 0; i < m; i++)
	pr[i] = vec[i];

    *array = tmp;
}

// ============================================================================
// Copy a dense CVector from TOAST to MATLAB format

void CopyVector (mxArray **array, const CVector &vec, VectorOrientation vo)
{
    int i, m = vec.Dim();
    mxArray *tmp;

    // note: these seem transposed, but produce correct result
    if (vo==ROWVEC) tmp = mxCreateDoubleMatrix (1, m, mxCOMPLEX);
    else            tmp = mxCreateDoubleMatrix (m, 1, mxCOMPLEX);
    double *pr = mxGetPr (tmp);
    double *pi = mxGetPi (tmp);

    for (i = 0; i < m; i++) {
        pr[i] = vec[i].real();
	pi[i] = vec[i].imag();
    }

    *array = tmp;
}

// ============================================================================
// Copy a dense IVector from TOAST to MATLAB format

void CopyVector (mxArray **array, const IVector &vec, VectorOrientation vo)
{
    int i, m = vec.Dim();

    const mwSize dim = (mwSize)vec.Dim();
    mxArray *tmp;

    if (vo==ROWVEC) tmp = mxCreateNumericMatrix (1,dim, mxINT32_CLASS, mxREAL);
    else            tmp = mxCreateNumericMatrix (dim,1, mxINT32_CLASS, mxREAL);
    int *pr = (int*)mxGetData (tmp);

    for (i = 0; i < m; i++)
	pr[i] = vec[i];

    *array = tmp;
}

// ============================================================================
// Copy a dense matrix from TOAST to MATLAB format

void CopyMatrix (mxArray **array, const RDenseMatrix &mat)
{
    int m = mat.nRows();
    int n = mat.nCols();
    int i, j, idx;

    mxArray *tmp = mxCreateDoubleMatrix (m, n, mxREAL);
    double *pr = mxGetPr (tmp);

    for (j = idx = 0; j < n; j++)
	for (i = 0; i < m; i++)
	    pr[idx++] = mat(i,j);

    *array = tmp;
}

// ============================================================================
// Copy a complex dense matrix from TOAST to MATLAB format

void CopyMatrix (mxArray **array, const CDenseMatrix &mat)
{
    int m = mat.nRows();
    int n = mat.nCols();
    int i, j, idx;

    mxArray *tmp = mxCreateDoubleMatrix (m, n, mxCOMPLEX);
    double *pr = mxGetPr (tmp);
    double *pi = mxGetPi (tmp);

    for (j = idx = 0; j < n; j++)
      for (i = 0; i < m; i++) {
	pr[idx] = mat(i,j).real();
	pi[idx] = mat(i,j).imag();
	idx++;
      }

    *array = tmp;
}

// ============================================================================
// Copy a sparse matrix from TOAST to MATLAB format

void CopyMatrix (mxArray **array, const RCompRowMatrix &mat)
{
    int i, j, k;
    int m = mat.nRows();
    int n = mat.nCols();
    int *rowptr = mat.rowptr;
    int *colidx = mat.colidx;
    const double *pval = mat.ValPtr();
    int nz = rowptr[m];
    mwIndex *rcount = new mwIndex[n];
    mwIndex idx;

    mxArray *tmp = mxCreateSparse (m, n, nz, mxREAL);

    for (i = 0; i < n; i++) rcount[i] = 0;
    for (i = 0; i < nz; i++) rcount[colidx[i]]++;

    double  *pr = mxGetPr(tmp);
    mwIndex *ir = mxGetIr(tmp);
    mwIndex *jc = mxGetJc(tmp);

    jc[0] = 0;
    for (i = 0; i < n; i++) jc[i+1] = jc[i]+rcount[i];
    for (i = 0; i < n; i++) rcount[i] = 0;
    for (i = 0; i < m; i++)
	for (k = rowptr[i]; k < rowptr[i+1]; k++) {
	    j = colidx[k];
	    idx = jc[j]+rcount[j];
	    ir[idx] = i;
	    pr[idx] = pval[k];
	    rcount[j]++;
	}
    delete []rcount;
    *array = tmp;
}

// ============================================================================
// Copy a sparse matrix from TOAST to MATLAB format

void CopyMatrix (mxArray **array, const CCompRowMatrix &mat)
{
    int i, j, k;
    int m = mat.nRows();
    int n = mat.nCols();
    int *rowptr = mat.rowptr;
    int *colidx = mat.colidx;
    const std::complex<double> *pval = mat.ValPtr();
    int nz = rowptr[m];
    mwIndex *rcount = new mwIndex[n];
    mwIndex idx;

    mxArray *tmp = mxCreateSparse (m, n, nz, mxCOMPLEX);

    for (i = 0; i < n; i++) rcount[i] = 0;
    for (i = 0; i < nz; i++) rcount[colidx[i]]++;

    double  *pr = mxGetPr(tmp);
    double  *pi = mxGetPi(tmp);
    mwIndex *ir = mxGetIr(tmp);
    mwIndex *jc = mxGetJc(tmp);

    jc[0] = 0;
    for (i = 0; i < n; i++) jc[i+1] = jc[i]+rcount[i];
    for (i = 0; i < n; i++) rcount[i] = 0;

    for (i = 0; i < m; i++)
	for (k = rowptr[i]; k < rowptr[i+1]; k++) {
	    j = colidx[k];
	    idx = jc[j]+rcount[j];
	    ir[idx] = i;
	    pr[idx] = pval[k].real();
	    pi[idx] = pval[k].imag();
	    rcount[j]++;	
	}
    delete []rcount;

    *array = tmp;
}

// ============================================================================
// Copy a symmetric sparse matrix from TOAST to MATLAB format

void CopyMatrix (mxArray **array, const CSymCompRowMatrix &mat)
{
    int i, j, k, nz;
    int m = mat.nRows();
    int n = mat.nCols();
    int *rowptr = mat.rowptr;
    int *colidx = mat.colidx;
    const std::complex<double> *pval = mat.ValPtr();
    mwIndex *rcount = new mwIndex[n];
    mwIndex idx;

    // calculate the number of nonzero entries in the equivalent
    // non-symmetric sparse matrix

    for (j = 0; j < n; j++) rcount[j] = 0;
    for (i = nz = 0; i < m; i++) {
	for (k = rowptr[i]; k < rowptr[i+1]; k++) {
	    j = colidx[k]; 
	    nz++;
	    rcount[j]++;
	    if (j < i) { // off-diagonal element: transpose
		nz++;
		rcount[i]++;
	    }
	}
    }

    mxArray *tmp = mxCreateSparse (m, n, nz, mxCOMPLEX);

    double  *pr = mxGetPr(tmp);
    double  *pi = mxGetPi(tmp);
    mwIndex *ir = mxGetIr(tmp);
    mwIndex *jc = mxGetJc(tmp);

    jc[0] = 0;
    for (i = 0; i < n; i++) jc[i+1] = jc[i]+rcount[i];
    for (i = 0; i < n; i++) rcount[i] = 0;

    for (i = 0; i < m; i++) {
	for (k = rowptr[i]; k < rowptr[i+1]; k++) {
	    j = colidx[k];
	    idx = jc[j]+rcount[j];
	    ir[idx] = i;
	    pr[idx] = pval[k].real();
	    pi[idx] = pval[k].imag();
	    rcount[j]++;

	    if (j < i) { // off-diagonal element: transpose
		idx = jc[i]+rcount[i];
		ir[idx] = j;
		pr[idx] = pval[k].real();
		pi[idx] = pval[k].imag();
		rcount[i]++;
	    }
	}
    }
    delete []rcount;
    *array = tmp;
}

// ============================================================================
// Transpose and copy a dense matrix from TOAST to MATLAB format

void CopyTMatrix (mxArray **array, const RDenseMatrix &mat)
{
    int n = mat.nRows();
    int m = mat.nCols();
    int i, j, idx;

    mxArray *tmp = mxCreateDoubleMatrix (m, n, mxREAL);
    double *pr = mxGetPr (tmp);

    for (j = idx = 0; j < n; j++)
	for (i = 0; i < m; i++)
	    pr[idx++] = mat(j,i);

    *array = tmp;
}

// ============================================================================
// Transpose and copy a sparse matrix from TOAST to MATLAB format

void CopyTMatrix (mxArray **array, const CCompRowMatrix &mat)
{
    int i;
    int m = mat.nCols();
    int n = mat.nRows();

    int *rowptr = mat.rowptr;
    int *colidx = mat.colidx;
    const std::complex<double> *pval = mat.ValPtr();
    int nz = rowptr[n];

    mxArray *tmp = mxCreateSparse (m, n, nz, mxCOMPLEX);
    double  *pr = mxGetPr (tmp);
    double  *pi = mxGetPi (tmp);
    mwIndex *ir = mxGetIr (tmp);
    mwIndex *jc = mxGetJc (tmp);

    for (i = 0; i < nz; i++) {
        pr[i] = pval[i].real();
	pi[i] = pval[i].imag();
	ir[i] = colidx[i];
    }
    for (i = 0; i <= n; i++)
	jc[i] = rowptr[i];

    *array = tmp;
}


// ============================================================================
// ============================================================================
// PART 2: MATLAB -> TOAST conversions
// ============================================================================
// ============================================================================

// ============================================================================
// Copy a RVector from MATLAB to TOAST format

void CopyVector (RVector &vec, const mxArray *array)
{
    mwIndex m = mxGetM(array);
    mwIndex n = mxGetN(array);
    mwIndex d = m*n;

    vec.New((int)d);
    double *val = vec.data_buffer();

    switch (mxGetClassID (array)) {
    case mxDOUBLE_CLASS: {
	double *pr = mxGetPr (array);
	memcpy (val, pr, d*sizeof(double));
        } break;
    case mxUINT32_CLASS:
    case mxINT32_CLASS: {
	mwIndex i;
	int *pr = (int*)mxGetData (array);
	for (i = 0; i < d; i++) val[i] = pr[i];
        } break;
    case mxUINT8_CLASS:
    case mxINT8_CLASS: {
	mwIndex i;
	char *pr = (char*)mxGetData (array);
	for (i = 0; i < d; i++) val[i] = pr[i];
        } break;
    default:
	cerr << mxGetClassName(array) << " not supported" << endl;
	break;
    }
}

// ============================================================================
// Copy a CVector from MATLAB to TOAST format

void CopyVector (CVector &vec, const mxArray *array)
{
    mwIndex m = mxGetM(array);
    mwIndex n = mxGetN(array);
    mwIndex d = m*n;
    double *pr = mxGetPr (array);
    double *pi = mxGetPi (array);
 
    vec.New((int)d);
    std::complex<double> *val = vec.data_buffer();
   
    for (mwIndex i = 0; i < d; i++)
        val[i] = std::complex<double> (pr[i], pi[i]);
}





// ============================================================================
// Copy a dense matrix from MATLAB to TOAST format

void CopyMatrix (RDenseMatrix &mat, const mxArray *array)
{
    mwIndex i, j;
    mwIndex m = mxGetM(array);
    mwIndex n = mxGetN(array);
    double *pr = mxGetPr (array);

    mat.New ((int)m,(int)n);
    double *val = mat.ValPtr();
    for (j = 0; j < n; j++)
	for (i = 0; i < m; i++)
	    val[i*n+j] = *pr++;
}

// ============================================================================
// Copy a dense integer matrix from MATLAB to TOAST format

void CopyMatrix (IDenseMatrix &mat, const mxArray *array)
{
    mwIndex i, j;
    mwIndex m = mxGetM(array);
    mwIndex n = mxGetN(array);
    double *pr = mxGetPr (array);

    mat.New ((int)m,(int)n);
    int *val = mat.ValPtr();
    for (j = 0; j < n; j++)
	for (i = 0; i < m; i++)
	    val[i*n+j] = (int)*pr++;
}

// ============================================================================
// Copy a dense matrix from MATLAB to TOAST format

void CopyTMatrix (RDenseMatrix &mat, const mxArray *array)
{
    mwIndex i, j;
    mwIndex m = mxGetM(array);
    mwIndex n = mxGetN(array);
    double *pr = mxGetPr (array);

    mat.New ((int)n,(int)m);
    double *val = mat.ValPtr();
    for (j = 0; j < n; j++)
	for (i = 0; i < m; i++)
	    val[j*m+i] = *pr++;
}

// ============================================================================
// Transpose and copy a sparse matrix from MATLAB to TOAST format

void CopyTMatrix (CCompRowMatrix &mat, const mxArray *array)
{
    mwIndex dim = mxGetNumberOfDimensions (array);
    if (dim > 2) mexErrMsgTxt ("CopyTMatrix: 2-D matrix expected");

    mwIndex m   = mxGetDimensions (array)[0];
    mwIndex n   = mxGetDimensions (array)[1];
    mwIndex nz  = mxGetNzmax (array);
    double *pr  = mxGetPr (array);
    double *pi  = mxGetPi (array);
    int *rowptr, *colidx;

#ifdef INDEX64
    mwIndex i;
    mwIndex *ir = mxGetIr (array);
    mwIndex *jc = mxGetJc (array);
    rowptr = new int[n+1];
    colidx = new int[nz];
    for (i = 0; i <= n; i++) rowptr[i] = (int)jc[i];
    for (i = 0; i < nz; i++) colidx[i] = (int)ir[i];
#else
    rowptr = mxGetJc (array);
    colidx = mxGetIr (array);
#endif

    std::complex<double> *val = new std::complex<double>[nz];
    for (mwIndex i = 0; i < nz; i++)
        val[i] = std::complex<double> (pr[i], pi[i]);

    mat.New ((int)n, (int)m);
    mat.Initialise (rowptr, colidx, val);

    delete []val;
#ifdef INDEX64
    delete []rowptr;
    delete []colidx;
#endif
}

// ============================================================================
// Transpose and copy a sparse matrix from MATLAB to TOAST format

void CopyTMatrix (RCompRowMatrix &mat, const mxArray *array)
{
    mwIndex dim = mxGetNumberOfDimensions (array);
    if (dim > 2) mexErrMsgTxt ("CopyTMatrix: 2-D matrix expected");

    mwIndex m   = mxGetDimensions (array)[0];
    mwIndex n   = mxGetDimensions (array)[1];
    mwIndex nz  = mxGetNzmax (array);
    double *pr  = mxGetPr (array);
    int *rowptr, *colidx;

#ifdef INDEX64
    mwIndex i;
    mwIndex *ir = mxGetIr (array);
    mwIndex *jc = mxGetJc (array);
    rowptr = new int[n+1];
    colidx = new int[nz];
    for (i = 0; i <= n; i++) rowptr[i] = (int)jc[i];
    for (i = 0; i < nz; i++) colidx[i] = (int)ir[i];
#else
    rowptr = mxGetJc (array);
    colidx = mxGetIr (array);
#endif

    double *val = new double[nz];
    for (mwIndex i = 0; i < nz; i++) {
	val[i] = pr[i];
    }

    mat.New ((int)n, (int)m);
    mat.Initialise (rowptr, colidx, val);

    delete []val;
#ifdef INDEX64
    delete []rowptr;
    delete []colidx;
#endif
}

// ============================================================================
// ============================================================================
// Assertion functions

void dAssert (bool cond, char *msg) {
    mxAssert (cond, msg);
}

void xAssert (bool cond, char *msg) {
    if (!cond) mexErrMsgTxt (msg);
}


// ============================================================================
// ============================================================================
// Misc. mex utility functions

// ============================================================================
// waitbar functions

static struct {
    mwSize ndim, dims[2];
    int flag;
    int range, d;
    mxArray *pmx1[1];
    mxArray *pmx2[2];
    mxArray *pmx3[1];
    mxArray *pmx4[1];
    double *pmx2_0;
    double *pmx3_0;
} waitbar = {0, {0,0}, 0, 0, 0};

void mxOpenWaitbar (char *s, int range)
{
    if (waitbar.flag < 1) {
	waitbar.flag++;
	waitbar.range = range;
	waitbar.d = std::max(1,waitbar.range/50); // granularity
	waitbar.ndim = 2; waitbar.dims[0] = 1; waitbar.dims[1] = 1;
	waitbar.pmx2[0] = mxCreateNumericArray(waitbar.ndim,waitbar.dims,mxUINT16_CLASS,mxREAL);
	waitbar.pmx2_0 = mxGetPr(waitbar.pmx2[0]);
	*waitbar.pmx2_0 = 0;
	waitbar.ndim = 2; waitbar.dims[0] = 1; waitbar.dims[1] = 1;
	waitbar.pmx3[0] = mxCreateNumericArray(waitbar.ndim,waitbar.dims,mxDOUBLE_CLASS,mxREAL);
	waitbar.pmx3_0 = mxGetPr(waitbar.pmx3[0]);
	waitbar.pmx2[1] = mxCreateString(s);
	mexCallMATLAB(1,waitbar.pmx1,2,waitbar.pmx2,"waitbar");
    }
}

void mxUpdateWaitbar (int i)
{
    double x;
    if (waitbar.flag > 0) {
	if ((i%waitbar.d) == 0) {
	    x = (double)(i+1)/waitbar.range;
	    *waitbar.pmx3_0 = x;
	    mexCallMATLAB(0,waitbar.pmx4,1,waitbar.pmx3,"waitbar");
	}
    }
}

void mxCloseWaitbar ()
{
    if (waitbar.flag > 0) {
	mexCallMATLAB (0, waitbar.pmx4, 1, waitbar.pmx1, "close");
	waitbar.flag = 0;
    }
}
