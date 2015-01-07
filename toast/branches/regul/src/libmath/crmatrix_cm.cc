// ==========================================================================
// Module mathlib
// File crmatrix_cm.cc
// Definition of class SCCompRowMatrixMixed ("complex mixed")
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include <iostream>
#include <iomanip>
#include <string.h>
#include "mathlib.h"
#include "slu_zdefs.h"
#include "supermatrix.h"
//#include "zsp_defs.h"

using namespace std;

SCCompRowMatrixMixed::SCCompRowMatrixMixed ()
  : TCompRowMatrix<std::complex<float> > ()
{}

SCCompRowMatrixMixed::SCCompRowMatrixMixed (int rows, int cols)
  : TCompRowMatrix<std::complex<float> > (rows, cols)
{}

SCCompRowMatrixMixed::SCCompRowMatrixMixed (int rows, int cols,
    const idxtype *_rowptr, const idxtype *_colidx,
    const std::complex<float> *data)
  : TCompRowMatrix<std::complex<float> > (rows, cols, _rowptr, _colidx, data)
{}

SCCompRowMatrixMixed::SCCompRowMatrixMixed (const SCCompRowMatrixMixed &m)
  : TCompRowMatrix<std::complex<float> > (m)
{}

SCCompRowMatrixMixed::~SCCompRowMatrixMixed ()
{}

int SCCompRowMatrixMixed::SparseRow (int r, idxtype *ci,
    std::complex<double> *rv) const
{
    int i, r0 = rowptr[r], nz = rowptr[r+1]-r0;
    for (i = 0; i < nz; i++) {
        ci[i] = colidx[r0+i];
	rv[i] = (std::complex<double>)val[r0+i];
    }
    return nz;
}

void SCCompRowMatrixMixed::Ax (const CVector &x, CVector &b) const
{
    dASSERT(x.Dim() == cols,
	"Parameter 1 invalid size (expected %d, actual %d)", cols, x.Dim());
    if (b.Dim() != rows) b.New(rows);

    int r, i, i2;
    std::complex<double> br;

    for (r = i = 0; r < rows;) {
	i2 = rowptr[r+1];
	for (br = std::complex<double>(0,0); i < i2; i++)
	    br += (std::complex<double>)val[i] * x[colidx[i]];
	b[r++] = br;
    }
}

void SCCompRowMatrixMixed::Ax (const CVector &x, CVector &b,
    int r1, int r2) const
{
    dASSERT(x.Dim() == cols,
        "Parameter 1 invalid size (expected %d, actual %d)", cols, x.Dim());

    idxtype r, i2;
    idxtype i = rowptr[r1];
    idxtype *pcolidx = colidx+i;
    std::complex<double> br;
    std::complex<float> *pval = val+i;

    if (b.Dim() != rows) b.New (rows);

    for (r = r1; r < r2;) {
	i2 = rowptr[r+1];
	for (br = std::complex<double>(0,0); i < i2; i++)
	    br += (std::complex<double>)(*pval++) * x[*pcolidx++];
	b[r++] = br;
    }
}

void SCCompRowMatrixMixed::ATx (const CVector &x, CVector &b) const
{
    dASSERT(x.Dim() == rows, "Invalid size - vector x");
    dASSERT(b.Dim() == cols, "Invalid size - vector b");

    if (!col_access) SetColAccess();
    int i, c;
    for (c = 0; c < cols; c++) b[c] = 0;
    for (i = c = 0; i < nval; i++) {
        while (colptr[c+1] <= i) c++;
	b[c] += conj((std::complex<double>)val[vofs[i]]) * x[rowidx[i]];
    }
}

// ==========================================================================
// Generalised minimal residuals (GMRES)
// ==========================================================================

int GMRES (const SCCompRowMatrixMixed &A, const CVector &b, CVector &x,
    double &elim, SCPreconditionerMixed *precon, int restart, int maxit,
    void (*clbk)(void*))
{
#ifdef VERBOSE_GMRES
    cerr << "gmres: Enter (tol=" << elim << ")" << endl;
#endif

    int N=b.Dim();
    int MAX_GMRES_STEPS=restart;
    int MAX_CYCLE = (maxit ? maxit : N+1);
    int j, k, l, cycle, n;
    double norm_b,norm_v,norm_x, tmp;
    std::complex<double> h1, h2, r, sum;
    TVector<std::complex<double> > y(MAX_GMRES_STEPS);
    TVector<std::complex<double> > s(MAX_GMRES_STEPS+1);
    TVector<std::complex<double> > c1(MAX_GMRES_STEPS);
    TVector<std::complex<double> > c2(MAX_GMRES_STEPS);
    TDenseMatrix<std::complex<double> > H(MAX_GMRES_STEPS,MAX_GMRES_STEPS+1);
    TVector<std::complex<double> > *v =
        new TVector<std::complex<double> >[MAX_GMRES_STEPS+1];
    for (j = 0; j <= MAX_GMRES_STEPS; j++)
        v[j].New (N);
 
    //----------------------------------------------------------------
    /* algorithm is prepaired to loop over #MAX_CYCLE iterations */
    cycle = j = 0;
    norm_b = sqrt(l2normsq(b));
    norm_x = sqrt(l2normsq(x));
 
    //cerr << endl;
    //cerr << "gmres: l2_norm(b) = " << norm_b;
    //cerr << "       l2_norm(x) = " << norm_x << endl;
 
    //----------------------------------------------------------------
 
    //M.S.: commented out
    //diag_inv = precond_setup(AC,b);
 
    /* v[0] = r = b - A * x(i) */
    v[0]  = (norm_x < elim*N) ? b : b - A*x;

    // MS: inserted
    if (precon) precon->Apply (v[0], v[0]);

    // for (n = 0; n < N; n++) v[0][n] = ctmp1[n];
 
    //----------------------------------------------------------------
 
    //cerr << "gmres: restart    = " << restart;
    //cerr << "       elim       = " << elim << endl;
 
    //----------------------------------------------------------------
 
    while ( ((norm_v=sqrt(l2normsq(v[0]))) > elim*norm_b)
            && (cycle < MAX_CYCLE)
	    )
      {
 
        // cerr << "c = " << cycle << ":: " << norm_v/norm_b << flush;
 
        //------------------------------------------------------------
	
	//M.S.: commented out
	//precon->Apply (v[0], v[0]);

        // for (n = 0; n < N; n++) v[0][n] = ctmp2[n];
	
        s[0] = norm_v;
	
        for (n = 0; n < N; n++) {
	    v[0][n] = v[0][n]/s[0];
        }
 
        //------------------------------------------------------------
	/* GMRES iteration */
        j = -1;
        cycle++;
        while ((norm_x=norm(s[j+1])) > elim*norm_b
                && (j < (MAX_GMRES_STEPS-1))
                && (cycle < MAX_CYCLE)
		)
	  {
	    j++;

	    // cerr << "  " << j << ": " << (norm_x/norm_b) << flush;
 
	    //---------------------------------------------------------
 
	    /* w = v[j+1] */
 
	    // for (n=0; n < N; n++) ctmp1[n] = v[j][n];

 	    TVector<std::complex<double> > &vj1 = v[j+1];
	    if (precon) precon->Apply (A*v[j], vj1);
	    else        vj1 = A*v[j];

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 0; k <= j; k++) {
	        TVector<std::complex<double> > &vk = v[k];
                sum = std::complex<double>(0,0);
		for (n = 0; n < N; n++) sum += vj1[n] * conj(vk[n]);
		H(j,k) = sum;  // H(j,k) = C-innerproduct(v[j+1],v[k]);
		for (n = 0; n < N; n++) vj1[n] -= sum * vk[n];
	    }
 
	    // for (n = 0; n < N; n++) ctmp1[n] = v[j+1][n];
 
	    tmp = sqrt(l2normsq(vj1));
	    H(j,j+1) = tmp;
	    tmp = 1.0/tmp;

	    for (n = 0; n < N; n++) vj1[n] *= tmp;
 
	    /* Q-R alg */
	    for (k = 0; k < j; k++) {
                h1 = c1[k]*H(j,k) + c2[k]*H(j,k+1);
		h2 = -c2[k]*H(j,k) + c1[k]*H(j,k+1);
		H(j,k) = h1;
		H(j,k+1) = h2;
	    }
	    r = sqrt(H(j,j)*H(j,j) + H(j,j+1)*H(j,j+1));
	    c1[j] = H(j,j)/r;
	    c2[j] = H(j,j+1)/r;
	    H(j,j) = r;
	    H(j,j+1) = 0.0;
	    s[j+1] = -c2[j]*s[j];
	    s[j] = s[j]*c1[j];
	  }
 
        //------------------------------------------------------------
 
#ifdef VERBOSE_GMRES
        cerr << "gmres: " << j+1 << ": " << (norm_x/norm_b) << endl;
#endif

        //------------------------------------------------------------
 
        /* Solving of the system of equations : H y = s */
        for (l = j; l >= 0; l--) {
	    y[l] = s[l]/H(l,l);
	    for (k = (l-1); k >= 0; k--) {
	        s[k] = s[k] - H(l,k)*y[l];
	    }
        }
 
        //------------------------------------------------------------
 
        /* updating solution */
        for (l = j; l >= 0; l--)
            for (n = 0; n < N; n++) x[n] = x[n] + y[l]*v[l][n];
 
        // cerr << "l2norm(x) = " << sqrt(l2normsq(x)) << endl;
 
        //------------------------------------------------------------
 
        /* computing new residual */
        v[0] = b - A*x;
	if (precon) precon->Apply (v[0], v[0]);

#ifdef VERBOSE_GMRES
        cerr << "gmres: norm(r) = " << sqrt(l2normsq(v[0])/norm_b);
        cerr << "  (" << cycle << ")" << endl;
#endif

	// pass the current iteration solution to the callback function
	if (clbk) clbk((void*)&x);
      }
 
    //----------------------------------------------------------------
 
    // for (j = 0; j <= MAX_GMRES_STEPS; j++) delete [] v[j];
 
    delete [] v;
 
    return cycle;
}
