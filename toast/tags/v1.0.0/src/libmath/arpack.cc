// Interface to FORTRAN routines in ARPACK
// Eigenpair calculations

#include "arpack.h"

using namespace std;

extern "C" {
  int dsaupd_ (int *ido, char *bmat, int *n, char *which, int *nev,
	       double *tol, double *resid, int *ncv, double *v,
	       int *ldv, int *iparam, int *ipntr, double *workd,
	       double *workl, int *lworkl, int *info);
  int dseupd_ (int *rvec, char *howmny, int *select, double *d, double *z,
	       int *ldv, double *sigma, char *bmat, int *n, char *which,
	       int *nev, double *tol, double *resid, int *ncv, double *v,
	       int *ldz, int *iparam, int *ipntr, double *workd,
	       double *workl, int *lworkl, int *info);
}

int ARPACK_Eigenpair (TMatrix<double> *A, int neig,
		      TVector<double> *eigval,
		      TDenseMatrix<double> *eigvec,
		      TVector<double> *residual)
{
    int i;
    int ido = 0;                  // revrse communication flag
    char bmat = 'I';              // standard eigenvalue problem A*x = lambda*x
    int n = A->nRows();           // dimension of eigenproblem
    char which[2] = {'S','A'};    // request smallest algebraic eigenvalues
    int nev = neig;               // requested number of eigenpairs
    //double tol = 0.0;             // use default tolerance
    double tol = 1e-10;
    int ncv = (12*nev)/5;         // problem dependent - CHECK THIS!
    if (ncv > n) ncv = n;
    int ldv = n;                  // leading dimension of v
    int lworkl = ncv*ncv + 8*ncv; // work array dimension
    int info = 0;                 // use random initial residual vector
    int iparam[11];               // input flags
    int ipntr[11];                // output list of starting pointers
    
    // set input flags
    iparam[0] = 1;               // ishift
    iparam[2] = 300;             // maxiter
    iparam[3] = 1;               // blocksize (must be 1)
    iparam[6] = 1;               // mode: A x = lambda x
    
    // allocate dynamic data
    double *resid = new double[n];      // residuals
    double *v = new double[n*ncv];      // Lanczos basis vectors
    double *workd = new double[3*n];    // work array
    double *workl = new double[lworkl]; // work array


    // MAIN LOOP (reverse communication)

    do {

        // *************************************************************
        dsaupd_ (&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
		 iparam, ipntr, workd, workl, &lworkl, &info);
	// *************************************************************

	if (ido == 1 || ido == -1) { // Matrix * vector operation

	    double *u = workd + (ipntr[0]-1);
	    double *v = workd + (ipntr[1]-1);
	    for (i = 0; i < n; i++)
	        v[i] = A->RowMult (i, u);

	}

    } while (ido == 1 || ido == -1);

    if (info < 0) {  // failed

        cerr << "Error in dsaupd, flag = " << info << endl;

    } else {  // post-processing to extract eigenpairs

        int rvec = 1;               // generate eigenvectors
	char howmny = 'A';          // compute all 'nev' Ritz vectors
	int *select = new int[ncv]; // workspace
	static double sigma = 0.0;  // dummy

	if (eigval->Dim() != nev)
	    eigval->New(nev);
	double *d = eigval->data_buffer();

	if (eigvec->nRows() != nev || eigvec->nCols() != n)
	    eigvec->New(nev, n);
	double *z = eigvec->data_buffer();

	// *************************************************************
	dseupd_ (&rvec, &howmny, select, d, z, &ldv, &sigma, &bmat, &n, which,
		 &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd,
		 workl, &lworkl, &info);
	// *************************************************************

	if (info != 0) { // failed
	    
	    cerr << "Error in dseupd, flag = " << info << endl;

	} else {

	    // need to compute residuals here

	}

	delete []select;
    }

    // cleanup
    delete []resid;
    delete []v;
    delete []workd;
    delete []workl;

    return info;
}
