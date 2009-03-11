// Interface to BLZPACK routine blzdrd:
// eigenvalue and eigenvector calculation
// Requires A symmetric

#include "eigpair.h"

extern "C" {
    int blzdrd_(int*, double*, double*, int*, double*, double*, int*, int*, 
		double*, double*);
}

template<class MT>
void CholeskySolve (const TCompRowMatrix<MT> &L, const TVector<MT> &d,
    const TVector<MT> &b, TVector<MT> &x);

template<class MT>
void CholeskySolve (const TCompRowMatrix<MT> &L, const TVector<MT> &d,
    const TDenseMatrix<MT> &BT, TDenseMatrix<MT> &XT, int n);

// ===========================================================================

int Eigenpair (TMatrix<double> *A, int neig,
	       TVector<double> *eigval,
	       TDenseMatrix<double> *eigvec,
	       TVector<double> *residual)
{
    dASSERT(A->nRows() == A->nCols(), Square matrix required);
    // strictly we should also test for symmetry

    int n = A->nRows();

    // ni: number of active rows of U, V and X on process i
    // sequential version: ni = n
    int ni = n;

    // leig: maximum number of eigenpairs.
    // 1st dimension of EIG, 2nd dimension of X
    int leig = (2*neig > n ? n : 2*neig);

    // nreig: number of required eigenpairs
    //int nreig = neig;
    int nreig = neig;

    // nvbin: number of vectors in a block
    //int nvbin = (nreig < 3 ? nreig : 3);
    int nvbin = 6;

    // nlsin: maximum number of steps per run
    // 0 = auto
    int nlsin = 2*neig;
    if (nlsin > n) nlsin = n;

    // nsvin: number of starting vectors provided
    int nsvin = 0;

    // nepin: number of eigenpairs provided
    int nepin = 0;

    // gnrzd: problem type flag
    // 0 = standard eigenproblem
    // 1 = generalised
    // 2 = buckling
    static int gnrzd = 0;

    // slice: spectrum slicing flag
    // 0 = slicing off
    // 1 = slicing on
    static int slice = 0;

    // purfy: eigenvectors purification flag
    // 0 = purification off
    // 1 = purification on
    static int purfy = 0;

    // lprnt: level of printing
    // 0 = nothing is printed
    // 1 = headers and exit messages
    // 2 = plus eigenvalues
    // 3 = plus statistics
    // 4 = plus history
    // 5 = plus eigenvectors
    // 6 = plus step history
    static int lprnt = 5;

    // lfile: file unit for output
    static int lfile = 6;

    // lcomm: MPI version: communicator
    //        PVM version: number of PEs
    //        sequential version: not used
    static int lcomm = 0;

    int k1 = nvbin*nlsin;
    if (!k1) k1 = (n < 180 ? n : 180);
    int k2 = (nvbin ? nvbin : 3);
    int k4 = (leig < n ? leig : n);
    int k3 = 484 + k1 * (13 + k1*2 + k2 + (k2+2 > 18 ? k2+2 : 18)) +
         k2*k2*3 + k4*2;

    // listor: dimension of istor
    int listor = 123+k1*12;
    
    // lrstor: dimension of rstor
    int lrstor = (gnrzd ? ni * (k2*4 + k1*2 + k4) + k3 :
		          ni * (k2*4 + k1) + k3);

    // eigl: inferior limit for eigenvalues
    double eigl = 0.0;

    // eigr: superior limit for eigenvalues
    double eigr = 0.0;

    // thrsh: threshold for convergence
    //double thrsh = 0.0;
    double thrsh = 1e-12;

    // allocate buffers =====================================================

    int    *istor = new int[listor+17];
    double *rstor = new double[lrstor+5];

    int ncu = (nvbin > 3 ? nvbin : 3);
    double *u   = new double[ni*ncu];
    double *v   = new double[ni*ncu];
    double *x   = new double[ni*leig];
    double *eig = new double[leig*2];

    // set input flags ======================================================

    istor[0]  = ni;
    istor[1]  = ni;
    istor[2]  = nreig;
    istor[3]  = leig;
    istor[4]  = nvbin;
    istor[5]  = nlsin;
    istor[6]  = nsvin;
    istor[7]  = nepin;
    istor[8]  = gnrzd;
    istor[9]  = slice;
    istor[10] = purfy;
    istor[11] = lprnt;
    istor[12] = lfile;
    istor[13] = lcomm;
    istor[14] = listor;

    rstor[0] = eigl;
    rstor[1] = eigr;
    rstor[2] = thrsh;
    rstor[3] = (double)lrstor;

    int nvopu, lflag = 0;
    int nneig;     // only used for generalised problem
    double sigma;  // only used for generalised problem

    do { // begin Lanczos loop

        // ****************************************************************
        blzdrd_(istor, rstor, &sigma, &nneig, u, v, &lflag, &nvopu, eig, x);
	// ****************************************************************

	if (lflag == 1) {  // compute V = AU
	    int i, j;
	    double *pu, *pv;
	    for (j = 0; j < nvopu; j++) {
		pu = u + j*ni;
	        pv = v + j*ni;
		for (i = 0; i < ni; i++)
		    pv[i] = A->RowMult (i, pu);
	    }
	}

    } while (lflag == 1);  // end Lanczos loop

    if (!lflag) {
        if (eigval->Dim() != neig) eigval->New(neig);
	memcpy (eigval->data_buffer(), eig, neig*sizeof(double));
	if (residual) {
	    if (residual->Dim() != neig) residual->New(neig);
	    memcpy (residual->data_buffer(), eig+leig, neig*sizeof(double));
	}
	if (eigvec->nRows() != neig || eigvec->nCols() != n)
	    eigvec->New(neig, n);
	memcpy (eigvec->data_buffer(), x, neig*n*sizeof(double));
    }

    // cleanup
    delete []istor;
    delete []rstor;
    delete []u;
    delete []v;
    delete []x;
    delete []eig;

    return lflag; // 0 = ok
}


// ===========================================================================

int Eigenpair_low (TCompRowMatrix<double> *L,
		   TVector<double> *d,
		   int neig, int &nconv,
		   TVector<double> *eigval,
		   TDenseMatrix<double> *eigvec,
		   TVector<double> *residual)
{
    int n = d->Dim();

    // ni: number of active rows of U, V and X on process i
    // sequential version: ni = n
    int ni = n;

    // leig: maximum number of eigenpairs.
    // 1st dimension of EIG, 2nd dimension of X
    int leig = (2*neig > n ? n : 2*neig);

    // nreig: number of required eigenpairs
    //int nreig = neig;
    int nreig = neig;

    // nvbin: number of vectors in a block
    int nvbin = (nreig < 3 ? nreig : 3);
    //int nvbin = 6;

    // nlsin: maximum number of steps per run
    // 0 = auto
    int nlsin = 3*neig;
    if (nlsin > n) nlsin = n;

    // nsvin: number of starting vectors provided
    static int nsvin = 0;

    // nepin: number of eigenpairs provided
    static int nepin = 0;

    // gnrzd: problem type flag
    // 0 = standard eigenproblem
    // 1 = generalised
    // 2 = buckling
    static int gnrzd = 0;

    // slice: spectrum slicing flag
    // 0 = slicing off
    // 1 = slicing on
    static int slice = 0;

    // purfy: eigenvectors purification flag
    // 0 = purification off
    // 1 = purification on
    static int purfy = 0;

    // lprnt: level of printing
    // 0 = nothing is printed
    // 1 = headers and exit messages
    // 2 = plus eigenvalues
    // 3 = plus statistics
    // 4 = plus history
    // 5 = plus eigenvectors
    // 6 = plus step history
    static int lprnt = 0;

    // lfile: file unit for output
    static int lfile = 6;

    // lcomm: MPI version: communicator
    //        PVM version: number of PEs
    //        sequential version: not used
    static int lcomm = 0;

    // workspace sizes
    int k1 = nvbin*nlsin;
    if (!k1) k1 = (n < 180 ? n : 180);
    int k2 = (nvbin ? nvbin : 3);
    int k4 = (leig < n ? leig : n);
    int k3 = 484 + k1 * (13 + k1*2 + k2 + (k2+2 > 18 ? k2+2 : 18)) +
         k2*k2*3 + k4*2;

    // listor: dimension of istor
    int listor = 123+k1*12;
    
    // lrstor: dimension of rstor
    int lrstor = (gnrzd ? ni * (k2*4 + k1*2 + k4) + k3 :
		          ni * (k2*4 + k1) + k3);

    // eigl: inferior limit for eigenvalues
    static double eigl = 0.0;

    // eigr: superior limit for eigenvalues
    static double eigr = 0.0;

    // thrsh: threshold for convergence
    //double thrsh = 0.0;
    double thrsh = 1e-10;

    // allocate buffers =====================================================

    int    *istor = new int[listor+17];
    double *rstor = new double[lrstor+5];

    int ncu = (nvbin > 3 ? nvbin : 3);
    double *x   = new double[ni*leig];
    double *eig = new double[leig*2];
    RDenseMatrix U(ncu,ni);
    RDenseMatrix V(ncu,ni);
    // direct access to matrix data for FORTRAN interface
    double *u = U.data_buffer();
    double *v = V.data_buffer();

    // set input flags ======================================================

    istor[0]  = ni;
    istor[1]  = ni;
    istor[2]  = nreig;
    istor[3]  = leig;
    istor[4]  = nvbin;
    istor[5]  = nlsin;
    istor[6]  = nsvin;
    istor[7]  = nepin;
    istor[8]  = gnrzd;
    istor[9]  = slice;
    istor[10] = purfy;
    istor[11] = lprnt;
    istor[12] = lfile;
    istor[13] = lcomm;
    istor[14] = listor;

    rstor[0] = eigl;
    rstor[1] = eigr;
    rstor[2] = thrsh;
    rstor[3] = (double)lrstor;

    int nvopu, lflag = 0;
    int nneig;     // only used for generalised problem
    double sigma;  // only used for generalised problem

    do { // begin Lanczos loop

        // ****************************************************************
        blzdrd_(istor, rstor, &sigma, &nneig, u, v, &lflag, &nvopu, eig, x);
	// ****************************************************************

	if (lflag == 1) {  // solve AV = U for V
	    CholeskySolve (*L, *d, U, V, nvopu);
	}

    } while (lflag == 1);  // end Lanczos loop

    if (!lflag) {
        nconv = istor[66]; // number of converged eigenpairs
        if (eigval->Dim() != nconv) eigval->New(nconv);
	for (int i = 0; i < nconv; i++)
	    (*eigval)[i] = (eig[i] ? 1.0/eig[i] : 0.0);
	if (residual) {
	    if (residual->Dim() != nconv) residual->New(nconv);
	    memcpy (residual->data_buffer(), eig+leig, nconv*sizeof(double));
	}
	if (eigvec->nRows() != nconv || eigvec->nCols() != n)
	    eigvec->New(nconv, n);
	memcpy (eigvec->data_buffer(), x, nconv*n*sizeof(double));
    }

    // cleanup
    delete []istor;
    delete []rstor;
    delete []x;
    delete []eig;

    return lflag; // 0 = ok
}
