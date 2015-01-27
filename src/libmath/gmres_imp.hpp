#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"
#include "timing.h"
#include <cmath>

using namespace std;

#define VERBOSE_GMRES_LOCAL

//int itcount = 0;
// global counter for GMRES iterations

#define VERBOSE_GMRES
// Uncomment this for verbose debugging output to stderr

extern double relax;
 
//CVector diag_inv;
//inline CVector precond_setup (const CGenericSparseMatrix&, const CVector&);
//inline CVector precond (const CGenericSparseMatrix&, const CVector&);
RVector ATA_lambdaI_x (const TMatrix<double> &A, double lambda,
    const RVector &x);

/*
    GMRES(restart) algorithm according to
    
    R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
    V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Solvers
    SIAM, Philadelphia, 1994
    
*/

// =====================================================================
// =====================================================================

template<class MT>
int gmres (int restart, const TMatrix<MT> &A, const TVector<MT> &b,
    TVector<MT> &x, TPreconditioner<MT> *precon, double &elim, int maxit,
    void (*clbk)(void*))
{
    if (toastVerbosity > 1) {
        cout << "gmres: Enter (tol=" << elim << ")" << endl << endl;
	tic();
    }

    int N=b.Dim();
    int MAX_GMRES_STEPS=restart;
    int MAX_CYCLE = (maxit ? maxit : N+1);
    int j, k, l, cycle, n;
    double norm_b,norm_v,norm_x, tmp;
    MT h1, h2, r, sum;
    TVector<MT> y(MAX_GMRES_STEPS);
    TVector<MT> s(MAX_GMRES_STEPS+1);
    TVector<MT> c1(MAX_GMRES_STEPS);
    TVector<MT> c2(MAX_GMRES_STEPS);
    TDenseMatrix<MT> H(MAX_GMRES_STEPS,MAX_GMRES_STEPS+1);
    TVector<MT> *v = new TVector<MT>[MAX_GMRES_STEPS+1];
    for (j = 0; j <= MAX_GMRES_STEPS; j++)
        v[j].New (N);
 
    //----------------------------------------------------------------
    /* algorithm is prepaired to loop over #MAX_CYCLE iterations */
    cycle = j = 0;
    TVector<MT> Mb(b);
    if (precon) precon->Apply (Mb, Mb); // MS 080310
    norm_b = sqrt(l2normsq(Mb));
    norm_x = sqrt(l2normsq(x));
 
    //cerr << endl;
    //cerr << "gmres: l2_norm(b) = " << norm_b;
    //cerr << "       l2_norm(x) = " << norm_x << endl;
 
    //----------------------------------------------------------------
 
    //M.S.: commented out
    //diag_inv = precond_setup(AC,b);
 
    /* v[0] = r = b - A * x(i) */
    //v[0]  = (norm_x < elim*N) ? b : b - A*x;

    if (norm_x < elim*N)
	v[0] = Mb;
    else {
	v[0] = A*x;
	if (precon) precon->Apply (v[0],v[0]);
	v[0] = Mb-v[0];
    }
	
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
          while ((norm_x=std::abs(s[j+1])) > elim*norm_b
                && (j < (MAX_GMRES_STEPS-1))
                && (cycle < MAX_CYCLE)
		)
	  {
	    j++;

	    if (toastVerbosity > 1)
	        cout << "\033[1Agmres(" << cycle << ',' << j << "): res="
		     << (norm_x/norm_b) << "    " << endl;

	    //---------------------------------------------------------
 
	    /* w = v[j+1] */
 
	    // for (n=0; n < N; n++) ctmp1[n] = v[j][n];

 	    TVector<MT> &vj1 = v[j+1];
	    if (precon) precon->Apply (A*v[j], vj1);
	    else        vj1 = A*v[j];

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 0; k <= j; k++) {
	        TVector<MT> &vk = v[k];
                sum = (MT)0;
     
// SP 24.01.15: See mathdef.h for details of the following OS X specifics
		for (n = 0; n < N; n++) sum += vj1[n] * toast::conj(vk[n]);
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
 
	if (toastVerbosity > 1)
	    cout << "\033[1Agmres(" << cycle << ',' << j+1 << "): res="
		 << (norm_x/norm_b) << "    " << endl;

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
	v[0] = A*x;
	if (precon) precon->Apply (v[0],v[0]);
	v[0] = Mb-v[0];

	// pass the current iteration solution to the callback function
	if (clbk) clbk((void*)&x);
    }
 
    //----------------------------------------------------------------
 
    // for (j = 0; j <= MAX_GMRES_STEPS; j++) delete [] v[j];
 
    delete [] v;
 
    if (toastVerbosity > 1)
        cout << "\033[1Agmres: Exit  (res=" << (norm_x/norm_b)
	     << ", it=" << cycle
	     << ", cpu=" << toc() << ")    " << endl;

    elim = norm_v/norm_b;
    return cycle;
}

// =====================================================================
// =====================================================================

template<class MT>
int gmres (int restart, TVector<MT> (*Av_clbk)(const TVector<MT> &v,
    void *context), const TVector<MT> &b, TVector<MT> &x,
    TPreconditioner<MT> *precon, double &elim, int maxit, void *context)
{
    // This "matrix-less" version of gmres can be used whenever matrix A is
    // not available in explicit form. The user provides a callback function
    // (Av_clbk) which is called whenever the product of the matrix with a
    // vector v is required.

    if (toastVerbosity > 1) {
        cout << "gmres: Enter (tol=" << elim << ")" << endl << endl;
	tic();
    }

    int N=b.Dim();
    int MAX_GMRES_STEPS=restart;
    int MAX_CYCLE = (maxit ? maxit : N+1);
    int j, k, l, cycle, n;
    double norm_b,norm_v,norm_x, tmp;
    MT h1, h2, r, sum;
    TVector<MT> y(MAX_GMRES_STEPS);
    TVector<MT> s(MAX_GMRES_STEPS+1);
    TVector<MT> c1(MAX_GMRES_STEPS);
    TVector<MT> c2(MAX_GMRES_STEPS);
    TDenseMatrix<MT> H(MAX_GMRES_STEPS,MAX_GMRES_STEPS+1);
    TVector<MT> *v = new TVector<MT>[MAX_GMRES_STEPS+1];
    for (j = 0; j <= MAX_GMRES_STEPS; j++)
        v[j].New (N);
 
    //----------------------------------------------------------------
    /* algorithm is prepaired to loop over #MAX_CYCLE iterations */
    cycle = j = 0;
    TVector<MT> Mb(b);
    if (precon) precon->Apply (Mb, Mb); // MS 080308
    norm_b = sqrt(l2normsq(Mb));
    norm_x = sqrt(l2normsq(x));
 
    //cerr << endl;
    //cerr << "gmres: l2_norm(b) = " << norm_b;
    //cerr << "       l2_norm(x) = " << norm_x << endl;
 
    //----------------------------------------------------------------
 
    //M.S.: commented out
    //diag_inv = precond_setup(AC,b);
 
    /* v[0] = r = b - A * x(i) */
    //v[0]  = (norm_x < elim*N) ? b : b - Av_clbk (x, context);

    if (norm_x < elim*N)
	v[0] = Mb;
    else {
	v[0] = Av_clbk (x, context);
	if (precon) precon->Apply (v[0], v[0]);
	v[0] = Mb-v[0];
    }

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
        while ((norm_x=std::abs(s[j+1])) > elim*norm_b
		   && (j < (MAX_GMRES_STEPS-1))
		   && (cycle < MAX_CYCLE)
		   )
		{
		    j++;
		    
		    if (toastVerbosity > 1) {
		      cout << "\033[1Agmres(" << cycle << ',' << j << "): res="
			   << (norm_x/norm_b) << "    " << endl;
		      cout << "gmres cycle "<< cycle << " iter " << j+1
			   << ", res = " << (norm_x/norm_b) << "    " << endl;
		    }
		    //---------------------------------------------------------

		    /* w = v[j+1] */

		    // for (n=0; n < N; n++) ctmp1[n] = v[j][n];

		    TVector<MT> &vj1 = v[j+1];
		    if (precon) precon->Apply (Av_clbk (v[j], context), vj1);
		    else        vj1 = Av_clbk(v[j], context);

		    /* modified Gram-Schmidt orthogonalization */
		    for (k = 0; k <= j; k++) {
			TVector<MT> &vk = v[k];
			sum = (MT)0;
                for (n = 0; n < N; n++) sum += vj1[n] * toast::conj(vk[n]);
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
	    
	    if (toastVerbosity > 1)
	        cout << "\033[1Agmres(" << cycle << ',' << j+1 << "): res="
		     << (norm_x/norm_b) << "    " << endl;
	    
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
	    //v[0] = b - Av_clbk(x, context);
	    //if (precon) precon->Apply (v[0], v[0]);
	    v[0] = Av_clbk(x, context);
	    if (precon) precon->Apply (v[0], v[0]);
	    v[0] = Mb-v[0];
	}
    
    //----------------------------------------------------------------
    
    // for (j = 0; j <= MAX_GMRES_STEPS; j++) delete [] v[j];
    
    delete [] v;
    
    if (toastVerbosity > 1) {
      cout << "\033[1Agmres: Exit  (res=" << (norm_x/norm_b)
	   << ", it=" << cycle
	   << ", cpu=" << toc() << ")    " << endl;
    }

    elim = norm_v/norm_b;
    return cycle;
}


// =====================================================================
// =====================================================================

#ifdef UNDEF
inline CVector precond_setup (const CGenericSparseMatrix& AC, const CVector& b)
{
    int j, k, N = b.Dim();
    toast::complex tmp1, tmp2;
    CVector diag(N);
 
    for (j = 0; j < N; j++) {
        diag[j] = AC.Get(j,j);
    }
 
    for (j = 0; j < N; j++) {
        if ( diag[j] == toast::complex(0.,0.) )
	    continue;
 
	diag[j] = toast::complex(1.,0.)/diag[j];
 
	for (k = (j+1); k < N; k++) {
	  if (
	      ((tmp1=AC.Get(j,k)) == toast::complex(0.,0.)) ||
	      ((tmp2=AC.Get(k,j)) == toast::complex(0.,0.))
	      )
	    continue;
 
	  diag[k] -= tmp1*tmp2*diag[j];
	}
    }
 
    return diag;
}

inline CVector precond (const CGenericSparseMatrix& AC,const CVector& y)
{
    int j, k, N = y.Dim();
 
    // return y;
 
    /**************************
 
    //
    // Jacobi
    //
    
    CVector x(N);
 
    for (j = 0; j < N; j++) {
    if ( AC.Get(j,j)==complex(0.,0.) )
    continue;
 
    x[j] = y[j]/AC.Get(j,j);
    }
  
    return x;
 
    ***************************/
    /**************************
 
    //
    // SOR
    //
 
    complex sum, tmp;
    CVector x_new(N);
    for (j = 0; j < N; j++)
    {
      sum = 0.0;
 
      for (k = 0; k < N; k++)
      {
          if ( (tmp=AC.Get(j,k))==complex(0.0,0.0) )
             continue;
 
          if ( k <= (j-1) )
          {
              sum += tmp*x_new[k];
          }
 
          if ( k >= (j+1) )
          {
             sum += tmp*y[k];
          }
      }
 
      sum = (y[j] - sum)/AC.Get(j,j);
 
      x_new[j] = y[j] + complex(relax,0.)*(sum - y[j]);
      // x_new[j] = complex(relax,0.)*sum;
   }
 
   return x_new;
 
***************************/
/**************************/
 
   //
   // D-ILU
   //
 
   toast::complex sum, tmp;
   CVector z(N), x_new(N);
 
   for (j = 0; j < N; j++)
   {
      sum = toast::complex(0.,0.);
 
      for (k = 0; k < j; k++)
      {
         if ((tmp=AC.Get(j,k)) == toast::complex(0.,0.))
            continue;
 
         sum += tmp*z[k];
      }
 
      z[j] = diag_inv[j]*(y[j] - sum);
   }
 
   for (j = (N-1); j >= 0; j--)
   {
      sum = toast::complex(0.,0.);
      for (k = (j+1); k < N; k++)
      {
         if ((tmp=AC.Get(j,k)) == toast::complex(0.,0.))
            continue;
 
         sum += tmp*x_new[k];
      }
 
      x_new[j] = z[j] - diag_inv[j]*sum;
   }
 
   return x_new;
 
/***************************/
}
#endif

// ==========================================================================
// instantiations

#ifdef UNDEF // NEED_EXPLICIT_INSTANTIATION

template int gmres (int restart, const TMatrix<float> &A,
    const TVector<float> &b, TVector<float> &x,
    TPreconditioner<float> *precon, double &elim, int maxit,
    void (*clbk)(void*));

template int gmres (int restart, const TMatrix<double> &A,
    const TVector<double> &b, TVector<double> &x,
    TPreconditioner<double> *precon, double &elim, int maxit,
    void (*clbk)(void*));

template int gmres (int restart, const TMatrix<toast::complex> &A,
    const TVector<toast::complex> &b, TVector<toast::complex> &x,
    TPreconditioner<toast::complex> *precon, double &elim, int maxit,
    void (*clbk)(void*));

template int gmres (int restart, const TMatrix<scomplex> &A,
    const TVector<scomplex> &b, TVector<scomplex> &x,
    TPreconditioner<scomplex> *precon, double &elim, int maxit,
    void (*clbk)(void*));

template int gmres (int restart,
    TVector<double> (*Av_clbk)(const TVector<double> &v, void *context),
    const TVector<double> &b, TVector<double> &x,
    TPreconditioner<double> *precon, double &elim, int maxit, void *context);

template int gmres (int restart,
    TVector<toast::complex> (*Av_clbk)(const TVector<toast::complex> &v,
    void *context), const TVector<toast::complex> &b,
    TVector<toast::complex> &x, TPreconditioner<toast::complex> *precon,
    double &elim, int maxit, void *context);

#endif // NEED_EXPLICIT_INSTANTIATION
