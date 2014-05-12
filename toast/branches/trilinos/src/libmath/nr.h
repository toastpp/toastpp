// -*-C++-*-
// ============================================================================
// TOAST V.14 (1999)  M. Schweiger and S.R. Arridge
// Routines derived (but not necessarily identical) from NR
// ============================================================================

// ============================================================================
// Contents of nrutil.h
// ============================================================================

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include "mathdef.h"

#ifdef UNDEF
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
// replaced by sqr (mathdef.h)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
// replaced by sqr (mathdef.h)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))
// replaced by max (mathdef.h)

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))
// replaced by min (mathdef.h)

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
         (maxarg1) : (maxarg2))
// replaced by max (mathdef.h)

static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
         (minarg1) : (minarg2))
// replaced by min (mathdef.h)

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))
// replaced by max (mathdef.h)

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
       (lminarg1) : (lminarg2))
// replaced by min (mathdef.h)

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
       (imaxarg1) : (imaxarg2))
// replaced by max (mathdef.h)

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
// replaced by min (mathdef.h)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
// replaced by sign (mathdef.h)

#endif // UNDEF

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(const char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
double *vector();
double **matrix();
double **submatrix();
double **convert_matrix();
double ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */

// ============================================================================
// End of nrutil.h
// ============================================================================

double gammln (double xx);
// Returns the logarithm of the gamma function, ln (Gamma (xx))

double poidev (double xm);
// Returns as a floating point number an integer value that is a random
// deviate drawn from a Poisson distribution of mean xm, using drand48
// as a source of uniform random deviates
// NR p. 294

double bnldev (double pp, int n);
// Returns as a floating point number an integer value that is a random
// deviate drawn from a binomial distribution on n trials each of
// probability pp, using drand48 as a source of uniform random deviates
// NR p. 295

double gasdev (double sigma);
// Returns as a floating point number an integer value that is a random
// derivate drawn from a Gaussian distribution ...

double bessi0 (double x);
// Returns the modified Bessel function I_0(x) for any real x
// NR p. 237

MATHLIB double bessi (int n, double x);
// Returns the modified Bessel function I_n(x) for any real x and n >= 2
// NR p. 239

double erff (double x);
// Returns the error function erf(x)
// NR p. 220

void fit (double x[], double y[], int ndata, double *a, double *b,
    double *siga, double *sigb, double *chi2);

void four1(double data[], int nn, int isign);
// Replaces data[] by its discrete Fourier transform, if isign is input as 1,
// or replaces data[] by nn times its inverse Fourier transform, if isign is
// input as -1. data is a complex array o length nn, or equivalently a real
// array of length 2*nn. nn must be an integer power of 2
// NR p. 507

double gammp (double a, double x);
// Returns the incomplete gamma function P(a,x)
// NR p. 218

void gcf (double *gammcf, double a, double x, double *gln);
// Returns the incomplete gamma function Q(a,x) evaluated by its continued
// fraction representation as gammcf. Also returns ln Gamma(a) as gln
// NR p. 219

void gser (double *gamser, double a, double x, double *gln);
// Returns the incomplete gamma function P(a,x) evaluated by its series
// representation as gamser. Also returns ln Gamma(a) as gln
// NR p. 218

void lnsrch (int n, double xold[], double fold, double g[], double p[],
    double x[], double *f, double stpmax, int *check,
    double (*func)(double[]));
// Given an n-dimensional point xold[1..n], the value of the function and
// gradient there, fold, and g[1..n], and a direction p[1..n], finds a new
// point x[1..n] along the direction p from xold where the function func
// has decreased "sufficiently". The new function value is returned in f.
// stpmax is an input quantity that limits the length of the steps so that
// you do not try to evaluate the function in regions where it is undefined
// or subject to overflow. p is usually the Newton direction. The output
// quantity check is false (0) on a normal exit. It is true (1) when x is
// too close to xold. In a minimization algorithm, this usually signals
// convergence and can be ignored. However, in a zero-finding algorithm the
// calling program should check whether the convergence is spurious. Some
// "difficult" problems may require double precision in this routine.
// NR p. 385

void dfpmin (double p[], int n, double gtol, int *iter, double *fret,
    double(*func)(double[]), void (*dfunc)(double[], double[]));
// Given a starting point p[1..n] that is a vector of length n, the Broyden-
// Fletcher-Goldfarb-Shanno variant of Davidon-Fletcher-Powell minimization
// is performed on a function func, using its gradient as calculated by a
// routine dfunc. The convergence requirement on zeroing the gradient is
// input as gtol. Returned quantities are p[1..n] (the location of the
// minimum), iter (the number of iterations that were performed), and fret
// (the minimum value of the function). The routine lnsrch is called to
// perform approximate line minimizations.

MATHLIB void fourn (double data[], unsigned long nn[], int ndim, int isign);
// Replaces data by its ndim-dimensional discrete Fourier transform, if
// isign is input as 1. nn[0..ndim-1] is an integer array containing the
// lenghts of each dimension (number of complex values), which must be
// powers of 2. data is a real array of length twice the product of these
// lenghts, in which the data are stored as in a multidimensional complex
// array: real and imaginary parts of each element are in consecutive
// locations, and the rightmost index of the array increases most rapidly
// as one proceeds along data. For a two-dimensional array, this is
// equivalent to storing the array by rows. If isisgn is input as -1, data
// is replaced by its inverse transform times the product of the lenghts
// of all dimensions.

void spline (double x[], double y[], int n, double yp1, double ypn,
    double y2[]);
// Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function,
// i.e. y_i = f(x_i), with x_0 < x_1 < ... < x_{N-1}, and given values
// yp1 and ypn for the first derivative of the interpolating function at
// points 0 and n-1, respectively, this routine returns an array y2[0..n-1]
// that contains the second derivatives of the interpolating function at
// the tabulated points x_i. If yp1 and/or ypn are equal to 1e30 or larger,
// the routine is signaled to set the corresponding boundary condition for
// a natural spline, with zero second derivative on that boundary.

void splint (double xa[], double ya[], double y2a[], int n, double x,
    double *y);
// Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function
// (with the xa_i's in order), and given the array y2a[0..n-1], which is
// the output from spline above, and given a value of x, this routine
// returns a cubic-spline interpolated value y.
