// ============================================================================
// TOAST V.14 (1999)  M. Schweiger and S.R. Arridge
// Routines derived (but not necessarily identical from NR
// ============================================================================

#define MATHLIB_IMPLEMENTATION

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "mathlib.h"

using namespace std;

// ============================================================================
// Contents of nrutil.c
// ============================================================================

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(const char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;

	/* allocate array of pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#else /* ANSI */
/* traditional - K&R */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl,long nrh,long ncl,long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl,long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;

	/* allocate array of pointers to rows */
	m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double*)));
	if (!m)	nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((unsigned int)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#endif /* ANSI */

// ============================================================================
// End of nrutil.c
// ============================================================================

// ============================================================================

double gammln (double xx)
{
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
			    24.01409824083091, -1.231739572450155,
			    0.1208650973866179e-2, -0.5395239384953e-5};
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x+0.5)* log (tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j]/++y;
    return -tmp + log (2.5066282746310005*ser/x);
}

// ============================================================================

double poidev (double xm)
{
    static double sq, alxm, g, oldm = -1.0;
    double em, t, y;

    if (xm < 12.0) {
        if (xm != oldm) {
	    oldm = xm;
	    g = exp (-xm);
	}
	em = -1;
	t = 1.0;
	do {
	    ++em;
	    t *= drand48();
	} while (t > g);
    } else {
        if (xm != oldm) {
	    oldm = xm;
	    sq = sqrt (2.0*xm);
	    alxm = log (xm);
	    g = xm * alxm - gammln (xm+1.0);
	}
	do {
	    do {
	        y = tan (Pi * drand48());
		em = sq * y + xm;
	    } while (em < 0.0);
	    em = floor (em);
	    t = 0.9 * (1.0+y*y) * exp (em*alxm - gammln (em+1.0) -g);
	} while (drand48() > t);
    }
    return em;
}

// ============================================================================

double bnldev (double pp, int n)
{
    int j;
    static int nold = -1;
    double am, em, g, angle, p, bnl, sq, t, y;
    static double pold = -1.0, pc, plog, pclog, en, oldg;

    p = (pp <= 0.5 ? pp : 1.0-pp);
    am = n*p;
    if (n < 25) {
        bnl = 0.0;
	for (j = 1; j <= n; j++)
	    if (drand48() < p) ++bnl;
    } else if (am < 1.0) {
        g = exp (-am);
	t = 1.0;
	for (j = 0; j <= n; j++) {
	    t *= drand48();
	    if (t < g) break;
	}
	bnl = (j <= n ? j : n);
    } else {
        if (n != nold) {
	    en = n;
	    oldg = gammln (en+1.0);
	    nold = n;
	}
	if (p != pold) {
	    pc = 1.0-p;
	    plog = log (p);
	    pclog = log (pc);
	    pold = p;
	}
	sq = sqrt (2.0*am*pc);
	do {
	    do {
	        angle = Pi*drand48();
		y = tan (angle);
		em = sq*y + am;
	    } while (em < 0.0 || em >= (en+1.0));
	    em = floor (em);
	    t = 1.2*sq*(1.0+y*y) * exp (oldg-gammln(em+1.0)
		-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
	} while (drand48() > t);
	bnl = em;
    }
    if (p != pp) bnl = n-bnl;
    return bnl;
}

// ============================================================================

double gasdev (double sigma)
{
    static int iset = 0;
    static double gset;
    float fac, r, v1, v2;

    if (iset == 0) {
	do {
	    v1 = (float)(2.0*drand48()-1.0);
	    v2 = (float)(2.0*drand48()-1.0);
	    r = v1*v1+v2*v2;
	} while (r >= 1.0);
	fac = (float)(sigma*sqrt(-2.0*log(r)/r));
	gset = v1*fac;
	iset = 1;
	return v2*fac;
    } else {
	iset = 0;
	return gset;
    }
}

// ============================================================================

double bessi0 (double x)
{
    double ax, ans, y;

    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
	y *= y;
	ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +
		    y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
    } else {
        y = 3.75/ax;
	ans=(exp(ax)/sqrt(ax)) * (0.39894228 + y*(0.1328592e-1 +
	     y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 +
             y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 +
	     y*0.392377e-2))))))));
    }
    return ans;
}

// ============================================================================

double bessi (int n, double x)
{
    const double ACC = 40.0;
    const double BIGNO = 1e10;
    const double BIGNI = 1e-10;

    double tox, bi, bip, bim, ans;
    int j;

    if (n < 2) {
        cerr << "Index n less than 2 in bessi" << endl;
	exit (1);
    }
    if (x == 0.0)
        return 0.0;
    else {
        tox = 2.0/fabs(x);
	bip = ans = 0.0;
	bi = 1.0;
	for (j = 2*(n+(int)sqrt(ACC*n)); j > 0; j--) {
	    bim = bip + j*tox * bi;
	    bip = bi;
	    bi = bim;
	    if (fabs(bi) > BIGNO) {
	        ans *= BIGNI;
		bi *= BIGNI;
		bip *= BIGNI;
	    }
	    if (j == n) ans = bip;
	}
	ans *= bessi0 (x) / bi;
	return x < 0.0 && (n & 1) ? -ans : ans;
    }
}

// ============================================================================

double erff (double x)
{
    return (x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x));
}

// ============================================================================

void fit (double x[], double y[], int ndata, double *a, double *b,
    double *siga, double *sigb, double *chi2)
{
    int i;
    double t, sxoss, sx=0.0, sy=0.0, st2=0.0, ss, sigdat;

    *b = 0.0;
    for (i=1; i<=ndata; i++) {
        sx += x[i];
	sy += y[i];
    }
    ss = ndata;
    sxoss = sx/ss;
    for (i=1; i<=ndata; i++) {
        t = x[i] - sxoss;
	st2 += t*t;
	*b += t*y[i];
    }
    *b /= st2;
    *a = (sy-sx*(*b))/ss;
    *siga = sqrt ((1.0+sx*sx/(ss*st2))/ss);
    *sigb = sqrt (1.0/st2);
    *chi2 = 0.0;
    for (i=1; i<=ndata; i++) {
        double tmp = y[i]-(*a)-(*b)*x[i];
	*chi2 += tmp*tmp;
    }
    sigdat = sqrt ((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
}

// ============================================================================

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], int nn, int isign)
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
	    SWAP(data[j],data[i]);
	    SWAP(data[j+1],data[i+1]);
	}
	m=n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=2*mmax;
	theta=6.28318530717959/(isign*mmax);
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0;
	wi=0.0;
	for (m=1;m<mmax;m+=2) {
	    for (i=m;i<=n;i+=istep) {
	        j=i+mmax;
		tempr=wr*data[j]-wi*data[j+1];
		tempi=wr*data[j+1]+wi*data[j];
		data[j]=data[i]-tempr;
		data[j+1]=data[i+1]-tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr=(wtemp=wr)*wpr-wi*wpi+wr;
	    wi=wi*wpr+wtemp*wpi+wi;
	}
	mmax=istep;
    }
}

#undef SWAP

// ============================================================================

double gammp (double a, double x)
{
    double gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0) nrerror ("Invalid arguments in routine gammp");
    if (x < (a+1.0)) {
        gser (&gamser, a, x, &gln);
	return gamser;
    } else {
        gcf (&gammcf, a, x, &gln);
	return 1.0-gammcf;
    }
}

// ============================================================================

void gcf (double *gammcf, double a, double x, double *gln)
{
    const int ITMAX = 100;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;

    int i;
    double an, b, c, d, del, h;

    *gln = gammln(a);
    b = x+1.0-a;
    c = 1.0/FPMIN;
    d = 1.0/b;
    h = d;
    for (i = 1; i <= ITMAX; i++) {
        an = -i*(i-a);
	b += 2.0;
	d = an*d + b;
	if (fabs (d) < FPMIN) d = FPMIN;
	c = b + an/c;
	if (fabs (c) < FPMIN) c = FPMIN;
	d = 1.0/d;
	del = d*c;
	h *= del;
	if (fabs (del-1.0) < EPS) break;
    }
    if (i > ITMAX) nrerror ("a too large, ITMAX too small in gcf");
    *gammcf = exp (-x+a*log(x)-(*gln))*h;
}

// ============================================================================

void gser (double *gamser, double a, double x, double *gln)
{
    const int ITMAX = 100;
    const double EPS = 3.0e-7;

    int n;
    double sum, del, ap;

    *gln = gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) nrerror ("x less than 0 in routine gser");
	*gamser = 0.0;
	return;
    } else {
        ap = a;
	del = sum = 1.0/a;
	for (n = 1; n <= ITMAX; n++) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs (del) < fabs (sum) * EPS) {
	        *gamser = sum * exp(-x+a*log(x)-(*gln));
		return;
	    }
	}
	nrerror ("a too large, ITMAX too small in routine gser");
	return;
    }
}

// ============================================================================

#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch (int n, double xold[], double fold, double g[], double p[],
    double x[], double *f, double stpmax, int *check, double (*func)(double[]))
{
    int i;
    double a, alam, alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope,
        sum, temp, test, tmplam;

    *check = 0;
    for (sum = 0.0, i = 1; i <= n; i++) sum += p[i]*p[i];
    sum = sqrt(sum);
    if (sum > stpmax)
        for (i = 1; i <= n; i++) p[i] *= stpmax/sum;
    for (slope = 0.0, i = 1; i <= n; i++)
        slope += g[i]*p[i];
    test = 0.0;
    for (i = 1; i <= n; i++) {
        temp = fabs (p[i])/max (fabs(xold[i]),1.0);
	if (temp > test) test = temp;
    }
    alamin = TOLX/test;
    alam = 1.0;
    for (;;) {
        for (i = 1; i <= n; i++) x[i] = xold[i]+alam*p[i];
	*f = (*func)(x);
	if (alam < alamin) {
	    for (i = 1; i <= n; i++) x[i] = xold[i];
	    *check = 1;
	    return;
	} else if (*f <= fold+ALF*alam*slope) return;
	else {
	    if (alam == 1.0)
	        tmplam = -slope/(2.0*(*f-fold-slope));
	    else {
	        rhs1 = *f-fold-alam*slope;
		rhs2 = f2-fold2-alam2*slope;
		a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
		b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/
		    (alam-alam2);
		if (a == 0.0) tmplam = -slope/(2.0*b);
		else {
		    disc = b*b-3.0*a*slope;
		    if (disc < 0.0) cerr << "Roundoff problem in lnsrch\n";
		    else tmplam = (-b+sqrt(disc))/(3.0*a);
		}
		if (tmplam > 0.5*alam)
		    tmplam = 0.5*alam;
	    }
	}
	alam2 = alam;
	f2 = *f;
	fold2 = fold;
	alam = max (tmplam,0.1*alam);
    }
}

#undef ALF
#undef TOLX

// ============================================================================

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);

void dfpmin (double p[], int n, double gtol, int *iter, double *fret,
    double(*func)(double[]), void (*dfunc)(double[], double[]))
{
    int check, i, its, j;
    double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, temp, test;
    double *dg, *g, *hdg, **hessin, *pnew, *xi;

    dg = vector (1,n);
    g = vector (1,n);
    hdg = vector (1,n);
    hessin = matrix (1,n,1,n);
    pnew = vector (1,n);
    xi = vector (1,n);
    fp = (*func)(p);
    cout << "initial objective = " << fp << endl;
    (*dfunc)(p,g);
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) hessin[i][j] = 0.0;
	hessin[i][i] = 1.0;
	xi[i] = -g[i];
	sum += p[i]*p[i];
    }
    stpmax = STPMX * max (sqrt(sum),(double)n);

    for (its = 1; its <= ITMAX; its++) {
        cout << "Starting iteration " << its << endl;

        *iter = its;
	lnsrch (n, p, fp, g, xi, pnew, fret, stpmax, &check, func);
	// The new function evaluation occurs in lnsrch; save the function
	// value in fp for the next line search. It is usually safe to
	// ignore the value of check.
	fp = *fret;
	cout << "objective = " << fp << endl;
	for (i = 1; i <= n; i++) {
	    xi[i] = pnew[i]-p[i];
	    p[i] = pnew[i];
	}
	test = 0.0;
	for (i = 1; i <= n; i++) {
	    temp = fabs(xi[i])/max (fabs(p[i]),1.0);
	    if (temp > test) test = temp;
	}
	if (test < TOLX) {
	    FREEALL
	    return;
	}
	for (i = 1; i <= n; i++) dg[i] = g[i];
	(*dfunc)(p,g);
	test = 0.0;
	den = max (*fret,1.0);
	for (i = 1; i <= n; i++) {
	    temp = fabs(g[i]) * max(fabs(p[i]),1.0)/den;
	    if (temp > test) test = temp;
	}
	if (test < gtol) {
	    FREEALL;
	    return;
	}
	for (i = 1; i <= n; i++) dg[i] = g[i]-dg[i];
	for (i = 1; i <= n; i++) {
	    hdg[i] = 0.0;
	    for (j = 1; j <= n; j++) hdg[i] += hessin[i][j]*dg[j];
	}
	fac = fae = sumdg = sumxi = 0.0;
	for (i = 1; i <= n; i++) {
	    fac += dg[i] * xi[i];
	    fae += dg[i] * hdg[i];
	    sumdg += sqr (dg[i]);
	    sumxi += sqr (xi[i]);
	}
	if (fac > sqrt (EPS*sumdg*sumxi)) {
	    fac = 1.0/fac;
	    fad = 1.0/fae;
	    for (i = 1; i <= n; i++) dg[i] = fac*xi[i]-fad*hdg[i];
	    for (i = 1; i <= n; i++) {
	        for (j = 1; j <= n; j++) {
		    hessin[i][j] += fac*xi[i]*xi[j]
		      -fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
		    hessin[j][i] = hessin[i][j];
		}
	    }
	}
	for (i = 1; i <= n; i++) {
	    xi[i] = 0.0;
	    for (j = 1; j <= n; j++) xi[i] -= hessin[i][j]*g[j];
	}
    }
    cerr << "too many iterations in dfpmin" << endl;
    FREEALL
}     

#undef ITMAX
#undef EPS

// ============================================================================

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn (double data[], unsigned long nn[], int ndim, int isign)
{
    int idim;
    unsigned long i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
    unsigned long ibit, k1, k2, n, nprev, nrem, ntot;
    double tempi, tempr;
    double theta, wi, wpi, wpr, wr, wtemp;

    for (ntot = 1, idim = 1; idim <= ndim; idim++)
        ntot *= nn[idim-1];
    nprev = 1;
    for (idim = ndim; idim >= 1; idim--) {
        n = nn[idim-1];
	nrem = ntot/(n*nprev);
	ip1 = nprev << 1;
	ip2 = ip1 * n;
	ip3 = ip2 * nrem;
	i2rev = 1;
	for (i2 = 1; i2 <= ip2; i2 += ip1) {
	    if (i2 < i2rev) {
	        for (i1 = i2; i1 <= i2+ip1-2; i1 += 2) {
		    for (i3 = i1; i3 <= ip3; i3 += ip2) {
		        i3rev = i2rev+i3-i2;
			SWAP(data[i3-1],data[i3rev-1]);
			SWAP(data[i3], data[i3rev]);
		    }
		}
	    }
	    ibit = ip2 >> 1;
	    while (ibit >= ip1 && i2rev > ibit) {
	        i2rev -= ibit;
		ibit >>= 1;
	    }
	    i2rev += ibit;
	}
	ifp1 = ip1;
	while (ifp1 < ip2) {
	    ifp2 = ifp1 << 1;
	    theta = isign * 6.28318530717959/(ifp2/ip1);
	    wtemp = sin(0.5*theta);
	    wpr = -2.0*wtemp*wtemp;
	    wpi = sin(theta);
	    wr = 1.0;
	    wi = 0.0;
	    for (i3 = 1; i3 <= ifp1; i3 += ip1) {
	        for (i1 = i3; i1 <= i3+ip1-2; i1 += 2) {
		    for (i2 = i1; i2 <= ip3; i2 += ifp2) {
		        k1 = i2;
			k2 = k1+ifp1;
			tempr = wr*data[k2-1]-wi*data[k2];
			tempi = wr*data[k2]+wi*data[k2-1];
			data[k2-1] = data[k1-1] - tempr;
			data[k2] = data[k1] - tempi;
			data[k1-1] += tempr;
			data[k1] += tempi;
		    }
		}
		wr = (wtemp=wr)*wpr-wi*wpi+wr;
		wi = wi*wpr+wtemp*wpi+wi;
	    }
	    ifp1 = ifp2;
	}
	nprev *= n;
    }
}

void spline (double x[], double y[], int n, double yp1, double ypn,
    double y2[])
{
    int i, k;
    double p, qn, sig, un, *u;
    u = new double[n-1];
    
    if (yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
	u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for (i = 1; i < n-1; i++) {
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
	p = sig*y2[i-1]+2.0;
	y2[i] = (sig-1.0)/p;
	u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
	u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
        qn = un = 0.0;
    else {
        qn = 0.5;
	un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k = n-2; k >= 0; k--)
        y2[k] = y2[k] * y2[k+1]+u[k];

    delete []u;
}

void splint (double xa[], double ya[], double y2a[], int n, double x,
    double *y)
{
    int klo, khi, k;
    double h, b, a;

    klo = 1;
    khi = n;
    while (khi-klo > 1) {
        k = (khi+klo) >> 1;
	if (xa[k-1] > x) khi = k;
	else klo = k;
    }
    h = xa[khi-1]-xa[klo-1];
    if (h == 0.0) cerr << "Bad xa input to routine splint" << endl;
    a = (xa[khi-1]-x)/h;
    b = (x-xa[klo-1])/h;
    *y = a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+
				  (b*b*b-b)*y2a[khi-1])*(h*h)/6.0;
}
