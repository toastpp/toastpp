/* ========================================================================== */
/* === mwmilupackmetise mexFunction ============================================ */
/* ========================================================================== */

/*
    Usage:

    Return METIS Nested Dissection by edges reordering combined
    with maximum weight matching
    
    Example:

    % for initializing parameters
    [pl,pr,Dl,Dr] = mwmilupackmetise(A);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	October 06, 2010. ILUPACK V2.3.  

    Notice:

	Copyright (c) 2010 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://ilupack.tu-bs.de/
*/

/* ========================================================================== */
/* === Include files and prototypes ========================================= */
/* ========================================================================== */

#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include <ilupack.h>

#define MAX_FIELDS 100

/* ========================================================================== */
/* === mexFunction ========================================================== */
/* ========================================================================== */

void mexFunction
(
    /* === Parameters ======================================================= */

    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs []	/* right-hand side matrices */
)
{
    Dmat A;
    DILUPACKparam param;
    mxArray *A_input;
    integer *p, *invq, nB=0;
    double  *prowscale, *pcolscale;
    int     ierr, i,j,k,l;
    size_t  mrows, ncols; 
    mwSize  nnz;
    double  *pr, *D, *A_a;
    mwIndex *A_ia, *A_ja;

    if (nrhs != 1)
       mexErrMsgTxt("One input argument required.");
    else if (nlhs!=4)
       mexErrMsgTxt("Four output arguments are required.");
    else if (!mxIsNumeric(prhs[0]))
       mexErrMsgTxt("Input must be a matrix.");
    else if (!mxIsSparse(prhs[0]))
       mexErrMsgTxt("Matrix must be sparse.");

    /* The first input must be a square matrix.*/
    A_input=(mxArray *)prhs[0];
    mrows = mxGetM(A_input);
    ncols = mxGetN(A_input);
    nnz = mxGetNzmax(A_input);
    if (mrows!=ncols) {
      mexErrMsgTxt("Input must be a square matrix.");
    }
    A_ja         = (mwIndex *)mxGetIr(A_input);
    A_ia         = (mwIndex *)mxGetJc(A_input) ;
    A_a          = (double *) mxGetPr(A_input) ;


    A.nc=A.nr=mrows;
    A.ia=(integer *)MAlloc((size_t)(A.nc+1)*sizeof(integer),"mwmilupackmetise");
    A.ja=(integer *)MAlloc((size_t)nnz     *sizeof(integer),"mwmilupackmetise");
    A.a =(double *) MAlloc((size_t)nnz     *sizeof(double), "mwmilupackmetise");
    A.ia[0]=1;
    for (i = 0 ; i < ncols ; i++) {
        A.ia[i+1]=A.ia[i];
	for (j = A_ia[i] ; j < A_ia[i+1] ; j++) {
  	    /* a_ik */
	    k=A_ja[j];
	    l=A.ia[i+1]-1;
	    A.ja[l]=k+1;
	    A.a[l]=A_a[j];
	    A.ia[i+1]=l+2;
	}
    }

    /* init parameters to their default values */
    DGNLAMGinit(&A,&param);

    p   =(integer *) MAlloc((size_t)A.nc*sizeof(integer),"mwmilupackmetise");
    invq=(integer *) MAlloc((size_t)A.nc*sizeof(integer),"mwmilupackmetise");
    prowscale=(double *) MAlloc((size_t)A.nc*sizeof(double),"mwmilupackmetise");
    pcolscale=(double *) MAlloc((size_t)A.nc*sizeof(double),"mwmilupackmetise");

#ifdef _MC64_MATCHING_    
    ierr=DGNLperm_mc64_metis_e(A, prowscale, pcolscale, p, invq, &nB, &param);
#elif defined _PARDISO_MATCHING_
    ierr=DGNLperm_mwm_metis_e(A, prowscale, pcolscale, p, invq, &nB, &param);
#else /* MUMPS matching */
    ierr=DGNLperm_matching_metis_e(A, prowscale, pcolscale, p, invq, &nB, &param);
#endif



    /* Create output vector */
    nlhs=4;

    plhs[0] =mxCreateDoubleMatrix(1,mrows, mxREAL);
    pr = (double *) mxGetPr(plhs[0]);
    for (i=0;i<mrows; i++)
        pr[i]=p[i];


    plhs[1] =mxCreateDoubleMatrix(1,mrows, mxREAL);
    pr = (double *) mxGetPr(plhs[1]);
    for (i=0;i<mrows; i++)
        pr[invq[i]-1]=i+1;


    plhs[2] =mxCreateDoubleMatrix(mrows,1, mxREAL);
    D = (double *) mxGetPr(plhs[2]);
    for (i=0;i<mrows; i++)
        D[i]=prowscale[i];


    plhs[3] =mxCreateDoubleMatrix(mrows,1, mxREAL);
    D = (double *) mxGetPr(plhs[3]);
    for (i=0;i<mrows; i++)
        D[i]=pcolscale[i];

    free(A.ia);
    free(A.ja);
    free(A.a);

    free(p);
    free(invq);
    free(prowscale);
    free(pcolscale);
    return;
}

