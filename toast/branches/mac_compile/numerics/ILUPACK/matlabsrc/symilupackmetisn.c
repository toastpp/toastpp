/* ========================================================================== */
/* === symilupackmetisn mexFunction ========================================= */
/* ========================================================================== */

/*
    Usage:

    Return METIS multilevel reordering by nodes
    
    Example:

    % for initializing parameters
    p = symilupackmetisn(A);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	March 09, 2008. ILUPACK V2.2.  

    Acknowledgements:

	This work was supported from 2002 to 2007 by the DFG research center
        MATHEON "Mathematics for key technologies"

    Notice:

	Copyright (c) 2008 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://www.math.tu-berlin.de/ilupack/
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
    double  *prowscale=NULL, *pcolscale=NULL;
    int     ierr, i,j,k,l;
    size_t  mrows, ncols; 
    mwSize  nnz;
    double  *pr;
    mwIndex *A_ia, *A_ja;

    if (nrhs != 1)
       mexErrMsgTxt("One input argument required.");
    else if (nlhs > 1)
       mexErrMsgTxt("Too many output arguments.");
    else if (!mxIsNumeric(prhs[0]))
       mexErrMsgTxt("First input must be a matrix.");

    /* The first input must be a square matrix.*/
    A_input=(mxArray *)prhs[0];
    mrows = mxGetM(A_input);
    ncols = mxGetN(A_input);
    nnz = mxGetNzmax(A_input);
    if (mrows!=ncols) {
      mexErrMsgTxt("First input must be a square matrix.");
    }
    A_ja         = (mwIndex *)mxGetIr(A_input);
    A_ia         = (mwIndex *)mxGetJc(A_input) ;


    A.nc=A.nr=mrows;
    A.ia=(integer *) MAlloc((size_t)(A.nc+1)*sizeof(integer),"symilupackmetisn");
    A.ja=(integer *) MAlloc((size_t)nnz     *sizeof(integer),"symilupackmetisn");
    A.a=NULL;
    A.ia[0]=1;
    for (i = 0 ; i < ncols ; i++) {
        A.ia[i+1]=A.ia[i];
	for (j = A_ia[i] ; j < A_ia[i+1] ; j++) {
	    k=A_ja[j];
	    if (k>=i) {
	       l=A.ia[i+1]-1;
	       A.ja[l]=k+1;
	       A.ia[i+1]=l+2;
	    }
	}
    }

    /* init parameters to their default values */
    DSPDAMGinit(&A,&param);
    param.ipar[7]=0;
    param.ipar[8]=0;

    p   =(integer *) MAlloc((size_t)A.nc*sizeof(integer),"symilupackmetisn");
    invq=(integer *) MAlloc((size_t)A.nc*sizeof(integer),"symilupackmetisn");
    ierr=DGNLperm_metis_n(A, prowscale, pcolscale, p, invq, &nB, &param);


    /* Create output vector */
    nlhs=1;
    plhs[0] =mxCreateDoubleMatrix(1,mrows, mxREAL);
    pr = (double *) mxGetPr(plhs[0]);
    for (i=0;i<mrows; i++)
        pr[i]=p[i];

    free(p);
    free(invq);
    free(A.ia);
    free(A.ja);
    return;
}

