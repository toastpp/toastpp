/* ========================================================================== */
/* === AMGsol mexFunction =================================================== */
/* ========================================================================== */

/*
    Usage:

    solve Ax=b approximately by solving the system with one step of the ILUPACK
    preconditioner PREC
    
    Example:

    % for initializing parameters
    x=ZHPDilupacksol(PREC,b);


    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	October 16, 2006. ILUPACK V2.2

    Acknowledgements:

	This work was supported from 2002 to 2007 by the DFG research center
        MATHEON "Mathematics for key technologies"

    Notice:

	Copyright (c) 2006 by TU Braunschweig.  All Rights Reserved.

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

#define MAX(A,B) (((A)>=(B))?(A):(B))

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
    ZAMGlevelmat *PRE;
    CAMGlevelmat *SPRE;
    ZILUPACKparam *param;
    integer n;

    const char **fnames;

    mxArray    *PRE_input, *b_input, *tmp, *fout;
    int        i,j,k,l,m, ifield,nfields;
    char       *pdata;
    double     *pr, *pi;
    doublecomplex *dbuff, *sol, *rhs;

    if (nrhs != 2)
       mexErrMsgTxt("Two input arguments required.");
    else if (nlhs !=1)
       mexErrMsgTxt("One output argument is required.");



    /* import pointer to the preconditioner */
    PRE_input = (mxArray*) prhs [0] ;
    /* get number of levels of input preconditioner structure `PREC' */
    /* nlev=mxGetN(PRE_input); */


    nfields = mxGetNumberOfFields(PRE_input);
    /* allocate memory  for storing pointers */
    fnames = mxCalloc((size_t)nfields, (size_t)sizeof(*fnames));
    for (ifield = 0; ifield < nfields; ifield++) {
        fnames[ifield] = mxGetFieldNameByNumber(PRE_input,ifield);
	/* check whether `PREC.ptr' exists */
	if (!strcmp("ptr",fnames[ifield])) {
	   /* field `ptr' */
	   tmp = mxGetFieldByNumber(PRE_input,0,ifield);
	   pdata = mxGetData(tmp);
	   memcpy(&PRE, pdata, (size_t)sizeof(size_t));
	}
	else if (!strcmp("param",fnames[ifield])) {
	   /* field `param' */
	   tmp = mxGetFieldByNumber(PRE_input,0,ifield);
	   pdata = mxGetData(tmp);
	   memcpy(&param, pdata, (size_t)sizeof(size_t));
	}
    }
    mxFree(fnames);



    /* copy right hand side `b' */
    b_input = (mxArray *) prhs [1] ;
    pr=mxGetPr(b_input);
    n=mxGetM(b_input);

    /* make sure that enough buffer memory is available */
    param->ndbuff=MAX(param->ndbuff,5*(size_t)n);
    param->dbuff=(doublecomplex*)ReAlloc(param->dbuff,
					 (size_t)param->ndbuff*sizeof(doublecomplex),
					 "ZHPDilupacksol:rhs");

    rhs=param->dbuff;
    if (!mxIsComplex(b_input)) {
       for (i=0; i<n; i++) {
	   rhs[i].r=pr[i];
	   rhs[i].i=0;;
       }
    }
    else {
       pi=mxGetPi(b_input);
       for (i=0; i<n; i++) {
	   rhs[i].r=pr[i];
	   rhs[i].i=pi[i];
       }
    }


    /* rescale right hand side */
    if (PRE->issingle) {
       SPRE=(CAMGlevelmat *)PRE;
       for (i=0; i <n; i++) {
	   rhs[i].r=rhs[i].r*(double)SPRE->rowscal[i].r;
	   rhs[i].i=rhs[i].i*(double)SPRE->rowscal[i].r;
       }
    }
    else {
       for (i=0; i <n; i++) {
	   rhs[i].r=rhs[i].r*PRE->rowscal[i].r;
	   rhs[i].i=rhs[i].i*PRE->rowscal[i].r;
       }
    }

    /* provide memory for the approximate solution */
    sol=param->dbuff+n;
    /* 3n spaces as buffer */
    dbuff=param->dbuff+2*n;
  
    ZHPDAMGsol_internal(PRE, param, rhs, sol, dbuff);

    /* Create a struct matrices for output */
    nlhs=1;

    plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)1, mxCOMPLEX);
    pr=mxGetPr(plhs[0]);
    pi=mxGetPi(plhs[0]);
    /* rescale approximate solution */
    if (PRE->issingle) {
       for (i=0; i <n; i++) {
	   pr[i]=sol[i].r*(double)SPRE->colscal[i].r;
	   pi[i]=sol[i].i*(double)SPRE->colscal[i].r;
       }
    }
    else {
       for (i=0; i <n; i++) {
	   pr[i]=sol[i].r*PRE->colscal[i].r;
	   pi[i]=sol[i].i*PRE->colscal[i].r;
       }
    }

    return;
}

