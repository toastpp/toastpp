/* ========================================================================== */
/* === AMGconvert mexFunction =============================================== */
/* ========================================================================== */

/*
    Usage:

    Return the structure 'options' and solution x of Ax=b based on ILUPACK V2.2
    
    Example:

    % for initializing parameters
    PREC = DSYMSPDilupackconvert(PREC);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	April 22, 2009. ILUPACK V2.3

    Notice:

	Copyright (c) 2009 by TU Braunschweig.  All Rights Reserved.

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
    DAMGlevelmat *PRE;
    DILUPACKparam *param;
    integer n;

    const char **fnames;
    const char *pnames[]= {"n","nB", "L","D","U", "E","F", "rowscal","colscal", "p","invq",
			   "param","ptr", "isreal","isdefinite","issymmetric","ishermitian"};


    const int  *dims;
    mxArray    *tmp, *fout;
    char       *pdata, *input_buf, *output_buf;
    int        mrows, ncols, ifield, jstruct, *classIDflags, buflen, status;
    int        NStructElems, nfields, ndim,nnz, ierr,i,j,k,l,m, *irs, *jcs;
    size_t     sizebuf;
    double     dbuf, *A_valuesR, *convert, *sr, *pr, *sol, *rhs;
    mxArray    *A_input , *b_input, *x0_input;
    mxArray    *PRE_input;
    int *A_ja ;                 /* row indices of input matrix A */
    int *A_ia ;                 /* column pointers of input matrix A */
    

    if (nrhs != 1)
       mexErrMsgTxt("one input argument required.");
    else if (nlhs !=1)
       mexErrMsgTxt("Too many output arguments.");


    /* import pointer to the preconditioner */
    PRE_input = (mxArray*) prhs [0] ;
    /* get number of levels of input preconditioner structure `PREC' */
    /* nlev=mxGetN(PRE_input); */

    nfields = mxGetNumberOfFields(PRE_input);
    /* allocate memory  for storing pointers */
    fnames = mxCalloc(nfields, sizeof(*fnames));
    for (ifield = 0; ifield < nfields; ifield++) {
        fnames[ifield] = mxGetFieldNameByNumber(PRE_input,ifield);
	/* check whether `PREC.ptr' exists */
	if (!strcmp("ptr",fnames[ifield])) {
	   /* field `ptr' */
	   tmp = mxGetFieldByNumber(PRE_input,0,ifield);
	   pdata = mxGetData(tmp);
	   memcpy(&PRE, pdata, sizeof(size_t));
	}
	else if (!strcmp("param",fnames[ifield])) {
	   /* field `param' */
	   tmp = mxGetFieldByNumber(PRE_input,0,ifield);
	   pdata = mxGetData(tmp);
	   memcpy(&param, pdata, sizeof(size_t));
	}
    }
    mxFree(fnames);

    DSYMSPDAMGconvert(PRE);

    /* Create a struct matrices for output */
    nlhs=1;

    plhs[0] = PRE_input; /*mxCreateStructMatrix(1, nlev, 17, pnames);*/

    return;
}

