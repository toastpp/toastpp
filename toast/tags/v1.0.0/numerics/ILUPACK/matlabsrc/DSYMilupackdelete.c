/* ========================================================================== */
/* === AMGdelete mexFunction ================================================ */
/* ========================================================================== */

/*
    Usage:

    release memory of the preconditioner obtained by ILUPACK V2.2
    
    Example:

    % for discarding the preconditioner
    DSYMilupackdelete(PREC);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	March 24, 2006. ILUPACK V2.2.  

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
    DAMGlevelmat *PRE;
    SAMGlevelmat *SPRE;
    DILUPACKparam *param;

    const char **fnames;       /* pointers to field names */
    mxArray    *tmp;
    char       *pdata, *input_buf;
    int        mrows, ncols, ifield, *classIDflags, buflen, status;
    int        NStructElems, nfields, ndim;
    size_t     sizebuf;
    double     dbuf;
    mxArray    *PRE_input, *PRE_output;


    if (nrhs != 1)
       mexErrMsgTxt("One input arguments required.");
    else if (nlhs > 0)
       mexErrMsgTxt("No output arguments.");
    else if (!mxIsStruct(prhs[0]))
       mexErrMsgTxt("Input must be a structure.");

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
	else if (!strcmp("n",fnames[ifield])) {
	   /* get size of the initial system */
	   tmp = mxGetFieldByNumber(PRE_input,0,ifield);
	   A.nr=A.nc=*mxGetPr(tmp);
	   A.ia=A.ja=NULL;
	   A.a=NULL;
	}
    }
    mxFree(fnames);


    /* finally release memory of the preconditioner */
    if (!(param->ipar[16]&DISCARD_MATRIX)) {
       if (PRE->issingle) {
	  SPRE=(SAMGlevelmat *)PRE;
	  free(SPRE->A.ia);
	  free(SPRE->A.ja);
	  free(SPRE->A.a);
       }
       else {
	  free(PRE->A.ia);
	  free(PRE->A.ja);
	  free(PRE->A.a);
       }
    }

    DSYMAMGdelete(&A,PRE,param);

    return;

}

