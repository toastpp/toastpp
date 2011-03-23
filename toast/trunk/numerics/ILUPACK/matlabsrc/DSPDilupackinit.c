/* ========================================================================== */
/* === AMGinit mexFunction ================================================== */
/* ========================================================================== */

/*
    Usage:

    Return the structure 'option' for ILUPACK V2.3
    
    Example:

    % for initializing parameters
    options = DSPDilupackinit(A,options);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	January 22, 2009. ILUPACK V2.3.  

    Notice:

	Copyright (c) 2006 by TU Braunschweig.  All Rights Reserved.

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

    const char **fnames;       /* pointers to field names */
    const mwSize  *dims;
    mxClassID  *classIDflags;
    mwSize     buflen, mrows, ncols, ndim;
    mxArray    *tmp, *fout, *A_input, *options_input, *options_output;
    char       *pdata, *input_buf, *output_buf;
    int        ifield, status, nfields;
    size_t     sizebuf;
    double     dbuf;


    if (nrhs != 2)
       mexErrMsgTxt("Two input arguments required.");
    else if (nlhs > 1)
       mexErrMsgTxt("Too many output arguments.");
    else if (!mxIsStruct(prhs[1]))
       mexErrMsgTxt("Second input must be a structure.");
    else if (!mxIsNumeric(prhs[0]))
       mexErrMsgTxt("First input must be a matrix.");

    /* The first input must be a square matrix.*/
    A_input=(mxArray *)prhs[0];
    mrows = mxGetM(A_input);
    ncols = mxGetN(A_input);
    if (mrows!=ncols) {
      mexErrMsgTxt("First input must be a square matrix.");
    }

    A.nc=A.nr=mrows;
    A.ia=A.ja=NULL;
    A.a=NULL;
    DSPDAMGinit(&A,&param);

    /* Get second input arguments */
    options_input=(mxArray *)prhs[1];
    nfields = mxGetNumberOfFields(options_input);

    /* Allocate memory  for storing pointers */
    fnames = mxCalloc((size_t)nfields, (size_t)sizeof(*fnames));

    /* Allocate memory  for storing classIDflags */
    classIDflags = (mxClassID *) mxCalloc((size_t)nfields, (size_t)sizeof(mxClassID));
    
    /* Get field name pointers */
    for (ifield = 0; ifield < nfields; ifield++) {
        fnames[ifield] = mxGetFieldNameByNumber(options_input,ifield);
    }

    /* Create a 1x1 struct matrix for output */
    nlhs=1;
    plhs[0] = mxCreateStructMatrix((mwSize)1,(mwSize)1, nfields, fnames);
    options_output=plhs[0];

    /* copy data */
    for (ifield = 0; ifield < nfields; ifield++) {
	tmp = mxGetFieldByNumber(options_input,(mwIndex)0,ifield);
	classIDflags[ifield] = mxGetClassID(tmp); 

	ndim = mxGetNumberOfDimensions(tmp);
	dims = mxGetDimensions(tmp);

	/* Create string/numeric array */
	if (classIDflags[ifield] == mxCHAR_CLASS) {
 	   /* Get the length of the input string. */
	   buflen = (mxGetM(tmp) * mxGetN(tmp)) + 1;

	   if (!strcmp("amg",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.amg)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param.amg);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("presmoother",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.presmoother)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param.presmoother);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("postsmoother",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.postsmoother)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param.postsmoother);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("typecoarse",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.typecoarse)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param.typecoarse);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("typetv",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.typetv)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param.typetv);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("FCpart",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.FCpart)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param.FCpart);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("solver",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.solver)+1, (size_t)sizeof(char));
	      strcpy(output_buf, param.solver);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("ordering",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param.ordering)+1, (size_t)sizeof(char));
	      strcpy(output_buf, param.ordering);
	      fout = mxCreateString(output_buf);
	   }
	   else {
	      /* Allocate memory for input and output strings. */
	      input_buf  = (char *) mxCalloc((size_t)buflen, (size_t)sizeof(char));
	      output_buf = (char *) mxCalloc((size_t)buflen, (size_t)sizeof(char));
	      
	      /* Copy the string data from tmp into a C string 
		 input_buf. */
	      status = mxGetString(tmp, input_buf, buflen);
	      
	      sizebuf = (size_t)buflen*sizeof(char);
	      memcpy(output_buf, input_buf, sizebuf);
	      fout = mxCreateString(output_buf);
	   }
	} 
	else {
	   fout = mxCreateNumericArray(ndim, dims, 
				       classIDflags[ifield], mxREAL);
	   pdata = mxGetData(fout);

	   sizebuf = mxGetElementSize(tmp);
	   if (!strcmp("elbow",fnames[ifield])) {
	      dbuf=param.elbow;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("lfilS",fnames[ifield])) {
	      dbuf=param.lfilS;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("lfil",fnames[ifield])) {
	      dbuf=param.lfil;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("maxit",fnames[ifield])) {
	      dbuf=param.maxit;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("droptolS",fnames[ifield])) {
	      dbuf=param.droptolS;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("droptol",fnames[ifield])) {
	      dbuf=param.droptol;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("condest",fnames[ifield])) {
	      dbuf=param.condest;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("restol",fnames[ifield])) {
	      dbuf=param.restol;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("npresmoothing",fnames[ifield])) {
	      dbuf=param.npresmoothing;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("npostmoothing",fnames[ifield])) {
	      dbuf=param.npostsmoothing;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("ncoarse",fnames[ifield])) {
	      dbuf=param.ncoarse;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("matching",fnames[ifield])) {
	      dbuf=param.matching;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("nrestart",fnames[ifield])) {
	      dbuf=param.nrestart;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("damping",fnames[ifield])) {
	      dbuf=param.damping;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      dbuf=param.mixedprecision;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else {
	      memcpy(pdata, mxGetData(tmp), sizebuf);
	   }
	}


	/* Set each field in output structure */
	mxSetFieldByNumber(options_output, (mwIndex)0, ifield, fout);
    }

    mxFree(fnames);
    return;
}

