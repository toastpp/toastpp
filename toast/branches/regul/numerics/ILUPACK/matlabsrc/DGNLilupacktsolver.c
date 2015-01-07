/* ========================================================================== */
/* === AMGsolver mexFunction ================================================ */
/* ========================================================================== */

/*
    Usage:

    Return the structure 'options' and solution x of Ax=b based on ILUPACK V2.2
    
    Example:

    % for initializing parameters
    [x, options] = DGNLilupacksolver(A,PREC,options, b,x0);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	October 10, 2006. ILUPACK V2.2.

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

    const mwSize  *dims;
    mxClassID  *classIDflags;
    mxArray    *tmp, *fout, *A_input , *b_input, *x0_input, *options_input, 
               *PRE_input, *options_output, *x_output;
    char       *pdata, *input_buf, *output_buf;
    mwSize     nnz, ndim,  buflen;
    mwIndex    *irs, *jcs;
    int        ifield, status, nfields,  ierr,i,j,k,l,m;
    size_t     mrows, ncols, sizebuf;
    double     dbuf, *A_valuesR, *convert, *sr, *pr, *sol, *rhs;
    mwIndex    *A_ja,                 /* row indices of input matrix A */
               *A_ia;                 /* column pointers of input matrix A */
    

    if (nrhs != 5)
       mexErrMsgTxt("five input arguments required.");
    else if (nlhs !=2)
       mexErrMsgTxt("Too many output arguments.");
    else if (!mxIsStruct(prhs[2]))
       mexErrMsgTxt("Third input must be a structure.");
    else if (!mxIsNumeric(prhs[0]))
       mexErrMsgTxt("First input must be a matrix.");

    /* The first input must be a square matrix.*/
    A_input = (mxArray *) prhs [0] ;
    /* get size of input matrix A */
    mrows = mxGetM (A_input) ;
    ncols = mxGetN (A_input) ;
    nnz = mxGetNzmax(A_input);
    if (mrows!=ncols) {
       mexErrMsgTxt("First input must be a square matrix.");
    }
    if (!mxIsSparse (A_input))
    {
        mexErrMsgTxt ("ILUPACK: input matrix must be in sparse format.") ;
    }





    /* copy input matrix to sparse row format */
    A.nc=A.nr=mrows;
    A.ia=(integer *) MAlloc((size_t)(A.nc+1)*sizeof(integer),"DGNLilupacksolver");
    A.ja=(integer *) MAlloc((size_t)nnz     *sizeof(integer),"DGNLilupacksolver");
    A. a=(double  *) MAlloc((size_t)nnz     *sizeof(double), "DGNLilupacksolver");

    A_ja         = (mwIndex *) mxGetIr (A_input) ;
    A_ia         = (mwIndex *) mxGetJc (A_input) ;
    A_valuesR    = (double *) mxGetPr(A_input);

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */

    /*
    for (i = 0 ; i < ncols ; i++)
      for (j = A_ia[i] ; j < A_ia[i+1] ; j++)
	printf("i=%d ja=%d  A.real=%e\n", i+1,  A_ja[j]+1, A_valuesR[j]);
    */


    for (i=0; i<=A.nr; i++)
        A.ia[i]=0;
    /* remember that MATLAB uses storage by columns and NOT by rows! */
    for (i=0; i<A.nr; i++) {
	for (j=A_ia[i]; j<A_ia[i+1]; j++) {
	    k=A_ja[j];
	    A.ia[k+1]++;
	}
    }
    /* now we know how many entries are located in every row */

    /* switch to pointer structure */
    for (i=0; i<A.nr; i++) 
        A.ia[i+1]+=A.ia[i];

    for (i=0; i<ncols; i++) {
        for (j=A_ia[i]; j<A_ia[i+1]; j++) {
	    /* row index l in C-notation */
	    l=A_ja[j];
	    /* where does row l currently start */
	    k=A.ia[l];
	    /* column index will be i in FORTRAN notation */
	    A.ja[k]=i+1;
	    A.a [k++]=A_valuesR[j];
	    A.ia[l]=k; 
	}
    }
    /* switch to FORTRAN style */
    for (i=A.nr; i>0; i--) 
        A.ia[i]=A.ia[i-1]+1;
    A.ia[0]=1;

    /*
    for (i = 0 ; i < A.nr ; i++)
      for (j = A.ia[i]-1 ; j < A.ia[i+1]-1 ; j++)
	  printf("i=%d ja=%d  A.real=%e\n", i+1,  A.ja[j], A.a[j]);
    */



    /* import pointer to the preconditioner */
    PRE_input = (mxArray*) prhs [1] ;
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


    /* rescale input matrix */
    /* obsolete
    for (i=0; i <A.nr; i++) {
	for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	    A.a[j]*=PRE->rowscal[i]*PRE->colscal[A.ja[j]-1];
	}
    }
    */



    /* Get third input argument `options' */
    options_input=(mxArray *) prhs[2];
    nfields = mxGetNumberOfFields(options_input);

    /* Allocate memory  for storing classIDflags */
    classIDflags = (mxClassID *) mxCalloc((size_t)nfields+1, (size_t)sizeof(mxClassID));
    
    /* allocate memory  for storing pointers */
    fnames = mxCalloc((size_t)nfields+1, (size_t)sizeof(*fnames));

    /* Get field name pointers */
    j=-1;
    for (ifield = 0; ifield < nfields; ifield++) {
        fnames[ifield] = mxGetFieldNameByNumber(options_input,ifield);
	/* check whether `options.niter' already exists */
	if (!strcmp("niter",fnames[ifield]))
	   j=ifield;
    }
    if (j==-1)
       fnames[nfields]="niter";
    /* mexPrintf("search for niter completed\n"); fflush(stdout); */


    /* import data */
    for (ifield = 0; ifield < nfields; ifield++) {
        /* mexPrintf("%2d\n",ifield+1); fflush(stdout); */
	tmp = mxGetFieldByNumber(options_input,0,ifield);
	classIDflags[ifield] = mxGetClassID(tmp); 

	ndim = mxGetNumberOfDimensions(tmp);
	dims = mxGetDimensions(tmp);

	/* Create string/numeric array */
	if (classIDflags[ifield] == mxCHAR_CLASS) {
	   /* Get the length of the input string. */
	   buflen = (mxGetM(tmp) * mxGetN(tmp)) + 1;

	   /* Allocate memory for input and output strings. */
	   input_buf = (char *) mxCalloc((size_t)buflen, (size_t)sizeof(char));

	   /* Copy the string data from tmp into a C string 
	      input_buf. */
	   status = mxGetString(tmp, input_buf, buflen);
	   
	   if (!strcmp("amg",fnames[ifield])) {
              if (strcmp(param->amg,input_buf)) {
		 param->amg=(char *)MAlloc((size_t)buflen*sizeof(char),"ilupacksolver");
		 strcpy(param->amg,input_buf);
	      }
	   }
	   else if (!strcmp("presmoother",fnames[ifield])) {
              if (strcmp(param->presmoother,input_buf)) {
		 param->presmoother=(char *)MAlloc((size_t)buflen*sizeof(char),
						   "ilupacksolver");
		 strcpy(param->presmoother,input_buf);
	      }
	   }
	   else if (!strcmp("postsmoother",fnames[ifield])) {
              if (strcmp(param->postsmoother,input_buf)) {
		 param->postsmoother=(char *)MAlloc((size_t)buflen*sizeof(char),
						    "ilupacksolver");
		 strcpy(param->postsmoother,input_buf);
	      }
	   }
	   else if (!strcmp("typecoarse",fnames[ifield])) {
              if (strcmp(param->typecoarse,input_buf)) {
		 param->typecoarse=(char *)MAlloc((size_t)buflen*sizeof(char),
						  "ilupacksolver");
		 strcpy(param->typecoarse,input_buf);
	      }
	   }
	   else if (!strcmp("typetv",fnames[ifield])) {
              if (strcmp(param->typetv,input_buf)) {
		 param->typetv=(char *)MAlloc((size_t)buflen*sizeof(char),
					      "ilupacksolver");
		 strcpy(param->typetv,input_buf);
	      }
	   }
	   else if (!strcmp("FCpart",fnames[ifield])) {
              if (strcmp(param->FCpart,input_buf)) {
		 param->FCpart=(char *)MAlloc((size_t)buflen*sizeof(char),
					      "ilupacksolver");
		 strcpy(param->FCpart,input_buf);
	      }
	   }
	   else if (!strcmp("solver",fnames[ifield])) {
              if (strcmp(param->solver,input_buf)) {
		 param->solver=(char *)MAlloc((size_t)buflen*sizeof(char),
					      "ilupacksolver");
		 strcpy(param->solver,input_buf);
	      }
	   }
	   else if (!strcmp("ordering",fnames[ifield])) {
              if (strcmp(param->ordering,input_buf)) {
	         param->ordering=(char *)MAlloc((size_t)buflen*sizeof(char),
						"ilupacksolver");
		 strcpy(param->ordering,input_buf);
	      }
	   }
	   else {
	      /* mexPrintf("%s ignored\n",fnames[ifield]);fflush(stdout); */
	   }
	} 
	else {
	   if (!strcmp("elbow",fnames[ifield])) {
	      param->elbow=*mxGetPr(tmp);
	   }
	   else if (!strcmp("lfilS",fnames[ifield])) {
	      param->lfilS=*mxGetPr(tmp);
	   }
	   else if (!strcmp("lfil",fnames[ifield])) {
	      param->lfil=*mxGetPr(tmp);
	   }
	   else if (!strcmp("maxit",fnames[ifield])) {
	      param->maxit=*mxGetPr(tmp);
	   }
	   else if (!strcmp("droptolS",fnames[ifield])) {
	      param->droptolS=*mxGetPr(tmp);
	   }
	   else if (!strcmp("droptol",fnames[ifield])) {
	      param->droptol=*mxGetPr(tmp);
	   }
	   else if (!strcmp("condest",fnames[ifield])) {
	      param->condest=*mxGetPr(tmp);
	   }
	   else if (!strcmp("restol",fnames[ifield])) {
	      param->restol=*mxGetPr(tmp);
	   }
	   else if (!strcmp("npresmoothing",fnames[ifield])) {
	      param->npresmoothing=*mxGetPr(tmp);
	   }
	   else if (!strcmp("npostmoothing",fnames[ifield])) {
	      param->npostsmoothing=*mxGetPr(tmp);
	   }
	   else if (!strcmp("ncoarse",fnames[ifield])) {
	      param->ncoarse=*mxGetPr(tmp);
	   }
	   else if (!strcmp("matching",fnames[ifield])) {
	      param->matching=*mxGetPr(tmp);
	   }
	   else if (!strcmp("nrestart",fnames[ifield])) {
	      param->nrestart=*mxGetPr(tmp);
	   }
	   else if (!strcmp("damping",fnames[ifield])) {
	      param->damping=*mxGetPr(tmp);
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      param->mixedprecision=*mxGetPr(tmp);
	   }
	   else {
	     /* mexPrintf("%s ignored\n",fnames[ifield]);fflush(stdout); */
	   }
	}
    }


    /* copy right hand side `b' */
    b_input = (mxArray *) prhs [3] ;
    /* get size of input matrix A */
    rhs=(double*) MAlloc((size_t)A.nr*sizeof(double),"DGNLilupacksolver:rhs");
    pr=mxGetPr(b_input);
    memcpy(rhs,pr,(size_t)A.nr*sizeof(double));





    /* copy initial solution `x0' */
    x0_input = (mxArray *) prhs [4] ;
    /* numerical solution */
    sol=(double *)MAlloc((size_t)A.nr*sizeof(double),"DGNLilupacksolver:sol");
    pr=mxGetPr(x0_input);
    memcpy(sol,pr,(size_t)A.nr*sizeof(double));


    /* set bit 2, transpose */
    param->ipar[4]|=4;
    /* set bit 3, conjugate */
    param->ipar[4]|=0;
    ierr=DGNLAMGsolver(&A, PRE, param, rhs, sol);


    /* Create a struct matrices for output */
    nlhs=2;
    if (j==-1)
       plhs[1] = mxCreateStructMatrix((mwSize)1, (mwSize)1, nfields+1, fnames);
    else
       plhs[1] = mxCreateStructMatrix((mwSize)1, (mwSize)1, nfields, fnames);
    if (plhs[1]==NULL)
       mexErrMsgTxt("Could not create structure mxArray\n");
    options_output=plhs[1];

    /* export data */
    for (ifield = 0; ifield<nfields; ifield++) {
	tmp = mxGetFieldByNumber(options_input,0,ifield);
	classIDflags[ifield] = mxGetClassID(tmp); 

	ndim = mxGetNumberOfDimensions(tmp);
	dims = mxGetDimensions(tmp);

	/* Create string/numeric array */
	if (classIDflags[ifield] == mxCHAR_CLASS) {
	   if (!strcmp("amg",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->amg)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param->amg);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("presmoother",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->presmoother)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param->presmoother);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("postsmoother",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->postsmoother)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param->postsmoother);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("typecoarse",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->typecoarse)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param->typecoarse);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("typetv",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->typetv)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param->typetv);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("FCpart",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->FCpart)+1, (size_t)sizeof(char));
	      strcpy(output_buf,param->FCpart);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("solver",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->solver)+1, (size_t)sizeof(char));
	      strcpy(output_buf, param->solver);
	      fout = mxCreateString(output_buf);
	   }
	   else if (!strcmp("ordering",fnames[ifield])) {
	      output_buf = (char *) mxCalloc((size_t)strlen(param->ordering)+1, (size_t)sizeof(char));
	      strcpy(output_buf, param->ordering);
	      fout = mxCreateString(output_buf);
	   }
	   else {
	      /* Get the length of the input string. */
	      buflen = (mxGetM(tmp) * mxGetN(tmp)) + 1;

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
	      dbuf=param->elbow;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("lfilS",fnames[ifield])) {
	      dbuf=param->lfilS;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("lfil",fnames[ifield])) {
	      dbuf=param->lfil;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("maxit",fnames[ifield])) {
	      dbuf=param->maxit;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("droptolS",fnames[ifield])) {
	      dbuf=param->droptolS;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("droptol",fnames[ifield])) {
	      dbuf=param->droptol;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("condest",fnames[ifield])) {
	      dbuf=param->condest;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("restol",fnames[ifield])) {
	      dbuf=param->restol;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("npresmoothing",fnames[ifield])) {
	      dbuf=param->npresmoothing;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("npostmoothing",fnames[ifield])) {
	      dbuf=param->npostsmoothing;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("ncoarse",fnames[ifield])) {
	      dbuf=param->ncoarse;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("matching",fnames[ifield])) {
	      dbuf=param->matching;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("nrestart",fnames[ifield])) {
	      dbuf=param->nrestart;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("damping",fnames[ifield])) {
	      dbuf=param->damping;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      dbuf=param->mixedprecision;
	      memcpy(pdata, &dbuf, sizebuf);
	   }
	   else {
	      memcpy(pdata, mxGetData(tmp), sizebuf);
	   }
	}


	/* Set each field in output structure */
	mxSetFieldByNumber(options_output, (mwIndex)0, ifield, fout);
    }

    /* store number of iteration steps */
    if (j==-1)
       ifield=nfields;
    else
       ifield=j;
    fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
    pr=mxGetPr(fout);
    
    *pr=param->ipar[26];
    
    /* set each field in output structure */
    mxSetFieldByNumber(options_output, (mwIndex)0, ifield, fout);      

    mxFree(fnames);
    mxFree(classIDflags);


    plhs[0] = mxCreateDoubleMatrix((mwSize)A.nr, (mwSize)1, mxREAL);
    x_output =plhs[0];
    pr=mxGetPr(x_output);
    memcpy(pr,sol,(size_t)A.nr*sizeof(double));  

    /* release right hand side */
    free(rhs);

    /* release solution */
    free(sol);

    /* release input matrix */
    free(A.ia);
    free(A.ja);
    free(A.a);
    
    

    switch (ierr) {
    case  0: /* perfect! */
      break;
    case -1: /* too many iteration steps */
      mexPrintf("!!! ILUPACK Warning !!!\n");
      mexPrintf("number of iteration steps exceeded its limit.\nEither increase `options.maxit'\n or recompute ILUPACK preconditioner using a smaller `options.droptol'");
      break;
    case -2: /* weird, should not occur */
      mexErrMsgTxt("not enough workspace provided.");
      plhs[0]=NULL;
      break;
    case -3: /* breakdown */
      mexErrMsgTxt("iterative solver breaks down.\nMost likely you need to recompute ILUPACK preconditioner using a smaller `options.droptol'");
      plhs[0]=NULL;
      break;
    default: /* zero pivot encountered at step number ierr */
      mexPrintf("iterative solver exited with error code %d",ierr);
      mexErrMsgTxt(".");
      plhs[0]=NULL;
      break;
    } /* end switch */

    return;
}

