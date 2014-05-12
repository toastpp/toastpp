/* ========================================================================== */
/* === AMGfactor mexFunction ================================================ */
/* ========================================================================== */

/*
    Usage:

    Return the structure 'options' and preconditioner 'PREC' for ILUPACK V2.3
    
    Example:

    % for initializing parameters
    [PREC, options, rcomflag,S,tv] = ZHPDilupackfactor(A,options,PRE,tv);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        January 23, 2009. ILUPACK V2.3.

    Acknowledgements:

	This work was supported from 2002 to 2007 by the DFG research center
        MATHEON "Mathematics for key technologies"

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
    Zmat A;
    ZAMGlevelmat *PRE, *current;
    CAMGlevelmat *SPRE, *scurrent;
    ZILUPACKparam *param;
    integer n, nnzU;
    int  tv_exists, tv_field;

    const char **fnames;
    const char *pnames[]= {"n","nB", "L","D","U", "E","F", "rowscal","colscal", "p","invq",
			   "param","ptr", "isreal","isdefinite","issymmetric","ishermitian", "issingle", "A_H", "errorL", "errorU", "errorS"};


    const mwSize  *dims;
    mxClassID  *classIDflags; 
    mxArray    *tmp, *fout, *PRE_input, *tv_input, *A_input, *options_input, 
               *PRE_output, *options_output, *S_output, *tv_output;
    char       *pdata, *input_buf, *output_buf;
    mwSize     ndim,nnz,buflen; 
    int        ifield, status, nfields, ierr,i,j,k,l,m;
    integer    *ibuff;
    mwIndex    jstruct, *irs, *jcs;
    const mwSize  mydims[]={1,1};
    size_t     sizebuf, mrows, ncols;
    double     dbuf, *A_valuesR,  *A_valuesI, *convert, *sr, *si, *pr, *pi, det, detr,deti;
    doublecomplex *p;
    complex       *sp;
    mwIndex *A_ja,                 /* row indices of input matrix A */
            *A_ia;                 /* column pointers of input matrix A */
    doublecomplex *mytv;
    

    if (nrhs!=2 && nrhs!=4)
       mexErrMsgTxt("Two/four input arguments required.");
    else if (nlhs!=2 && nlhs!=5 )
       mexErrMsgTxt("wrong number of output arguments.");
    else if (!mxIsStruct(prhs[1]))
       mexErrMsgTxt("Second input must be a structure.");
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
    A.ia=(integer *)       MAlloc((size_t)(A.nc+1)*sizeof(integer),"ZHPDilupackfactor");
    A.ja=(integer *)       MAlloc((size_t)nnz     *sizeof(integer),"ZHPDilupackfactor");
    A.a =(doublecomplex *) MAlloc((size_t)nnz     *sizeof(doublecomplex),"ZHPDilupackfactor");


    A_ja         = (mwIndex *) mxGetIr (A_input) ;
    A_ia         = (mwIndex *) mxGetJc (A_input) ;
    A_valuesR    = (double *)  mxGetPr(A_input);
    A_valuesI    = (double *)  mxGetPi(A_input);

    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */

    /*    
    for (i = 0 ; i < ncols ; i++)
      for (j = A_ia[i] ; j < A_ia[i+1] ; j++)
	printf("i=%d ja=%d  A.real=%e A.imag=%e\n", i+1,  A_ja[j]+1, A_valuesR[j], A_valuesI[j]);
    */

    A.ia[0]=1;
    for (i = 0 ; i < ncols ; i++) {
        A.ia[i+1]=A.ia[i];
	for (j = A_ia[i] ; j < A_ia[i+1] ; j++) {
	    k=A_ja[j];
	    if (k>=i) {
	       l=A.ia[i+1]-1;
	       A.ja[l]=k+1;
	       A.a [l].r= A_valuesR[j];
	       A.a [l].i=-A_valuesI[j]; /* conjugate since MATLAB uses storage by columns, but ILUPACK
					   requires storage by rows. For Hermitian matrices we simply
					   have to switch the sign of the imaginary part */
	       A.ia[i+1]=l+2;
	    }
	}
    }

    /*
    for (i = 0 ; i < A.nr ; i++)
      for (j = A.ia[i]-1 ; j < A.ia[i+1]-1 ; j++)
	  printf("i=%d ja=%d  A.real=%e  A.imag=%e\n", i+1,  A.ja[j], A.a[j].r, A.a[j].i);
    */

    param=(ZILUPACKparam *)MAlloc((size_t)sizeof(ZILUPACKparam),"ZHPDilupackfactor:param");
    ZHPDAMGinit(&A,param);

    /* Get second input arguments */
    options_input=(mxArray *)prhs[1];
    nfields = mxGetNumberOfFields(options_input);

    /* Allocate memory  for storing classIDflags */
    classIDflags = (mxClassID *) mxCalloc((size_t)nfields, (size_t)sizeof(mxClassID));
    
    /* allocate memory  for storing pointers */
    fnames = mxCalloc((size_t)nfields, (size_t)sizeof(*fnames));

    /* Get field name pointers */
    for (ifield = 0; ifield < nfields; ifield++) {
        fnames[ifield] = mxGetFieldNameByNumber(options_input,ifield);
    }

    /* import data */
    tv_exists=0;
    tv_field=-1;
    for (ifield = 0; ifield < nfields; ifield++) {
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
		 param->amg=(char *)MAlloc((size_t)buflen*sizeof(char),
					   "ilupackfactor");
		 strcpy(param->amg,input_buf);
	      }
	   }
	   else if (!strcmp("presmoother",fnames[ifield])) {
              if (strcmp(param->presmoother,input_buf)) {
		 param->presmoother=(char *)MAlloc((size_t)buflen*sizeof(char),
						   "ilupackfactor");
		 strcpy(param->presmoother,input_buf);
	      }
	   }
	   else if (!strcmp("postsmoother",fnames[ifield])) {
              if (strcmp(param->postsmoother,input_buf)) {
		 param->postsmoother=(char *)MAlloc((size_t)buflen*sizeof(char),
						    "ilupackfactor");
		 strcpy(param->postsmoother,input_buf);
	      }
	   }
	   else if (!strcmp("typecoarse",fnames[ifield])) {
              if (strcmp(param->typecoarse,input_buf)) {
		 param->typecoarse=(char *)MAlloc((size_t)buflen*sizeof(char),
						  "ilupackfactor");
		 strcpy(param->typecoarse,input_buf);
	      }
	   }
	   else if (!strcmp("typetv",fnames[ifield])) {
              if (strcmp(param->typetv,input_buf)) {
		 param->typetv=(char *)MAlloc((size_t)buflen*sizeof(char),
					      "ilupackfactor");
		 strcpy(param->typetv,input_buf);
	      }
	      /* 'static' or 'dynamic' test vector */
	      if (strcmp("none",input_buf)) tv_exists=-1;
	   }
	   else if (!strcmp("FCpart",fnames[ifield])) {
              if (strcmp(param->FCpart,input_buf)) {
		 param->FCpart=(char *)MAlloc((size_t)buflen*sizeof(char),
					      "ilupackfactor");
		 strcpy(param->FCpart,input_buf);
	      }
	   }
	   else if (!strcmp("solver",fnames[ifield])) {
              if (strcmp(param->solver,input_buf)) {
		 param->solver=(char *)MAlloc((size_t)buflen*sizeof(char),
					      "ilupackfactor");
		 strcpy(param->solver,input_buf);
	      }
	   }
	   else if (!strcmp("ordering",fnames[ifield])) {
              if (strcmp(param->ordering,input_buf)) {
	         param->ordering=(char *)MAlloc((size_t)buflen*sizeof(char),
						"ilupackfactor");
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
	      param->damping.r=*mxGetPr(tmp);
	      if (mxIsComplex(tmp))
		 param->damping.i=*mxGetPi(tmp);
	      else
		 param->damping.i=0;
	   }
	   else if (!strcmp("tv",fnames[ifield])) {
	      tv_field=ifield;
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      param->mixedprecision=*mxGetPr(tmp);
	   }
	   else {
	     /* mexPrintf("%s ignored\n",fnames[ifield]);fflush(stdout); */
	   }
	}
    }
    if (param->droptolS>0.125*param->droptol) {
       mexPrintf("!!! ILUPACK Warning !!!\n");
       mexPrintf("`param.droptolS' is recommended to be one order of magnitude less than `param.droptol'\n");
    } 

    if (tv_exists && tv_field>=0 && 
        (nrhs==2   || (nrhs==4 && mxIsNumeric(prhs[2])))) {
       tmp=mxGetFieldByNumber(options_input,0,tv_field);
       pr=mxGetPr(tmp);
       mytv=(doublecomplex *)MAlloc((size_t)A.nr*sizeof(doublecomplex),"ZHPDilupackfactor");
       param->tv=mytv;
       if (!mxIsComplex(tmp)) {
	  for (i=0; i<A.nr; i++) {
	      param->tv[i].r=pr[i];
	      param->tv[i].i=0;
	  }
       }
       else {
	  pi=mxGetPi(tmp);
	  for (i=0; i<A.nr; i++) {
	      param->tv[i].r=pr[i];
	      param->tv[i].i=pi[i];
	  }
       } 
    }
    mxFree(fnames);


    /*mexPrintf("factor the matrix\n");fflush(stdout);*/
    PRE=(ZAMGlevelmat *)MAlloc((size_t)sizeof(ZAMGlevelmat),"ZHPDilupackfactor");

#ifdef PRINT_INFO
    mexPrintf("ZHPDilupackfactor: factorize matrix\n");fflush(stdout);
#endif
    if (nrhs==4) {
       
       /* at least the second time that we call this routine */
       if (!mxIsNumeric(prhs[2])) {
	  /* import pointer to the preconditioner */
	  PRE_input = (mxArray*) prhs [2] ;
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


	  /* import pointer to the test vector */
	  tv_input = (mxArray*) prhs [3] ;
	  /* size of the test vector */
	  j=mxGetM(tv_input);
	  pr=mxGetPr(tv_input);

	  if (!mxIsComplex(tv_input)) {
	     for (i=0; i<j; i++) {
	         param->tv[i].r=pr[i];
		 param->tv[i].i=0;
	     }
	  }
	  else {
	     pi=mxGetPi(tv_input);
	     for (i=0; i<j; i++) {
	         param->tv[i].r=pr[i];
		 param->tv[i].i=pi[i];
	     }
	  }

	  /* remap A if necessary */
	  if (A.nr==param->mstack[0].nr) {
	    param->mstack[0].ia=A.ia;
	    param->mstack[0].ja=A.ja;
	    param->mstack[0].a =A.a;
	    /* mexPrintf("original matrix remapped\n");fflush(stdout); */
	  }
       }
    }

    ierr=ZHPDAMGfactor(&A, PRE, param);
    if (tv_exists && tv_field>=0 && 
        (nrhs==2   || (nrhs==4 && mxIsNumeric(prhs[2])))) 
       free(mytv);
    /*mexPrintf("matrix factored\n");fflush(stdout);*/

#ifdef PRINT_INFO
    mexPrintf("ZHPDilupackfactor: matrix factored\n");fflush(stdout);
    if (param->rcomflag!=0) {
       mexPrintf("ZHPDilupackfactor: request for reverse communication\n");fflush(stdout);
    }
#endif

    if (nlhs==5 && ierr==0) {
       plhs[2]=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
       pr = mxGetPr(plhs[2]);
       *pr=param->rcomflag;
	
       /* mexPrintf("rcomflag set\n");fflush(stdout); */

       if (param->rcomflag!=0) {
	  nnz=param->A.ia[param->A.nr]-1;
	  /* mexPrintf("extract sparse matrix %ld,%ld,%ld\n",param->A.nr,param->A.nr,nnz);fflush(stdout); */
	  plhs[3]=mxCreateSparse((mwSize)param->A.nr,(mwSize)param->A.nc, (mwSize)nnz, mxCOMPLEX);
	  S_output=plhs[3];

	  sr  = (double *)  mxGetPr(S_output);
	  si  = (double *)  mxGetPi(S_output);
	  irs = (mwIndex *) mxGetIr(S_output);
	  jcs = (mwIndex *) mxGetJc(S_output);

	  k=0;
	  for (i=0; i<param->A.nr; i++) {
	      jcs[i]=k;
	      /* strict lower triangular part */
	      for (j=param->A.ia[i]-1; j<param->A.ia[i+1]-1; j++) {
		  irs[k] =param->A.ja[j]-1;
		  sr[k]  =param->A.a[j].r;
		  si[k++]=param->A.a[j].i;
	      }
	  } /* end for i */
	  jcs[i]=k;

	  /* mexPrintf("extract test vector\n");fflush(stdout);*/
	  /* output test vector */
	  plhs[4]=mxCreateDoubleMatrix((mwSize)A.nr,(mwSize)1, mxCOMPLEX);
	  tv_output=plhs[4];
	  pr = mxGetPr(tv_output);
	  pi = mxGetPi(tv_output);
	  for (i=0; i<A.nr; i++) {
	      *pr++=param->tv[i].r;
	      *pi++=param->tv[i].i;
	  }
       }
       else {
	  plhs[3]=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	  plhs[4]=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
       }
    }


    if (param->ipar[16]&DISCARD_MATRIX) {
       free(A.ia);
       free(A.ja);
       free(A.a);
    }
    else {
       if (PRE->issingle) {
	  SPRE=(CAMGlevelmat *)PRE;
	  SPRE->A.ia=A.ia;
	  SPRE->A.ja=A.ja;
	  SPRE->A.a =(complex *)A.a;
       }
       else {
	  PRE->A.ia=A.ia;
	  PRE->A.ja=A.ja;
	  PRE->A.a =A.a;
       }
    }
    
    if (ierr) {
       nlhs=0;
       /* finally release memory of the preconditioner */
       free(PRE);
       free(param);
    }

    switch (ierr) {
    case  0: /* perfect! */
      break;
    case -1: /* Error. input matrix may be wrong.
		(The elimination process has generated a
		row in L or U whose length is .gt.  n.) */
      mexErrMsgTxt("ILUPACK error, data may be wrong.");
      break;
    case -2: /* The matrix L overflows the array alu */
      mexErrMsgTxt("memory overflow, please increase `options.elbow' and retry");
      break;
    case -3: /* The matrix U overflows the array alu */
      mexErrMsgTxt("memory overflow, please increase `options.elbow' and retry");
      break;
    case -4: /* Illegal value for lfil */
      mexErrMsgTxt("Illegal value for `options.lfil'\n");
      break;
    case -5: /* zero row encountered */
      mexErrMsgTxt("zero row encountered, please reduce `options.droptol'\n");
      break;
    case -6: /* zero column encountered */
      mexErrMsgTxt("zero column encountered, please reduce `options.droptol'\n");
      break;
    case -7: /* buffers too small */
      mexErrMsgTxt("memory overflow, please increase `options.elbow' and retry");
      break;
    default: /* zero pivot encountered at step number ierr */
      mexErrMsgTxt("zero pivot encountered, please reduce `options.droptol'\n");
      break;
    } /* end switch */
    if (ierr) {
      plhs[0]=NULL;
      return;
    }


    /* mexPrintf("export parameters\n");fflush(stdout);*/
    /* prepare a struct matrices for output */
    nfields = mxGetNumberOfFields(options_input);

    /* allocate memory  for storing pointers */
    fnames = mxCalloc((size_t)nfields, (size_t)sizeof(*fnames));
    /* Get field name pointers */
    for (ifield=0; ifield<nfields; ifield++) {
        fnames[ifield] = mxGetFieldNameByNumber(options_input,ifield);
    }

    plhs[1] = mxCreateStructMatrix((mwSize)1, (mwSize)1, nfields, fnames);
    if (plhs[1]==NULL)
       mexErrMsgTxt("Could not create structure mxArray");
    options_output=plhs[1];

    /* export data */
    for (ifield = 0; ifield<nfields; ifield++) {
      /*mexPrintf("%2d\n",ifield+1);fflush(stdout);*/
	tmp = mxGetFieldByNumber(options_input,0,ifield);
	classIDflags[ifield] = mxGetClassID(tmp); 

	ndim = mxGetNumberOfDimensions(tmp);
	dims = mxGetDimensions(tmp);

	/* Create string/numeric array */
	if (classIDflags[ifield] == mxCHAR_CLASS) {
	   if (!strcmp("amg",fnames[ifield])) {
	      output_buf = (char *)mxCalloc((size_t)strlen(param->amg)+1, (size_t)sizeof(char));
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
	      
	      sizebuf = buflen*sizeof(char);
	      memcpy(output_buf, input_buf, (size_t)sizebuf);
	      fout = mxCreateString(output_buf);
	   }
	} 
	else {
	   /* real case */
	   if (mxGetPi(tmp)==NULL && strcmp("damping",fnames[ifield]))
	      fout = mxCreateNumericArray((mwSize)ndim, dims, 
					  classIDflags[ifield], mxREAL);
	   else { /* complex case */
	      fout = mxCreateNumericArray((mwSize)ndim, dims, 
					  classIDflags[ifield], mxCOMPLEX);
	   }
	   pdata = mxGetData(fout);

	   sizebuf = mxGetElementSize(tmp);
	   if (!strcmp("elbow",fnames[ifield])) {
	      dbuf=param->elbow;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("lfilS",fnames[ifield])) {
	      dbuf=param->lfilS;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("lfil",fnames[ifield])) {
	      dbuf=param->lfil;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("maxit",fnames[ifield])) {
	      dbuf=param->maxit;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("droptolS",fnames[ifield])) {
	      dbuf=param->droptolS;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("droptol",fnames[ifield])) {
	      dbuf=param->droptol;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("condest",fnames[ifield])) {
	      dbuf=param->condest;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("restol",fnames[ifield])) {
	      dbuf=param->restol;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("npresmoothing",fnames[ifield])) {
	      dbuf=param->npresmoothing;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("npostmoothing",fnames[ifield])) {
	      dbuf=param->npostsmoothing;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("ncoarse",fnames[ifield])) {
	      dbuf=param->ncoarse;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("matching",fnames[ifield])) {
	      dbuf=param->matching;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("nrestart",fnames[ifield])) {
	      dbuf=param->nrestart;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("damping",fnames[ifield])) {
	      pr=mxGetPr(fout);
	      pi=mxGetPi(fout);
	      *pr=param->damping.r;
	      *pi=param->damping.i;
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      dbuf=param->mixedprecision;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else {
	      memcpy(pdata, mxGetData(tmp), (size_t)sizebuf);
	   }
	}


	/* Set each field in output structure */
	mxSetFieldByNumber(options_output, 0, ifield, fout);
    }
    /*mexPrintf("parameters exported\n");fflush(stdout);*/
    mxFree(fnames);
    mxFree(classIDflags);

    /*mexPrintf("export preconditioner\n");fflush(stdout);*/
    plhs[0] = mxCreateStructMatrix((mwSize)1, (mwSize)PRE->nlev, 22, pnames);
    if (plhs[0]==NULL)
      mexErrMsgTxt("Could not create structure mxArray\n");
    PRE_output=plhs[0];



    current=PRE;
    if (PRE->issingle) {
       SPRE=(CAMGlevelmat *)PRE;
       scurrent=SPRE;
    }
    n=A.nr;
    ibuff  =(integer *)MAlloc((size_t)n*sizeof(integer),"ZHPDilupackfactor:ibuff");
    convert=(double *) MAlloc((size_t)n*sizeof(double), "ZHPDilupackfactor:convert");
    for (jstruct = 0; jstruct < PRE->nlev; jstruct++) {
      
        /*  1. save level size to field `n' */
        ifield=0;

	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);
	
	if (PRE->issingle)
	   *pr=scurrent->n;
	else
	   *pr=current->n;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwSize)jstruct, ifield, fout);



	/*  2. save leading block size to field `nB' */
	++ifield;
	/*mexPrintf("%2d,%2d\n",jstruct+1,ifield);fflush(stdout);*/
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	if (PRE->issingle)
	   *pr=scurrent->nB;
	else
	   *pr=current->nB;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/*  3. field `L' */
	++ifield;
	if (param->rcomflag==0) {
	   /* switched to full-matrix processing */
	   if (jstruct==PRE->nlev-1 && 
	       ((PRE->issingle)?scurrent->LU.ja:current->LU.ja)==NULL) {

	      if (PRE->issingle) {
		 fout=mxCreateDoubleMatrix((mwSize)scurrent->nB,(mwSize)scurrent->nB, mxCOMPLEX);
		 for (i=0; i<scurrent->nB; i++) 
		     ibuff[i]=i*scurrent->nB-(i*(i-1))/2;
	      }
	      else { /* !PRE->issingle */
		 fout=mxCreateDoubleMatrix((mwSize)current->nB,(mwSize)current->nB, mxCOMPLEX);
		 for (i=0; i<current->nB; i++) 
		     ibuff[i]=i*current->nB-(i*(i-1))/2;
	      } /* end if-else PRE->issingle */
	   
	      sr=mxGetPr(fout);
	      si=mxGetPi(fout);
	      if (PRE->issingle)
		 sp=scurrent->LU.a;
	      else
		 p =current->LU.a;
	      for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
#ifndef USE_LAPACK_DRIVER
		  l=((PRE->issingle)?scurrent->LU.ia[i]:current->LU.ia[i])-1;
#else
		  l=i;
#endif
		  for (j=0; j<i; j++) {
#ifndef USE_LAPACK_DRIVER
		      k=((PRE->issingle)?scurrent->LU.ia[j]:current->LU.ia[j])-1;
#else
		      k=j;
#endif
		      m=k-l;
		      if (m>0)
			 m=ibuff[l]+m;
		      else
			 m=ibuff[k]-m;
		      if (PRE->issingle) {
			 *sr++=sp[m].r;
			 *si++=sp[m].i;
		      }
		      else {
			 *sr++=p[m].r;
			 *si++=p[m].i;
		      }
		  }
#ifndef USE_LAPACK_DRIVER
		  if (PRE->issingle) {
		     det=1.0/(sp[ibuff[l]].r*sp[ibuff[l]].r+sp[ibuff[l]].i*sp[ibuff[l]].i);
		     *sr++= det*sp[ibuff[l]].r;
		     *si++=-det*sp[ibuff[l]].i; /* should be zero */
		  }
		  else {
		     det=1.0/(p[ibuff[l]].r*p[ibuff[l]].r+p[ibuff[l]].i*p[ibuff[l]].i);
		     *sr++= det*p[ibuff[l]].r;
		     *si++=-det*p[ibuff[l]].i; /* should be zero */
		  }
#else
		  if (PRE->issingle) {
		     *sr++=sp[ibuff[l]].r;
		     *si++=sp[ibuff[l]].i; /* should be zero */
		  }
		  else {
		     *sr++=p[ibuff[l]].r;
		     *si++=p[ibuff[l]].i; /* should be zero */
		  }
#endif
		  for (j=i+1; j<((PRE->issingle)?scurrent->nB:current->nB); j++) {
		      *sr++=0;
		      *si++=0;
		  }
	      }	   

	      /* set each field in output structure */
	      mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	   }
	   else {
	      if (PRE->issingle)
	         nnz=scurrent->LU.ja[scurrent->nB-1]-scurrent->LU.ja[0]+scurrent->nB;
	      else
		 nnz=current->LU.ja[current->nB-1]-current->LU.ja[0]+current->nB;
	      if (param->flags&COARSE_REDUCE) {
		 if (PRE->issingle)
		    fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)nnz, mxCOMPLEX);
		 else
		    fout=mxCreateSparse((mwSize)current->nB, (mwSize)current->nB,  (mwSize)nnz, mxCOMPLEX);
	      }
	      else {
		 if (PRE->issingle)
		    fout=mxCreateSparse((mwSize)scurrent->n, (mwSize)scurrent->nB, (mwSize)nnz, mxCOMPLEX);
		 else
		    fout=mxCreateSparse((mwSize)current->n, (mwSize)current->nB, (mwSize)nnz, mxCOMPLEX);
	      }
	      
	      sr =(double *) mxGetPr(fout);
	      si =(double *) mxGetPi(fout);
	      irs=(mwIndex *)mxGetIr(fout);
	      jcs=(mwIndex *)mxGetJc(fout);

	      if (PRE->issingle) {
		 nnzU=scurrent->LU.ja[scurrent->nB];
		 scurrent->LU.ja[scurrent->nB]=scurrent->LU.nnz+1;
	      }
	      else {
		 nnzU=current->LU.ja[current->nB];
		 current->LU.ja[current->nB]=current->LU.nnz+1;
	      }

	      k=0;
	      if (PRE->issingle) {
		 for (i=0; i<scurrent->nB; i++) {
		     /* extract diagonal entry */
		     jcs[i]=k;
		     irs[k]=i;
		     det=1.0/(scurrent->LU.a[i].r*scurrent->LU.a[i].r+scurrent->LU.a[i].i*scurrent->LU.a[i].i);
		     sr[k]  = det*scurrent->LU.a[i].r;
		     si[k++]=-det*scurrent->LU.a[i].i; /* should be zero */
		     
		     for (j=scurrent->LU.ja[i]-1; j<scurrent->LU.ja[i+1]-1; j++) {
		         irs[k] = scurrent->LU.ja[j]-1;
			 sr[k]  = scurrent->LU.a[j].r;
			 si[k++]=-scurrent->LU.a[j].i; /* conjugate U! */
		     }
		 }
		 scurrent->LU.ja[scurrent->nB]=nnzU;
	      }
	      else { /* !PRE->issingle */
		 for (i=0; i<current->nB; i++) {
		     /* extract diagonal entry */
		     jcs[i]=k;
		     irs[k]=i;
		     det=1.0/(current->LU.a[i].r*current->LU.a[i].r+current->LU.a[i].i*current->LU.a[i].i);
		     sr[k]  = det*current->LU.a[i].r;
		     si[k++]=-det*current->LU.a[i].i; /* should be zero */
		     
		     for (j=current->LU.ja[i]-1; j<current->LU.ja[i+1]-1; j++) {
		         irs[k] = current->LU.ja[j]-1;
			 sr[k]  = current->LU.a[j].r;
			 si[k++]=-current->LU.a[j].i; /* conjugate U! */
		     }
		 }
		 current->LU.ja[current->nB]=nnzU;
	      } /* end if-else PRE->issingle */
	      jcs[i]=k;
	      
	      /* set each field in output structure */
	      mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	   }
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}



	/*  4. field `D' */
	++ifield;
	if (param->rcomflag==0) {
	   if (PRE->issingle)
	      fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)scurrent->nB, mxCOMPLEX);
	   else
	      fout=mxCreateSparse((mwSize)current->nB, (mwSize)current->nB,  (mwSize)current->nB,  mxCOMPLEX);

	   sr =(double *) mxGetPr(fout);
	   si =(double *) mxGetPi(fout);
	   irs=(mwIndex *)mxGetIr(fout);
	   jcs=(mwIndex *)mxGetJc(fout);
	   
	   for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
	       jcs[i]=i;
	       irs[i]=i;
	   }
	   jcs[i]=i;
	   
	   /* switched to full-matrix processing */
	   if (jstruct==PRE->nlev-1 && 
	       ((PRE->issingle)?scurrent->LU.ja:current->LU.ja)==NULL) {
	     
  	      if (PRE->issingle)
	         sp=scurrent->LU.a;
	      else
	         p=current->LU.a;
	      for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
#ifndef USE_LAPACK_DRIVER
		  if (PRE->issingle) {
		     l=scurrent->LU.ia[i]-1;
		     det=1.0/(sp[ibuff[l]].r*sp[ibuff[l]].r+sp[ibuff[l]].i*sp[ibuff[l]].i);
		     sr[i]= det*sp[ibuff[l]].r;
		     si[i]=-det*sp[ibuff[l]].i; /* should be zero */
		  }
		  else {
		     l=current->LU.ia[i]-1;
		     det=1.0/(p[ibuff[l]].r*p[ibuff[l]].r+p[ibuff[l]].i*p[ibuff[l]].i);
		     sr[i]= det*p[ibuff[l]].r;
		     si[i]=-det*p[ibuff[l]].i; /* should be zero */
		  }
#else
		  l=i;
		  if (PRE->issingle) {
		     sr[i]=sp[ibuff[l]].r;
		     si[i]=sp[ibuff[l]].i; /* should be zero */
		  }
		  else {
		     sr[i]=p[ibuff[l]].r;
		     si[i]=p[ibuff[l]].i; /* should be zero */
		  }
#endif
	      }	   
	   }
	   else {
	      for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
#ifndef USE_LAPACK_DRIVER
		  if (PRE->issingle) {
		     det=1.0/(scurrent->LU.a[i].r*scurrent->LU.a[i].r+scurrent->LU.a[i].i*scurrent->LU.a[i].i);
		     sr[i]= det*scurrent->LU.a[i].r;
		     si[i]=-det*scurrent->LU.a[i].i; /* should be zero */
		  }
		  else { /* !PRE->issingle */
		     det=1.0/(current->LU.a[i].r*current->LU.a[i].r+current->LU.a[i].i*current->LU.a[i].i);
		     sr[i]= det*current->LU.a[i].r;
		     si[i]=-det*current->LU.a[i].i; /* should be zero */
		  } /* end if-else PRE->issingle */
#else
		  if (PRE->issingle) {
		     sr[i]=scurrent->LU.a[i].r;
		     si[i]=scurrent->LU.a[i].i; /* should be zero */
		  }
		  else { /* !PRE->issingle */
		     sr[i]=current->LU.a[i].r;
		     si[i]=current->LU.a[i].i; /* should be zero */
		  } /* end if-else PRE->issingle */
#endif
	      }
	   }
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	   
	   
	/*  5. upper triangular matrix `U' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/*  6. field `E' E=F^T */
	++ifield;
	if (param->rcomflag==0) {
	   if (jstruct<PRE->nlev-1){ 
	
	      if (param->flags&COARSE_REDUCE) {
		 if (PRE->issingle) {
		    nnz=scurrent->F.ia[scurrent->nB]-1;
		    fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)scurrent->nB, (mwSize)nnz, mxCOMPLEX);
		 }
		 else { /* !PRE->issingle */
		    nnz=current->F.ia[current->nB]-1;
		    fout=mxCreateSparse((mwSize)n-current->nB,(mwSize)current->nB, (mwSize)nnz, mxCOMPLEX);
		 } /* end if-else PRE->issingle */

		 sr =(double *) mxGetPr(fout);
		 si =(double *) mxGetPi(fout);
		 irs=(mwIndex *)mxGetIr(fout);
		 jcs=(mwIndex *)mxGetJc(fout);
		 
		 k=0;
		 if (PRE->issingle) {
		    for (i=0; i<scurrent->nB; i++) {
		        jcs[i]=k;
			for (j=scurrent->F.ia[i]-1; j<scurrent->F.ia[i+1]-1; j++) {
		            irs[k] = scurrent->F.ja[j]-1;
			    sr[k]  = scurrent->F.a[j].r;
			    si[k++]=-scurrent->F.a[j].i; /* E=F' */
			}
		    }
		 }
		 else { /* !PRE->issingle */
		    for (i=0; i<current->nB; i++) {
		        jcs[i]=k;
			for (j=current->F.ia[i]-1; j<current->F.ia[i+1]-1; j++) {
		            irs[k] = current->F.ja[j]-1;
			    sr[k]  = current->F.a[j].r;
			    si[k++]=-current->F.a[j].i; /* E=F' */
			}
		    }
		 } /* end if-else PRE->issingle */
		 jcs[i]=k;
	      }
	      else {
		 nnz=0;
		 if (PRE->issingle)
		    fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)scurrent->nB, (mwSize)nnz, mxCOMPLEX);
		 else
		    fout=mxCreateSparse((mwSize)n-current->nB, (mwSize)current->nB,  (mwSize)nnz, mxCOMPLEX);
	      }

	      /* set each field in output structure */
	      mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	   }
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	
	
	/*  7. field `F' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/*  8. field `rowscal' */
	++ifield;
	
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxCOMPLEX);
	sr=mxGetPr(fout);
	si=mxGetPi(fout);
	
	if (param->rcomflag==0) {
	   if (PRE->issingle) {
	      for (i=0; i<n; i++) {
		  sr[i]=scurrent->rowscal[i].r;
		  si[i]=scurrent->rowscal[i].i;
	      }
	   }
	   else { /* !PRE->issingle */
	      for (i=0; i<n; i++) {
		  sr[i]=current->rowscal[i].r;
		  si[i]=current->rowscal[i].i;
	      }
	   } /* end if-else PRE->issingle */
	}
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	
	
	
	/*  9. field `colscal' */
	++ifield;
	
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxCOMPLEX);
	sr=mxGetPr(fout);
	si=mxGetPi(fout);
	
	if (param->rcomflag==0) {
	   if (PRE->issingle) {
	      for (i=0; i<n; i++) {
		  sr[i]=scurrent->colscal[i].r;
		  si[i]=scurrent->colscal[i].i;
	      }
	   }
	   else { /* !PRE->issingle */
	      for (i=0; i<n; i++) {
		  sr[i]=current->colscal[i].r;
		  si[i]=current->colscal[i].i;
	      }
	   } /* end if-else PRE->issingle */
	}
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	
	
	
	/* 10. field `p' */
	++ifield;
	if (param->rcomflag==0) {
	   fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);
	   pdata = mxGetData(fout);
	
	   for (i=0;i<n; i++)
	       convert[i]=((PRE->issingle)?scurrent->p[i]:current->p[i]);
	   memcpy(pdata, convert, (size_t)n*sizeof(double));
	
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	
	
	/* 11. field `invq' */
	++ifield;
	if (param->rcomflag==0) {	
	   fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);
	   pdata = mxGetData(fout);
	
	   for (i=0;i<n; i++)
	       convert[i]=((PRE->issingle)?scurrent->invq[i]:current->invq[i]);
	   memcpy(pdata, convert, (size_t)n*sizeof(double));
	   
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}


	/* 12. field `param' */
	++ifield;	
	fout=mxCreateNumericArray((mwSize)1,mydims, mxUINT64_CLASS, mxREAL);
	pdata = mxGetData(fout);
	
	memcpy(pdata, &param, (size_t)sizeof(size_t));
	
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* 13. field `ptr' */
	++ifield;
	
	fout=mxCreateNumericArray((mwSize)1,mydims, mxUINT64_CLASS, mxREAL);
	pdata = mxGetData(fout);
	
	memcpy(pdata, &PRE, (size_t)sizeof(size_t));
	
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* 14. save non-real property to field `isreal' */
	++ifield;

	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=0;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* 15. save positive definite property to field `isdefinite' */
	++ifield;

	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=1;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* 16. save non-symmetry property to field `issymmetric' */
	++ifield;

	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=0;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* 17. save Hermitian property to field `ishermitian' */
	++ifield;

	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=1;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	


	/* 18. save trivial Hermitian property to field `issingle' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	if (PRE->issingle)
	   *pr=scurrent->issingle;
	else
	   *pr=current->issingle;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);


	
	/* 19. save coarse grid system `A_H' */
	++ifield;
	if (jstruct>=PRE->nlev-1) {
	   fout=mxCreateSparse((mwSize)0,(mwSize)0,(mwSize)0, mxCOMPLEX);
	}
	else if (param->ipar[16]&DISCARD_MATRIX) {
 	   if (PRE->issingle)
	      fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)n-scurrent->nB, (mwSize)0, mxCOMPLEX);
	   else
	      fout=mxCreateSparse((mwSize)n-current->nB, (mwSize)n-current->nB,  (mwSize)0, mxCOMPLEX);
	}
	else {
	   /* switched to full-matrix processing */
	   if (jstruct==PRE->nlev-2 && 
	       ((PRE->issingle)?scurrent->next->LU.ja:current->next->LU.ja)==NULL) 
	      fout=mxCreateSparse((mwSize)0,(mwSize)0,(mwSize)0, mxCOMPLEX);
	   else {
	      if (PRE->issingle) {
		 nnz=scurrent->next->A.ia[scurrent->next->A.nr]-1;
		 fout=mxCreateSparse((mwSize)scurrent->next-A.nr, (mwSize)scurrent->next->A.nc, (mwSize)nnz, mxCOMPLEX);
	      }
	      else { /* !PRE->issingle */
		 nnz=current->next->A.ia[current->next->A.nr]-1;

		 /*
		   mexPrintf("level %d, %d,%d,%d\n",jstruct+1,current->next->A.nr,current->next->A.nc,nnz);fflush(stdout);
		   fout=mxCreateSparse(0,0,0, mxCOMPLEX);
		 */
	      
		 fout=mxCreateSparse((mwSize)current->next-A.nr, (mwSize)current->next->A.nc, (mwSize)nnz, mxCOMPLEX);
	      } /* end if-else PRE->issingle */
  
	      sr =(double *) mxGetPr(fout);
	      si =(double *) mxGetPi(fout);
	      irs=(mwIndex *)mxGetIr(fout);
	      jcs=(mwIndex *)mxGetJc(fout);
	      
	      k=0;
	      if (PRE->issingle) {
		 for (i=0; i<scurrent->next->A.nr; i++) {
		     jcs[i]=k;
		     for (j=scurrent->next->A.ia[i]-1; j<scurrent->next->A.ia[i+1]-1; j++) {
		         irs[k] =scurrent->next->A.ja[j]-1;
			 sr[k]  =scurrent->next->A.a[j].r;
			 si[k++]=scurrent->next->A.a[j].i;
		     }
		 }
	      }
	      else { /* !PRE->issingle */
		 for (i=0; i<current->next->A.nr; i++) {
		     jcs[i]=k;
		     for (j=current->next->A.ia[i]-1; j<current->next->A.ia[i+1]-1; j++) {
		         irs[k] =current->next->A.ja[j]-1;
			 sr[k]  =current->next->A.a[j].r;
			 si[k++]=current->next->A.a[j].i;
		     }
		 }
	      } /* end if-else PRE->issingle */
	      jcs[i]=k;
	   }
	}
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);

	
	
	/* 20. save error in L to field `errorL' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);
	if (PRE->issingle)
	   *pr=scurrent->errorL;     
	else
	   *pr=current->errorL;     
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);

	/* 21. save error in U to field `errorU' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);
	if (PRE->issingle)
	   *pr=scurrent->errorU;
	else
	   *pr=current->errorU;
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);

	/* 22. save error in S to field `errorS' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);
	if (PRE->issingle)
	   *pr=scurrent->errorS;
	else
	   *pr=current->errorS;
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);


	if (PRE->issingle) {
	   n-=scurrent->nB;
	   scurrent=scurrent->next;
	}
	else {
	   n-=current->nB;
	   current=current->next;
	}
    }
    /*mexPrintf("preconditioner exported\n");fflush(stdout);*/
    free(ibuff);
    free(convert);
    /*mexPrintf("memory released\n");fflush(stdout);*/
    
    return;
}

