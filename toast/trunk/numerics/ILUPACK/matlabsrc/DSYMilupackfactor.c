/* ========================================================================== */
/* === AMGfactor mexFunction ================================================ */
/* ========================================================================== */

/*
    Usage:

    Return the structure 'options' and preconditioner 'PREC' for ILUPACK V2.3
    
    Example:

    % for initializing parameters
    [PREC, options, rcomflag,S,tv] = DSYMilupackfactor(A,options,PRE,tv);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	March 17, 2009. ILUPACK V2.3.  

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
#include <lapack.h>

#define MAX_FIELDS 100
#define MAX(A,B) (((A)>=(B))?(A):(B))
/* #define PRINT_INFO */

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
    Dmat A, M;
    DAMGlevelmat *PRE, *current;
    SAMGlevelmat *SPRE, *scurrent;
    DILUPACKparam *param;
    integer n, nnzU, len, k_old, *stack, *Ibuff, ii,blocksize;
    int tv_exists, tv_field, ind_exists, ind_field, ind_shiftmatrix=-1;

    const char **fnames;
    const char *pnames[]= {"n","nB", "L","D","U", "E","F", "rowscal","colscal", "p","invq",
			   "param","ptr", "isreal","isdefinite","issymmetric","ishermitian", "issingle", "A_H", "errorL", "errorU", "errorS", "isblock"};


    const mwSize  *dims;
    const mwSize  mydims[]={1,1};
    mxClassID  *classIDflags;
    mxArray    *tmp, *fout, *PRE_input, *tv_input, *A_input, *options_input,
               *PRE_output, *options_output, *S_output, *tv_output;
    char       *pdata, *input_buf, *output_buf;
    mwSize     nnz, ndim, buflen;
    mwIndex    jstruct, *irs, *jcs;
    int        ifield, status, nfields, ierr,i,j,k,l,m;
    size_t     mrows, ncols, sizebuf;
    double     dbuf, *A_valuesR, *convert, *sr, *pr, det, *Dbuff;
    mwIndex    *A_ja,                 /* row indices of input matrix A */
               *A_ia;                 /* column pointers of input matrix A */
    

    if (nrhs!=2 && nrhs!=4)
       mexErrMsgTxt("Two/four input arguments required.");
    else if (nlhs!=2 && nlhs!=5)
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
    A.ia=(integer *) MAlloc((size_t)(A.nc+1)*sizeof(integer),"DSYMilupackfactor");
    A.ja=(integer *) MAlloc((size_t)nnz     *sizeof(integer),"DSYMilupackfactor");
    A. a=(double  *) MAlloc((size_t)nnz     *sizeof(double), "DSYMilupackfactor");
    stack=(integer *) MAlloc((size_t)A.nc   *sizeof(integer), "DSYMilupackfactor");


    A_ja         = (mwIndex *)mxGetIr(A_input);
    A_ia         = (mwIndex *)mxGetJc(A_input);
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

    A.ia[0]=1;
    for (i = 0 ; i < ncols ; i++) {
        A.ia[i+1]=A.ia[i];
	for (j = A_ia[i] ; j < A_ia[i+1] ; j++) {
	    k=A_ja[j];
	    if (k>=i) {
	       l=A.ia[i+1]-1;
	       A.ja[l]=k+1;
	       A.a [l]=A_valuesR[j];
	       A.ia[i+1]=l+2;
	    }
	}
    }

    /*
    for (i = 0 ; i < A.nr ; i++)
      for (j = A.ia[i]-1 ; j < A.ia[i+1]-1 ; j++)
	  printf("i=%d ja=%d  A.real=%e\n", i+1,  A.ja[j], A.a[j]);
    */

    param=(DILUPACKparam *)MAlloc((size_t)sizeof(DILUPACKparam),"DSYMilupackfactor:param");
    DSYMAMGinit(&A,param);


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
    ind_exists=0;
    ind_field=-1;
    for (ifield = 0; ifield < nfields; ifield++) {
	tmp = mxGetFieldByNumber(options_input,0,ifield);
	classIDflags[ifield] = mxGetClassID(tmp); 

	ndim = mxGetNumberOfDimensions(tmp);
	dims = mxGetDimensions(tmp);

	/* Create string/numeric array */
	if (classIDflags[ifield] == mxCHAR_CLASS) {
	   /* Get the length of the input string. */
	   buflen=(mxGetM(tmp)*mxGetN(tmp))+1;

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
	      param->damping=*mxGetPr(tmp);
	   }
	   else if (!strcmp("tv",fnames[ifield])) {
	      tv_field=ifield;
	   }
	   else if (!strcmp("ind",fnames[ifield])) {
	      ind_exists=-1;
	      ind_field=ifield;
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      param->mixedprecision=*mxGetPr(tmp);
	   }
	   else if (!strcmp("shift0",fnames[ifield])) {
	      param->shift0=*mxGetPr(tmp);
	      /* mexPrintf("shift0=%8.1le\n",param->shift0);*/
	      fflush(stdout);
	   }
	   else if (!strcmp("shiftmax",fnames[ifield])) {
	      param->shiftmax=*mxGetPr(tmp);
	      /* mexPrintf("shiftmax=%8.1le\n",param->shiftmax);*/
	      fflush(stdout);
	   }
	   else if (!strcmp("shifts",fnames[ifield])) {
	      mrows=mxGetM(tmp);
	      ncols=mxGetN(tmp);
	      param->nshifts=MAX(mrows,ncols);
	      pr=mxGetPr(tmp);
	      /* mexPrintf("shifts (%d,%d):\n",mrows,ncols); */
	      param->shifts=(doubleprecision *)MAlloc(param->nshifts*sizeof(doubleprecision),
						      "DSYMilupackfactor");
	      for (i=0; i<param->nshifts; i++) {
		  param->shifts[i]=pr[i];
		  /* mexPrintf("%8.1le",param->shifts[i]); */
	      }
	      /* mexPrintf("\n");
		 fflush(stdout);
	      */
	   }
	   else if (!strcmp("shiftmatrix",fnames[ifield])) {
	      ind_shiftmatrix=ifield;
	      mrows=mxGetM(tmp);
	      ncols=mxGetN(tmp);
	      nnz  =mxGetNzmax(tmp);
	      if (mrows!=ncols) {
		 mexErrMsgTxt("Shift matrix must be a square matrix.");
	      }
	      if (!mxIsSparse (tmp)) {
		  mexErrMsgTxt ("ILUPACK: shift matrix must be in sparse format.") ;
	      }

	      /* copy input matrix to sparse row format */
	      M.nc=M.nr=mrows;
	      M.ia=(integer *) MAlloc((size_t)(M.nc+1)*sizeof(integer),"DSYMilupackfactor");
	      M.ja=(integer *) MAlloc((size_t)nnz     *sizeof(integer),"DSYMilupackfactor");
	      M. a=(double  *) MAlloc((size_t)nnz     *sizeof(double), "DSYMilupackfactor");

	      A_ja     =(mwIndex *)mxGetIr(tmp);
	      A_ia     =(mwIndex *)mxGetJc(tmp);
	      A_valuesR=(double *) mxGetPr(tmp);

	      /* -------------------------------------------------------------------- */
	      /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	      /*     notation.                                                        */
	      /* -------------------------------------------------------------------- */
    
	      /*
		for (i = 0 ; i < ncols ; i++)
		for (j = A_ia[i] ; j < A_ia[i+1] ; j++)
		printf("i=%d ja=%d  A.real=%e\n", i+1,  A_ja[j]+1, A_valuesR[j]);
	      */

	      M.ia[0]=1;
	      for (i=0; i<ncols; i++) {
		  M.ia[i+1]=M.ia[i];
		  for (j=A_ia[i]; j<A_ia[i+1]; j++) {
		      k=A_ja[j];
		      if (k>=i) {
			 l=M.ia[i+1]-1;
			 M.ja[l]=k+1;
			 M.a [l]=A_valuesR[j];
			 M.ia[i+1]=l+2;
		      }
		  }

		  /*
		  mexPrintf("%8d\n",i+1);
		  for (j=M.ia[i]; j<M.ia[i+1]; j++) {
		      mexPrintf("%8d",M.ja[j-1]);
		  }
		  mexPrintf("\n");
		  for (j=M.ia[i]; j<M.ia[i+1]; j++) {
		      mexPrintf("%8.1le",M.a[j-1]);
		  }
		  mexPrintf("\n");
		  */
	      }
	      /* fflush(stdout); */
	      param->shiftmatrix=&M;
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
       param->tv=mxGetPr(tmp);
    }
    if (ind_exists && ind_field>=0) {
       tmp=mxGetFieldByNumber(options_input,0,ind_field);
       param->ind=(integer *)MAlloc(A.nc*sizeof(integer),"DSYMilupackfactor");
       pr=mxGetPr(tmp);
       for (i=0;i<A.nc; i++)
	   param->ind[i]=pr[i];
    }
    mxFree(fnames);


    PRE =(DAMGlevelmat *)MAlloc((size_t)sizeof(DAMGlevelmat),"DSYMilupackfactor");

#ifdef PRINT_INFO
    mexPrintf("DSYMilupackfactor: factorize matrix\n");fflush(stdout);
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
	  tv_input=(mxArray*)prhs[3] ;
	  /* size of the test vector */
	  j=mxGetM(tv_input);
	  pr=mxGetPr(tv_input);
	  memcpy(param->tv,pr,j*sizeof(double));

	  /* remap A if necessary */
	  if (A.nr==param->mstack[0].nr) {
	    param->mstack[0].ia=A.ia;
	    param->mstack[0].ja=A.ja;
	    param->mstack[0].a =A.a;
	    /* mexPrintf("original matrix remapped\n");fflush(stdout); */
	  }
       }
    }

#ifdef PRINT_INFO
    mexPrintf("DSYMilupackfactor: factor matrix \n");fflush(stdout);
#endif

    ierr=DSYMAMGfactor(&A, PRE, param);
#ifdef PRINT_INFO
    mexPrintf("DSYMilupackfactor: matrix factored\n");fflush(stdout);
    if (param->rcomflag!=0) {
       mexPrintf("DSYMilupackfactor: request for reverse communication\n");fflush(stdout);
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
	  plhs[3]=mxCreateSparse((mwSize)param->A.nr,(mwSize)param->A.nc, (mwSize)nnz, mxREAL);
	  S_output=plhs[3];

	  sr  = (double *) mxGetPr(S_output);
	  irs = (mwIndex *)mxGetIr(S_output);
	  jcs = (mwIndex *)mxGetJc(S_output);

	  k=0;
	  for (i=0; i<param->A.nr; i++) {
	      jcs[i]=k;
	      /* strict lower triangular part */
	      for (j=param->A.ia[i]-1; j<param->A.ia[i+1]-1; j++) {
		  irs[k] =param->A.ja[j]-1;
		  sr[k++]=param->A.a[j];
	      }
	  } /* end for i */
	  jcs[i]=k;

	  /* mexPrintf("extract test vector\n");fflush(stdout);*/
	  /* output test vector */
	  plhs[4]=mxCreateDoubleMatrix((mwSize)A.nr,(mwSize)1, mxREAL);
	  tv_output=plhs[4];
	  pr = mxGetPr(tv_output);
	  for (i=0; i<A.nr; i++) 
	      *pr++=param->tv[i];
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
       free(param->ind);
    }
    else {
       if (PRE->issingle) {
	  SPRE=(SAMGlevelmat *)PRE;
	  SPRE->A.ia=A.ia;
	  SPRE->A.ja=A.ja;
	  SPRE->A.a =(real *)A.a;
       }
       else {
	  PRE->A.ia=A.ia;
	  PRE->A.ja=A.ja;
	  PRE->A.a=A.a;
       }
    }
    if (param->shiftmatrix!=NULL) {
       free(M.ia);
       free(M.ja);
       free(M.a);
    }
    if (param->shifts!=NULL) {
       free(param->shifts);
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
    for (ifield=0; ifield<nfields; ifield++) {
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
	      buflen=(mxGetM(tmp)*mxGetN(tmp))+1;

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
	else if (!strcmp("shiftmatrix",fnames[ifield])) {
	   mrows=mxGetM(tmp);
	   ncols=mxGetN(tmp);
	   nnz  =mxGetNzmax(tmp);

	   A_ja     =(mwIndex *)mxGetIr(tmp);
	   A_ia     =(mwIndex *)mxGetJc(tmp);
	   A_valuesR=(double *) mxGetPr(tmp);

	   fout=mxCreateSparse((mwSize)mrows,(mwSize)ncols, (mwSize)nnz, mxREAL);

	   irs=(mwIndex *)mxGetIr(fout);
	   jcs=(mwIndex *)mxGetJc(fout);
	   sr =(double *) mxGetPr(fout);
   
	   for (i=0; i<=ncols; i++)
	       jcs[i]=A_ia[i];
	   for (i=0; i<nnz; i++)
	       irs[i]=A_ja[i];
	   for (i=0; i<nnz; i++)
	       sr[i]=A_valuesR[i];
	}
	else {
	   fout = mxCreateNumericArray((mwSize)ndim, dims, 
				       classIDflags[ifield], mxREAL);
	   pdata = mxGetData(fout);

	   sizebuf = mxGetElementSize(tmp);
	   if (!strcmp("elbow",fnames[ifield])) {
	      dbuf=param->elbow;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("lfilS",fnames[ifield])) {
	      dbuf=param->lfilS;
	      memcpy(pdata, &dbuf, sizebuf);
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
	      dbuf=param->damping;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("mixedprecision",fnames[ifield])) {
	      dbuf=param->mixedprecision;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("shift0",fnames[ifield])) {
	      dbuf=param->shift0;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("shiftmax",fnames[ifield])) {
	      dbuf=param->shiftmax;
	      memcpy(pdata, &dbuf, (size_t)sizebuf);
	   }
	   else if (!strcmp("shifts",fnames[ifield])) {
	      mrows=mxGetM(tmp);
	      ncols=mxGetN(tmp);
	      pr=mxGetPr(tmp);
	      sr=mxGetPr(fout);
	      for (i=0; i<mrows*ncols; i++)
		  sr[i]=pr[i];
	   }
	   else {
	      memcpy(pdata, mxGetData(tmp), (size_t)sizebuf);
	   }
	}


	/* Set each field in output structure */
	mxSetFieldByNumber(options_output, (mwIndex)0, ifield, fout);
    }
    mxFree(fnames);
    mxFree(classIDflags);



    /* create preconditioner for output */
    plhs[0] = mxCreateStructMatrix((mwSize)1, (mwSize)PRE->nlev, 23, pnames);
    if (plhs[0]==NULL)
       mexErrMsgTxt("Could not create structure mxArray\n");
    PRE_output=plhs[0];


    current=PRE;
    if (PRE->issingle) {
       SPRE=(SAMGlevelmat *)PRE;
       scurrent=SPRE;
    }
    n=A.nr;
    convert=(double*)MAlloc((size_t)n*sizeof(double),"DSYMilupackfactor:convert");
    for (jstruct = 0; jstruct<PRE->nlev; jstruct++) {

        /*  1. save level size to field `n' */
        ifield=0;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr=mxGetPr(fout);
	
	if (PRE->issingle) 
	   *pr=scurrent->n;
	else
	   *pr=current->n;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/*  2. save leading block size to field `nB' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr=mxGetPr(fout);

	if (PRE->issingle) 
	   *pr=scurrent->nB;
	else
	   *pr=current->nB;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);


	
	/*  3. lower triangular matrix field `L=U^T' */
	++ifield;
	if (param->rcomflag==0) {
	   /* switched to full-matrix / PARDISO processing */
	   if (jstruct==PRE->nlev-1 && 
	       ((PRE->issingle)?scurrent->LU.ja:current->LU.ja)==NULL) {

	      /* PARDISO case */
	      if (((PRE->issingle)?(void *)scurrent->LU.a:(void *)current->LU.a)==NULL) {
		 
		 if (PRE->issingle)
		    fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)scurrent->nB, mxREAL);
		 else
		    fout=mxCreateSparse((mwSize)current->nB,(mwSize)current->nB, (mwSize)current->nB, mxREAL);

		 sr =(double *) mxGetPr(fout);
		 irs=(mwIndex *)mxGetIr(fout);
		 jcs=(mwIndex *)mxGetJc(fout);
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
		     sr[i]=1.0;
		     jcs[i]=i;
		     irs[i]=i;
		 }
		 jcs[i]=i;
	      }
	      else {
		 if (PRE->issingle)
		    fout=mxCreateDoubleMatrix((mwSize)scurrent->nB,(mwSize)scurrent->nB, mxREAL);
		 else
		    fout=mxCreateDoubleMatrix((mwSize)current->nB,(mwSize)current->nB, mxREAL);

		 /* complete data using the return values from LAPACK's
		    Bunch-Kaufman pivoting */
		 pr  = mxGetPr(fout);
		 for (i=0; i<((PRE->issingle)?scurrent->nB*scurrent->nB:current->nB*current->nB); i++) {
		     pr[i]=0.0;
		 }
		 i=0;
		 j=0;
		 k=0;
		 while (i<((PRE->issingle)?scurrent->nB:current->nB)) {
		       if (((PRE->issingle)?scurrent->LU.ia[i]:current->LU.ia[i])>0) {
			  j+=i;
			  pr[j++]=1.0;
			  l=((PRE->issingle)?scurrent->nB:current->nB)-i-1;
			  k++;
			  while (l--) {
			        pr[j++]=((PRE->issingle)?scurrent->LU.a[k++]:current->LU.a[k++]);
			  }
			  i++;
		       }
		       else {
			  j+=i;
			  pr[j++]=1.0;
			  pr[j++]=0.0;
			  l=((PRE->issingle)?scurrent->nB:current->nB)-i-2;
			  k+=2;
			  while (l--) {
			        pr[j++]=((PRE->issingle)?scurrent->LU.a[k++]:current->LU.a[k++]);
			  }
			  j+=i+1;
			  pr[j++]=1.0;
			  l=((PRE->issingle)?scurrent->nB:current->nB)-i-2;
			  k++;
			  while (l--) {
			        pr[j++]=((PRE->issingle)?scurrent->LU.a[k++]:current->LU.a[k++]);
			  }
			  i+=2;
		       }
		 }
		 /* reorder L-factor */
		 i=0;
		 while (i<((PRE->issingle)?scurrent->nB:current->nB)) {
		       k=((PRE->issingle)?scurrent->LU.ia[i]:current->LU.ia[i]);
		       if (k>0) {
			  if (k!=i+1) {
			     if (PRE->issingle) {
			        for (j=0; j<i; j++) {
				    dbuf=pr[j*scurrent->nB+i];
				    pr[j*scurrent->nB+i]  =pr[j*scurrent->nB+k-1];
				    pr[j*scurrent->nB+k-1]=dbuf;
				}
			     }
			     else {
			        for (j=0; j<i; j++) {
				    dbuf=pr[j*current->nB+i];
				    pr[j*current->nB+i]  =pr[j*current->nB+k-1];
				    pr[j*current->nB+k-1]=dbuf;
				}
			     }
			  }
			  i++;
		       }
		       else {
			  if (PRE->issingle) {
			     for (j=0; j<i; j++) {
			         dbuf=pr[j*scurrent->nB+i+1];
				 pr[j*scurrent->nB+i+1]=pr[j*scurrent->nB-k-1];
				 pr[j*scurrent->nB-k-1]=dbuf;
			     }
			  }
			  else {
			     for (j=0; j<i; j++) {
			         dbuf=pr[j*current->nB+i+1];
				 pr[j*current->nB+i+1]=pr[j*current->nB-k-1];
				 pr[j*current->nB-k-1]=dbuf;
			     }
			  }
			  i+=2;
		       } /* end else */
		 } /* end while */
	      } /* end if-else (current->LU.a=NULL) */
	      /* set each field in output structure */
	      mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	   }
	   else { /* sparse case */
	      if (PRE->isblock) {
		 if (PRE->issingle)
		    nnz=scurrent->LU.nnz;
		 else
		    nnz=current->LU.nnz;
		 if (param->flags&COARSE_REDUCE) {
		    if (PRE->issingle)
		       fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)nnz, mxREAL);
		    else
		       fout=mxCreateSparse((mwSize)current->nB, (mwSize)current->nB,  (mwSize)nnz, mxREAL);
		 }
		 else {
		    if (PRE->issingle)
		       fout=mxCreateSparse((mwSize)scurrent->n, (mwSize)scurrent->nB, (mwSize)nnz-scurrent->nB, mxREAL);
		    else
		       fout=mxCreateSparse((mwSize)current->n,  (mwSize)current->nB,  (mwSize)nnz-current->nB, mxREAL);
		 }
		 /* mexPrintf("nnz=%d\n",nnz-current->nB); */

		 sr  = (double *) mxGetPr(fout);
		 irs = (mwIndex *)mxGetIr(fout);
		 jcs = (mwIndex *)mxGetJc(fout);
		 
		 if (PRE->issingle) {
		    nnzU=scurrent->LU.ja[scurrent->nB];
		    scurrent->LU.ja[scurrent->nB]=scurrent->LU.nnz+1;
		 }
		 else {
		    nnzU=current->LU.ja[current->nB];
		    current->LU.ja[current->nB]=current->LU.nnz+1;
		 }

		 m=(PRE->issingle)?SPRE->maxblocksize:PRE->maxblocksize;
		 Dbuff=(double *)MAlloc((size_t)2*m*m*sizeof(double),"DSYMilupackfactor");
		 Ibuff=(integer *)MAlloc((size_t)m*sizeof(integer),"DSYMilupackfactor");
		 
		 m=0;
		 /* mexPrintf("nB=%d\n",((PRE->issingle)?scurrent->nB:current->nB)); */
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); ) {
		     blocksize=((PRE->issingle)?scurrent->blocksize[i]:current->blocksize[i]);
		     ii=((PRE->issingle)?
			 scurrent->LU.ja[i]-1:current->LU.ja[i]-1);
		     /* copy block */
		     for (k=0; k<blocksize; k++) {
		         for (l=0; l<blocksize; l++) {
			     Dbuff[l+k*blocksize]=(PRE->issingle)?
			       scurrent->LU.a[ii]:current->LU.a[ii];
			     ii++;
			 }
		     }
		     /* factorize diagonal block */
		     k=blocksize*blocksize;
		     dsytrf("u",&blocksize,Dbuff,&blocksize,
			    Ibuff,Dbuff+k,&k,&l,1);
		     if (l<0) {
		        mexPrintf("_SYTRF: argument %d has an illegal value",-l);
			return ;
		     }
		     if (l>0) {
		        mexPrintf("_SYTRF: matrix is singular");
			return ;
		     }
		     /* invert diagonal block */
		     k=blocksize*blocksize;
		     dsytri("u",&blocksize,Dbuff,&blocksize,
			    Ibuff,Dbuff+k,&l,1);
		     if (l<0) {
		        mexPrintf("_SYTRI: argument %d has an illegal value",-l);
			return ;
		     }
		     /* complete upper triangular part */
		     for (k=0; k<blocksize; k++) {
		         for (l=k+1; l<blocksize; l++) {
			     Dbuff[l+k*blocksize]=Dbuff[k+l*blocksize];
			 }
		     }
		     /*
		     for (k=0; k<blocksize; k++) {
		         for (l=0; l<blocksize; l++) {
			   mexPrintf("%16.9le",Dbuff[l+k*blocksize]);
			 }
			 mexPrintf("\n");
		     }
		     */

		     /* transfer column i+l */
		     for (l=0; l<blocksize; l++) {
		         jcs[i+l]=m;
			 /* mexPrintf("jcs[%d]=%d\n",i+l,m);fflush(stdout); */
			 /* transfer column i+l of the diagonal block */
		         for (k=0; k<blocksize; k++) {
			     /* row i+k */
			     irs[m]=i+k;
			     sr[m++]=Dbuff[l+k*blocksize];
			     /* mexPrintf("%d,%d,%8.1le\n",irs[m-1],i+l+1,sr[m-1]); */
			 }
			 k=(PRE->issingle)?
			   scurrent->LU.ja[i]-1:current->LU.ja[i]-1;
			 ii=k+l+blocksize*blocksize;
			 k+=blocksize;
			 /* mexPrintf("%d,%d\n",k,(PRE->issingle)? 
				  scurrent->LU.ja[i+1]-1:
				  current->LU.ja[i+1]-1); */
		         for (; k<((PRE->issingle)?
				  scurrent->LU.ja[i+1]-1:
				   current->LU.ja[i+1]-1); k++) {
			     /* row i+k */
			     irs[m]=(PRE->issingle)?
			            scurrent->LU.ja[k]-1:
			            current->LU.ja[k]-1;
			     sr[m++]=(PRE->issingle)?
			            scurrent->LU.a[ii]:
			            current->LU.a[ii];
			     /* mexPrintf("%d,%d,%d,%8.1le\n",k,irs[m-1],i+l+1,sr[m-1]); */
			     ii+=blocksize;
			 }
		     }

		     jcs[i+blocksize]=m;
		     /* mexPrintf("jcs[%d]=%d\n",i+blocksize,m);fflush(stdout); */
		     i+=blocksize;
		 } /* end for i */
		 free(Dbuff);
		 free(Ibuff);

		 if (PRE->issingle) 
		    scurrent->LU.ja[scurrent->nB]=nnzU;
		 else
		    current->LU.ja[current->nB]=nnzU;

	      } /* end if PRE->isblock */
	      else { /* scalar version */
		 if (PRE->issingle)
		    nnz=scurrent->LU.nnz+1-scurrent->LU.ja[0]+3*scurrent->nB-2;
		 else
		    nnz=current->LU.nnz+1-current->LU.ja[0]+3*current->nB-2;

		 if (param->flags&COARSE_REDUCE) {
		    if (PRE->issingle)
		       fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)nnz, mxREAL);
		    else
		       fout=mxCreateSparse((mwSize)current->nB, (mwSize)current->nB,  (mwSize)nnz, mxREAL);
		 }
		 else {
		    if (PRE->issingle)
		       fout=mxCreateSparse((mwSize)scurrent->n, (mwSize)scurrent->nB, (mwSize)nnz, mxREAL);
		    else
		       fout=mxCreateSparse((mwSize)current->n,  (mwSize)current->nB,  (mwSize)nnz, mxREAL);
		 }

		 sr  = (double *) mxGetPr(fout);
		 irs = (mwIndex *)mxGetIr(fout);
		 jcs = (mwIndex *)mxGetJc(fout);
		 
		 if (PRE->issingle) {
		    nnzU=scurrent->LU.ja[scurrent->nB];
		    scurrent->LU.ja[scurrent->nB]=scurrent->LU.nnz+1;
		 }
		 else {
		    nnzU=current->LU.ja[current->nB];
		    current->LU.ja[current->nB]=current->LU.nnz+1;
		 }
		 k=0;
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); ) {
		     /* extract (block) diagonal entry */
		     if (((PRE->issingle)?scurrent->LU.ja[scurrent->nB+1+i]:current->LU.ja[current->nB+1+i])==0) {
		        /* column i */
		        jcs[i]=k;

			/* diagonal entry */
			irs[k]=i;
			sr[k++]=1.0/((PRE->issingle)?scurrent->LU.a[i]:current->LU.a[i]);
			
			/* strict lower triangular part */
			k_old=k;
			if (PRE->issingle) {
			   for (j=scurrent->LU.ja[i]-1; j<scurrent->LU.ja[i+1]-1; j++) {
			       irs[k] =scurrent->LU.ja[j]-1;
			       sr[k++]=scurrent->LU.a[j];
			   }
			}
			else {
			   for (j=current->LU.ja[i]-1; j<current->LU.ja[i+1]-1; j++) {
			       irs[k] =current->LU.ja[j]-1;
			       sr[k++]=current->LU.a[j];
			   }
			}
			len=k-k_old;
			/* sort entries in increasing order */
			Dqsort(sr+k_old,(integer *)irs+k_old,stack,&len);
			i++;
		     }
		     else {
		        if (PRE->issingle)
			   det=1.0/(scurrent->LU.a[i]*scurrent->LU.a[i+1]
				   -scurrent->LU.a[scurrent->nB+1+i]*scurrent->LU.a[scurrent->nB+1+i]);
			else
			   det=1.0/(current->LU.a[i]*current->LU.a[i+1]
				   -current->LU.a[current->nB+1+i]*current->LU.a[current->nB+1+i]);
			
			/* column i */
			jcs[i]=k;
			
			/* diagonal entry */
			irs[k]=i;
			if (PRE->issingle)
			   sr[k++]=det*scurrent->LU.a[i+1];
			else
			   sr[k++]=det*current->LU.a[i+1];
			
			/* sub-diagonal entry */
			irs[k]=i+1;
			if (PRE->issingle)
			   sr[k++]=-det*scurrent->LU.a[scurrent->nB+1+i];
			else
			   sr[k++]=-det*current->LU.a[current->nB+1+i];

			/* strict lower triangular part */
			k_old=k;
			if (PRE->issingle) {
		           l=scurrent->LU.ja[i]-1;
			   for (j=l; j<scurrent->LU.ja[i+1]-1; j++) {
			       irs[k] =scurrent->LU.ja[j]-1;
			       sr[k++]=scurrent->LU.a[l];
			       l+=2;
			   }
			}
			else {
			   l=current->LU.ja[i]-1;
			   for (j=l; j<current->LU.ja[i+1]-1; j++) {
			       irs[k] =current->LU.ja[j]-1;
			       sr[k++]=current->LU.a[l];
			       l+=2;
			   }
			}
			len=k-k_old;
			/* sort entries in increasing order */
			Dqsort(sr+k_old,(integer *)irs+k_old,stack,&len);
		     
			
			/* column i+1 */
			jcs[i+1]=k;
			
			/* super-diagonal entry */
			irs[k]=i;
			if (PRE->issingle) 
			   sr[k++]=-det*scurrent->LU.a[scurrent->nB+1+i];
			else
			   sr[k++]=-det*current->LU.a[current->nB+1+i];
			
			/* diagonal entry */
			irs[k]=i+1;
			if (PRE->issingle) 
			   sr[k++]=det*scurrent->LU.a[i];
			else
			   sr[k++]=det*current->LU.a[i];
			
			k_old=k;
			if (PRE->issingle) {
			   l=scurrent->LU.ja[i]-1;
			   /* strict lower triangular part */
			   for (j=l; j<scurrent->LU.ja[i+1]-1; j++) {
			       irs[k] =scurrent->LU.ja[j]-1;
			       sr[k++]=scurrent->LU.a[l+1];
			       l+=2;
			   }
			}
			else {
			   l=current->LU.ja[i]-1;
			   /* strict lower triangular part */
			   for (j=l; j<current->LU.ja[i+1]-1; j++) {
			       irs[k] =current->LU.ja[j]-1;
			       sr[k++]=current->LU.a[l+1];
			       l+=2;
			   }
			}
			len=k-k_old;
			/* sort entries in increasing order */
			Dqsort(sr+k_old,(integer *)irs+k_old,stack,&len);	     
			i+=2;
		     }
		 } /* end for i */
		 jcs[i]=k;
		 if (PRE->issingle) 
		    scurrent->LU.ja[scurrent->nB]=nnzU;
		 else
		    current->LU.ja[current->nB]=nnzU;
	      } /* end if-else isblock */
	
	      /* set each field in output structure */
	      mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	   } /* end if-else dense/sparse case */
	} /* end if rcomflag==0 */
	else {
	  fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	  mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}


	/* mexPrintf("4. field `D'\n");fflush(stdout); */
	/*  4. diagonal matrix field `D' */
	++ifield;
	if (param->rcomflag==0) {
	   if (PRE->issingle) {
	      if (scurrent->isblock)
		 fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)SPRE->maxblocksize*scurrent->nB, mxREAL);
	      else
		 fout=mxCreateSparse((mwSize)scurrent->nB,(mwSize)scurrent->nB, (mwSize)2*scurrent->nB, mxREAL);
	   }
	   else {
	      if (current->isblock) {
		 fout=mxCreateSparse((mwSize)current->nB, (mwSize)current->nB,  (mwSize)PRE->maxblocksize*current->nB, mxREAL);
		 /* mexPrintf("nnz=%d\n",PRE->maxblocksize*current->nB); */
	      }
	      else
	         fout=mxCreateSparse((mwSize)current->nB, (mwSize)current->nB,  (mwSize)2*current->nB, mxREAL);
	   }

	   sr =(double *) mxGetPr(fout);
	   irs=(mwIndex *)mxGetIr(fout);
	   jcs=(mwIndex *)mxGetJc(fout);
   
	   /* switched to full-matrix / PARDISO processing */
	   if (jstruct==PRE->nlev-1 && ((PRE->issingle)?scurrent->LU.ja:current->LU.ja)==NULL) {
  	      /* PARDISO case */
	      if (((PRE->issingle)?(void *)scurrent->LU.a:(void *)current->LU.a)==NULL) {
		 /* nothing */
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
		     sr[i]=1.0;
		     jcs[i]=i;
		     irs[i]=i;
		 }
		 jcs[i]=i;
	      }
	      else {/* complete data using the return values from LAPACK's
		       Bunch-Kaufman pivoting */
		k=0;
		j=0;
		for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); ) {
	            /* extract (block) diagonal entry */
		    if (((PRE->issingle)?scurrent->LU.ia[i]:current->LU.ia[i])>0) {
		       /* column i */
		       jcs[i]=k;

		       /* diagonal entry */
		       irs[k]=i;
		       sr[k++]=((PRE->issingle)?scurrent->LU.a[j]:current->LU.a[j]);
		       
		       j+=((PRE->issingle)?scurrent->nB:current->nB)-i;
		       i++;
		    }
		    else {

		       /* column i */
		       jcs[i]=k;

		       /* diagonal entry */
		       irs[k]=i;
		       sr[k++]=((PRE->issingle)?scurrent->LU.a[j]:current->LU.a[j]);
		       
		       /* sub-diagonal entry */
		       irs[k]=i+1;
		       sr[k++]=((PRE->issingle)?scurrent->LU.a[j+1]:current->LU.a[j+1]);
		       
		       /* column i+1 */
		       jcs[i+1]=k;
		       
		       /* super-diagonal entry */
		       irs[k]=i;
		       sr[k++]=((PRE->issingle)?scurrent->LU.a[j+1]:current->LU.a[j+1]);
		       
		       /* diagonal entry */
		       irs[k]=i+1;
		       if (PRE->issingle) {
			  sr[k++]=scurrent->LU.a[j+scurrent->nB-i];
			  j+=2*(scurrent->nB-i)-1;
		       }
		       else {
			  sr[k++]=current->LU.a[j+current->nB-i];
			  j+=2*(current->nB-i)-1;
		       }
		       i+=2;
		    }
		}
		jcs[i]=k;
	      } /* if-else (current->LU.a=NULL) */
	   } /* end-if (jstruct==PRE->nlev-1 && current->LU.ja==NULL) */
	   else { /* sparse case */
	      if (PRE->isblock) {
		 m=((PRE->issingle)?SPRE->maxblocksize:PRE->maxblocksize);
		 Dbuff=(double *)MAlloc((size_t)2*m*m*sizeof(double),"DSYMilupackfactor");
		 Ibuff=(integer *)MAlloc((size_t)m*sizeof(integer),"DSYMilupackfactor");
		 
		 m=0;
		 /* mexPrintf("nB=%d\n",current->nB); */
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); ) {
		     blocksize=((PRE->issingle)?scurrent->blocksize[i]:current->blocksize[i]);
		     ii=((PRE->issingle)?
			scurrent->LU.ja[i]-1:current->LU.ja[i]-1);
		     /* copy block */
		     for (k=0; k<blocksize; k++) {
		         for (l=0; l<blocksize; l++) {
			     Dbuff[l+k*blocksize]=(PRE->issingle)?
			       scurrent->LU.a[ii]:current->LU.a[ii];
			     ii++;
			 }
		     }
		     /* factorize diagonal block */
		     k=blocksize*blocksize;
		     dsytrf("u",&blocksize,Dbuff,&blocksize,
			    Ibuff,Dbuff+k,&k,&l,1);
		     if (l<0) {
		        mexPrintf("_SYTRF: argument %d has an illegal value",-l);
			return ;
		     }
		     if (l>0) {
		        mexPrintf("_SYTRF: matrix is singular");
			return ;
		     }
		     /* invert diagonal block */
		     k=blocksize*blocksize;
		     dsytri("u",&blocksize,Dbuff,&blocksize,
			    Ibuff,Dbuff+k,&l,1);
		     if (l<0) {
		        mexPrintf("_SYTRI: argument %d has an illegal value",-l);
			return ;
		     }
		     for (k=0; k<blocksize; k++) {
		         for (l=k+1; l<blocksize; l++) {
			     Dbuff[l+k*blocksize]=Dbuff[k+l*blocksize];
			 }
		     }
		     for (k=0; k<blocksize; k++) {
		         /* column i+k */
		         jcs[i+k]=m;
			 /* mexPrintf("jcs[%d]=%d\n",i+k,m); */
		         for (l=0; l<blocksize; l++) {
			     /* row i+l */
			     irs[m]=i+l;
			     sr[m++]=Dbuff[l+k*blocksize];
			     /* mexPrintf("%d,%d,%8.1le\n",i+k,irs[m-1],sr[m-1]); */
			 }
		     }
		     jcs[i+k]=m;
		     /* mexPrintf("jcs[%d]=%d\n",i+k,m); */
		     i+=blocksize;
		 } /* end for i */
		 free(Dbuff);
		 free(Ibuff);
	      }
	      else { /* non-block version */
		 k=0;
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); ) {
		     /* extract (block) diagonal entry */
		     if (((PRE->issingle)?scurrent->LU.ja[scurrent->nB+1+i]:current->LU.ja[current->nB+1+i])==0) {
		        /* column i */
		        jcs[i]=k;

			/* diagonal entry */
			irs[k]=i;
			sr[k++]=1.0/((PRE->issingle)?scurrent->LU.a[i]:current->LU.a[i]);
			
			i++;
		     }
		     else { /* 2x2 block */
		        if (PRE->issingle)
			   det=1.0/(scurrent->LU.a[i]*scurrent->LU.a[i+1]
				   -scurrent->LU.a[scurrent->nB+1+i]*scurrent->LU.a[scurrent->nB+1+i]);
			else
			   det=1.0/(current->LU.a[i]*current->LU.a[i+1]
				   -current->LU.a[current->nB+1+i]*current->LU.a[current->nB+1+i]);
		     
			/* column i */
			jcs[i]=k;
			
			/* diagonal entry */
			irs[k]=i;
			if (PRE->issingle)
			   sr[k++]=det*scurrent->LU.a[i+1];
			else
			   sr[k++]=det*current->LU.a[i+1];
		     
			/* sub-diagonal entry */
			irs[k]=i+1;
			if (PRE->issingle)
			   sr[k++]=-det*scurrent->LU.a[scurrent->nB+1+i];
			else
			   sr[k++]=-det*current->LU.a[current->nB+1+i];
			
			/* column i+1 */
			jcs[i+1]=k;
			
			/* super-diagonal entry */
			irs[k]=i;
			if (PRE->issingle)
			   sr[k++]=-det*scurrent->LU.a[scurrent->nB+1+i];
			else
			   sr[k++]=-det*current->LU.a[current->nB+1+i];
		     
			/* diagonal entry */
			irs[k]=i+1;
			if (PRE->issingle)
			   sr[k++]=det*scurrent->LU.a[i];
			else
			   sr[k++]=det*current->LU.a[i];
			
			i+=2;
		     } /* end-else */
		 } /* end for i */
		 jcs[i]=k;
	      } /* end if-else PRE->isblock */
	   } /* end if-else */
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}


	/* mexPrintf("5. field `U'\n");fflush(stdout); */
	/*  5. field `U' */
	++ifield;
	/* switched to full-matrix processing */
	fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* mexPrintf("6. field `E'\n");fflush(stdout); */
	/*  6. save level size to field `E' */
	++ifield;
	if (param->rcomflag==0) {
	   if (jstruct<PRE->nlev-1){ 
	
	      if (param->flags&COARSE_REDUCE) {
		 if (PRE->issingle) {
		    nnz=scurrent->F.ia[scurrent->nB]-1;
		    fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)scurrent->nB, nnz, mxREAL);
		 }
		 else {
		    nnz=current->F.ia[current->nB]-1;
		    fout=mxCreateSparse((mwSize)n-current->nB,(mwSize)current->nB, nnz, mxREAL);
		 }

		 sr =(double *) mxGetPr(fout);
		 irs=(mwIndex *)mxGetIr(fout);
		 jcs=(mwIndex *)mxGetJc(fout);
		 
		 k=0;
		 for (i=0; i<((PRE->issingle)?scurrent->nB:current->nB); i++) {
		     jcs[i]=k;
		     k_old=k;
		     if (PRE->issingle) {
		        for (j=scurrent->F.ia[i]-1; j<scurrent->F.ia[i+1]-1; j++) {
			    irs[k] =scurrent->F.ja[j]-1;
			    sr[k++]=scurrent->F.a[j];
			}
		     }
		     else {
		        for (j=current->F.ia[i]-1; j<current->F.ia[i+1]-1; j++) {
			    irs[k] =current->F.ja[j]-1;
			    sr[k++]=current->F.a[j];
			}
		     }
		     len=k-k_old;
		     /* sort entries in increasing order */
		     Dqsort(sr+k_old,(integer *)irs+k_old,stack,&len);	     
		 }
		 jcs[i]=k;
	      }
	      else {
		 nnz=0;
		 if (PRE->issingle) 
		    fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)scurrent->nB, (mwSize)nnz, mxREAL);
		 else
		    fout=mxCreateSparse((mwSize)n-current->nB, (mwSize)current->nB,  (mwSize)nnz, mxREAL);
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


	
	/* mexPrintf("7. field `F'\n");fflush(stdout); */
	/*  7. field `F' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	


	/* mexPrintf("8. field `rowscal'\n");fflush(stdout); */
	/*  8. row scaling field `rowscal' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);

	if (param->rcomflag==0) {
	   if (PRE->issingle) {
	      pr=mxGetPr(fout);
	      for (i=0; i<n; i++)
		  pr[i]=(double)scurrent->rowscal[i];
	   }
	   else {
	      pdata=mxGetData(fout);
	      memcpy(pdata, current->rowscal, (size_t)n*sizeof(double));
	   }
	}
	
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	


	/* mexPrintf("9. field `colscal'\n");fflush(stdout); */
	/*  9. column scaling field `colscal' */
	++ifield;	
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);
	
	if (param->rcomflag==0) {
	   if (PRE->issingle) {
	      pr=mxGetPr(fout);
	      for (i=0; i<n; i++)
		  pr[i]=(double)scurrent->colscal[i];
	   }
	   else {
	      pdata=mxGetData(fout);
	      memcpy(pdata, current->colscal, (size_t)n*sizeof(double));
	   }
	}

	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* mexPrintf("10. field `p'\n");fflush(stdout); */
	/* 10. field `p' */
	++ifield;	
	if (param->rcomflag==0) {
	   fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);
	   pdata = mxGetData(fout);
	
	   if (PRE->issingle)
	      for (i=0;i<n; i++)
		  convert[i]=scurrent->p[i];
	   else
	      for (i=0;i<n; i++)
		  convert[i]=current->p[i];
	   /* switched to full-matrix processing, adapt permutation */
	   if (jstruct==PRE->nlev-1 && 
	       ((PRE->issingle)?scurrent->LU.ja:current->LU.ja)==NULL && 
	       ((PRE->issingle)?(void *)scurrent->LU.a:(void *)current->LU.a )!=NULL) {

	     /* complete data using the return values from LAPACK's
		Bunch-Kaufman pivoting */
	     i=0;
	     while (i<n) {
	           /* 2x2 pivot */
	           j=((PRE->issingle)?scurrent->LU.ia[i]:current->LU.ia[i]);
		   if (j<0) {
		      k=convert[i+1];
		      convert[i+1]=convert[-j-1];
		      convert[-j-1]=k;
		      i+=2;
		   }
		   else { /* 1x1 pivot */
		      if (j!=i+1) {
			 k=convert[i];
			 convert[i]=convert[j-1];
			 convert[j-1]=k;
		      }
		      i++;
		   }
	     } 
	   }
	   
	   memcpy(pdata, convert, (size_t)n*sizeof(double));
	
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}

	
	
	/* mexPrintf("11. field `invq'\n");fflush(stdout); */
	/* 11. field `invq' */
	++ifield;	
	if (param->rcomflag==0) {
	   fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)n, mxREAL);
	   pr = mxGetPr(fout);
	   
	   for (i=0;i<n; i++) {
	       k=convert[i];
	       pr[k-1]=i+1;
	   }
	   /* memcpy(pdata, convert, (size_t)n*sizeof(double)); */
	
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	else {
	   fout=mxCreateDoubleMatrix((mwSize)0,(mwSize)0, mxREAL);
	   /* set each field in output structure */
	   mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	}
	
		
	
	/* mexPrintf("12. field `param'\n");fflush(stdout); */
	/* 12. field `param' */
	++ifield;	
	fout=mxCreateNumericArray((mwSize)1,mydims, mxUINT64_CLASS, mxREAL);
	pdata = mxGetData(fout);
	
	memcpy(pdata, &param, (size_t)sizeof(size_t));
	
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	
	
	
	
	/* mexPrintf("13. field `ptr'\n");fflush(stdout); */
	/* 13. field `ptr' */
	++ifield;	
	fout=mxCreateNumericArray((mwSize)1,mydims, mxUINT64_CLASS, mxREAL);
	pdata = mxGetData(fout);
	
	memcpy(pdata, &PRE, (size_t)sizeof(size_t));
	
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);
	


	/* mexPrintf("14. field `isreal'\n");fflush(stdout);*/
	/* 14. save real property to field `isreal' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=1;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* mexPrintf("15. field `isdefinite'\n");fflush(stdout); */
	/* 15. save non-positive definite property to field `isdefinite' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=0;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* mexPrintf("16. field `issymmetric'\n");fflush(stdout); */
	/* 16. save symmetry property to field `issymmetric' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=1;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* mexPrintf("17. field `ishermitian'\n");fflush(stdout); */
	/* 17. save trivial Hermitian property to field `ishermitian' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=1;
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);


	
	/* mexPrintf("18. field `issingle'\n");fflush(stdout); */
	/* 18. save trivial Hermitian property to field `issingle' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	*pr=((PRE->issingle)?scurrent->issingle:current->issingle);
       
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);


	
	/* mexPrintf("19. field `A_H'\n");fflush(stdout); */
	/* 19. save coarse grid system `A_H' */
	++ifield;
	if (jstruct>=PRE->nlev-1) {
	   if (PRE->issingle)
	      fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)n-scurrent->nB,(mwSize)0, mxREAL);
	   else
	      fout=mxCreateSparse((mwSize)n-current->nB, (mwSize)n-current->nB, (mwSize)0, mxREAL);
	}
	else if (param->ipar[16]&DISCARD_MATRIX) {
	     if (PRE->issingle)
	        fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)n-scurrent->nB,(mwSize)0, mxREAL);
	     else
	        fout=mxCreateSparse((mwSize)n-current->nB, (mwSize)n-current->nB, (mwSize)0, mxREAL);
	}
	else {
	   /* switched to full-matrix processing */
	   if (jstruct==PRE->nlev-2 && 
	       ((PRE->issingle)?scurrent->next->LU.ja:current->next->LU.ja)==NULL) {
	      if (PRE->issingle)
	         fout=mxCreateSparse((mwSize)n-scurrent->nB,(mwSize)n-scurrent->nB,(mwSize)0, mxREAL);
	      else
	         fout=mxCreateSparse((mwSize)n-current->nB, (mwSize)n-current->nB, (mwSize)0, mxREAL);
	   }
	   else {
	      if (PRE->issingle)
		 nnz=scurrent->next->A.ia[scurrent->next->A.nr]-1;
	      else
		 nnz=current->next->A.ia[current->next->A.nr]-1;

	      /*
	      mexPrintf("level %d, %d,%d,%d\n",jstruct+1,current->next->A.nr,current->next->A.nc,nnz);fflush(stdout);
	      fout=mxCreateSparse(0,0,0, mxREAL);
	      */
	      
	      if (PRE->issingle)
		 fout=mxCreateSparse((mwSize)scurrent->next->A.nr,(mwSize)scurrent->next->A.nc,(mwSize)nnz, mxREAL);
	      else
		 fout=mxCreateSparse((mwSize)current->next->A.nr, (mwSize)current->next->A.nc, (mwSize)nnz, mxREAL);
  
	      sr =(double *) mxGetPr(fout);
	      irs=(mwIndex *)mxGetIr(fout);
	      jcs=(mwIndex *)mxGetJc(fout);
	      
	      k=0;
	      for (i=0; i<((PRE->issingle)?scurrent->next->A.nr:current->next->A.nr); i++) {
		  jcs[i]=k;
		  k_old=k;
		  if (PRE->issingle) {
		     for (j=scurrent->next->A.ia[i]-1; j<scurrent->next->A.ia[i+1]-1; j++) {
		         irs[k] =scurrent->next->A.ja[j]-1;
			 sr[k++]=scurrent->next->A.a[j];
		     }
		  }
		  else {
		     for (j=current->next->A.ia[i]-1; j<current->next->A.ia[i+1]-1; j++) {
		         irs[k] =current->next->A.ja[j]-1;
			 sr[k++]=current->next->A.a[j];
		     }
		  }
		  len=k-k_old;
		  /* sort entries in increasing order */
		  Dqsort(sr+k_old,(integer *)irs+k_old,stack,&len);	     
	      }
	      jcs[i]=k;
	   }
	}
	/* set each field in output structure */
	mxSetFieldByNumber(PRE_output, (mwIndex)jstruct, ifield, fout);



	/* mexPrintf("20. field `errorL'\n");fflush(stdout); */
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

	/* mexPrintf("21. field `errorU'\n");fflush(stdout); */
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

	/* mexPrintf("22. field `errorS'\n");fflush(stdout); */
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


	/* mexPrintf("23. field `isblock'\n");fflush(stdout); */
	/* 23. save trivial block property to field `isblock' */
	++ifield;
	fout=mxCreateDoubleMatrix((mwSize)1,(mwSize)1, mxREAL);
	pr = mxGetPr(fout);

	if (PRE->issingle)
	   *pr=scurrent->isblock;
	else
	   *pr=current->isblock;
       
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

    free(convert);
    free(stack);

#ifdef PRINT_INFO
    mexPrintf("DSYMilupackfactor: memory released\n");fflush(stdout);
#endif
    
    return;
}

