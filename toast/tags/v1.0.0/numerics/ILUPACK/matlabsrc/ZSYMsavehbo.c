/* ========================================================================== */
/* === savehbo mexFunction ================================================== */
/* ========================================================================== */

/*
    Usage:

    saves matrix A and optionally b, x0 and x to a file in Harwell-Boeing
    format

    Example:

    % for initializing parameters
    ZSYMsavehbo(filename, A,b,x0,x);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	March 05, 2008. ILUPACK V2.2.  

    Acknowledgements:

	This work was partially supported by the DFG research center
        MATHEON "Mathematics for key technologies", 2002-2006

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

    char       *pdata, fname[205], fname2[8],tname[4];
    int        i,j,k,l, status;
    integer    m, ncolumns;
    mwSize     nnz, buflen;
    size_t     mrows, ncols, sizebuf;
    double     *A_valuesR, *A_valuesI, *pr, *pi;
    doublecomplex *sol, *rhs, *x0;
    mxArray    *f_input, *A_input , *b_input, *x0_input, *x_input;
    mwIndex    *A_ja,                 /* row indices of input matrix A */
               *A_ia;                 /* column pointers of input matrix A */
    


    if (nrhs<2)
       mexErrMsgTxt("At least two input arguments required.");
    else if (nlhs>5)
       mexErrMsgTxt("At most five input arguments are allowed.");
    else if (mxGetClassID(prhs[0])!=mxCHAR_CLASS)
       mexErrMsgTxt("First input must be a string.");
    else if (!mxIsNumeric(prhs[1]))
       mexErrMsgTxt("Second input must be a matrix.");

    if (nrhs>2) {
       if (!mxIsNumeric(prhs[2]))
	  mexErrMsgTxt("Third input must be a matrix or vector.");
       if (mxIsSparse(prhs[2]))
	  mexErrMsgTxt("b must be dense.");
    }
    if (nrhs>3) {
       if (!mxIsNumeric(prhs[3]))
	  mexErrMsgTxt("Fourth input must be a matrix or vector.");
       if (mxIsSparse(prhs[3]))
	  mexErrMsgTxt("x0 must be dense.");
    }
    if (nrhs>4) {
       if (!mxIsNumeric(prhs[4]))
	  mexErrMsgTxt("Fifth input must be a matrix or vector.");
       if (mxIsSparse(prhs[4]))
	  mexErrMsgTxt("x0 must be dense.");
    }


    /* get filename */
    f_input = (mxArray *) prhs[0];
    mrows = mxGetM (f_input) ;
    ncols = mxGetN (f_input) ;
    /* Get the length of the input string. */
    buflen = (mrows*ncols) + 1;
    if (buflen>200)
       buflen=200;

    /* Allocate memory for input and output strings. */
    pdata = (char *) mxCalloc((size_t)buflen, (size_t)sizeof(char));

    /* Copy the string data from tmp into a C string pdata */
    status = mxGetString(f_input, pdata, buflen);

    sizebuf = buflen*sizeof(char);
    memcpy(fname, pdata, (size_t)sizebuf);
    memcpy(fname2, pdata, (size_t)((sizebuf<8*sizeof(char))?sizebuf:8*sizeof(char)));

    /* skip '\0' character */
    buflen--;
    i=buflen;
    while (i<8)
          fname2[i++]=' ';
    fname[buflen++]='.';
    fname[buflen++]='c';
    fname[buflen++]='s';
    fname[buflen++]='a';
    fname[buflen]  ='\0';
    /* printf("\n%s\n",fname); */


    A_input = (mxArray *) prhs [1] ;
    /* get size of input matrix A */
    mrows = mxGetM (A_input) ;
    ncols = mxGetN (A_input) ;
    nnz = mxGetNzmax(A_input);
    if (!mxIsSparse (A_input)) {
       mexErrMsgTxt ("ILUPACK: input matrix must be in sparse format.") ;
    }


    /* copy input matrix to sparse row format */
    A.nc=ncols;
    A.nr=mrows;
    A.ia=(integer *)      MAlloc((size_t)(A.nc+1)*sizeof(integer),       "ZSYMilupackfactor");
    A.ja=(integer *)      MAlloc((size_t)nnz     *sizeof(integer),       "ZSYMilupackfactor");
    A. a=(doublecomplex*) MAlloc((size_t)nnz     *sizeof(doublecomplex), "ZSYMilupackfactor");

    A_ja         = (mwIndex *)mxGetIr(A_input) ;
    A_ia         = (mwIndex *)mxGetJc(A_input) ;
    A_valuesR    = (double *) mxGetPr(A_input);
    A_valuesI    = (double *) mxGetPi(A_input);

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
    for (i=0; i<ncols; i++) {
        A.ia[i+1]=A.ia[i];
	for (j=A_ia[i]; j<A_ia[i+1]; j++) {
	    k=A_ja[j];
	    if (k>=i) {
	       l=A.ia[i+1]-1;
	       A.ja[l]=k+1;
	       A.a [l].r=A_valuesR[j];
	       A.a [l].i=A_valuesI[j];
	       A.ia[i+1]=l+2;
	    }
	}
    }

    /*
    for (i = 0 ; i < A.nr ; i++)
      for (j = A.ia[i]-1 ; j < A.ia[i+1]-1 ; j++)
	  printf("i=%d ja=%d  A.real=%e\n", i+1,  A.ja[j], A.a[j]);
    */


    /* copy right hand side `b' */
    if (nrhs>2) {
       b_input = (mxArray *) prhs [2] ;
       k=mxGetM(b_input);
       l=mxGetN(b_input);
       if (k!=mrows)
	 mexErrMsgTxt("b must have same number of rows as A");
       mrows=k;
       ncols=l;
       i=mrows*ncols;
    }
    else {
       i=0;
       ncols=0;
    }

    /* copy initial guess `x0' */
    if (nrhs>3) {
       x0_input = (mxArray *) prhs [3] ;
       k=mxGetM(x0_input);
       l=mxGetN(x0_input);
       if (k!=mrows)
	 mexErrMsgTxt("x0 must have same number of rows as A and b");
       if (l!=ncols)
	 mexErrMsgTxt("x0 must have same number of columns as b");
       i+=mrows*ncols;
    }


    /* copy initial guess `x' */
    if (nrhs>4) {
       x_input = (mxArray *) prhs [4] ;
       k=mxGetM(x_input);
       l=mxGetN(x_input);
       if (k!=mrows)
	 mexErrMsgTxt("x must have same number of rows as A and b and x0");
       if (l!=ncols)
	 mexErrMsgTxt("x must have same number of columns as b and x0");
       i+=mrows*ncols;
    }

    rhs=(doublecomplex*)MAlloc(i*sizeof(doublecomplex),"ZSYMsavehbo:rhs");
    sol=rhs;
    tname[0]=' ';
    tname[1]=' ';
    tname[2]=' ';
    if (nrhs>2) {
      pr=mxGetPr(b_input);
      if (mxIsComplex(b_input)) 
	 pi=mxGetPr(b_input);
      for (j=0; j<mrows*ncols; j++) {
	  sol->r=*pr++;
	  if (mxIsComplex(b_input)) 
	     sol->i=*pi++;
	  else
 	     sol->i=0;
	  sol++;
      }
      tname[0]='F';
    }

    if (nrhs>3) {
      pr=mxGetPr(x0_input);
      if (mxIsComplex(x0_input)) 
	 pi=mxGetPr(x0_input);
      for (j=0; j<mrows*ncols; j++) {
	  sol->r=*pr++;
	  if (mxIsComplex(x0_input)) 
	     sol->i=*pi++;
	  else
 	     sol->i=0;
	  sol++;
      }
      tname[1]='G';
    }

    if (nrhs>4) {
      pr=mxGetPr(x_input);
      if (mxIsComplex(x_input)) 
	 pi=mxGetPr(x_input);
      for (j=0; j<mrows*ncols; j++) {
	  sol->r=*pr++;
	  if (mxIsComplex(x_input)) 
	     sol->i=*pi++;
	  else
 	     sol->i=0;
	  sol++;
      }
      tname[2]='X';
    }

    nnz=A.ia[A.nr]-1;
    ncolumns=ncols;
    m=nnz;
    Zwritemtc(fname,A.a,A.ja,A.ia,rhs,&ncolumns,
	      tname,&(A.nr),&(A.nc),&m,
	      "ILUPACK toolbox for MATLAB. Export routine to Harwell-Boeing format     ",fname2,"csa",strlen(fname),3,72,8,3);

    free(A.ia);
    free(A.ja);
    free(A.a);
    if (i>0)
       free(rhs);
    return;
}

