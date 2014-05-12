/* ========================================================================== */
/* === loadhbo mexFunction ================================================== */
/* ========================================================================== */

/*
    Usage:

    saves matrix A and optionally b, x0 and x to a file in Harwell-Boeing
    format

    Example:

    % load Harwell-Boeing matrix
    [A,rhs,rhstyp]=Dloadhbo(filename);



    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

	January 14, 2007. ILUPACK V2.2.  

    Acknowledgements:

	This work was supported from 2002-2006 by the DFG research center
        MATHEON "Mathematics for key technologies" in Berlin

    Notice:

	Copyright (c) 2007 by TU Braunschweig.  All Rights Reserved.

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
#include <stdio.h>


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

    char       *fname, rhstyp[4], title[81], key[9], type[4];
    int        i,j,k,l, status, mynrhs;
    integer    nrhsix, m,*rhsptr=NULL, *rhsind=NULL, ierr, nr,nc,nz, nmax,nzmax,job,ncolumns;
    mwSize     nnz, buflen;
    mwIndex    *irs, *jcs;
    size_t     mrows, ncols, sizebuf;
    double     *pr, *sol, *rhs=NULL, *x0, *rhsval=NULL, *sr;
    mxArray    *fout, *f_input;
    FILE *fp;


    if (nrhs!=1)
       mexErrMsgTxt("Only one input argument required.");
    else if (nlhs!=3)
       mexErrMsgTxt("Three output arguments are required.");
    else if (mxGetClassID(prhs[0])!=mxCHAR_CLASS)
       mexErrMsgTxt("Input must be a string.");


    /* get filename */
    f_input = (mxArray *) prhs[0];
    mrows = mxGetM (f_input) ;
    ncols = mxGetN (f_input) ;
    /* Get the length of the input string. */
    buflen = (mrows*ncols) + 1;

    /* Allocate memory for input and output strings. */
    fname = (char *)mxCalloc((size_t)buflen, (size_t)sizeof(char));

    /* Copy the string data from tmp into a C string pdata */
    status = mxGetString(f_input, fname, buflen);


    ncolumns = 0;
    nmax=0;
    nzmax=0;
    job=0;
    nr=nc=nz=0;
    nrhsix=0;

    if ((fp=fopen(fname,"r"))==NULL) {
       mexPrintf(" file %s",fname);
       mexErrMsgTxt(" not found");
       return;
    }
    fclose(fp);

    rhstyp[0]=' ';
    rhstyp[1]=' ';
    rhstyp[2]=' ';


    A.nc=0;
    Dreadmtc(&nmax,&nzmax,&job,fname,A.a,A.ja,A.ia,
	     rhs,&ncolumns,rhstyp,&nr,&nc,&nz,title,key,type,
	     &nrhsix,rhsptr,rhsind,rhsval,&ierr,strlen(fname),2,72,8,3);
    ncols=ncolumns;


    if (ierr) {
        mexPrintf(" ierr = %d\n",ierr);
        mexErrMsgTxt(" error in reading the matrix, stop.\n");
	switch(ierr) {
	case 1:
	  mexErrMsgTxt("too many columns\n");
	  break;  
	case 2:
	  mexErrMsgTxt("too many nonzeros\n");
	  break;  
	case 3:
	  mexErrMsgTxt("too many columns and nonzeros\n");
	  break;  
	case 4:
	  mexErrMsgTxt("right hand side has incompatible type\n");
	  break;  
	case 5:
	  mexErrMsgTxt("too many right hand side entries\n");
	  break;  
	case 6:
	  mexErrMsgTxt("wrong type (real/complex)\n");
	  break;  
	}
        exit(ierr);
    }

    if (ncols>0) {
       m=1;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') {
	 m++;
       }
       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	 m++;
       }
    }
    else
       m=0;
    m*=nr*ncols;

    rhsptr=NULL;
    rhsind=NULL;
    rhsval=NULL;
    if (ncols!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
       rhsptr=(integer *)  MAlloc((size_t)(ncols+1)*sizeof(integer),"Dloadhbo:rhsptr");
       rhsind=(integer *)  MAlloc((size_t)nrhsix*sizeof(integer),   "Dloadhbo:rhsind");
       rhsval=(double *)   MAlloc((size_t)nrhsix*sizeof(double),    "Dloadhbo:rhsval");
    }
    

    A.ia=(integer *)MAlloc((size_t)(nc+1)*sizeof(integer),"Dloadhbo:A.ia");
    A.ja=(integer *)MAlloc((size_t)nz*sizeof(integer),    "Dloadhbo:A.ja");
    A.a =(double *) MAlloc((size_t)nz*sizeof(double),     "Dloadhbo:A.a");
    A.nr=nr;
    A.nc=nc;
    rhs =(double *)MAlloc((size_t)m*sizeof(double),  "Dloadhbo:rhs");
    /* advance pointer to reserve space when uncompressing the right hand side */
    if (ncols!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       rhs+=nr*ncols;
    
    job = 3;
    nmax = nc;
    nzmax= nz;

    if (ncols!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       m-=nr*ncols;
    Dreadmtc(&nmax,&nzmax,&job,fname,A.a,A.ja,A.ia,
	     rhs,&m,rhstyp,&nr,&nc,&nz,title,key,type,
	     &nrhsix,rhsptr,rhsind,rhsval,&ierr,strlen(fname),2,72,8,3);


    if (ncols!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       m+=nr*ncols;



    if (ierr) {
        mexPrintf(" ierr = %d\n",ierr);
        mexErrMsgTxt(" error in reading the matrix, stop.\n");
	return;
    }


    /* set output parameters */
    nlhs=3;
    fout =mxCreateSparse(nr,nc, nz, mxREAL);
    plhs[0]=fout;

    sr  = (double *)  mxGetPr(fout);
    irs = (mwIndex *) mxGetIr(fout);
    jcs = (mwIndex *) mxGetJc(fout);

    jcs[0]=0;
    for (i=0; i<nc; i++) {
        jcs[i+1]=A.ia[i+1]-1;
	for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	    irs[j]=A.ja[j]-1;
	    sr[j] =A.a[j];
	}
    }


    if (ncols>0) {
       m=1;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') {
	 m++;
       }
       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	 m++;
       }
    }
    else
       m=0;


    fout=mxCreateDoubleMatrix((mwSize)nr,(mwSize)m*ncols, mxREAL);
    plhs[1]=fout;

    if (ncols!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
    }
    else {
	pr = mxGetPr(fout);
	memcpy(pr,rhs,(size_t)nr*m*ncols*sizeof(double));
    }

    rhstyp[3]='\0';
    plhs[2] = mxCreateString(rhstyp);



    free(A.a);
    free(A.ia);
    free(A.ja);
    if (ncols>0)
      free(rhs);

    if (ncols!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
      free(rhsptr);
      free(rhsind);
      free(rhsval);
    }

    return;
}

