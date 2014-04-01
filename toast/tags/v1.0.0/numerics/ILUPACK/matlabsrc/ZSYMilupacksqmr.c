/* ========================================================================== */
/* === DSYMilupacksqmr mexFunction ========================================== */
/* ========================================================================== */

/*
    Usage:

    simplified QMR for real symmetric matrices and real symmetric preconditioning
    
    Example:

    % for initializing parameters
    [x,src,control,iter,relres]=DSYMilupacksqmr(b,x,drain,tol,maxit,control,nrmA,typeres)

    b       right hand side
    x      initial guess
    drain   returns the result of drain=A*src of matrix vector-multiplication
            or the result of preconditioning drain=P*src, depending on the
	    value of control
    tol     tolerance for termination of the iterative process
    maxit   maximum number of iteration steps
    control parameter for reverse communication
    nrmA    estimate for the 2-norm of A
    typeres type of residual 
            1 ||Ax_k-b|| <= tol ||Ax_0-b||
            2 ||Ax_k-b|| <= tol ||b||
            3 ||Ax_k-b|| <= tol (||A||+||x_k||+||b||)

    x       approximate solution of Ax=b
    src     source vector for matrix-vector multiplication or preconditioning 
    control parameter for reverse communication
    iter    current number of iterative steps
    relres  residual from the current step

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        May 02, 2008. ILUPACK V2.2.  

    Acknowledgements:

	This work was supported from 2002 to 2007 by the DFG research center
        MATHEON "Mathematics for key technologies" at TU Berlin

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
    mxArray *plhs[],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs[]	/* right-hand side matrices */
)
{
    mxArray *b_input;
    integer i, m;
    size_t  n; 
    double  *pr, *pi;
    doublecomplex *rhs, *sol;

    static integer ipar[16];
    static double  fpar[16];
    static doublecomplex *workspace;
    static int init=0;

    if (nrhs != 8)
       mexErrMsgTxt("Eight input argument required.");
    else if (nlhs!=5)
       mexErrMsgTxt("Five output arguments are required.");
    else if (!mxIsNumeric(prhs[0]))
       mexErrMsgTxt("First input must be a matrix.");



    /* First input parameter: right hand side */
    b_input = (mxArray *) prhs[0] ;
    pr=mxGetPr(b_input);
    n=mxGetM(b_input);
    rhs=(doublecomplex *) MAlloc((size_t)n*sizeof(doublecomplex), "DSYMilupacksqmr");
    if (!mxIsComplex(b_input)) {
       for (i=0; i<n; i++) {
	   rhs[i].r=pr[i];
	   rhs[i].i=0;
       }
    }
    else {
       pi=mxGetPi(b_input);
       for (i=0; i<n; i++) {
	   rhs[i].r=pr[i];
	   rhs[i].i=pi[i];
       }
    }

    /* Second input parameter: initial guess */
    b_input = (mxArray *) prhs[1] ;
    pr=mxGetPr(b_input);
    sol=(doublecomplex *) MAlloc((size_t)n*sizeof(doublecomplex), "DSYMilupacksqmr");
    if (!mxIsComplex(b_input)) {
       for (i=0; i<n; i++) {
	   sol[i].r=pr[i];
	   sol[i].i=0;
       }
    }
    else {
       pi=mxGetPi(b_input);
       for (i=0; i<n; i++) {
	   sol[i].r=pr[i];
	   sol[i].i=pi[i];
       }
    }

    /* Third input parameter: drain, result from matrix-vector multiplication or preconditioning */
    /* we must have called the routine earlier */
    if (init) {
       b_input = (mxArray *) prhs[2] ;
       pr=mxGetPr(b_input);
       workspace--;
       if (!mxIsComplex(b_input)) {
	  for (i=0; i<n; i++) {
	      workspace[ipar[8]+i].r=pr[i];
	      workspace[ipar[8]+i].i=0;
	  }
       }
       else {
	  pi=mxGetPi(b_input);
	  for (i=0; i<n; i++) {
	      workspace[ipar[8]+i].r=pr[i];
	      workspace[ipar[8]+i].i=pi[i];
	  }
       }
       workspace++;
    }

    
    /* Fourth input parameter: tol */
    b_input = (mxArray *) prhs[3] ;
    pr=mxGetPr(b_input);
    fpar[0]=*pr;

    /* Fifth input parameter: maxit */
    b_input = (mxArray *) prhs[4] ;
    pr=mxGetPr(b_input);
    ipar[5]=*pr;

    /* Sixth input parameter: control */
    if (init) {
       b_input = (mxArray *) prhs[5] ;
       pr=mxGetPr(b_input);
       if (*pr==1) 
	  ipar[0]=1;
       else if (*pr==2)
	  ipar[0]=5;
    }


    /* so far use 'right' preconditioning */
    ipar[1]=2;

    /* size of our workspace */
    ipar[3]=6*n;
    /* not referenced */
    ipar[4]=0;

    /* no absolute tolerance */
    fpar[1]=0.0;

    /* init parameters for sqmr when called for the first time */
    if (!init) {
       init=1;
       /* init solver */
       ipar[0]=0;
       /* init with zeros */
       ipar[6]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=
	 ipar[14]=ipar[15]=0;
       fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=
	 fpar[11]=fpar[12]=fpar[13]=fpar[14]=fpar[15]=0.0;

       /* provide workspace */
       workspace=(doublecomplex *) MAlloc(6*(size_t)n*sizeof(doublecomplex),"DSYMilupacksqmr");
    }


    /* Seventh input parameter: nrmA */
    b_input = (mxArray *) prhs[6] ;
    pr=mxGetPr(b_input);
    fpar[11]=*pr;

    /* Eighth input parameter: typeres */
    b_input = (mxArray *) prhs[7] ;
    pr=mxGetPr(b_input);
    /* use backward error as stopping criterion */
    ipar[2]=*pr;

    m=n;
    ZSYMqmr(&m,rhs,sol,ipar,fpar,workspace);



    /* Create output vector */
    nlhs=5;

    /* first left hand side: x */
    plhs[0] =mxCreateDoubleMatrix(n,1, mxCOMPLEX);
    pr = (double *) mxGetPr(plhs[0]);
    pi = (double *) mxGetPi(plhs[0]);
    for (i=0;i<n; i++) {
        pr[i]=sol[i].r;
        pi[i]=sol[i].i;
    }

    /* second left hand side: src for matrix-vector multiplication or preconditioning */
    plhs[1] =mxCreateDoubleMatrix(n,1, mxCOMPLEX);
    pr = (double *) mxGetPr(plhs[1]);
    pi = (double *) mxGetPi(plhs[1]);
    workspace--;
    for (i=0; i<n; i++) {
        pr[i]=workspace[ipar[7]+i].r;
        pi[i]=workspace[ipar[7]+i].i;
    }
    workspace++;


    /* third left hand side: control parameter */
    plhs[2] =mxCreateDoubleMatrix(1,1, mxREAL);
    pr = (double *) mxGetPr(plhs[2]);
    if (ipar[0]<=1)
      *pr=ipar[0];
    else if (ipar[0]==5)
      *pr=2;
    else if (ipar[0]==-1)
      *pr=-1;
    else if (ipar[0]==-2)
      *pr=-2;
    else if (ipar[0]>0)
      *pr=ipar[0];
    else 
      mexErrMsgTxt("undefined sqmr error code");

    /* fourth left hand side: iter */
    plhs[3] =mxCreateDoubleMatrix(1,1, mxREAL);
    pr = (double *) mxGetPr(plhs[3]);
    *pr=ipar[6];

    /* fifth left hand side: relres */
    plhs[4] =mxCreateDoubleMatrix(1,1, mxREAL);
    pr = (double *) mxGetPr(plhs[4]);
    *pr=fpar[4];


    free(sol);
    free(rhs);

    /* if we finally decide to quit, the we release the workspace */
    if (ipar[0]!=1 && ipar[0]!=5) {
       init=0;
       free(workspace);
    }

    return;
}

