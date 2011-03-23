#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <sparspak.h>
#include <ilupack.h>

#include <ilupackmacros.h>

#define MAX_LINE        255
#define STDERR          stdout
#define STDOUT          stdout
#define PRINT_INFO
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))

// maximum number of iterations independent on n
#define MAX_IT          500

// measure for terminating iterative solver
#define RESTOL_FUNC(A)  sqrt(A)
//#define RESTOL_FUNC(A)  (A)

// use SPD scaling
#define USE_SCALE


// use an iterative solver from SPARSKIT defined by variable SOLVER
#define USE_SPARSKIT







int main(int argc, char **argv)
{
    /* SOLVER choice:   1  pcg
                        2  cgnr
                        3  bcg
                        4  dbcg
                        5  bcgstab
                        6  tfqmr
                        7  fom
                        8  gmres
                        9  fgmres
                       10  dqgmres */
    integer SOLVER=1; /* gmres */

    CSRMAT A, ilutmat;

    FILE *fp, *fo; 
    char rhstyp[3], title[72], key[8], type[3], fname[100], foname[100];
    char line[MAX_LINE], *tmpstring, timeout[7], residout[7];
    integer  i,j,k,m,fnamelen,n,nc,nz,nrhs,tmp0,tmp,tmp2,tmp3,ierr,ipar[16],flag,
         *invperm, *buff, *ws, *symmmd,
         nrhsix, *rhsptr, *rhsind;
    REALS eps,fpar[16], DROP_TOL;
    FLOAT *rhs,*sol,*w, *dbuff, *scale, *rhsval;
    float  systime, time_start,   time_stop,   secnds, secndsprep;
    integer ELBOW, param;
    size_t nibuff, ndbuff;

    /* the last argument passed serves as file name */
    if (argc!=4) {
      printf("usage '%s <drop tol.> <elbow space> <matrix>'\n",argv[0]);
       exit(0);
    }
    i=0;
    while (argv[argc-1][i]!='\0')
    {
          fname[i]=argv[argc-1][i];
          i++;
    } /* end while */
    fname[i]='\0';
    fnamelen=i;
    while (i>=0 && fname[i]!='/')
          i--;
    i++;
    j=0;
    while (i<fnamelen && fname[i]!='.')
          foname[j++]=fname[i++];
    while (j<16)
          foname[j++]=' ';
    foname[j]='\0';

    ELBOW=atoi(argv[argc-2]);
    DROP_TOL=atof(argv[argc-3]);


/*----------------------------------------------------------------------
|     Read a Harwell-Boeing matrix.
|     Use readmtc first time to determine sizes of arrays.
|     Read in values on the second call.
|---------------------------------------------------------------------*/

    nrhs = 0;
    tmp0 = 0;

    if ((fp=fopen(fname,"r"))==NULL) {
        fprintf(STDERR," file %s not found\n",fname);
        exit(0);
    }
    fclose(fp);

    READMTC(&tmp0,&tmp0,&tmp0,fname,A.a,A.ja,A.ia,
	    rhs,&nrhs,rhstyp,&n,&nc,&nz,title,key,type,
	    &nrhsix,rhsptr,rhsind,rhsval,&ierr,fnamelen,2,72,8,3);
    if (ierr) {
        fprintf(STDERR," ierr = %d\n",ierr);
        fprintf(STDERR," error in reading the matrix, stop.\n");
	switch(ierr) {
	case 1:
	  fprintf(STDERR,"too many columns\n");
	  break;  
	case 2:
	  fprintf(STDERR,"too many nonzeros\n");
	  break;  
	case 3:
	  fprintf(STDERR,"too many columns and nonzeros\n");
	  break;  
	case 4:
	  fprintf(STDERR,"right hand side has incompatible type\n");
	  break;  
	case 5:
	  fprintf(STDERR,"too many right hand side entries\n");
	  break;  
	case 6:
	  fprintf(STDERR,"wrong type (real/complex)\n");
	  break;  
	}
        exit(ierr);
    }
    printf("Matrix: %s: size (%d,%d), nnz=%d(%4.1lfav.)\n", fname, n,nc,
	   nz,((double)nz)/n);

    if (fname[fnamelen-1]=='5')
      fo = fopen("out_mc64","aw");
    else
      fo = fopen("out_normal","aw");

    fprintf(fo,"%s|%7.1e|       |",foname,DROP_TOL);

    m=1;
    if (nrhs>0) {
       printf ("Number of right hand sides supplied: %d \n", nrhs) ;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') {
	 m++;
	 printf ("Initial solution(s) offered\n") ;
       }
       else
	 printf ("\n") ;

       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	 m++;
	 printf ("Exact solution(s) provided\n") ;
       }
       else
	 printf ("\n") ;
    }
    else 
      printf("\n\n\n");
    printf("\n");

    rhsptr=NULL;
    rhsind=NULL;
    rhsval=NULL;
    if (rhstyp[0]=='M' || rhstyp[0]=='m') {
       rhsptr=(integer *)  MALLOC((size_t)(nrhs+1)*sizeof(integer),"main:rhsptr");
       rhsind=(integer *)  MALLOC((size_t)nrhsix*sizeof(integer),  "main:rhsind");
       rhsval=(FLOAT *)MALLOC((size_t)nrhsix*sizeof(FLOAT),"main:rhsval");
       // no dense right hand side
       m--;
       m*=n*MAX(1,nrhs);
       // in any case we need at least one vector for the r.h.s.
       m+=n;
    }
    else
       m*=n*MAX(1,nrhs);
    A.ia=(integer *)  MALLOC((size_t)(n+1)*sizeof(integer),"main:A.ia");
    A.ja=(integer *)  MALLOC((size_t)nz*sizeof(integer),   "main:A.ja");
    A.a =(FLOAT *)MALLOC((size_t)nz*sizeof(FLOAT), "main:A.a");
    A.nr=n;
    A.nc=n;
    rhs       =(FLOAT *) MALLOC((size_t)m*sizeof(FLOAT),  "main:rhs");
    // advance pointer to reserve space when uncompressing the right hand side
    if (rhstyp[0]=='M' || rhstyp[0]=='m')
       rhs+=n;
    
    sol       =(FLOAT *) MALLOC((size_t)n*sizeof(FLOAT),  "main:sol");
    dbuff     =(FLOAT *) MALLOC(3*(size_t)n*sizeof(FLOAT),"main:dbuff");
    ndbuff=3*(size_t)n;

    tmp = 3;
    tmp2 = n;
    tmp3 = nz;

    if (rhstyp[0]=='M' || rhstyp[0]=='m')
       m-=n;
    READMTC(&tmp2,&tmp3,&tmp,fname,A.a,A.ja,A.ia,
	    rhs,&m,rhstyp,&n,&nc,&nz,title,key,type,
	    &nrhsix,rhsptr,rhsind,rhsval,&ierr,fnamelen,2,72,8,3);
    if (rhstyp[0]=='M' || rhstyp[0]=='m')
       m+=n;



    if (ierr) {
        fprintf(STDERR," ierr = %d\n",ierr);
        fprintf(STDERR," error in reading the matrix, stop.\n");
	fprintf(fo,"out of memory\n");fclose(fo); 
        exit(ierr);
    }
    /*
    for (i=0; i<n;i++) {
      printf("%4d:\n",i+1);
      for (j=A.ia[i];j<A.ia[i+1]; j++)
	printf("%16d",A.ja[j-1]);
      printf("\n");
      fflush(stdout);
      for (j=A.ia[i];j<A.ia[i+1]; j++) 
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	printf("%16.1e",A.a[j-1]);
#else
	printf("%8.1e%8.1e",A.a[j-1].r,A.a[j-1].i);
#endif
      printf("\n");
      fflush(stdout);
    }// end for i
    */

    // release part of rhs that may store the uncompressed rhs
    if (rhstyp[0]=='M' || rhstyp[0]=='m')
       rhs-=n;

    // if no right hand side is provided, then set sol=1 and rhs=A*sol
    if (nrhs==0) { 
       for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   sol[i]=1;
#else
	   sol[i].r=1;
	   sol[i].i=0;
#endif
       } // end for i
       HERMATVEC(A,sol,rhs);
    }
    else {
       if (rhstyp[0]=='M' || rhstyp[0]=='m') {
	  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	      rhs[i]=0;
#else
	      rhs[i].r=rhs[i].i=0;
#endif
	  }// end for i

	  // uncompress rhs
	  for (i=0; i<rhsptr[1]-rhsptr[0]; i++) {
	      j=rhsind[i]-1;
	      rhs[j]=rhsval[i];
	  } // end for i
       } // end if
    } // end if-else
    
    
    evaluate_time(&time_start,&systime);
    scale=(FLOAT *)MALLOC((size_t)n*sizeof(FLOAT),"mainildlc:scale");
    for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        scale[i]=1.0;
#else
        scale[i].r=1.0;
        scale[i].i=0.0;
#endif
    } // end for i

#ifdef USE_SCALE
    i=0;
    SPDSCALE(&(A.nc),A.a, A.ja, A.ia, scale, &i);
    
    if(i!=0) {
      fprintf(STDERR, "SPDSCALE: a zero row/column...\n" );
      return;
    }
#endif /* USE_SCALE */

    // rescale right hand side
    for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        rhs[i]*=scale[i];
#else
	rhs[i].r*=scale[i].r;
	rhs[i].i*=scale[i].r;
#endif
    } // end for i


    /*
    for (i=0; i<n; i++) {
        printf("%4d:\n",i+1);
        for (j=A.ia[i]; j<A.ia[i+1]; j++) {
	  printf("%8d",A.ja[j-1]);
	} // end for j
	printf("\n");fflush(stdout);
        for (j=A.ia[i]; j<A.ia[i+1]; j++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    printf("%16.1e",A.a[j-1]);
#else
	    printf("%8.1e%8.1e",A.a[j-1].r,A.a[j-1].i);
#endif
	} // end for j
	printf("\n");fflush(stdout);
    } // end for i
    */


    fprintf(fo,"   s/   s|");
    evaluate_time(&time_stop,&systime);
    secndsprep=time_stop-time_start;
    printf("time: %8.1e [sec]\n",(double)secndsprep);



    /* ------------------------------------------------------------------- */
    /* -----------------------   Start MAIN part   ----------------------- */
    /* ------------------------------------------------------------------- */

    fprintf(fo,"IB-ILDLC  |");










    /* set up paramaters for ILUT */
    /* n       = integer. The row dimension of the matrix A. The matrix      */
    /* a,ja,ia = matrix stored in Compressed Sparse Row format.              */
    /* lfil    = integer. The fill-in parameter. Each row of L and each row
                 of U will have a maximum of lfil elements (excluding the 
		 diagonal element). lfil must be .ge. 0.
		 ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
		 EARLIER VERSIONS.                                           */
    ipar[0]=n+1;
    /* droptol = real*8. Sets the threshold for dropping small terms in the
                 factorization. See below for details on dropping strategy.  */
    fpar[0]=DROP_TOL;
    /* iwk     = integer. The lengths of arrays alu and jlu. If the arrays
                 are not big enough to store the ILU factorizations, ilut
		 will stop with an error message.                            */
    ipar[1]=ELBOW*nz;
    /* On return:
       ===========

       alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
                 the L and U factors together. The diagonal (stored in
		 alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
		 contains the i-th row of L (excluding the diagonal entry=1)
		 followed by the i-th row of U.                              */
    ilutmat.ja=(integer *)  MALLOC((size_t)ipar[1]*sizeof(integer),  "mainildlc:ilutmat.ja");
    ilutmat.a =(FLOAT *)MALLOC((size_t)ipar[1]*sizeof(FLOAT),"mainildlc:ilutmat.a");
    /* ju      = integer array of length n containing the pointers to
                 the beginning of each row of U in the matrix alu,jlu.       */
    ilutmat.ia=(integer *)  MALLOC((size_t)n*sizeof(integer),"mainildlc:ilutmat.ia");   
    /* ierr    = integer. Error message with the following meaning.
                 ierr  = 0    --> successful return.
		 ierr .gt. 0  --> zero pivot encountered at step number ierr.
		 ierr  = -1   --> Error. input matrix may be wrong.
                                  (The elimination process has generated a
				  row in L or U whose length is .gt.  n.)
                 ierr  = -2   --> The matrix L overflows the array alu.
		 ierr  = -3   --> The matrix U overflows the array alu.
		 ierr  = -4   --> Illegal value for lfil.
		 ierr  = -5   --> zero row encountered.                      */
    ipar[2]=0;
    /* work arrays:
       =============
       jw      = integer work array of length 2*n.
       w       = real work array of length n                                 */
    w =(FLOAT *)MALLOC(3*(size_t)n*sizeof(FLOAT),"mainildlc:w");
    ws=(integer *)  MALLOC(7*(size_t)n*sizeof(integer),"mainildlc:ws");

    // perform ILUT decomposition
    // DROP_INVERSE|TISMENETSKY_SC|NO_SHIFT|REPEAT_FACT|IMPROVED_ESTIMATE
    // param=DROP_INVERSE|TISMENETSKY_SC|IMPROVED_ESTIMATE|DIAGONAL_COMPENSATION;
    param=DROP_INVERSE|IMPROVED_ESTIMATE|ENSURE_SPD;

    printf("ILDLC PARAMETERS:\n");
    printf("   droptol=%8.1e\n",DROP_TOL);
    printf("   lfil   =%d\n",  ipar[0]);
    printf("   elbow space factor=%4d\n",  ELBOW);
    if (param & DROP_INVERSE)
       printf("   inverse-based dropping\n");
    else
       printf("   dual threshold dropping\n");
    if (param & TISMENETSKY_SC)
       printf("   Tismenetsky-like Schur complement\n");
    else
       printf("   simple Schur complement\n");
    if (param & NO_SHIFT)
       printf("   zero pivots are kept\n");
    else
       printf("   zero pivots are shifted away\n");
    if (param & REPEAT_FACT)
       printf("   second pass\n");
    else
       printf("   first pass\n");
    if (param & IMPROVED_ESTIMATE)
       printf("   improved estimate of the inverse\n");
    else
       printf("   simple estimate of the inverse\n");
    if (param & DIAGONAL_COMPENSATION)
       printf("   diagonal compensation\n");
    else
       printf("   diagonal entries unmodified\n");
    if (param & ENSURE_SPD)
       printf("   SPD property preserved\n");
    else
       printf("   no attempt to preserve SPD\n");
    fflush(stdout);

    evaluate_time(&time_start,&systime);
    ILDLC(&n,A.a,A.ja,A.ia,ipar,fpar,&param,
	 ilutmat.a,ilutmat.ja,ipar+1,w,ws,ipar+2);
    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
    fprintf(fo,"   |");
    free(w);
    free(ws);

    /*
    printf("diagonal entries\n");
    for (i=0; i<n; i++)
      printf("%8.1e",ilutmat.a[i]);
      // printf("%8.1e%8.1e",ilutmat.a[i].r,ilutmat.a[i].i);
    printf("\n");fflush(stdout);
    printf("U entries\n");
    for (i=0; i<n; i++) {
      printf("%8d:\n",i+1);
      for (j=ilutmat.ja[i]; j<ilutmat.ja[i+1]; j++)
	printf("%16d",ilutmat.ja[j-1]);
      printf("\n");
      for (j=ilutmat.ja[i]; j<ilutmat.ja[i+1]; j++)
	printf("%8.1e",ilutmat.a[j-1]);
        // printf("%8.1e%8.1e",ilutmat.a[j-1].r,ilutmat.a[j-1].i);
      printf("\n");
    }
    fflush(stdout);
    */

    nc=0;
    for (i=0; i<n; i++)
        nc+=ilutmat.ja[i+1]-ilutmat.ja[i];
    printf("\nfill-in U: %5d(%4.1fav.), totally %5d(%4.1fav.)\n",
	   nc,(1.0*nc)/n,nc+n,(1.0*(nc+n))/n);
    printf("fill-in factor: %5.1f\n",(1.0*(nc+n))/nz);
    fprintf(fo,"%5.1f|",(1.0*(tmp3+nc+n))/nz);
    printf("time: %8.1e [sec]\n\n",(double)secnds);
    fprintf(fo,"%7.1le|",(double)secnds+secndsprep);
    switch (ipar[2]) {
           case  0: /* perfect! */
	            break;
           case -1: /* Error. input matrix may be wrong.
                       (The elimination process has generated a
			row in L or U whose length is .gt.  n.) */
	            printf("Error. input matrix may be wrong\n");
		    break;
           case -2: /* The matrix L overflows the array alu */
	            printf("The matrix L overflows the array alu\n");
		    break;
           case -3: /* The matrix U overflows the array alu */
	            printf("The matrix U overflows the array alu\n");
		    break;
           case -4: /* Illegal value for lfil */
	            printf("Illegal value for lfil\n");
		    break;
           case -5: /* zero row encountered */
	            printf("zero row encountered\n");
		    break;
           default: /* zero pivot encountered at step number ierr */
	            printf("zero pivot encountered at step number %d\n",ipar[2]);
		    break;
    } // end switch


    if (ipar[2]!=0) {
	printf("Iterative solver(s) cannot be applied\n\n");
	fprintf(fo,"Iterative solver(s) cannot be applied\n");
	fflush(fo);
	exit(0);
    } // end if

    /* --------------------   apply preconditioner   --------------------- */
    // initial solution
    if (rhstyp[1]=='G' || rhstyp[1]=='g') {
       j=nrhs*n;
       if (rhstyp[0]=='M' || rhstyp[0]=='m')
	  j=n;
       for (i=0; i<n; i++,j++)
	   sol[i]=rhs[j];
    }
    else
       for (i=0; i<n; i++)
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   sol[i]=0;
#else
           sol[i].r=sol[i].i=0;
#endif
    
    // rescale initial solution
    for (i=0; i<n; i++) {  
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        sol[i]/=scale[i];
#else
	sol[i].r/=scale[i].r;
	sol[i].i/=scale[i].r;
#endif
    } // end for i


    // init solver
    ipar[0]=0;
    // right preconditioning
    ipar[1]=2;
    if (ipar[1]==0)
      printf("no preconditioning\n");
    else if (ipar[1]==1)
      printf("left preconditioning\n");
    else if (ipar[1]==2)
      printf("right preconditioning\n");
    else if (ipar[1]==3)
      printf("both sides preconditioning\n");
       
    // stopping criteria, energy norm
    ipar[2]=3;
    // number of restart length for GMRES,FOM,...
    ipar[4]=30;
    // maximum number of iteration steps
    ipar[5]=1.1*n+10;
    ipar[5]=MIN(ipar[5],MAX_IT);
    // init with zeros
    ipar[6]=ipar[7]=ipar[8]=ipar[9]=ipar[10]=ipar[11]=ipar[12]=ipar[13]=0;
    // relative error tolerance
    fpar[0]=1e-8;
    // absolute error tolerance
    fpar[1]=0;
    // init with zeros
    fpar[2]=fpar[3]=fpar[4]=fpar[5]=fpar[6]=fpar[7]=fpar[8]=fpar[9]=fpar[10]=0;
       
    // allocate memory for iterative solver depending on our choice */
    i=ipar[4];
    switch (SOLVER) {
    case  1:   // pcg
      i=5*n;
      printf("use PCG as iterative solver\n");
      break;
    case  2:   // cgnr 
      i=5*n;
      printf("use CGNR as iterative solver\n");
      break;
    case  3:   // bcg
      printf("use BCG as iterative solver\n");
      i=7*n;
      break;
    case  4:   // dbcg
      i=11*n;
      printf("use DBCG as iterative solver\n");
      break;
    case  5:   // bcgstab
      i=8*n;
      printf("use BCGSTAB as iterative solver\n");
      break;
    case  6:   // tfqmr
      i=11*n;
      printf("use TFQMR as iterative solver\n");
      break;
    case  7:   // fom
      i=(n+3)*(i+2)+(i+1)*i/2;
      printf("use FOM(%d) as iterative solver\n",ipar[4]);
      break;
    case  8:   // gmres
      i=(n+3)*(i+2)+(i+1)*i/2;
      printf("use GMRES(%d) as iterative solver\n",ipar[4]);
      break;
    case  9:   // fgmres
      i=2*n*(i+1)+(i+1)*i/2;
      printf("use FGMRES(%d) as iterative solver\n",ipar[4]);
      break;
    case 10:   // dqgmres
      i=n+(i+1)*(2*n+4);
      printf("use DQGMRES(%d) as iterative solver\n",ipar[4]);
      break;
    } // end switch
    w=(FLOAT *)MALLOC((size_t)i*sizeof(FLOAT),"mainildlc:w");
    ipar[3]=i;
       
    // start iterative solver
    evaluate_time(&time_start,&systime);
    flag=-1;
    while (flag) {
          // which solver do we use?
	  switch (SOLVER) {
	  case  1:   // pcg
	    PCG(&n,rhs,sol,ipar,fpar,w);
	    break;
	    /*
	     case  2:   // cgnr
	       cgnr(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  3:   // bcg
	       bcg(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  4:   // dbcg
	       dbcg(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  5:   // bcgstab
	       bcgstab(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  6:   // tfqmr
	       tfqmr(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  7:   // fom 
	       fom(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  8:   // gmres
	       GMRES(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case  9:   // fgmres
	       fgmres(&n,rhs,sol,ipar,fpar,w);
	       break;
	     case 10:   // dqgmres
	       dqgmres(&n,rhs,sol,ipar,fpar,w);
	       break;
	    */
	  } // end switch
	  // what is needed for the solver?
	  w--;
	  switch (ipar[0]) {
	  case  1:   // apply matrix vector multiplication
	    HERMATVEC(A,w+ipar[7],w+ipar[8]);
	    break;
	  case  2:   // apply transposed matrix vector multiplication
	    HERMATVEC(A,w+ipar[7],w+ipar[8]);
	    break;
	  case  3:   // (no) left preconditioner solver
	    i=1;
	    if (ipar[1]==1)
	       ILDLCSOL(&n, w+ipar[7],w+ipar[8], ilutmat.a,ilutmat.ja);
	    else
	       COPY(&n,w+ipar[7],&i,w+ipar[8],&i);
	    break;
	  case  4:   // (no) left preconditioner transposed solver
	    i=1;
	    if (ipar[1]==1)
	       ILDLCSOL(&n, w+ipar[7],w+ipar[8], ilutmat.a,ilutmat.ja);
	    else
	       COPY(&n,w+ipar[7],&i,w+ipar[8],&i);
	    break;
	  case  5:   // (no) right preconditioner solver
	    i=1;
	    if (ipar[1]==2)
	       ILDLCSOL(&n, w+ipar[7],w+ipar[8], ilutmat.a,ilutmat.ja);
	    else
	       COPY(&n,w+ipar[7],&i,w+ipar[8],&i);
	    break;
	  case  6:   // (no) right preconditioner transposed solver
	    i=1;
	    if (ipar[1]==2)
	       ILDLCSOL(&n, w+ipar[7],w+ipar[8], ilutmat.a,ilutmat.ja);
	    else
	       COPY(&n,w+ipar[7],&i,w+ipar[8],&i);
	    break;
	  case 10:   /* call my own stopping test routine */
	    break;
	  default:   /* the iterative solver terminated with code=ipar[0] */
	    flag=0;
	    break;
	  } // end switch
	  w++;
    } // end while

    // rescale solution
    for (i=0; i<n; i++) {  
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        sol[i]*=scale[i];
#else
	sol[i].r*=scale[i].r;
	sol[i].i*=scale[i].r;
#endif
    } // end for i
    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
    // why did the iterative solver stop?
    switch (ipar[0]) {
    case  0:  // everything is fine
      break;
    case -1:  // too many iterations
      printf("number of iteration steps exceeds its limit\n");
      break;
    case -2:  // not enough work space
      printf("not enough work space provided\n");
      break;
    case -3:  // not enough work space
      printf("algorithm breaks down ");
      // Arnoldi-type solvers
      if (SOLVER>=7) {
	switch (ipar[11]) {
	case -1:  
	  printf("due to zero input vector");
	  break;
	case -2:  
	  printf("since input vector contains abnormal numbers");
	  break;
	case -3:  
	  printf("since input vector is a linear combination of others");
	  break;
	case -4:  
	  printf("since triangular system in GMRES/FOM/etc. has null rank");
	  break;
	} // end switch
      } // end if
      printf("\n");
      break;
    default:  // unknown
      printf("solver exited with error code %d\n",ipar[0]);
    } // end switch

    printf("number of iteration steps: %d\n",ipar[6]);
    printf("time: %8.1le [sec]\n\n", (double)secnds);
    if (ipar[6]>=MAX_IT) {
       fprintf(fo,"  inf  |");
       fprintf(fo," inf\n");
    } // end if
    else 
       if (ipar[0]!=0) {
	  fprintf(fo,"  NaN  |"); 
	  fprintf(fo," NaN\n");
       } // end if 
       else {
	  fprintf(fo,"%7.1le|",(double)secnds);
	  fprintf(fo," %3d\n",ipar[6]);
       } // if-else

    fflush(fo);
    
    printf("residual norms:\n");
    printf("initial: %8.1e\n",fpar[2]);
    printf("target:  %8.1e\n",fpar[3]);
    printf("current: %8.1e\n",fpar[5]);

    if (nrhs==0) {
       for (i=0; i<n; i++) 
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   sol[i]-=1;
#else
	   sol[i].r-=1;
#endif
       i=1;
       if ((rhstyp[2]=='X' || rhstyp[2]=='x') || (rhstyp[0]!='M' && rhstyp[0]!='m'))
	  printf("rel. error in the solution: %8.1le\n\n",NRM(&n,sol,&i)/sqrt(1.0*n));
       else printf("\n");
    } // end if
    else 
       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	 j=nrhs*n;
	 if (rhstyp[0]=='M' || rhstyp[0]=='m')
	    j=n;
	 if (rhstyp[1]=='G' || rhstyp[1]=='g') 
	    j+=nrhs*n;
	 for (i=0; i<n; i++,j++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	     sol[i]-=rhs[j];
#else
	     sol[i].r-=rhs[j].r;
	     sol[i].i-=rhs[j].i;
#endif
	 } // end for i
	 j-=n;
	 i=1;
	 printf("rel. error in the solution: %8.1le\n\n",NRM(&n,sol,&i)/NRM(&n,rhs+j,&i));
       } // end if

    free(w);



    /* ------------------------------------------------------------------- */
    /* ------------------------   END MAIN part   ------------------------ */
    /* ------------------------------------------------------------------- */
} // end main

