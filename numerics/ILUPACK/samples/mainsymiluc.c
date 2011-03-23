#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <ilupack.h>

#include <ilupackmacros.h>

#include "symswitches.c"

#define MAX_LINE        255
#define STDERR          stdout
#define STDOUT          stdout
//#define PRINT_INFO
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))

// maximum number of iterations independent on n
#define MAX_IT          500

// measure for terminating iterative solver
#define RESTOL_FUNC(A)  sqrt(A)
//#define RESTOL_FUNC(A)  (A)


// all the following reorderings are use in combination with symmetric
// maximum weight matching obtained by MC64 

// reorder the system according to the symmetric minimum degree ordering
//#define SMC64_MMD

// reorder the system according to MeTiS multilevel nested dissection ordering
//#define SMC64_METIS_E
//#define SMC64_METIS_N

// reorder the system according to the reverse Cuthill-McKee ordering
//#define SMC64_RCM

// reorder system using approximate minimum fill by Patrick Amestoy, 
// Tim Davis and Iain Duff
//#define SMC64_AMF


// reorder system using approximate minimum degree by Patrick Amestoy
// Tim Davis and Iain Duff
#define SMC64_AMD


// all the following reorderings are use in combination with symmetric
// maximum weight matching obtained by MPS as part of PARDISO

// reorder the system according to the symmetric minimum degree ordering
//#define SMWM_MMD

// reorder the system according to MeTiS multilevel nested dissection ordering
//#define SMWM_METIS_E
//#define SMWM_METIS_N

// reorder the system according to the reverse Cuthill-McKee ordering
//#define SMWM_RCM

// reorder system using approximate minimum fill by Patrick Amestoy, 
// Tim Davis and Iain Duff
//#define SMWM_AMF


// reorder system using approximate minimum degree by Patrick Amestoy
// Tim Davis and Iain Duff
#define SMWM_AMD



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
    integer SOLVER=4; /* gmres */

    CSRMAT A, symilumat;
    AMGLEVELMAT PRE, *next;

    FILE *fp, *fo; 
    char rhstyp[3], title[72], key[8], type[3], fname[100], foname[100];
    char line[MAX_LINE], *tmpstring, timeout[7], residout[7];
    integer  i,j,k,l,m,fnamelen,n,nc,nz,nrhs,tmp0,tmp,tmp2,tmp3,ierr,ipar[16],flag,
      *invperm, *buff, *ws, *symmmd, nlev=1, nnzU,
         nrhsix, *rhsptr, *rhsind, mynrhs=2;
    REALS eps,fpar[16], DROP_TOL, droptols[3], restol, condest, CONDEST, val,vb;
    FLOAT *rhs,*sol,*w, *dbuff, *scale, *rhsval;
    float  systime, time_start,   time_stop,   secnds, secndsprep,
           timeAx_start, timeAx_stop, secndsAx, sumtime;
    integer ELBOW, *p, *invq, nB, flags, iwk, elbow,  max_it, *ibuff, sumit;
    ILUPACKPARAM param;
    size_t nibuff=0, ndbuff=0;

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


    /* read in a symmetric matrix in Harwell-Boeing format. 
       By definition Harwell-Boeing format stores a matrix in compressed sparse
       COLUMN format. In the symmetric case this is equivalent to compressed 
       sparse ROW format which ILUPACK uses. 
       Note that only the upper triangular part is stored by rows.
                                 / 3.5 -1.0  0 \
       A symm. pos. def. matrix  |-1.0  2.0  0 | is stored as follows
                                 \  0    0  1.5/ 

       A.ia:   1  3  4  5           pointer to the start of every compressed 
                                    upper triangular row plus pointer to the
				    first space behind the compressed rows

       A.ja:   1    2    2    3     nonzero column indices

       A.a:   3.5 -1.0  2.0  1.5    nonzero numerical values

       The read part finally yields the following data structures
        -  A:  matrix in compressed sparse row format
	    o  A.nr, A.nc: number of rows and columns of A
            o  A.ia:  pointer array
            o  A.ja:  nonzero column index array 
	    o  A.a:   nonzero numerical values
	-  rhs:  right hand side(s) and additional data like exact solution
	         or initial guess
	-  n:  same as A.nr,A.nc
	-  nz:  number of nonzero entries
     */
#include "spdreadmatrix.c"

    
    // allocate memory for the solution vector and some buffer
    sol  =(FLOAT *) MALLOC(mynrhs*(size_t)n*sizeof(FLOAT), "mainspd:sol");


    // set parameters to the default settings
    MYSYMAMGINIT(A, &param);

    
    evaluate_time(&time_start,&systime);
    p   =(integer *)MALLOC((size_t)n*sizeof(integer),"mainsymiluc:p");
    invq=(integer *)MALLOC((size_t)n*sizeof(integer),"mainsymiluc:invq");
    scale=(FLOAT *)MALLOC((size_t)n*sizeof(FLOAT),"mainsymiluc:scale");

    for (i=0; i<n; i++) {
        p[i]=invq[i]=i+1;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
        scale[i]=1.0;
#else
        scale[i].r=1.0;
        scale[i].i=0.0;
#endif
    } // end for i
    ierr=0;
    nB=n;



#if defined SMC64_AMD
    ierr=PERMSMC64AMD(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mc64ad/  |");
    printf("prescribe symmetric MC64 + approximate minimum degree ordering\n");
#elif defined SMC64_AMF
    ierr=PERMSMC64AMF(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mc64af/  |");
    printf("prescribe symmetric MC64 + approximate minimum fill ordering\n");
#elif defined SMC64_MMD
    ierr=PERMSMC64MMD(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mc64md/  |");
    printf("prescribe symmetric MC64 + minimum degree ordering\n");
#elif defined SMC64_METIS_E
    ierr=PERMSMC64METISE(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mc64me/  |");
    printf("prescribe symmetric MC64 + MeTiS multilevel nested dissection (by edges)\n");
#elif defined SMC64_METIS_N
    ierr=PERMSMC64METISN(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mc64mn/  |");
    printf("prescribe symmetric MC64 + MeTiS multilevel nested dissection (by nodes)\n");
#elif defined SMC64_RCM
    ierr=PERMSMC64RCM(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mc64rc/  |");
    printf("prescribe symmetric MC64 + reverse Cuthill-McKee ordering\n");
#elif defined SMWM_AMD
    ierr=PERMSMWMAMD(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mwad/    |");
    printf("prescribe symmetric Max. Weight Match. + approximate minimum degree ordering\n");
#elif defined SMWM_AMF
    ierr=PERMSMWMAMF(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mwaf/    |");
    printf("prescribe symmetric Max. Weight Match. + approximate minimum fill ordering\n");
#elif defined SMWM_MMD
    ierr=PERMSMWMMMD(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mwmd/    |");
    printf("prescribe symmetric Max. Weight Match. + minimum degree ordering\n");
#elif defined SMWM_METIS_E
    ierr=PERMSMWMMETISE(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mwme/    |");
    printf("prescribe symmetric Max. Weight Match. + MeTiS multilevel nested dissection (by edges)\n");
#elif defined SMWM_METIS_N
    ierr=PERMSMWMMETISN(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mwmn/    |");
    printf("prescribe symmetric Max. Weight Match. + MeTiS multilevel nested dissection (by nodes)\n");
#elif defined SMWM_RCM
    ierr=PERMSMWMRCMRCM(A,scale,scale,p,invq,&nB,&param);
    fprintf(fo,"mwrc/    |");
    printf("prescribe symmetric Max. Weight Match. + reverse Cuthill-McKee ordering\n");
#else
    fprintf(fo,"    /    |");
#endif


    evaluate_time(&time_stop,&systime);
    secndsprep=time_stop-time_start;
    printf("time preprocessing: %8.1e [sec]\n",(double)secndsprep);
    ILUPACK_secnds[7]=ILUPACK_secnds[8]=ILUPACK_secnds[0]=secndsprep;



    /* ------------------------------------------------------------------- */
    /* -----------------------   Start MAIN part   ----------------------- */
    /* ------------------------------------------------------------------- */

    fprintf(fo,"IB-SYMILUC|");



    MYSYMAMGGETPARAMS(param, &flags, &elbow, droptols, &condest,
		      &restol, &max_it);
    //flags|=DROP_INVERSE;

    // turn off aggressive dropping
    flags&=~AGGRESSIVE_DROPPING;

    // overwrite the default value for the residual norm
    // restol=1e-8;

    // rewrite the updated parameters
    MYSYMAMGSETPARAMS(A, &param, flags, elbow, droptols, condest,
		      restol, max_it);

    // analyze how the size of the buffer has to be chosen
    if (flags&TISMENETSKY_SC) {
        if (flags&MULTI_PILUC)
	   nibuff=MAX(nibuff,17*(size_t)n);
	nibuff=MAX(nibuff,16*(size_t)n);
    }
    else {
      if (flags&MULTI_PILUC)
          nibuff=MAX(nibuff,14*(size_t)n);
      nibuff=MAX(nibuff,13*(size_t)n); // recommended version!
    }

    if (flags&DROP_INVERSE) {
      // if a more reliable estimate for the inverse is required
      if (flags&IMPROVED_ESTIMATE || !(flags&(SIMPLE_SC|TISMENETSKY_SC)))
        ndbuff=MAX(ndbuff,4*(size_t)n); // recommended version!
      else 
        ndbuff=MAX(ndbuff,3*(size_t)n);
    }
    else {
      // currently the norm of the inverse factors are computed
      // even if classical dropping is desired
      // The comments in (m)piluc still refer to their predecessor
      // iluc, which allows classical dropping in a strict sense
      
      // if a more reliable estimate for the inverse is required
      if (flags&IMPROVED_ESTIMATE || !(flags&(SIMPLE_SC|TISMENETSKY_SC)))
        ndbuff=MAX(ndbuff,4*(size_t)n);
      else 
        ndbuff=MAX(ndbuff,3*(size_t)n);
    }

    dbuff=(FLOAT *)MALLOC((size_t)ndbuff*sizeof(FLOAT),"mainsymiluc:dbuff");
    ibuff=(integer *)  MALLOC((size_t)nibuff*sizeof(integer),  "mainsymiluc:ibuff");

    
    iwk=ELBOW*nz;
    symilumat.ja=(integer *)  MALLOC((size_t)iwk*sizeof(integer),  "mainsymiluc:symilumat.ja");
    symilumat.a =(FLOAT *)MALLOC((size_t)iwk*sizeof(FLOAT),"mainsymiluc:symilumat.a");


    printf("SYMILUC PARAMETERS:\n");
    printf("   droptol=%8.1e\n",DROP_TOL);
    printf("   elbow space factor=%4d\n",  ELBOW);
    if (flags & DROP_INVERSE)
       printf("   inverse-based dropping\n");
    else
       printf("   dual threshold dropping\n");
    if (flags & TISMENETSKY_SC)
       printf("   Tismenetsky-like Schur complement\n");
    else
       printf("   simple Schur complement\n");
    if (flags & NO_SHIFT)
       printf("   zero pivots are kept\n");
    else
       printf("   zero pivots are shifted away\n");
    if (flags & IMPROVED_ESTIMATE)
       printf("   improved estimate of the inverse\n");
    else
       printf("   simple estimate of the inverse\n");
    if (flags & DIAGONAL_COMPENSATION)
       printf("   diagonal compensation\n");
    else
       printf("   diagonal entries unmodified\n");
    fflush(stdout);

    evaluate_time(&time_start,&systime);
    MYSYMILUC(&n,A.a,A.ja,A.ia,&n,&DROP_TOL,&flags,
	      p,invq,symilumat.a,symilumat.ja,&iwk,
	      dbuff,ibuff,&ierr);

    // dummy preconditioner in order to be able to apply ILUPACK routines
    PRE.next=NULL;
    PRE.prev=NULL;

    PRE.rowscal=PRE.colscal=scale;
    PRE.p=p;
    PRE.invq=invq;

    PRE.n=n;
    PRE.nB=n;

    PRE.LU.nr=PRE.LU.nc=n;
    PRE.LU.ia=NULL;
    PRE.LU.ja=symilumat.ja;
    PRE.LU.a =symilumat.a;
    PRE.LUperm=NULL;
    PRE.E.nr=PRE.E.nc=0;
    PRE.F.nr=PRE.F.nc=0;

    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
   
    ILUPACK_secnds[3]=secnds;
    ILUPACK_secnds[7]+=secnds;
    ILUPACK_secnds[8]+=secnds;
    secnds+=secndsprep;

    switch (ierr)
    {
           case  0: /* perfect! */
	            break;
           case -1: /* Error. input matrix may be wrong.
                       (The elimination process has generated a
			row in L or U whose length is .gt.  n.) */
	            printf("Error. input matrix may be wrong at level %d\n",nlev);
		    fprintf(fo,"input matrix may be wrong\n");fclose(fo); 
		    break;
           case -2: /* The matrix L overflows the array alu */
	            printf("The matrix L overflows the array alu at level %d\n",nlev);
		    fprintf(fo,"out of memory\n");fclose(fo); 
		    break;
           case -3: /* The matrix U overflows the array alu */
	            printf("The matrix U overflows the array alu at level %d\n",nlev);
		    fprintf(fo,"out of memory\n");fclose(fo); 
		    break;
           case -4: /* Illegal value for lfil */
	            printf("Illegal value for lfil at level %d\n",nlev);
		    fprintf(fo,"Illegal value for lfil\n");fclose(fo); 
		    break;
           case -5: /* zero row encountered */
	            printf("zero row encountered at level %d\n",nlev);
		    fprintf(fo,"zero row encountered\n");fclose(fo); 
		    break;
           case -6: /* zero column encountered */
	            printf("zero column encountered at level %d\n",nlev);
		    fprintf(fo,"zero column encountered\n");fclose(fo); 
		    break;
           case -7: /* buffers too small */
	            printf("buffers are too small\n");
		    // This error message would not be necessary if AMGsetup is
		    // called with the correct values of nibuff, ndbuff
		    printf("increase buffer size to at least %ld (float), %ld (integer)\n",
			   ndbuff, nibuff);
		    fflush(stdout);
		    fprintf(fo,"buffers are too small\n");fclose(fo); 
           default: /* zero pivot encountered at step number ierr */
	            printf("zero pivot encountered at step number %d of level %d\n",ierr,nlev);
		    fprintf(fo,"zero pivot encountered\n");fclose(fo); 
		    break;
    } /* end switch */
    if (ierr) {
       fprintf(fo,"Iterative solver(s) cannot be applied\n");
       fflush(fo);
       exit(ierr);
    }


#ifdef PRINT_INFO
    next=&PRE;
    for (i=1; i<=nlev; i++) {
        // fill-in LU
        printf("level %3d, block size %7d\n",i,next->LU.nr); fflush(stdout);
	if (i<nlev || next->LU.ja!=NULL) {
	  printf("U-factor");
	  printf("\n");fflush(stdout);
	  for (l=0; l<next->LU.nr; ) {
	    if (next->LU.ja[next->LU.nr+1+l]==0){
	      for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		printf("%8d",next->LU.ja[j-1]);
	      }
	      printf("\n");fflush(stdout);
	      for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		printf("%8.1le",next->LU.a[j-1]);
	      }
	      l++;
	    }
	    else {
	      for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		printf("%8d",next->LU.ja[j-1]);
	      }
	      printf("\n");fflush(stdout);
	      for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		printf("%8.1le",next->LU.a[next->LU.ja[l]+2*(j-next->LU.ja[l])-1]);
	      }
	      printf("\n");fflush(stdout);
	      for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		printf("%8.1le",next->LU.a[next->LU.ja[l]+2*(j-next->LU.ja[l])]);
	      }
	      l+=2;
	    }
	    printf("\n");fflush(stdout);
	  }

	  printf("Block diagonal factor\n");
	  for (l=0; l<next->LU.nr;) {
	    if (next->LU.ja[next->LU.nr+1+l]==0){
	      printf("%8.1le",next->LU.a[l]);
	      l++;
	    }
	    else {
	      printf("%8.1le%8.1le",next->LU.a[l],next->LU.a[next->LU.nr+1+l]);
	      l+=2;
	    }
	  }
	  printf("\n");fflush(stdout);
	  for (l=0; l<next->LU.nr; ) {
	    if (next->LU.ja[next->LU.nr+1+l]==0) {
	      printf("        ");
	      l++;
	    }
	    else {
	      printf("%8.1le%8.1le",next->LU.a[next->LU.nr+1+l],next->LU.a[l+1]);
	      l+=2;
	    }
	  }
	  printf("\n");fflush(stdout);



	}
	if (i==nlev) {
	   if (next->LU.ja==NULL) {
	      printf("switched to full matrix processing\n");fflush(STDOUT);
	   }
	}

	if (i<nlev) {
	   // fill-in F
	   nnzU+=next->F.ia[next->F.nr]-1;
	   printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,next->LU.nr,next->F.nc);
	   printf("  local fill-in F %7d(%5.1lfav pr)\n",
		  next->F.ia[next->F.nr]-1,(1.0*(next->F.ia[next->F.nr]-1))/next->LU.nr);
	}
	next=next->next;
    }
#endif



    // print some statistics about the levels, their size and the 
    // computation time
#include "symprintperformance.c"



    /* some decisions about the right hand side, the exact solution and the 
       initial guess

       - if a right hand side is provided, then solve the system according
         to this given right hand side.
	 Otherwise set x=(1,...,1)^T and use b=Ax as right hand side
	 --> rhs

       - if an initial guess is provided, then use this as starting vector.
         Otherwise use x=0 as initial guess
	 --> sol
       
       - for statistics also compute the norm of the initial residual --> val
    */
    sumtime=0.0;
    sumit=0;
    for (l=0; l<mynrhs; l++) {
#include "syminitvectors.c"


        evaluate_time(&time_start,&systime);
	ierr=MYSYMAMGSOLVER(&A, &PRE, nlev, &param, rhs+A.nr*l, sol+A.nr*l);
	evaluate_time(&time_stop,&systime);
	secnds=time_stop-time_start;

	// why did the iterative solver stop?
	switch (ierr) {
	case  0:  // everything is fine
	      sumtime+=ILUPACK_secnds[5];
	      sumit+=param.ipar[26];
	      if (l==mynrhs-1) {
		 sumtime/=mynrhs;
		 sumit   =sumit/((double)mynrhs)+.5;
		 fprintf(fo,"%7.1le|",(double)sumtime);
		 fprintf(fo," %3d|",sumit);
		 fprintf(fo,"%7.1le\n",(double)sumtime+ILUPACK_secnds[7]);
	      }
	      break;
	case -1:  // too many iterations
	      printf("number of iteration steps exceeds its limit\n");
	      fprintf(fo,"  inf  |");
	      fprintf(fo," inf|");
	      fprintf(fo,"  inf  \n");
	      break;
	case -2:  /* not enough work space */
              printf("not enough work space provided\n");
	      break;
	case -3:  /* not enough work space */
	      printf("algorithm breaks down ");
	      printf("\n");
	      fprintf(fo,"  NaN  |");
	      fprintf(fo," NaN|");
	      fprintf(fo,"  NaN  \n");
	      break;
	default:  /* unknown */
              printf("solver exited with error code %d\n",ierr);
	} // end switch 
	fflush(fo);

	// stop if necessary
	if (ierr)
	   mynrhs=l;

	printf("%3d. right hand side\n",l+1);
	printf("number of iteration steps: %d\n",param.ipar[26]);
	printf("time: %8.1le [sec]\n",  (double)secnds);
	printf("      %8.1le [sec]\n\n",(double)ILUPACK_secnds[5]); 
       
	printf("time matrix-vector multiplication: %8.1le [sec]\n",ILUPACK_secnds[6]);

	printf("residual norms:\n");
	printf("initial: %8.1le\n",val);
	printf("target:  %8.1le\n",param.fpar[23]);


	/* some final statistics 
	   - about the true current residual of the computed solution and
	   - the relative error in the solution though the exact solution is known
	*/
#include "symfinalres.c"
    } // end for l
    fclose(fo);


    // outputfile name, add ".sol" extension
    i=fnamelen;
    j=-1;
    while (j) {
          i--;
	  if (i<0)
	    j=0;
	  if (fname[i]=='.')
	    j=0;
    }
    fname[++i]='s';
    fname[++i]='o';
    fname[++i]='l';
    i++;
    j=i;
    while (i<100)
      fname[i++]='\0';

    WRITEVECTORS(fname,sol,&(A.nr),&mynrhs,
		 "solution vectors computed by ILUPACK                                    ",
		 "sol     ",
		 j, 72, 8);

    /*
      // sample code for rereading the vectors from C
      READVECTORS(fname,sol,&(A.nr),&mynrhs,title,key,j, 72, 8);
      title[72]='\0';
      key[8]='\0';
    */

    
    // finally release memory of the preconditioner
    MYSYMAMGDELETE(A,PRE,nlev,&param);

    // release memory used for the input matrix
    free(A.ia);
    free(A.ja);
    free(A.a );
    free(rhs );
    free(dbuff);
    free(ibuff);
} // end mainsymiluc

