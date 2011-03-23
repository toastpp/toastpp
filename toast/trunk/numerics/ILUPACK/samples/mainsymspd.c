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





// reorder the system according to the symmetric minimum degree ordering
//#define MINIMUM_DEGREE

// reorder the system according to the nested dissection ordering
//#define NESTED_DISSECTION 

// reorder the system according to the reverse Cuthill-McKee ordering
//#define REVERSE_CM

// reorder the system according to some independent set ordering
//#define IND_SET

// reorder system using approximate minimum fill by Patrick Amestoy, 
// Tim Davis and Iain Duff
//#define AMF

// reorder system using approximate minimum degree by Patrick Amestoy
// Tim Davis and Iain Duff
//#define AMD

// reorder the columns and rows of a system differently using a new unsymmetric
// reordering strategy by Yousef Saad. This is the adapted symmetric version
// of it
//#define PP_PERM


// mixed strategies that finally switch to PP if necessary
//#define MMD_PP
//#define AMF_PP
//#define AMD_PP
//#define RCM_PP
//#define FC_PP
//#define METIS_E_PP
//#define METIS_N_PP

// alternative Symmetric Maximum Weight Matching interface provided by PARDISO
//#define SMC64_RCM_PP
//#define SMC64_MMD_PP
//#define SMC64_AMF_PP
#define SMC64_AMD_PP
//#define SMC64_METIS_E_PP
//#define SMC64_METIS_N_PP

//#define SMWM_RCM_PP
//#define SMWM_MMD_PP
//#define SMWM_AMF_PP
//#define SMWM_AMD_PP
//#define SMWM_METIS_E_PP
//#define SMWM_METIS_N_PP

//#define FC_SMWM_RCM
//#define FC_SMWM_AMD
//#define FC_SMWM_AMF
//#define FC_SMWM_MMD
//#define FC_SMWM_METIS_N
//#define FC_SMWM_METIS_E



// variant of PILUC that uses a repeated multiple factorization approach
//#define USE_MPILUC









int main(int argc, char **argv)
{
    /* SOLVER choice:   1  pcg
                        2  sbcg
                        3  bcg
                        4  sqmr
                        8  gmres
                        9  fgmres */
    integer SOLVER=4; /* sqmr */

    CSRMAT A, ilutmat;
    AMGLEVELMAT PRE, *next;
    integer nlev=0, nprev, nB;
    integer (*perm0)(),(*perm)(),(*permf)();

    FILE *fp, *fo; 
    char rhstyp[3], title[73], key[9], type[3], fname[100], foname[100];
    char line[MAX_LINE], *tmpstring, timeout[7], residout[7];
    integer  i,j,k,m,fnamelen,n,nc,nz,nrhs,mynrhs=2,tmp0,tmp,tmp2,tmp3,ierr,flag,
         *invperm, *buff, *ws, *symmmd,flags,elbow,max_it,ipar[20],
         nrhsix, *rhsptr, *rhsind, sumit;
    REALS eps, DROP_TOL, CONDEST,condest,droptols[2],restol, val,vb;
    FLOAT *rhs,*sol,*w, *scale, *rhsval, *dbuff;
    float  systime, time_start,   time_stop,   secnds, secndsprep,
           timeAx_start, timeAx_stop, secndsAx, sumtime;
    integer ELBOW, nnzU, l, nAx;
    ILUPACKPARAM param;
    size_t  nibuff, ndbuff;
     
    /* the last argument passed serves as file name */
    if (argc!=5) {
      printf("usage '%s <drop tol.> <bound for L^{-1},U^{-1}> <elbow space> <matrix>'\n",argv[0]);
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
    CONDEST=atof(argv[argc-3]);
    DROP_TOL=atof(argv[argc-4]);

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

    
    // if right hand sides are provided, then run AMGSOLVER for any of these
    // right hand sides. Otherwise use own set of right hand sides

    // allocate memory for the solution vector and some buffer
    sol  =(FLOAT *) MALLOC(mynrhs*(size_t)n*sizeof(FLOAT), "mainspd:sol");
    dbuff=(FLOAT *) MALLOC(3*(size_t)n*sizeof(FLOAT),      "mainspd:dbuff");
    ndbuff=3*(size_t)n;


    // set parameters to the default settings
    MYSYMAMGINIT(A, &param);
    param.dbuff=dbuff;
    param.ndbuff=ndbuff;
    
    // select reordering functions
    perm0=SYMPERMNULL;
    perm =SYMPERMNULL;
    permf=SYMPERMNULL;
#ifdef MINIMUM_DEGREE
    perm0=SYMPERMMMD;
    perm =SYMPERMMMD;
    permf=SYMPERMMMD;
    fprintf(fo,"mmds/mmds|");
#elif defined REVERSE_CM
    perm0=SYMPERMRCM;
    perm =SYMPERMRCM;
    permf=SYMPERMRCM;
    fprintf(fo,"rcms/rcms|");
#elif defined NESTED_DISSECTION
    perm0=SYMPERMND;
    perm =SYMPERMND;
    permf=SYMPERMND;
    fprintf(fo,"nds /nds |");
#elif defined IND_SET
    perm0=SYMPERMINDSET;
    perm =SYMPERMINDSET;
    permf=SYMPERMINDSET;
    fprintf(fo,"inds/inds|");
#elif defined AMF
    perm0=SYMPERMAMF;
    perm =SYMPERMAMF;
    permf=SYMPERMAMF;
    fprintf(fo,"amfs/amfs|");
#elif defined AMD
    perm0=SYMPERMAMD;
    perm =SYMPERMAMD;
    permf=SYMPERMAMD;
    fprintf(fo,"amds/amds|");
#elif defined PP_PERM
    perm0=SPDPERMPP;
    perm =SPDPERMPP;
    permf=SPDPERMPP;
    fprintf(fo,"PPs /PPs |");
#elif defined MMD_PP
    perm0=SYMPERMMMD;
    perm =SYMPERMMMD;
    permf=SPDPERMPP;
    fprintf(fo,"mmds/PPs |");
#elif defined AMF_PP
    perm0=SYMPERMAMF;
    perm =SYMPERMAMF;
    permf=SPDPERMPP;
    fprintf(fo,"amfs/PPs |");
#elif defined AMD_PP
    perm0=SYMPERMAMD;
    perm =SYMPERMAMD;
    permf=SPDPERMPP;
    fprintf(fo,"amds/PPs |");
#elif defined RCM_PP
    perm0=SYMPERMRCM;
    perm =SYMPERMRCM;
    permf=SPDPERMPP;
    fprintf(fo,"rcms/PPs |");
#elif defined FC_PP
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=SPDPERMPP;
    fprintf(fo,"FCs /PPs |");
#elif defined METIS_E_PP
    perm0=SYMPERMMETISE;
    perm =SYMPERMMETISE;
    permf=SPDPERMPP;
    fprintf(fo,"mes /PQs |");
#elif defined METIS_N_PP
    perm0=SYMPERMMETISN;
    perm =SYMPERMMETISN;
    permf=SPDPERMPP;
    fprintf(fo,"mns /PPs |");
#elif defined SMWM_MMD_PP
    perm0=PERMSMWMMMD;
    perm =PERMSMWMMMD;
    permf=SPDPERMPP;
    fprintf(fo,"mwmd/PPs |");
#elif defined SMWM_AMF_PP
    perm0=PERMSMWMAMF;
    perm =PERMSMWMAMF;
    permf=SPDPERMPP;
    fprintf(fo,"mwaf/PPs |");
#elif defined SMWM_AMD_PP
    perm0=PERMSMWMAMD;
    perm =PERMSMWMAMD;
    permf=SPDPERMPP;
    fprintf(fo,"mwad/PPs |");
#elif defined SMWM_RCM_PP
    perm0=PERMSMWMRCM;
    perm =PERMSMWMRCM;
    permf=SPDPERMPP;
    fprintf(fo,"mwrc/PPs |");
#elif defined SMWM_METIS_E_PP
    perm0=PERMSMWMMETISE;
    perm =PERMSMWMMETISE;
    permf=SPDPERMPP;
    fprintf(fo,"mwme/PQs |");
#elif defined SMWM_METIS_N_PP
    perm0=PERMSMWMMETISN;
    perm =PERMSMWMMETISN;
    permf=SPDPERMPP;
    fprintf(fo,"mwmn/PPs |");
#elif defined SMC64_MMD_PP
    perm0=PERMSMC64MMD;
    perm =PERMSMC64MMD;
    permf=SPDPERMPP;
    fprintf(fo,"mc64md/PP|");
#elif defined SMC64_AMF_PP
    perm0=PERMSMC64AMF;
    perm =PERMSMC64AMF;
    permf=SPDPERMPP;
    fprintf(fo,"mc64af/PP|");
#elif defined SMC64_AMD_PP
    perm0=PERMSMC64AMD;
    perm =PERMSMC64AMD;
    permf=SPDPERMPP;
    fprintf(fo,"mc64ad/PP|");
#elif defined SMC64_RCM_PP
    perm0=PERMSMC64RCM;
    perm =PERMSMC64RCM;
    permf=SPDPERMPP;
    fprintf(fo,"mc64rc/PP|");
#elif defined SMC64_METIS_E_PP
    perm0=PERMSMC64METISE;
    perm =PERMSMC64METISE;
    permf=SPDPERMPP;
    fprintf(fo,"mc64me/PP|");
#elif defined SMC64_METIS_N_PP
    perm0=PERMSMC64METISN;
    perm =PERMSMC64METISN;
    permf=SPDPERMPP;
    fprintf(fo,"mc64mn/PP|");
#elif defined FC_SMWM_METIS_E
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=PERMSMWMMETISE;
    fprintf(fo,"fc/mwme|");
#elif defined FC_SMWM_METIS_N
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=PERMSMWMMETISN;
    fprintf(fo,"fc/mwmn|");
#elif defined FC_SMWM_AMD
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=PERMSMWMAMD;
    fprintf(fo,"fc/mwad|");
#elif defined FC_SMWM_AMF
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=PERMSMWMAMF;
    fprintf(fo,"fc/mwaf|");
#elif defined FC_SMWM_MMD
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=PERMSMWMMMD;
    fprintf(fo,"fc/mwmd|");
#elif defined FC_SMWM_RCM
    perm0=SYMPERMFC;
    perm =SYMPERMFC;
    permf=PERMSMWMRCM;
    fprintf(fo,"fc/mwrc|");
#else
    fprintf(fo,"    /    |");
#endif



    // modify the default settings
    MYSYMAMGGETPARAMS(param, &flags, &elbow, droptols, &condest,
		      &restol, &max_it);

    // ONLY for mixed reordering strategies it is useful to
    // set the 'FINAL_PIVOTING' flag
#if !defined RCM_PP && !defined AMF_PP && !defined AMD_PP && !defined MMD_PP && !defined FC_PP && !defined METIS_E_PP && !defined METIS_N_PP && !defined FC_SMWM_METIS_N && !defined FC_SMWM_METIS_E && !defined FC_SMWM_AMD  && !defined FC_SMWM_AMF && !defined FC_SMWM_MMD && !defined FC_SMWM_RCM 
    flags&=~FINAL_PIVOTING;
#endif
    // change flags if mpiluc is desired
#ifdef USE_MPILUC
    flags|=MULTI_PILUC;
#endif

    // overwrite the default drop tolerances
    droptols[1]=DROP_TOL; 

    // choose the iterative solver, default: 8 (GMRES)
    param.ipar[5]=SOLVER;

    // overwrite the default value for elbow space
    elbow=ELBOW;

    // overwrite default values for condest
    condest=CONDEST;

    // overwrite the default value for the residual norm
    // restol=1e-8;

    // turn of scaling
    // param.ipar[7]&=~(3+24+192);
    // param.ipar[8]&=~((1+32+1024)*4);
    // param.ipar[8]|=((1+32+1024)*5);

    // use Standard Schur-complement
    //flags|=TISMENETSKY_SC;
    //flags&=~SIMPLE_SC;

    // turn off aggressive dropping
    flags&=~AGGRESSIVE_DROPPING;

    // do not discard matrix entries on an early stage
    // flags&=~DISCARD_MATRIX;

    // rewrite the updated parameters
    MYSYMAMGSETPARAMS(A, &param, flags, elbow, droptols, condest,
		      restol, max_it);

    // manually change the drop tolerance for the approximate Schur complement
    // param.fpar[8]=DROP_TOL;


    // print some messages that give information about flags and reorderings
#include "symmessages.c"
       
    evaluate_time(&time_start,&systime);
    ierr=MYSYMAMGFACTOR(&A, &PRE, &nlev, &param, perm0,perm, permf);


    // convert matrix to postive definite preconditioner
    SYMSPDCONVERT(&PRE, nlev);

    // update buffer size
    nibuff=param.nibuff;
    ndbuff=param.ndbuff;
    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
   

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
} /* end mainsym */








