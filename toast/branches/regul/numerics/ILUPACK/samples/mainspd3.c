#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <ilupack.h>

#include <ilupackmacros.h>




#define STDERR          stdout
#define STDOUT          stdout
//#define PRINT_INFO
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))




int main(int argc, char **argv)
{
    CSRMAT       A, ilutmat;
    FILE         *fp, *fo; 
    char         rhstyp[3],title[73],key[9],type[3], fname[100],foname[100],
                 extension[3];
    integer      i,j,k,l,m,fnamelen,n,nc,nz,nrhs,tmp0,tmp,tmp2,tmp3,ierr,
                 nrhsix, *rhsptr, *rhsind, nnzU, nAx=0, mynrhs=2, sumit;
    REALS        DROP_TOL, CONDEST,ELBOW, val,vb;
    FLOAT        *rhs,*sol, *scale, *rhsval;
    float        systime, time_start, time_stop, secnds, sumtime;
    AMGLEVELMAT  PRE, *next;
    ILUPACKPARAM param;
     
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

    ELBOW   =atof(argv[argc-2]);
    CONDEST =atof(argv[argc-3]);
    DROP_TOL=atof(argv[argc-4]);

    /* read in a symmetric positive definite matrix in Harwell-Boeing format. 
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


    // set parameters to the default settings
    SPDAMGINIT(&A, &param);
    
    // **********************************************************************
    // !!! new parameter selection scheme !!!
    // (the old getparams/setparams interface still works but it is much less
    //  convenient)
    // **********************************************************************

    // Here you can (if you like) change some of the default options
    
    //  1. reordering strategies
    // by default Approximate Minimum Degree by P. Amestoy, T. Davis and 
    // I. Duff is pre-selected. If you prefer a different ordering, change
    // the following field. Available are
    //   (a)  "amd"      (default)
    //   (b)  "metisn"   Metis Multilevel Nested Dissection by Nodes by 
    //                   G. Karypis and V. Kumar
    //   (c)  "metise"   Metis Multilevel Nested Dissection by Edges by 
    //                   G. Karypis and V. Kumar
    //   (d)  "rcm"      Reverse Cuthill-McKee from the SPARPAK package by
    //                   A. George and J. Liu
    //   (e)  "amf"      HALOAMD by P. Amestoy, T. Davis, I. Duff
    //   (d)  "mmd"      Minimum Degree from the SPARPAK package by
    //                   A. George and J. Liu
    //   (f)  "indset"   Independent set strategy from ARMS by Y. Saad
    //   (g)  "pq"       ddPQ strategy from ARMS by Y. Saad
    //   (h)  any other string (say "mydriver") that does not match  (a)-(g) 
    //        will force you to provide your own permutation drivers perm0,
    //        perm and permf when calling AMGfactor
    // param.ordering="metisn";

    //  2. drop tolerance for the LU factors
    // by default, 1e-2 is chosen. Here you can overwrite the default values
    // This sample program is designed to pass your own drop tolerance (called
    // (DROP_TOL) via the command line when being called
    param.droptol=DROP_TOL;

    //  3. drop tolerance for the approximate Schur complement
    // by default, 1e-2 is chosen. Here you can overwrite the default values
    // This sample program is designed to pass your own drop tolerance (called
    // (DROP_TOL) via the command line when being called. In this example, we
    // simply use the same drop tolerance as for the LU factors
    param.droptolS=DROP_TOL;

    //  4. norm of the inverse triangular factors
    // by default, 1e+2 is chosen. Here you can overwrite the default values.
    // This sample program is designed to pass your own inverse-bound (called
    // (CONDEST) via the command line when being called. 
    // As a rule of thumb, small CONDEST will allow more entries to be dropped
    // (which may accelerate the computation and save memory) but at the same 
    // time, more levels will be necessary (which in turn may slow down the 
    // computation and increase the memory). Typically, values  between 5 and 
    // 100 make sense. 
    // CONDEST=5 will make ILUPACK behave like AMG and select many coarse grid
    // nodes. If you have a PDE-based problem, this might be the right choice.
    // Otherwise, CONDEST=100 will safeguard the ILU computation and prevent
    // the norm of the inverse triangular factors from becoming too large.
    param.condest=CONDEST;

    //  5. residual tolerance for the iterative solver
    // The built-in iterative solver (here restarted GMRES) will use this
    // tolerance to terminate whenever the backward error is less than this
    // threshold. By default, sqrt(eps)~1e-8 is chosen for double precision, 
    // eps^(3/4)~1e-6 is chosen for single precision. 
    // param.restol=1e-12;

    //  6. maximum number of iteration steps
    // By default the iterative solver will quit once the number of iteration
    // steps exceeds MIN(1.1*n+10,500). 
    // param.maxit=1000;

    //  7. elbow space for the ILU
    // Here please pass an estimate how much memory you are willing to spend.
    // ILUPACK will try to keep the ILU inside the range you passed. The elbow
    // space is a real number measuring the number of nonzeros of the ILU 
    // relative to the fill of the original matrix. By default, 10 is chosen. 
    // Note however, if your estimate is too small, ILUPACK will adapt elbow 
    // and overwrite this parameter. As long as enough memory is available, the
    // ILU will be successfully computed.
    // This sample program is designed to pass your own elbow space estimate
    // (called ELBOW) via the command line when being called. 
    param.elbow=ELBOW;

    //  8. maximum number of nonzeros per column in L (resp. per row in U)
    // by default n+1 is chosen, i.e. this option is disabled. You can limit
    // the amount of memory by using some smaller value, e.g. A.ia[A.nc]-1
    // is the fill of A and ELBOW*(A.ia[A.nc]-1.0)/A.nc would restrict the
    // maximum number of fill to the average number of nonzeros of A per column
    // (or per row) times the ELBOW. Note however that this parameter cuts off
    // the fill in L and U by brute force. It recommended NOT to use it.
    // parameter.lfil=ELBOW*(A.ia[A.nc]-1.0)/A.nc;

    //  9. maximum number of nonzeros per row in S (approximate Schur 
    //     complement)
    // by default n+1 is chosen, i.e. this option is disabled. You can limit
    // the amount of memory by using some smaller value, e.g. A.ia[A.nc]-1
    // is the fill of A and ELBOW*(A.ia[A.nc]-1.0)/A.nc would restrict the
    // maximum number of fill to the average number of nonzeros of A per column
    // (or per row) times the ELBOW. Note however that this parameter cuts off
    // the fill in S by brute force. It strongly recommended NOT to use it.
    // parameter.lfilS=ELBOW*(A.ia[A.nc]-1.0)/A.nc;

    // 10. type of test vector
    // for some PDE-based problems it might be sensible to ensure that the ILU
    // is exact when being applied to some given test vector. By default this
    // option is disabled ("none"). If you want to use this feature you can
    // either use "static" to pass a fixed test vector to the ILU. Or you can
    // use any other string. In this case, using reverse communication
    // principle, on every level you need to pass a test vector to the ILU. The
    // ILU passes to you the current coarse grid system and an initial guess for
    // the test vector. You have to return your own test vector. On entry to the
    // first level, this initial guess is simply the test vector you prescribed.
    // On any subsequent level, this will be your old test vector restricted to 
    // the coarse grid.
    //param.typetv="static";

    // 11. test vector
    // If you decide to pass a test vector, then pass the associated pointer.
    // ILUPACK will make its own copy inside AMGfactor, and you can release the
    // memory if you like. In PDE-based applications, a typical guess is the 
    // vector with all ones.
    //param.flags&=~COARSE_REDUCE;
    //mytestvector=(FLOAT *)MALLOC(A.nc*sizeof(FLOAT),"mainsym");
    //for (i=0; i<A.nc; i++) mytestvector[i]=1.0;
    //param.tv=mytestvector;

    // 12. type of multigrid
    //   (a) "ilu"   (default) multilevel ILU
    //   (b) "amli"  on each coarse grid, an inner iteration is used based on
    //               flexible GMRES to solve the inner coarse grid system, 
    //               preconditioned by the associated ILU. Note that this 
    //               requires to maintain all coarse grid systems and increases
    //               the amount of memory.
    //   (c) "mg"    full multigrid with pre- and post smoothing, V-cycle,
    //               W-cycle or flexible cycle is chosen. Essentially, the
    //               multilevel ILU is used to define interpolation and 
    //               restriction operators as well as the coarse grid systems,
    //               while the other components are set up as in the usual
    //               multigrid framework. Not that the flexible cycle does not
    //               pre-select the number of coarse grid solves a priori (e.g.
    //               (1 for V-cycle, 2 for W-cycle), but on on each coarse grid,
    //               an inner iteration is used based on flexible GMRES to solve
    //               the inner coarse grid system, preconditioned by the 
    //               associated full multigrid solver. Note that this type of
    //               multigrid preconditioning requires to maintain all coarse
    //               grid systems and increases the amount of memory.
    // param.amg="amli";

    // 13. number of pre-smoothing steps
    // If classical multigrid is selected (param.amg="mg";), then here you can
    // set the number of pre-smoothing steps. default: 1
    // param.npresmoothing=1;

    // 14. number of coarse grid solves
    // Except for multilevel ILU (i.e. param.amg="amli"; or param.amg="mg";),
    // here you define how often the coarse grid solve is performed. By default,
    // only one coarse grid solve is used (V-cycle). The choice param.ncoarse=2;
    // would correspond to a W-cycle. Note however, if a negative value is
    // passed, a flexible solver is invoked, i.e. the number of coarse grid
    // solves varies from one grid to another and from one step to the next one.
    // param.ncoarse=-1;

    // 15. type of pre-smoother
    // if full multigrid is used (param.amg="mg";), then here you can choose 
    // between built-in smoothers or your own hand-made smoother.
    //   (a) "gsf"     (default) Gauss-Seidel forward
    //   (b) "gsb"     Gauss-Seidel backward
    //   (c) "j"       (damped) Jacobi
    //   (d) "ilu"     ILU on the fine grid system
    //   (e) any other string that does not match (a)-(d) will cause AMGsolver
    //       to use reverse communication principle in order to let you provide
    //       your own smoother. In that case ILUPACK will give you the matrix,
    //       the right hand side and an initial solution (typically 0). You
    //       have to override the initial solution 
    // param.presmoother="gsf";

    // 16. pre-selection of coarse grid nodes
    // In some PDE-based applications it might be useful to select some coarse
    // grid nodes in advance. Essentially this strategy uses a Ruge-Stueben-like
    // heuristic strategy. If a test vector is available, the coarsening 
    // strategy is applied to the matrix, which is diagonally scaled from the 
    // right with the test vector. 
    //  (a)  "none"  (default)  leave the coarsening process to ILUPACK, 
    //               inverse-based strategy will construc a coarse grid on its 
    //               own.
    //  (b)  "yes"   Some nodes are pre-selected as coarse grid nodes, ILUPACK
    //               might add some further nodes. 
    // param.FCpart="none";
    
    // 17. type of coarse grid system
    // By default the coarse grid system S is computed from A and the ILU in
    //                                 /L11   0\ /D11   0\ /U11   U12\ 
    // typical ILU manner, i.e. if A ~ |       |*|       |*|         |, then S 
    //                                 \L21   I/ \ 0    S/ \ 0     I /
    // is defined via S:= A22-L21*D11*U12.
    // Alternatively one could compute W21 ~ L21*L11^{-1}, Z12 ~ U11^{-1}*U12 
    //                                  /-Z12\
    // and define S via S:= [-W21  I]*A*|    |. This would refer to an AMG-like
    //                                  \  I /
    // strategy to compute a coarse grid system.
    // available are
    //    (a) "ilu"    (default)  ILU-type coarse grid system
    //    (b) "amg"    AMG-type coarse grid system
    // param.typecoarse="ilu";

    // indicate that we want to compute more than one ILU
    param.flags|=RE_FACTOR;


    // print some messages that give information about flags and reorderings
#include "spdmessages.c"

    evaluate_time(&time_start,&systime);
    ierr=SPDAMGFACTOR(&A, &PRE, &param);
    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
   

    switch (ierr) {
           case  0: /* perfect! */
	        break;
           case -1: /* Error. input matrix may be wrong.
                       (The elimination process has generated a
			row in L or U whose length is .gt.  n.) */
	        printf("Error. input matrix may be wrong at level %d\n",
		       PRE.nlev);
		fprintf(fo,"input matrix may be wrong\n");
		break;
           case -2: /* The matrix L overflows the array alu */
	        printf("The matrix L overflows the array alu at level %d\n",
		       PRE.nlev);
		fprintf(fo,"out of memory\n");
		break;
           case -3: /* The matrix U overflows the array alu */
	        printf("The matrix U overflows the array alu at level %d\n",
		       PRE.nlev);
		fprintf(fo,"out of memory\n");
		break;
           case -4: /* Illegal value for lfil */
	        printf("Illegal value for lfil at level %d\n",PRE.nlev);
		fprintf(fo,"Illegal value for lfil\n");
		break;
           case -5: /* zero row encountered */
	        printf("zero row encountered at level %d\n",PRE.nlev);
		fprintf(fo,"zero row encountered\n"); 
		break;
           case -6: /* zero column encountered */
	        printf("zero column encountered at level %d\n",PRE.nlev);
		fprintf(fo,"zero column encountered\n");
		break;
           case -7: /* buffers too small */
	        printf("buffers are too small\n");
		fprintf(fo,"buffers are too small\n");
		break;
           default: /* zero pivot encountered at step number ierr */
	        printf("zero pivot encountered at step number %d of level %d\n",
		       ierr,PRE.nlev);
		fprintf(fo,"zero pivot encountered\n");
    } /* end switch */
    if (ierr) {
       fprintf(fo,"Iterative solver(s) cannot be applied\n");
       fflush(fo);
       exit(ierr);
    }

    // print some statistics about the levels, their size and the 
    // computation time
#include "spdprintperformance.c"



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
#include "spdinitvectors.c"
         
        evaluate_time(&time_start,&systime);
	ierr=SPDAMGSOLVER(&A, &PRE, &param, rhs+A.nr*l, sol+A.nr*l);
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
#include "spdfinalres.c"
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
    extension[0]=fname[++i]; fname[i]='s';
    extension[1]=fname[++i]; fname[i]='o';
    extension[2]=fname[++i]; fname[i]='l';
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

    // outputfile name, rewrite ".rsa,..." extension
    i=fnamelen;
    j=-1;
    while (j) {
          i--;
	  if (i<0)
	    j=0;
	  if (fname[i]=='.')
	    j=0;
    }
    fname[++i]=extension[0];
    fname[++i]=extension[1];
    fname[++i]=extension[2];
    i++;
    j=i;
    while (i<100)
      fname[i++]='\0';


    // finally release memory of the preconditioner
    SPDAMGDELETE(&A,&PRE,&param);

    // release memory used for the input matrix
    free(A.ia);
    free(A.ja);
    free(A.a );
    free(rhs );
    free(sol );



    // ------------------------------------------------------------------------
    // read in another matrix, here: same matrix for the second time
    // note that we keep the RE_FACTOR flag to indicate that we want
    // to solve a third system
    // ------------------------------------------------------------------------

    
#include "spdreadmatrix.c"
    // if right hand sides are provided, then run AMGSOLVER for any of these
    // right hand sides. Otherwise use own set of right hand sides

    // allocate memory for the solution vector and some buffer
    sol  =(FLOAT *) MALLOC(mynrhs*(size_t)n*sizeof(FLOAT), "mainspd:sol");


    // print some messages that give information about flags and reorderings
#include "spdmessages.c"

    evaluate_time(&time_start,&systime);
    ierr=SPDAMGFACTOR(&A, &PRE, &param);
    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
   

    switch (ierr) {
           case  0: /* perfect! */
	        break;
           case -1: /* Error. input matrix may be wrong.
                       (The elimination process has generated a
			row in L or U whose length is .gt.  n.) */
	        printf("Error. input matrix may be wrong at level %d\n",
		       PRE.nlev);
		fprintf(fo,"input matrix may be wrong\n");
		break;
           case -2: /* The matrix L overflows the array alu */
	        printf("The matrix L overflows the array alu at level %d\n",
		       PRE.nlev);
		fprintf(fo,"out of memory\n");
		break;
           case -3: /* The matrix U overflows the array alu */
	        printf("The matrix U overflows the array alu at level %d\n",
		       PRE.nlev);
		fprintf(fo,"out of memory\n");
		break;
           case -4: /* Illegal value for lfil */
	        printf("Illegal value for lfil at level %d\n",PRE.nlev);
		fprintf(fo,"Illegal value for lfil\n"); 
		break;
           case -5: /* zero row encountered */
	        printf("zero row encountered at level %d\n",PRE.nlev);
		fprintf(fo,"zero row encountered\n"); 
		break;
           case -6: /* zero column encountered */
	        printf("zero column encountered at level %d\n",PRE.nlev);
		fprintf(fo,"zero column encountered\n");
		break;
           case -7: /* buffers too small */
	        printf("buffers are too small\n");
		fprintf(fo,"buffers are too small\n");
		break;
           default: /* zero pivot encountered at step number ierr */
	        printf("zero pivot encountered at step number %d of level %d\n",
		       ierr,PRE.nlev);
		fprintf(fo,"zero pivot encountered\n");
    } /* end switch */
    if (ierr) {
       fprintf(fo,"Iterative solver(s) cannot be applied\n");
       fflush(fo);
       exit(ierr);
    }

    // print some statistics about the levels, their size and the 
    // computation time
#include "spdprintperformance.c"



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
#include "spdinitvectors.c"
         
        evaluate_time(&time_start,&systime);
	ierr=SPDAMGSOLVER(&A, &PRE, &param, rhs+A.nr*l, sol+A.nr*l);
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
#include "spdfinalres.c"
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
    extension[0]=fname[++i]; fname[i]='s';
    extension[1]=fname[++i]; fname[i]='o';
    extension[2]=fname[++i]; fname[i]='l';
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

    // outputfile name, rewrite ".rsa,..." extension
    i=fnamelen;
    j=-1;
    while (j) {
          i--;
	  if (i<0)
	    j=0;
	  if (fname[i]=='.')
	    j=0;
    }
    fname[++i]=extension[0];
    fname[++i]=extension[1];
    fname[++i]=extension[2];
    i++;
    j=i;
    while (i<100)
      fname[i++]='\0';


    // finally release memory of the preconditioner
    SPDAMGDELETE(&A,&PRE,&param);

    // release memory used for the input matrix
    free(A.ia);
    free(A.ja);
    free(A.a );
    free(rhs );
    free(sol );






    // ------------------------------------------------------------------------
    // read in another matrix, here: same matrix for the third time
    // NOW note that we CLEAR the RE_FACTOR flag to indicate that this is our
    // final system.
    // ------------------------------------------------------------------------

    

#include "spdreadmatrix.c"
    // if right hand sides are provided, then run AMGSOLVER for any of these
    // right hand sides. Otherwise use own set of right hand sides

    // allocate memory for the solution vector and some buffer
    sol  =(FLOAT *) MALLOC(mynrhs*(size_t)n*sizeof(FLOAT), "mainspd:sol");


    // NOW indicate that we DO NOT want to compute more than one ILU
    param.flags&=~RE_FACTOR;


    // print some messages that give information about flags and reorderings
#include "spdmessages.c"

    evaluate_time(&time_start,&systime);
    ierr=SPDAMGFACTOR(&A, &PRE, &param);
    evaluate_time(&time_stop,&systime);
    secnds=time_stop-time_start;
   

    switch (ierr) {
           case  0: /* perfect! */
	        break;
           case -1: /* Error. input matrix may be wrong.
                       (The elimination process has generated a
			row in L or U whose length is .gt.  n.) */
	        printf("Error. input matrix may be wrong at level %d\n",
		       PRE.nlev);
		fprintf(fo,"input matrix may be wrong\n");
		break;
           case -2: /* The matrix L overflows the array alu */
	        printf("The matrix L overflows the array alu at level %d\n",
		       PRE.nlev);
		fprintf(fo,"out of memory\n");
		break;
           case -3: /* The matrix U overflows the array alu */
	        printf("The matrix U overflows the array alu at level %d\n",
		       PRE.nlev);
		fprintf(fo,"out of memory\n");
		break;
           case -4: /* Illegal value for lfil */
	        printf("Illegal value for lfil at level %d\n",PRE.nlev);
		fprintf(fo,"Illegal value for lfil\n"); 
		break;
           case -5: /* zero row encountered */
	        printf("zero row encountered at level %d\n",PRE.nlev);
		fprintf(fo,"zero row encountered\n"); 
		break;
           case -6: /* zero column encountered */
	        printf("zero column encountered at level %d\n",PRE.nlev);
		fprintf(fo,"zero column encountered\n"); 
		break;
           case -7: /* buffers too small */
	        printf("buffers are too small\n");
		fprintf(fo,"buffers are too small\n");
		break;
           default: /* zero pivot encountered at step number ierr */
	        printf("zero pivot encountered at step number %d of level %d\n",
		       ierr,PRE.nlev);
		fprintf(fo,"zero pivot encountered\n");
    } /* end switch */
    if (ierr) {
       fprintf(fo,"Iterative solver(s) cannot be applied\n");
       fflush(fo);
       exit(ierr);
    }

    // print some statistics about the levels, their size and the 
    // computation time
#include "spdprintperformance.c"



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
#include "spdinitvectors.c"
         
        evaluate_time(&time_start,&systime);
	ierr=SPDAMGSOLVER(&A, &PRE, &param, rhs+A.nr*l, sol+A.nr*l);
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
#include "spdfinalres.c"
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
    extension[0]=fname[++i]; fname[i]='s';
    extension[1]=fname[++i]; fname[i]='o';
    extension[2]=fname[++i]; fname[i]='l';
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

    // outputfile name, rewrite ".rsa,..." extension
    i=fnamelen;
    j=-1;
    while (j) {
          i--;
	  if (i<0)
	    j=0;
	  if (fname[i]=='.')
	    j=0;
    }
    fname[++i]=extension[0];
    fname[++i]=extension[1];
    fname[++i]=extension[2];
    i++;
    j=i;
    while (i<100)
      fname[i++]='\0';


    // finally release memory of the preconditioner
    SPDAMGDELETE(&A,&PRE,&param);

    // release memory used for the input matrix
    free(A.ia);
    free(A.ja);
    free(A.a );
    free(rhs );
    free(sol );

} /* end spdmain */








