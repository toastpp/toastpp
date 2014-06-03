#if !defined _DOUBLE_REAL_ &&   !defined _SINGLE_REAL_
#define _COMPLEX_SYMMETRIC_
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <ilupack.h>

#include <ilupackmacros.h>


#include "symswitches.c"



#define MAX_LINE         100
#define STDERR          stdout
#define STDOUT          stdout
//#define PRINT_INFO
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))



void laplacematrix2D(integer m, integer n, integer *ia, integer *ja, FLOAT *a, integer *nz)
{
  /* Local variables */
  integer y,ni,i,j,k,ElementsPut, nn[5];
  // 1/h^2
  REALS bb=((REALS) m)*((REALS) m);

  ElementsPut = 0;
  for (j=0; j<m; ++j) {
      for (i=0; i<m; ++i) {
	  y=i + j*m;
	  nn[0]=y;

	  if (i<m-1)
	    nn[1]=(i+1)%m+j*m;
	  else
	    nn[1]=0;
	  if (i>0) 
	    nn[2]=(m+i-1)%m+j*m;
	  else
	    nn[2]=0;
	  
	  if (j<m-1)
	    nn[3]=i+((j+1)%m)*m;
	  else
	    nn[3]=0;
	  if (j>0) 
	    nn[4]=i+((m+j-1)%m)*m;
	  else
	    nn[4]=0;

	  ++ElementsPut;
	  ia[y]=ElementsPut;
	  ja[ElementsPut-1]=nn[0]+1;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  a[ElementsPut-1]=4.0*bb;
#else
	  a[ElementsPut-1].r=4.0*bb;
	  a[ElementsPut-1].i=0.0;
#endif

	  for (ni=1; ni<=4; ++ni) {
	    if (nn[ni]>y) {
	      ++ElementsPut;
	      ja[ElementsPut-1]=nn[ni]+1;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	      a[ElementsPut-1]=-1.0*bb;
#else
	      a[ElementsPut-1].r=-1.0*bb;
	      a[ElementsPut-1].i= 0.0;
#endif
	    }
	  }
      }
  }
  ia[n]=ElementsPut+1;
  *nz = ElementsPut;
} /* laplacematrix2D */




void laplacematrix3D(integer m, integer n, integer *ia, integer *ja, FLOAT *a, integer *nz)
{
  /* Local variables */
  integer y,ni,i,j,k,ElementsPut, nn[7];
  // 1/h^2
  REALS bb=((REALS) m)*((REALS) m);

  ElementsPut = 0;
  for (k=0; k<=m-1; ++k) {
      for (j=0; j<=m-1; ++j) {
	    for (i=0; i<=m-1; ++i) {
		y=i + j*m + k*m*m;

		nn[0]=y;

		if (i<m-1)
		   nn[1]=(i+1)%m+j*m+k*m*m;
		else
		   nn[1]=0;
		if (i>0) 
		   nn[2]=(m+i-1)%m+j*m+k*m*m;
		else
		   nn[2]=0;

		if (j<m-1)
		   nn[3]=i+((j+1)%m)*m+k*m*m;
		else
		   nn[3]=0;
		if (j>0) 
		   nn[4]=i+((m+j-1)%m)*m+k*m*m;
		else
		   nn[4]=0;

		if (k<m-1)
		   nn[5]=i+j*m+((k+1)%m)*m*m;
		else
		   nn[5]=0;
		if (k>0) 
		   nn[6]=i+j*m+((m+k-1)%m)*m*m;
		else
		   nn[6]=0;

		++ElementsPut;
		ia[y]=ElementsPut;
		ja[ElementsPut-1]=nn[0]+1;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		a[ElementsPut-1]=6.0*bb;
#else
		a[ElementsPut-1].r=6.0*bb;
		a[ElementsPut-1].i=0.0;
#endif

		for (ni=1; ni<=6; ++ni) {
		    if (nn[ni]>y) {
		       ++ElementsPut;
		       ja[ElementsPut-1]=nn[ni]+1;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
		       a[ElementsPut-1]=-1.0*bb;
#else
		       a[ElementsPut-1].r=-1.0*bb;
		       a[ElementsPut-1].i= 0.0;
#endif
		    }
		}
	    }
      }
  }
  ia[n]=ElementsPut+1;
  *nz = ElementsPut;
} /* laplacematrix3D */






int main(int argc, char **argv)
{
    CSRMAT       A, ilutmat;
    FILE         *fp, *fo; 
    char         rhstyp[3], line[MAX_LINE];
    integer      i,j,k,l,m,n,nz,nrhs=0,mynrhs=2,tmp0,tmp,tmp2,tmp3,
                 ierr, sumit, dim, wave, nnzU, nAx;
    REALS        DROP_TOL, CONDEST, ELBOW, val,vb, shiftr,shifti;
    FLOAT        *rhs,*sol, *scale;
    float        systime, time_start, time_stop, secnds, sumtime;
    AMGLEVELMAT  PRE, *next;
    ILUPACKPARAM param;
     

    fp=fopen("helmholtz.dat","r");
    if (fp==NULL) {
       printf("file `helmholtz.dat' not found\n");
       exit(0);
    }
    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    DROP_TOL=atof(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    CONDEST=atof(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    ELBOW=atof(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    dim=atoi(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    m=atoi(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    wave=atoi(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    shiftr=atof(line);

    fgets(line, MAX_LINE, fp);
    fgets(line, MAX_LINE, fp);
    shifti=atof(line);

    fclose(fp);

    printf("settings:\n");
    printf("drop tolerance       %12.4e\n",DROP_TOL);
    printf("bound for ||L^{-1}|| %12.4e\n",CONDEST);
    printf("elbow space          %8d\n",ELBOW);
    if (dim==3) {
       printf("3D problem\n");
       printf("n=m^3, where n=%10d, m=%3d\n",m*m*m,m);
    }
    else {
       printf("2D problem\n");
       printf("n=m^2, where n=%10d, m=%3d\n",m*m,m);
    }
    printf("wave length          %8d\n",wave);
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
    printf("shift  %12.4e\n",shiftr);
#else
    printf("shift  %12.4e+i*%12.4e\n",shiftr,shifti);
#endif

    if (dim==3) {
       A.nc=A.nr=n=m*m*m;
       nz=4*A.nr;
       A.ia=(integer *)MALLOC((n+1)*sizeof(integer),"helholtz:ia");
       A.ja=(integer *)MALLOC(nz*   sizeof(integer),"helholtz:ja");
       A.a =(FLOAT *)  MALLOC(nz*   sizeof(FLOAT),  "helholtz:a");
       laplacematrix3D(m,n,A.ia,A.ja,A.a,&nz);

       ilutmat.nc=ilutmat.nr=n=m*m*m;
       nz=4*ilutmat.nr;
       ilutmat.ia=(integer *)MALLOC((n+1)*sizeof(integer),"helholtz:ia");
       ilutmat.ja=(integer *)MALLOC(nz*   sizeof(integer),"helholtz:ja");
       ilutmat.a =(FLOAT *)  MALLOC(nz*   sizeof(FLOAT),  "helholtz:a");
       laplacematrix3D(m,n,ilutmat.ia,ilutmat.ja,ilutmat.a,&nz);
    }
    else {
       A.nc=A.nr=n=m*m;
       nz=3*A.nr;
       A.ia=(integer *)MALLOC((n+1)*sizeof(integer),"helholtz:ia");
       A.ja=(integer *)MALLOC(nz*   sizeof(integer),"helholtz:ja");
       A.a =(FLOAT *)  MALLOC(nz*   sizeof(FLOAT),  "helholtz:a");
       laplacematrix2D(m,n,A.ia,A.ja,A.a,&nz);

       ilutmat.nc=ilutmat.nr=n=m*m;
       nz=3*ilutmat.nr;
       ilutmat.ia=(integer *)MALLOC((n+1)*sizeof(integer),"helholtz:ia");
       ilutmat.ja=(integer *)MALLOC(nz*   sizeof(integer),"helholtz:ja");
       ilutmat.a =(FLOAT *)  MALLOC(nz*   sizeof(FLOAT),  "helholtz:a");
       laplacematrix2D(m,n,ilutmat.ia,ilutmat.ja,ilutmat.a,&nz);
    }


    // shift operator and preconditioner
    for (i=0; i<n; i++) {
        for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	    k=A.ja[j]-1;
	    if (k==i) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	       A.a[j]-=((REALS)wave)*((REALS)wave);
	       ilutmat.a[j]+=shiftr*((REALS)wave)*((REALS)wave);
#else
	       A.a[j].r-=((REALS)wave)*((REALS)wave);
	       ilutmat.a[j].r+=shiftr*((REALS)wave)*((REALS)wave);
	       ilutmat.a[j].i+=shifti*((REALS)wave)*((REALS)wave);
#endif
	    }
	}
    }


    /*
    for (i=0; i<n; i++) {
      printf("%8d\n",i+1);
      for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	printf("%8d",A.ja[j]);
      }
      printf("\n");
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
      for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	printf("%8.1e",A.a[j]);
      }
      printf("\n");
#else
      for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	printf("%8.1e",ilutmat.a[j].r);
      }
      printf("\n");
      for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	printf("%8.1e",ilutmat.a[j].i);
      }
      printf("\n");
#endif
      fflush(stdout);
    }
    */

    // if right hand sides are provided, then run AMGSOLVER for any of these
    // right hand sides. Otherwise use own set of right hand sides

    // allocate memory for the solution vector and some buffer
    rhs  =(FLOAT *) MALLOC(mynrhs*(size_t)n*sizeof(FLOAT), "main:rhs");
    sol  =(FLOAT *) MALLOC(mynrhs*(size_t)n*sizeof(FLOAT), "main:sol");


    // set parameters to the default settings
    MYSYMAMGINIT(&A, &param);

    
    fo = fopen("out_normal","aw");
    fprintf(fo,"helmholtz");
    if (dim==3)
       fprintf(fo,"_3D    |");
    else
       fprintf(fo,"_2D    |");
    fprintf(fo,"%7.1e|",DROP_TOL);
    fprintf(fo,"%7.1e|",CONDEST);

    // **********************************************************************
    // !!! new parameter selection scheme !!!
    // (the old getparams/setparams interface still works but it is much less
    //  convenient)
    // **********************************************************************

    // Here you can (if you like) change some of the default options
    
    //  1. maximum weight matching 
    // by default matching is turned on and the drivers are properly 
    // pre-selected. If you do not like to have matching, then simply turn it
    // off
    // param.matching=1;

    //  2. reordering strategies
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

    //  3. drop tolerance for the LU factors
    // by default, 1e-2 is chosen. Here you can overwrite the default values
    // This sample program is designed to pass your own drop tolerance (called
    // (DROP_TOL) via the command line when being called
    param.droptol=DROP_TOL;

    //  4. drop tolerance for the approximate Schur complement
    // by default, 1e-2 is chosen. Here you can overwrite the default values
    // This sample program is designed to pass your own drop tolerance (called
    // (DROP_TOL) via the command line when being called. In this example, we
    // simply use the same drop tolerance as for the LU factors
    param.droptolS=DROP_TOL;

    //  5. norm of the inverse triangular factors
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

    //  6. residual tolerance for the iterative solver
    // The built-in iterative solver (here restarted GMRES) will use this
    // tolerance to terminate whenever the backward error is less than this
    // threshold. By default, sqrt(eps)~1e-8 is chosen for double precision, 
    // eps^(3/4)~1e-6 is chosen for single precision. 
    // param.restol=1e-8;

    //  7. maximum number of iteration steps
    // By default the iterative solver will quit once the number of iteration
    // steps exceeds MIN(1.1*n+10,500). 
    // param.maxit=1000;

    //  8. elbow space for the ILU
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

    //  9. maximum number of nonzeros per column in L (resp. per row in U)
    // by default n+1 is chosen, i.e. this option is disabled. You can limit
    // the amount of memory by using some smaller value, e.g. A.ia[A.nc]-1
    // is the fill of A and ELBOW*(A.ia[A.nc]-1.0)/A.nc would restrict the
    // maximum number of fill to the average number of nonzeros of A per column
    // (or per row) times the ELBOW. Note however that this parameter cuts off
    // the fill in L and U by brute force. It recommended NOT to use it.
    // parameter.lfil=ELBOW*(A.ia[A.nc]-1.0)/A.nc;

    // 10. maximum number of nonzeros per row in S (approximate Schur 
    //     complement)
    // by default n+1 is chosen, i.e. this option is disabled. You can limit
    // the amount of memory by using some smaller value, e.g. A.ia[A.nc]-1
    // is the fill of A and ELBOW*(A.ia[A.nc]-1.0)/A.nc would restrict the
    // maximum number of fill to the average number of nonzeros of A per column
    // (or per row) times the ELBOW. Note however that this parameter cuts off
    // the fill in S by brute force. It strongly recommended NOT to use it.
    // parameter.lfilS=ELBOW*(A.ia[A.nc]-1.0)/A.nc;

    // 11. type of test vector
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
    // param.typetv="static";

    // 12. test vector
    // If you decide to pass a test vector, then pass the associated pointer.
    // ILUPACK will make its own copy inside AMGfactor, and you can release the
    // memory if you like. In PDE-based applications, a typical guess is the 
    // vector with all ones.
    //param.flags&=~COARSE_REDUCE;
    //mytestvector=(FLOAT *)MALLOC(A.nc*sizeof(FLOAT),"mainsym");
    //for (i=0; i<A.nc; i++) mytestvector[i]=1.0;
    //param.tv=mytestvector;

    // 13. type of multigrid
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

    // 14. number of pre-smoothing steps
    // If classical multigrid is selected (param.amg="mg";), then here you can
    // set the number of pre-smoothing steps. default: 1
    // param.npresmoothing=1;

    // 15. number of coarse grid solves
    // Except for multilevel ILU (i.e. param.amg="amli"; or param.amg="mg";),
    // here you define how often the coarse grid solve is performed. By default,
    // only one coarse grid solve is used (V-cycle). The choice param.ncoarse=2;
    // would correspond to a W-cycle. Note however, if a negative value is
    // passed, a flexible solver is invoked, i.e. the number of coarse grid
    // solves varies from one grid to another and from one step to the next one.
    // param.ncoarse=-1;

    // 16. type of pre-smoother
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

    // 17. pre-selection of coarse grid nodes
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
    
    // 18. type of coarse grid system
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

    // print some messages that give information about flags and reorderings
#include "symmessages.c"
       
    evaluate_time(&time_start,&systime);
    ierr=MYSYMAMGFACTOR(&ilutmat, &PRE, &param);


    for (i=0; i<A.nc; i++) {
        for (j=A.ia[i]-1; j<A.ia[i+1]-1; j++) {
	    k=A.ja[j]-1;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    A.a[j]*=PRE.colscal[i]*PRE.colscal[k];
#else
	    A.a[j].r*=PRE.colscal[i].r*PRE.colscal[k].r;
	    A.a[j].i*=PRE.colscal[i].r*PRE.colscal[k].r;
#endif
	}
    }
    

    // update buffer size
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
    rhstyp[0]=' ';
    rhstyp[1]=' ';
    rhstyp[2]=' ';
    for (l=0; l<mynrhs; l++) {
        // remember that scaling is explicitly applied to A
        scale=PRE.colscal;

	// prescribe artificial solution
	for (i=0; i<n; i++)  {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    sol[i+A.nr*l]=1.0+i*l;
#else
	    sol[i+A.nr*l].r=1.0+i*(l-1);
	    sol[i+A.nr*l].i=0.1*(A.nr-i+l-1);
#endif
	} // end for i

	// construct associated right hand side rhs=A*sol
	// remember that A is already rescaled!!!
	// rescale solution
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	    sol[i+A.nr*l]/=scale[i];
#else
	    val=1.0/(scale[i].r*scale[i].r + scale[i].i*scale[i].i);
	    vb=sol[i+A.nr*l].r;
	    sol[i+A.nr*l].r=( vb*scale[i].r+sol[i+A.nr*l].i*scale[i].i)*val;
	    sol[i+A.nr*l].i=(-vb*scale[i].i+sol[i+A.nr*l].i*scale[i].r)*val;
#endif
	} // end for i
	MYSMATVEC(A,sol+A.nr*l,rhs+A.nr*l);

	// rescale right hand side and solution back to their original
	// counterparts associated with the input matrix
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	    sol[i+A.nr*l]*=scale[i];
	    rhs[i+A.nr*l]/=scale[i];
#else
	    vb=sol[i+A.nr*l].r;
	    sol[i+A.nr*l].r=vb*scale[i].r-sol[i+A.nr*l].i*scale[i].i;
	    sol[i+A.nr*l].i=vb*scale[i].i+sol[i+A.nr*l].i*scale[i].r;
	    val=1.0/(scale[i].r*scale[i].r + scale[i].i*scale[i].i);
	    vb=rhs[i+A.nr*l].r;
#ifdef _COMPLEX_SYMMETRIC_
	    // transposed scaling
	    rhs[i+A.nr*l].r=( vb*scale[i].r+rhs[i+A.nr*l].i*scale[i].i)*val;
	    rhs[i+A.nr*l].i=(-vb*scale[i].i+rhs[i+A.nr*l].i*scale[i].r)*val;
#else
	    // conjugate transposed scaling
	    rhs[i+A.nr*l].r=( vb*scale[i].r-rhs[i+A.nr*l].i*scale[i].i)*val;
	    rhs[i+A.nr*l].i=( vb*scale[i].i+rhs[i+A.nr*l].i*scale[i].r)*val;
#endif
#endif
	} // end for i
	// norm of the initial residual
	i=1;
	val=NRM(&n, rhs+A.nr*l, &i);

	// initial solution
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    sol[i+A.nr*l]=0;
#else
	    sol[i+A.nr*l].r=sol[i+A.nr*l].i=0;
#endif
	} 


        evaluate_time(&time_start,&systime);
	ierr=MYSYMAMGSOLVER(&A, &PRE, &param, rhs+A.nr*l, sol+A.nr*l);
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

	// -------   compute final residual   ------
	// compute final residual
	// remember that A is already rescaled!!!
	// rescale solution
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	    sol[i+A.nr*l]/=scale[i];
#else
	    val=1.0/(scale[i].r*scale[i].r + scale[i].i*scale[i].i);
	    vb=sol[i+A.nr*l].r;
	    sol[i+A.nr*l].r=( vb*scale[i].r+sol[i+A.nr*l].i*scale[i].i)*val;
	    sol[i+A.nr*l].i=(-vb*scale[i].i+sol[i+A.nr*l].i*scale[i].r)*val;
#endif
	} // end for i
	MYSMATVEC(A,sol+A.nr*l,param.dbuff);
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	    sol[i+A.nr*l]*=scale[i];
	    param.dbuff[i]/=scale[i];
#else
	    vb=sol[i+A.nr*l].r;
	    sol[i+A.nr*l].r=vb*scale[i].r-sol[i+A.nr*l].i*scale[i].i;
	    sol[i+A.nr*l].i=vb*scale[i].i+sol[i+A.nr*l].i*scale[i].r;
	    val=1.0/(scale[i].r*scale[i].r+scale[i].i*scale[i].i);
	    vb=param.dbuff[i].r;
	    param.dbuff[i].r=( vb*scale[i].r+param.dbuff[i].i*scale[i].i)*val;
	    param.dbuff[i].i=(-vb*scale[i].i+param.dbuff[i].i*scale[i].r)*val;
#endif
	} // end for i
	
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	    param.dbuff[i]-=rhs[i+A.nr*l];
#else
#ifdef _COMPLEX_SYMMETRIC_
	    // transposed scaling
	    param.dbuff[i].r-=rhs[i+A.nr*l].r;
	    param.dbuff[i].i-=rhs[i+A.nr*l].i;
#else
	    // conjugate transposed scaling
	    param.dbuff[i].r-= rhs[i+A.nr*l].r;
	    param.dbuff[i].i-=-rhs[i+A.nr*l].i;
#endif
#endif
	} // end for i
	i=1;
	val=NRM(&n, param.dbuff, &i);
	// -----------------------------------------
	printf("current: %8.1le\n",val);

	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    sol[i+A.nr*l]-=1.0+i*l;
#else
	    sol[i+A.nr*l].r-=1.0+i*(l-1);
	    sol[i+A.nr*l].i-=0.1*(A.nr-i+l-1);
#endif
	}
	i=1;
	val=NRM(&n,sol+A.nr*l,&i);
	for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	    sol[i+A.nr*l]+=1.0+i*l;
#else
	    sol[i+A.nr*l].r+=1.0+i*(l-1);
	    sol[i+A.nr*l].i+=0.1*(A.nr-i+l-1);
#endif
	}
	i=1;
	vb=NRM(&n,sol+A.nr*l,&i);
	printf("rel. error in the solution: %8.1le\n\n",val/vb);
    } // end for l
    fclose(fo);

    
    // finally release memory of the preconditioner
    MYSYMAMGDELETE(&A,&PRE,&param);
} /* end helmholtz */








