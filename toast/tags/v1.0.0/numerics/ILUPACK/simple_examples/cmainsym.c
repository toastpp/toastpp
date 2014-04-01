#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ilupack.h>



int main(int argc, char **argv)
{
    long int      ia[9] ={   1,   5,   8,  10,  12,  15,  17,  18,  19},
                  ja[18]={   1,        3,             6,   7,
                                  2,   3,        5,
                                       3,                       8,
                                            4,             7,
                                                 5,   6,   7,
                                                      6,        8,
                                                           7,
                                                                8};
    complex       a[18], rhs[8], sol[8];
    Cmat          A;
    long int      i,ierr;
    CAMGlevelmat  PRE;
    CILUPACKparam param;


    rhs[0].r=1.0; rhs[1].r=1.0; rhs[2].r=1.0; rhs[3].r=1.0;
    rhs[0].i=0.0; rhs[1].i=0.0; rhs[2].i=0.0; rhs[3].i=0.0;
    rhs[4].r=1.0; rhs[5].r=1.0; rhs[6].r=1.0; rhs[7].r=1.0;
    rhs[4].i=0.0; rhs[5].i=0.0; rhs[6].i=0.0; rhs[7].i=0.0;


    /* initialize matrix*/
    a[ 0].r= 7.0; a[ 1].r= 1.0; a[ 2].r= 2.0; a[3].r= 7.0;
    a[ 0].i= 1.0; a[ 1].i= 0.0; a[ 2].i=-1.0; a[3].i=-1.0;

    a[ 4].r=-4.0; a[ 5].r= 8.0; a[ 6].r= 2.0;
    a[ 4].i= 1.0; a[ 5].i= 0.0; a[ 6].i=-1.0;

    a[ 7].r= 1.0; a[ 8].r= 5.0;
    a[ 7].i= 0.0; a[ 8].i=-1.0;

    a[ 9].r= 7.0; a[10].r= 9.0;
    a[ 9].i= 0.0; a[10].i=-1.0;

    a[11].r= 5.0; a[12].r= 1.0; a[13].r= 5.0;
    a[11].i= 0.0; a[12].i=-1.0; a[13].i=-1.0;

    a[14].r= 0.0; a[15].r= 5.0;
    a[14].i= 0.0; a[15].i=-1.0;

    a[16].r=11.0;
    a[16].i= 0.0;

    a[17].r= 5.0;
    a[17].i= 0.0;

    A.nr=A.nc=8;
    A.nnz=18;
    A.ia=ia;
    A.ja=ja;
    A.a=a;

    /* ILUPACK uses compressed sparse ROW format. 
                                   / 3.5   -1.0+i  0 \
       A complex symmetric matrix  |-1.0+i  2.0-i  0 | is stored as follows
                                   \   0      0   1.5/ 

       A.ia:   1  3  4  5           pointer to the start of every compressed 
                                    row (upper triangular part) plus pointer
                                    to the first space behind the compressed
                                    rows

       A.ja:   1     2     2     3     nonzero column indices

       A.a:   3.5 -1.0+i 2.0-i  1.5    nonzero numerical values

       The read part finally yields the following data structures
        -  A:  matrix in compressed sparse row format
	    o  A.nr, A.nc: number of rows and columns of A
            o  A.ia:  pointer array
            o  A.ja:  nonzero column index array 
	    o  A.a:   nonzero numerical values
	-  rhs:  right hand side(s) and additional data like exact solution
	         or initial guess
	-  n:  same as A.nr,A.nc
	-  nnz:  number of nonzero entries
     */


    // set parameters to the default settings
    CSYMAMGinit(&A, &param);


    // Here you can (if you like) change some of the default options
    
    //  1. maximum weight matching 
    // by default matching is turned on and the drivers are properly 
    // pre-selected. If you do not like to have matching, then simply turn it
    // off
    param.matching=1;

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
    param.ordering="metisn";

    //  3. drop tolerance for the LU factors
    // by default, 1e-2 is chosen. Here you can overwrite the default values
    // This sample program is designed to pass your own drop tolerance (called
    // (DROP_TOL) via the command line when being called
    param.droptol=0.1;

    //  4. drop tolerance for the approximate Schur complement
    // by default, 1e-2 is chosen. Here you can overwrite the default values
    // This sample program is designed to pass your own drop tolerance (called
    // (DROP_TOL) via the command line when being called. In this example, we
    // use the drop tolerance as for the LU factors but multiplied by 0.1
    param.droptolS=0.1*param.droptol;

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
    param.condest=5;

    //  6. residual tolerance for the iterative solver
    // The built-in iterative solver (here restarted GMRES) will use this
    // tolerance to terminate whenever the backward error is less than this
    // threshold. By default, sqrt(eps)~1e-8 is chosen for double precision, 
    // eps^(3/4)~1e-6 is chosen for single precision. 
    // param.restol=1e-14;

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
    param.elbow=10;

    //  9. maximum number of nonzeros per column in L (resp. per row in U)
    // by default n+1 is chosen, i.e. this option is disabled. You can limit
    // the amount of memory by using some smaller value, e.g. A.ia[A.nc]-1
    // is the fill of A and ELBOW*(A.ia[A.nc]-1.0)/A.nc would restrict the
    // maximum number of fill to the average number of nonzeros of A per column
    // (or per row) times the ELBOW. Note however that this parameter cuts off
    // the fill in L and U by brute force. It recommended NOT to use it.
    // param.lfil=ELBOW*(A.ia[A.nc]-1.0)/A.nc;

    // 10. maximum number of nonzeros per row in S (approximate Schur 
    //     complement)
    // by default n+1 is chosen, i.e. this option is disabled. You can limit
    // the amount of memory by using some smaller value, e.g. A.ia[A.nc]-1
    // is the fill of A and ELBOW*(A.ia[A.nc]-1.0)/A.nc would restrict the
    // maximum number of fill to the average number of nonzeros of A per column
    // (or per row) times the ELBOW. Note however that this parameter cuts off
    // the fill in S by brute force. It strongly recommended NOT to use it.
    // param.lfilS=ELBOW*(A.ia[A.nc]-1.0)/A.nc;

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
    // for (i=0; i<A.nc; i++) sol[i]=1.0; param.tv=sol;

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

    // 15. number of post-smoothing steps
    // If classical multigrid is selected (param.amg="mg";), then here you can
    // set the number of post-smoothing steps. default: 1
    // param.npostsmoothing=1;

    // 16. number of coarse grid solves
    // Except for multilevel ILU (i.e. param.amg="amli"; or param.amg="mg";),
    // here you define how often the coarse grid solve is performed. By default,
    // only one coarse grid solve is used (V-cycle). The choice param.ncoarse=2;
    // would correspond to a W-cycle. Note however, if a negative value is
    // passed, a flexible solver is invoked, i.e. the number of coarse grid
    // solves varies from one grid to another and from one step to the next one.
    // param.ncoarse=-1;

    // 17. type of pre-smoother
    // if full multigrid is used (param.amg="mg";), then here you can choose 
    // between built-in smoothers or your own hand-made smoother.
    //   (a) "gsf"     (default) Gauss-Seidel forward
    //   (b) "gsb"     Gauss-Seidel backward
    //   (c) "j"       (damped) Jacobi
    //   (d) "ilu"     ILU on the fine grid system
    //   (e) any other string that does not match (a)-(d) will cause AMGsolver
    //       to use reverse communication principle in order to let you provide
    //       your own smoother. In that case ILUPACK will give you the matrix,
    //       the right hand side and an initial solution (typically 0). You have
    //       to override the initial solution 
    // param.presmoother="gsf";

    // 18. type of post-smoother
    // if full multigrid is used (param.amg="mg";), then here you can choose 
    // between built-in smoothers or your own hand-made smoother.
    //   (a) "gsb"     (default) Gauss-Seidel backward
    //   (b) "gsf"     Gauss-Seidel forward
    //   (c) "j"       (damped) Jacobi
    //   (d) "ilu"     ILU on the fine grid system
    //   (e) any other string that does not match (a)-(d) will cause AMGsolver
    //       to use reverse communication principle in order to let you provide
    //       your own smoother. In that case ILUPACK will give you the matrix,
    //       the right hand side and an initial solution (typically 0). You have
    //       to override the initial solution 
    // param.postsmoother="gsb";

    // 19. pre-selection of coarse grid nodes
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
    // param.FCpart="yes";
    
    // 20. type of coarse grid system
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

    // 21. number of steps before GMRES is restarted
    // The iterative solver uses restarted GMRES (resp. FGMRES). By default, 30
    // steps are computed, before the method is restarted. Note that a smaller
    // number reduces the memory, while a larger number can improve the 
    // convergence.
    // param.nrestart=30;

    // 22. require the computation of the preconditioner in
    //     single precision
    // param.mixedprecision=1;

    if (param.mixedprecision)
       param.solver="fgmres";

    ierr=CSYMAMGfactor(&A, &PRE, &param);

    switch (ierr)
    {
           case  0: /* perfect! */
	            printf("factorization successful with %d levels completed\n", 
			   PRE.nlev);
		    printf("final elbow space factor=%8.2f\n",param.elbow+0.005);
	            break;
           case -1: /* Error. input matrix may be wrong.
                       (The elimination process has generated a
			row in L or U whose length is .gt.  n.) */
	            printf("Error. input matrix may be wrong at level %d\n",
			   PRE.nlev);
		    break;
           case -2: /* The matrix L overflows the array alu */
	            printf("The matrix L overflows the array alu at level %d\n",
			   PRE.nlev);
		    break;
           case -3: /* The matrix U overflows the array alu */
	            printf("The matrix U overflows the array alu at level %d\n",
			   PRE.nlev);
		    break;
           case -4: /* Illegal value for lfil */
	            printf("Illegal value for lfil at level %d\n",PRE.nlev);
		    break;
           case -5: /* zero row encountered */
	            printf("zero row encountered at level %d\n",PRE.nlev);
		    break;
           case -6: /* zero column encountered */
	            printf("zero column encountered at level %d\n",PRE.nlev);
		    break;
           case -7: /* buffers too small */
	            printf("buffers are too small\n");
           default: /* zero pivot encountered at step number ierr */
	            printf("zero pivot encountered at step number %d of level %d\n",
			   ierr,PRE.nlev);
		    break;
    } /* end switch */
    if (ierr) {
       exit(ierr);
    }



    /* -------------------------------------------------------------------- */
    // demonstrate a simple forward/backward solve with the multilevel ILU
    CSYMAMGsol(&PRE, &param, rhs,sol);
    printf("approximate solution after one preconditioning step\n");
    for (i=0; i<A.nc; i++) {
        printf("%24.16le+%24.16lei\n",sol[i].r,sol[i].i);
    }
    printf("\n");
    /* -------------------------------------------------------------------- */



    /* -------------------------------------------------------------------- */
    // demonstrate iterative solver
    for (i=0; i<A.nc; i++) {
        sol[i].r=sol[i].i=0;
    }

    ierr=CSYMAMGsolver(&A, &PRE, &param, rhs, sol);
	
    // why did the iterative solver stop?
    switch (ierr) {
    case  0:  // everything is fine
              printf("iteration successful completed after %d steps\n",param.ipar[25]);
	      break;
    case -1:  // too many iterations
              printf("number of iteration steps exceeds its limit\n");
	      break;
    case -2:  /* not enough work space */
              printf("not enough work space provided\n");
	      break;
    case -3:  /* not enough work space */
	      printf("algorithm breaks down\n");
	      break;
    default:  /* unknown */
	      printf("solver exited with error code %d\n",ierr);
    } // end switch 

    printf("approximate solution after completing the iterative process\n");
    for (i=0; i<A.nc; i++) {
        printf("%24.16le+%24.16lei\n",sol[i].r,sol[i].i);
    }
    printf("\n");
    /* -------------------------------------------------------------------- */




    // finally release memory of the preconditioner
    CSYMAMGdelete(&A,&PRE,&param);
} // end main

