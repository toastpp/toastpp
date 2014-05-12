#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <sparspak.h>
#include <ilupack.h>

#include <ilupackmacros.h>
#include "symswitches.c"

//#define PRINT_INFO
#define STDERR          stdout
#define STDOUT          stdout
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))




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
    char *tmpstring, timeout[7], residout[7];
    integer  i,j,k,m,fnamelen,n,nc,nz,nrhs,mynrhs=1,tmp0,tmp,tmp2,tmp3,ierr,flag,
         *invperm, *buff, *ws, *symmmd,flags,elbow,max_it,ipar[20],
         nrhsix, *rhsptr, *rhsind, sumit, TYPEPERM;
    REALS eps, DROP_TOL, CONDEST,condest,droptols[2],restol, val,vb;
    FLOAT *rhs,*sol,*w, *scale, *rhsval, *dbuff;
    float  systime, time_start,   time_stop,   secnds, secndsprep,
           timeAx_start, timeAx_stop, secndsAx, sumtime;
    integer ELBOW, nnzU, l, nAx;
    ILUPACKPARAM param;
    long Fparam, FPREC;
    size_t nibuff, ndbuff;
     
    /* the last argument passed serves as file name */
    if (argc!=6) {
      printf("usage '%s <drop tol.> <bound for L^{-1},U^{-1}> <elbow space> <type permutation> <matrix>'\n",argv[0]);
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

    TYPEPERM=atoi(argv[argc-2]);
    ELBOW   =atoi(argv[argc-3]);
    CONDEST =atof(argv[argc-4]);
    DROP_TOL=atof(argv[argc-5]);

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
    fprintf(fo,"\n");
    fflush(fo);
    fclose(fo);

    /* Depending on whether a right hand side (and a solution) has been provided,
       generate an artificial solution vector in order to measure the quality of
       the iterative solver
    */

    // allocate memory for the solution vector
    sol  =(FLOAT *) MALLOC((size_t)n*sizeof(FLOAT), "mainspd:sol");

    // seek for an exact solution
    if (rhstyp[2]=='X' || rhstyp[2]=='x') {
       j=nrhs*A.nr;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') 
	 j*=2;
       if (rhstyp[0]=='M' || rhstyp[0]=='m')
	 j-=nrhs*A.nr;
       for (i=0; i<A.nr; i++,j++)
	   sol[i]=rhs[j];
    }

    // release part of rhs that may store the uncompressed rhs
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
       rhs-=A.nr;
    }

    // set up right hand side
    if (nrhs==0) {
      // prescribe artificial solution
      for (i=0; i<A.nr; i++)  {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	  sol[i]=1.0;
#else
	  sol[i].r=1.0+i*(-1);
  sol[i].i=0.1*(A.nr-i-1);
#endif
      } // end for i

      MYSMATVEC(A,sol,rhs);

    } // end if (nrhs==0)
    else {
       if (rhstyp[0]=='M' || rhstyp[0]=='m') {
	  for (i=0; i<A.nr; i++) {
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


    // initial solution
    if (rhstyp[1]=='G' || rhstyp[1]=='g') {
       j=nrhs*A.nr;
       if (rhstyp[0]=='M' || rhstyp[0]=='m')
	  j=A.nr;
       for (i=0; i<A.nr; i++,j++)
	   sol[i]=rhs[j];
    }
    else
       for (i=0; i<A.nr; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   sol[i]=0;
#else
	   sol[i].r=sol[i].i=0;
#endif
       }

    /*
      TYPEPERM   type of reordering scheme
                 0:          nothing
		 1:          symmetric minimum weight matching  + multiple minimum degree
		 2:          symmetric minimum weight matching  + approximate minimum fill
		 3:          symmetric minimum weight matching  + approximate minimum degree from UMFPACK 4.3
		 4:          symmetric minimum weight matching  + reverse Cuthill-McKee
		 5:          symmetric minimum weight matching  + MeTiS multilevel nested dissection by edges
		 6:          symmetric minimum weight matching  + MeTiS multilevel nested dissection by nodes
    */
    //ierr=SYMILUPACK(&n, A.ia, A.ja, A.a, rhs, sol, &DROP_TOL, &CONDEST, &ELBOW, &TYPEPERM);
    
    restol=1e-14;
    max_it=200;

    ierr=SYMILUPACKFAC(&Fparam, &FPREC, &n, A.ia, A.ja, A.a, rhs, sol,
		       &DROP_TOL, &CONDEST, &ELBOW, &TYPEPERM, &nlev, &restol,
		       &max_it);

    switch (ierr) {
    case  0: /* perfect! */
      break;
    case -1: /* Error. input matrix may be wrong.
		(the elimination process has generated a
		row in L or U whose length is .gt.  n.) */
      fprintf(stderr,"Error. input matrix may be wrong at level %d\n",nlev);
      break;
    case -2: /* The matrix L overflows the array alu */
      fprintf(stderr,"The matrix L overflows the array alu at level %d\n",nlev);
      break;
    case -3: /* The matrix U overflows the array alu */
      fprintf(stderr,"The matrix U overflows the array alu at level %d\n",nlev);
      break;
    case -4: /* Illegal value for lfil */
      fprintf(stderr,"Illegal value for lfil at level %d\n",nlev);
      break;
    case -5: /* zero row encountered */
      fprintf(stderr,"zero row encountered at level %d\n",nlev);
      break;
    case -6: /* zero column encountered */
      fprintf(stderr,"zero column encountered at level %d\n",nlev);
      break;
    case -7: /* buffers too small */
      fprintf(stderr,"buffers are too small\n");
      // This error message would not be necessary if AMGsetup is
      // called with the correct values of nibuff, ndbuff
    default: /* zero pivot encountered at step number ierr */
      fprintf(stderr,"zero pivot encountered at step number %d of level %d\n",ierr,nlev);
      break;
    } /* end switch */
    if (ierr) return (ierr);


    ierr=SYMILUPACKSOL(&Fparam, &FPREC, &n, A.ia, A.ja, A.a, rhs, sol,
		       &DROP_TOL, &CONDEST, &ELBOW, &TYPEPERM, &nlev, &restol,
		       &max_it);

    // why did the iterative solver stop?
    switch (ierr) {
    case  0:  // everything is fine
      break;
    case -1:  // too many iterations
      fprintf(stderr,"number of iteration steps exceeds its limit\n");
      break;
    case -2:  /* not enough work space */
      fprintf(stderr,"not enough work space provided\n");
      break;
    case -3:  /* not enough work space */
      fprintf(stderr,"algorithm breaks down\n");
      break;
    default:  /* unknown */
      fprintf(stderr, "solver exited with error code %d\n",ierr);
    } // end switch 
    if (ierr) return(ierr);

    // finally release memory of the preconditioner
    ierr=SYMILUPACKDEL(&Fparam, &FPREC, &n, A.ia, A.ja, A.a, rhs, sol,
		       &DROP_TOL, &CONDEST, &ELBOW, &TYPEPERM, &nlev, &restol,
		       &max_it);


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
   
    // release memory used for the input matrix
    free(A.ia);
    free(A.ja);
    free(A.a );
    free(rhs );
    free(sol );
} /* end mainsym */








