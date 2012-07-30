// ==========================================================================
// Module mathlib
// File precon.cc
// Definition of template classes
//     TPreconditioner
//     TPrecon_Null
//     TPrecon_Diag
// ==========================================================================

#define MATHLIB_IMPLEMENTATION
#include "mathlib.h"

#ifdef HAVE_ILU
#include <pthread.h>
#include "ilutoast.h"
#include <ilupack.h>

pthread_mutex_t sid_mutex = PTHREAD_MUTEX_INITIALIZER;

template<>
void TPrecon_ILU<toast::complex>::Reset (TCompRowMatrix<toast::complex> &_A, int matching, char *ordering, double droptol, int condest, int elbow)
{
    CreateZmat(_A, &A);

    ZGNLAMGinit(&A, &param);
    param.matching=matching;
    param.ordering = ordering;
    param.droptol = droptol;
    param.droptolS = 0.1*param.droptol;
    param.condest=condest;
    param.elbow=elbow;

    long int ierr;
    ierr=ZGNLAMGfactor(&A, &PRE, &param);
    
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

    rhs = new ilu_doublecomplex[A.nr];
    sol = new ilu_doublecomplex[A.nr];
}

template<>
void TPrecon_ILU<toast::complex>::Apply (const TVector<toast::complex> &rh, TVector<toast::complex> &s)
{
	pthread_mutex_lock(&sid_mutex);
	memcpy(rhs, rh.data_buffer(), A.nr*sizeof(ilu_doublecomplex));
	ZGNLAMGsol(&PRE, &param, rhs, sol);
	memcpy(s.data_buffer(), sol, A.nr*sizeof(ilu_doublecomplex));
	pthread_mutex_unlock(&sid_mutex);
}

#endif // HAVE_ILU

