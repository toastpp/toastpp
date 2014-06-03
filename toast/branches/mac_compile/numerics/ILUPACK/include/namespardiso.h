#ifndef _NAMES_PARDISO_H
#define _NAMES_PARDISO_H

#include "f2c.h"


/* on several architectures names of fortran routines are passed to C in 
   different ways. To cover this different architectures use in C only lower
   letters for the fortran names. Dependent on the switch you use they are
   replaced by the correct function name
*/

/* only use capital letters */
#if defined __CAPS__ && !defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define pardisoinit         PARDISOINIT
#define pardiso             PARDISO

/* no capital letters */
#elif defined __UNDERSCORE__ && !defined __CAPS__ && !defined __2UNDERSCORES__
#define pardisoinit         pardisoinit_
#define pardiso             pardiso_


/* both are defined */
#elif defined __CAPS__ && defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define pardisoinit         PARDISOINIT_
#define pardiso             PARDISO_


/* CAPS and 2 underscores are defined */
#elif defined __CAPS__ && defined __2UNDERSCORES__
#define pardisoinit         PARDISOINIT__
#define pardiso             PARDISO__


/* no capital letters but 2 underscores */
#elif defined __2UNDERSCORES__ && !defined __CAPS__
#define pardisoinit         pardisoinit__
#define pardiso             pardiso__


// no switch defined use lower case letters in FORTRAN
#else
#define pardisoinit         pardisoinit
#define pardiso             pardiso

#endif /* defined __CAPS__ && ... */


#endif /* _NAMES_ILU_PACK_H */
