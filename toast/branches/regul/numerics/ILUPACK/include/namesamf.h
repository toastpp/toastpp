#ifndef _NAMES_AMF_H
#define _NAMES_AMF_H

#include "f2c.h"


/* on several architectures names of fortran routines are passed to C in 
   different ways. To cover this different architectures use in C only lower
   letters for the fortran names. Dependent on the switch you use they are
   replaced by the correct function name
*/

/* only use capital letters */
#if defined __CAPS__ && !defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define r2hamdf4            R2HAMDF4


/* no capital letters */
#elif defined __UNDERSCORE__ && !defined __CAPS__ && !defined __2UNDERSCORES__
#define r2hamdf4            r2hamdf4_


/* both are defined */
#elif defined __CAPS__ && defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define r2hamdf4            R2HAMDF4_


/* CAPS and 2 underscores are defined */
#elif defined __CAPS__ && defined __2UNDERSCORES__
#define r2hamdf4            R2HAMDF4__


/* no capital letters but 2 underscores */
#elif defined __2UNDERSCORES__ && !defined __CAPS__
#define r2hamdf4            r2hamdf4__


// no switch defined use lower case letters in FORTRAN
#else
#define r2hamdf4            r2hamdf4

#endif /* defined __CAPS__ && ... */


#endif /* _NAMES_AMF_H */
