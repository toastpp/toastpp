#ifndef _NAMES_MUMPS_MATCHING_H
#define _NAMES_MUMPS_MATCHING_H

#include "f2c.h"

/* on several architectures names of fortran routines are passed to C in 
   different ways. To cover this different architectures use in C only lower
   letters for the fortran names. Dependent on the switch you use they are
   replaced by the correct function name
*/

/* only use capital letters */
#if defined __CAPS__ && !defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define smumps_match          SMUMPS_MATCH
#define dmumps_match          DMUMPS_MATCH
#define cmumps_match          CMUMPS_MATCH
#define zmumps_match          ZMUMPS_MATCH

/* no capital letters */
#elif defined __UNDERSCORE__ && !defined __CAPS__ && !defined __2UNDERSCORES__
#define smumps_match          smumps_match_
#define dmumps_match          dmumps_match_
#define cmumps_match          cmumps_match_
#define zmumps_match          zmumps_match_

/* both are defined */
#elif defined __CAPS__ && defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define smumps_match          SMUMPS_MATCH_
#define dmumps_match          DMUMPS_MATCH_
#define cmumps_match          CMUMPS_MATCH_
#define zmumps_match          ZMUMPS_MATCH_


/* CAPS and 2 underscores are defined */
#elif defined __CAPS__ && defined __2UNDERSCORES__
#define smumps_match          SMUMPS_MATCH__
#define dmumps_match          DMUMPS_MATCH__
#define cmumps_match          CMUMPS_MATCH__
#define zmumps_match          ZMUMPS_MATCH__


/* no capital letters but 2 underscores */
#elif defined __2UNDERSCORES__ && !defined __CAPS__
#define smumps_match          smumps_match__
#define dmumps_match          dmumps_match__
#define cmumps_match          cmumps_match__
#define zmumps_match          zmumps_match__

#else
#define smumps_match          smumps_match
#define dmumps_match          dmumps_match
#define cmumps_match          cmumps_match
#define zmumps_match          zmumps_match

#endif /* defined __CAPS__ && ... */

#endif /* _NAMES_MUMPS_MATCHING_H */


