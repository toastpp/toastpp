#ifndef _MUMPS_MATCHING_H
#define _MUMPS_MATCHING_H

/**************** include nameshsl.h if necessary *********/
#include "namesmumps_matching.h"

#include "long_integer.h"

void dmumps_match(integer *, integer *, integer *, 
		  integer *, integer *, doubleprecision *, 
		  doubleprecision *, doubleprecision *);
void smumps_match(integer *, integer *, integer *, 
		  integer *, integer *, real *, 
		  real *, real *);
void cmumps_match(integer *, integer *, integer *, 
		  integer *, integer *, complex *, 
		  real *, real *);
void zmumps_match(integer *, integer *, integer *, 
		  integer *, integer *, ilu_doublecomplex *, 
		  doubleprecision *, doubleprecision *);

#endif /* _MUMPS_MATCHING_H */
