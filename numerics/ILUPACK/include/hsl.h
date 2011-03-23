#ifndef _HSL_H
#define _HSL_H

/**************** include nameshsl.h if necessary *********/
#include "nameshsl.h"

#include "long_integer.h"

/* Common Block Declarations and Initialized Data */


void mc64ad(integer *, integer *, integer *, integer *, integer *, 
	    doubleprecision *, integer *, integer *, integer *, integer *, 
	    integer *, doubleprecision *, integer *, integer *);
void mc64as(integer *, integer *, integer *, integer *, integer *, 
	    real *, integer *, integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *);
void mc64az(integer *, integer *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *);
void mc64ac(integer *, integer *, integer *, integer *, integer *, 
	    complex *, integer *, integer *, integer *, integer *, 
	    integer *, complex *, integer *, integer *);

void mc64id(integer *);
void mc64is(integer *);
void mc64iz(integer *);
void mc64ic(integer *);

doubleprecision fd05ad(integer *); 
real            fd05as(integer *); 

#endif /* _HSL_H */
