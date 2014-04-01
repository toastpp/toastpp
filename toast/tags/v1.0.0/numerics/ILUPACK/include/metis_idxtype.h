#ifndef _METIS_IDXTYPE_H_
#define _METIS_IDXTYPE_H_
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * idxtype.h extracted from struct.h
 *
 * This file contains data structures for ILU routines.
 *
 * Started 9/26/95
 * George
 *
 * $Id: struct.h,v 1.1 1998/11/27 17:59:31 karypis Exp $
 */

#include "long_integer.h"

/* Undefine the following #define in order to use short integer as the idxtype */


#define IDXTYPE_INT

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef integer idxtype;
#else
typedef short idxtype;
#endif

#ifdef _LONG_INTEGER_
/* not quite correct but bypasses the << operator problems for long int */
/* #define MAXIDX	((long)16*(1<<8*sizeof(int)-2))*/
#define MAXIDX	((long)1<<(long)8*sizeof(idxtype)-2)
#else
#define MAXIDX	(1<<8*sizeof(idxtype)-2)
#endif 


#endif
