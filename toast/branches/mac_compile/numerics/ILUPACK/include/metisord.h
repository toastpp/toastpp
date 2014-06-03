#ifndef _METIS_ORD_H_
#define _METIS_ORD_H_

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1 1998/11/27 17:59:21 karypis Exp $
 */

#include "long_integer.h"


/* Undefine the following #define in order to use short integer as the idxtype */
#include "metis_idxtype.h"

void METIS_EdgeND(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 
void METIS_NodeND(integer *, idxtype *, idxtype *, integer *, integer *, idxtype *, idxtype *); 

#endif
