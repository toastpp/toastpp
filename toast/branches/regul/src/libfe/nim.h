// -*-C++-*-
// Utilities for nodal images

#ifndef __NIM_H
#define __NIM_H

#include "mathlib.h"
#include "felib.h"

FELIB int ReadNim (const char *name, int idx, RVector &img,
    char *meshname=NULL);

FELIB bool ReadNimAll (char *name, RDenseMatrix &img);

FELIB void WriteNim (const char *name, const char *meshname,
    const RVector &img);

FELIB int NimPhaseUnwrap (const Mesh *mesh, RVector &phase, Point seed);

#endif // !__NIM_H
