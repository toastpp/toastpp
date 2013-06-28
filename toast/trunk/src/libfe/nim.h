// -*-C++-*-
// Utilities for nodal images

#ifndef __NIM_H
#define __NIM_H

#include "mathlib.h"
#include "felib.h"

void NimPhaseUnwrap (const Mesh *mesh, RVector &phase, Point seed);

#endif // !__NIM_H
