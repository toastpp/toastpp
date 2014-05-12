#ifndef __SOURCE_H
#define __SOURCE_H

#include "mathlib.h"
#include "felib.h"

typedef enum {
    SRCMODE_NEUMANN,
    SRCMODE_ISOTROPIC
} SourceMode;

RVector QVec_Gaussian (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode);

RVector QVec_Cosine (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode);

RVector QVec_Point (const Mesh &mesh, const Point &cnt, SourceMode mode);

#endif // !__SOURCE_H
