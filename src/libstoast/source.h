#ifndef __SOURCE_H
#define __SOURCE_H

#include "stoastlib.h"
#include "toasttype.h"

STOASTLIB RVector QVec_Gaussian (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode);

STOASTLIB RVector QVec_Cosine (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode);

STOASTLIB RVector QVec_Point (const Mesh &mesh, const Point &cnt, SourceMode mode);

#endif // !__SOURCE_H
