#include <stdlib.h>
#include <math.h>
#include "mathdef.h"

// ==========================================================================
// function instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template double min (const double x, const double y);
template float min (const float x, const float y);
template int min (const int x, const int y);

template double max (const double x, const double y);
template float max (const float x, const float y);
template int max (const int x, const int y);

#endif // NEED_EXPLICIT_INSTANTIATION
