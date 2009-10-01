#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

using namespace std;
using namespace toast;

// ==========================================================================
// class Raster_Blob

Raster_Blob::Raster_Blob (const IVector &_bdim, const IVector &_gdim,
    Mesh *mesh, RDenseMatrix *bb)
    : Raster (_bdim, _gdim, mesh, bb)
{
}
