// -*-C++-*-
// ==========================================================================
// File supertoast_util.h
// General utility function for supertoast-style reconstruction applications
// ==========================================================================

#ifndef __SUPERTOAST_UTIL_H
#define __SUPERTOAST_UTIL_H

void ReadDataFile (char *fname, RVector &data);

void WriteJacobian (const RMatrix *J, const Raster &raster,
    const QMMesh &mesh);

#endif // !__SUPERTOAST_UTIL_H
