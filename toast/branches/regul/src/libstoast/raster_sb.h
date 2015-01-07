// -*-C++-*-

#ifndef __RASTER_SB_H
#define __RASTER_SB_H

/**
 * \file Implements class Raster_SplineBlob (radially symmetric basis
 *   functions with cubic spline profile)
 */

class STOASTLIB Raster_SplineBlob: public Raster_Blob {
public:
    Raster_SplineBlob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sup, RDenseMatrix *bb = 0);

protected:
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

private:
    double xi[5]; // cubic B-spline nodes
    double iomega[5];
    double scale;
};

#endif // !__RASTER_SB_H
