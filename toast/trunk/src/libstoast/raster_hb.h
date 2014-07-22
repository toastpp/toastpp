// -*-C++-*-

#ifndef __RASTER_HB_H
#define __RASTER_HB_H

/**
 * \file Implements class Raster_HanningBlob (radially symmetric basis
 *   functions with Hanning profile)
 *
 * b_a(r) = 0.5(1+cos(Pi*r/a))    (a=radius of support)
 * Refs:
 * KM Hanson, GW Wecksung, Appl Opt 24, 4028-4039 (1985)
 */

class STOASTLIB Raster_HanningBlob: public Raster_Blob {
public:
    Raster_HanningBlob (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        double _sup, RDenseMatrix *bb = 0);

protected:
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

private:
    double fac, a2, scale;
};

#endif // !__RASTER_HB_H
