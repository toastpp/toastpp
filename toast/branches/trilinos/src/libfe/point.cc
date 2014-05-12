// ==========================================================================
// Module libfe
// File point.cc
// Definition of class Point
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include "mathlib.h"
#include "felib.h"

Point Point::operator= (const Point& p)
{
    Copy (p);
    return *this;
}

Point Point::operator= (const RVector& v)
{
    Copy (v);
    return *this;
}

void Swap (Point &p1, Point &p2)
{
    dASSERT(p1.Dim() == p2.Dim(), "Point dimensions differ");
    double tmp;
    for (int d = p1.Dim()-1; d >= 0; d--)
        tmp = p1[d], p1[d] = p2[d], p2[d] = tmp;
}

double Point::Dist (const Point &p) const
{
    dASSERT(size == p.Dim(), "Points have different dimension.");
    double d2;
    int i;
    for (d2 = 0.0, i = 0; i < size; i++)
	d2 += (data[i]-p.data[i]) * (data[i]-p.data[i]);
    return sqrt (d2);
}

void Point::Rotate (const RDenseMatrix &rot, Point &rotpoint) const
{
    dASSERT (rot.nRows() == size && rot.nCols() == size,
	     "Wrong matrix dimension");
    dASSERT (rotpoint.Dim() == size, "Wrong point dimension");
    int r, c;

    for (r = 0; r < size; r++) {
        rotpoint[r] = 0.0;
        for (c = 0; c < size; c++)
	    rotpoint[r] += data[c] * rot.Get(r,c);
    }
}
