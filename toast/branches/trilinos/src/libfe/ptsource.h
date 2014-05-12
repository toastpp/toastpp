// ==========================================================================
// Module libfe
// File ptsource.h
// Declaration of classes PointSource
//
// Inheritance:
// ------------
// RVector ----> Point ----> PointSource
// ==========================================================================

#ifndef __PTSOURCE_H
#define __PTSOURCE_H

// ==========================================================================
// class PointSource

class PointSource: public Point {
public:
    PointSource ();
    PointSource (int dim): Point (dim) { strength = 0.0; }
    PointSource (const PointSource &ps);

    // friends
    friend std::ostream &operator<< (std::ostream &os, const PointSource &ps);
    friend std::istream &operator>> (std::istream &is, PointSource &ps);

    double strength;
};

// ==========================================================================
// friend prototypes

std::ostream &operator<< (std::ostream &os, const PointSource &ps);
std::istream &operator>> (std::istream &is, PointSource &ps);

#endif // !__PTSOURCE_H
