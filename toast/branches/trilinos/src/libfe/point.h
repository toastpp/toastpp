// ==========================================================================
// Module libfe
// File point.h
// Declaration of classes Point, Point2D, Point3D
//
// Inheritance:
// ------------
// RVector ----> Point ----> Point2D
//                     ----> Point3D
// ==========================================================================

#ifndef __POINT_H
#define __POINT_H

// ==========================================================================
// class Point

class FELIB Point: public RVector {
public:
    // constructors
    Point () : RVector () {}
    Point (int dim) : RVector (dim) {}
    Point (const Point &p) : RVector (p) {}
    Point (const RVector &v) : RVector (v) {}

    // assignment
    Point operator= (const Point& p);
    Point operator= (const RVector &v);

    friend void Swap (Point &p1, Point &p2);
    // Swap two points. Assumes same dimension

    double Dist (const Point& p) const;
    // distance between points

    // rotate point around origin, by multiplying with matrix `rot',
    // and return rotated point in `rotpoint'
    void Rotate (const RDenseMatrix &rot, Point &rotpoint) const;
};


// ==========================================================================
// class Point2D

class Point2D: public Point {
public:
    Point2D (): Point (2) {}
    Point2D (const Point2D &p): Point (p) {}
    Point2D (double x1, double x2): Point (2)
	{ data[0] = x1; data[1] = x2; }
};


// ==========================================================================
// class Point3D

class Point3D: public Point {
public:
    Point3D (): Point (3) {}
    Point3D (const Point3D &p): Point (p) {}
    Point3D (double x1, double x2, double x3): Point (3)
	{ data[0] = x1; data[1] = x2; data[2] = x3; }
};

#endif // !__POINT_H


// ==========================================================================
// friend prototypes

void Swap (Point &p1, Point &p2);
