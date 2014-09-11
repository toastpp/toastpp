// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe
// File: tet_qr
//
// Quadrature rules over the local tetrahedron, defined by
// xi >= 0 && eta >= 0 && zeta >= 0 && xi+eta+zeta <= 1
// with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include "felib.h"
#include "tet_qr.h"

FELIB int QRule_tet_1_1 (const double **wght, const Point **absc)
{
    // Region: Tetrahedron
    // Degree: 1
    // Points: 1
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-1-1s.html

    static const double w[1] = {
        0.16666666666666666
    };
    static const Point3D a[1] = {
        Point3D(0.25, 0.25, 0.25)
    };
    *wght = w;
    *absc = a;
    return 1;
}

FELIB int QRule_tet_2_4 (const double **wght, const Point **absc)
{
    // Region: Tetrahedron
    // Degree: 2
    // Points: 4
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-2-4bs.html

    static const double w[4] = {
        0.041666666666666666, 0.041666666666666666,
	0.041666666666666666, 0.041666666666666666
    };
    static const Point3D a[4] = {
        Point3D(0.13819660112501051, 0.13819660112501051, 0.5854101966249685),
	Point3D(0.13819660112501051, 0.5854101966249685, 0.13819660112501051),
	Point3D(0.5854101966249685, 0.13819660112501051, 0.13819660112501051),
	Point3D(0.13819660112501051, 0.13819660112501051, 0.13819660112501051)
    };
    *wght = w;
    *absc = a;
    return 4;
}

FELIB int QRule_tet_4_14 (const double **wght, const Point **absc)
{
    // Region: Tetrahedron
    // Degree: 4
    // Points: 14
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-4-14s.html

    static const double w[14] = {
        3.1746031746031746e-3, 3.1746031746031746e-3, 3.1746031746031746e-3,
	3.1746031746031746e-3, 3.1746031746031746e-3, 3.1746031746031746e-3,
	0.014764970790496785, 0.014764970790496785,
	0.014764970790496785, 0.014764970790496785,
	0.022139791114265119, 0.022139791114265119,
	0.022139791114265119, 0.022139791114265119
    };
    static const Point3D a[14] = {
        Point3D(0.5, 0.0, 0.0), Point3D(0.0, 0.5, 0.0),
	Point3D(0.0, 0.0, 0.5), Point3D(0.5, 0.5, 0.0),
	Point3D(0.5, 0.0, 0.5), Point3D(0.0, 0.5, 0.5),
	Point3D(0.10052676522520447, 0.10052676522520447, 0.6984197043243866),
	Point3D(0.10052676522520447, 0.6984197043243866, 0.10052676522520447),
	Point3D(0.6984197043243866, 0.10052676522520447, 0.10052676522520447),
	Point3D(0.10052676522520447, 0.10052676522520447, 0.10052676522520447),
	Point3D(0.31437287349319219, 0.31437287349319219, 0.0568813795204234),
	Point3D(0.31437287349319219, 0.0568813795204234, 0.31437287349319219),
	Point3D(0.0568813795204234, 0.31437287349319219, 0.31437287349319219),
	Point3D(0.31437287349319219, 0.31437287349319219, 0.31437287349319219)
    };
    *wght = w;
    *absc = a;
    return 14;
}

FELIB int QRule_tet_6_24 (const double **wght, const Point **absc)
{
    // Region: Tetrahedron
    // Degree: 6
    // Points: 24
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-6-24s.html

    static const double w[24] = {
        6.6537917096945820e-3, 6.6537917096945820e-3,
	6.6537917096945820e-3, 6.6537917096945820e-3,
	1.6795351758867738e-3, 1.6795351758867738e-3,
	1.6795351758867738e-3, 1.6795351758867738e-3,
	9.2261969239424536e-3, 9.2261969239424536e-3,
	9.2261969239424536e-3, 9.2261969239424536e-3,
	8.0357142857142857e-3, 8.0357142857142857e-3,
	8.0357142857142857e-3, 8.0357142857142857e-3,
	8.0357142857142857e-3, 8.0357142857142857e-3,
	8.0357142857142857e-3, 8.0357142857142857e-3,
	8.0357142857142857e-3, 8.0357142857142857e-3,
	8.0357142857142857e-3, 8.0357142857142857e-3
    };
    static const Point3D a[24] = {
        Point3D(0.21460287125915202, 0.21460287125915202, 0.3561913862225439),
	Point3D(0.21460287125915202, 0.3561913862225439, 0.21460287125915202),
	Point3D(0.3561913862225439, 0.21460287125915202, 0.21460287125915202),
	Point3D(0.21460287125915202, 0.21460287125915202, 0.21460287125915202),
	Point3D(0.040673958534611353,0.040673958534611353,0.87797812439616594),
	Point3D(0.040673958534611353,0.87797812439616594,0.040673958534611353),
	Point3D(0.87797812439616594,0.040673958534611353,0.040673958534611353),
	Point3D(0.040673958534611353,0.040673958534611353,0.040673958534611353),
	Point3D(0.32233789014227551, 0.32233789014227551, 0.0329863295731735),
	Point3D(0.32233789014227551, 0.0329863295731735, 0.32233789014227551),
	Point3D(0.0329863295731735, 0.32233789014227551, 0.32233789014227551),
	Point3D(0.32233789014227551, 0.32233789014227551, 0.32233789014227551),
	Point3D(0.063661001875017525, 0.26967233145831580, 0.6030056647916492),
	Point3D(0.063661001875017525, 0.6030056647916492, 0.26967233145831580),
	Point3D(0.26967233145831580, 0.063661001875017525, 0.6030056647916492),
	Point3D(0.26967233145831580, 0.6030056647916492, 0.063661001875017525),
	Point3D(0.6030056647916492, 0.063661001875017525, 0.26967233145831580),
	Point3D(0.6030056647916492, 0.26967233145831580, 0.063661001875017525),
	Point3D(0.063661001875017525,0.063661001875017525,0.6030056647916492),
	Point3D(0.063661001875017525,0.6030056647916492,0.063661001875017525),
	Point3D(0.6030056647916492,0.063661001875017525,0.063661001875017525),
	Point3D(0.063661001875017525,0.063661001875017525,0.26967233145831580),
	Point3D(0.063661001875017525,0.26967233145831580,0.063661001875017525),
	Point3D(0.26967233145831580,0.063661001875017525,0.063661001875017525)
    };
    *wght = w;
    *absc = a;
    return 24;
}
