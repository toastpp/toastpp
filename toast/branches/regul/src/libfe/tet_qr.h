// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe      File: tet_qr.h
//
// Quadrature rules over the local tetrahedron, defined by
// xi >= 0 && eta >= 0 && zeta >= 0 && xi+eta+zeta <= 1
// with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
// ==========================================================================

// All quadrature rules return the weights in array 'wght', the integration
// abscissae in array 'absc', and the number of points as the function
// return value

#ifndef __TET_QR_H
#define __TET_QR_H

#include <math.h>
#include "point.h"

FELIB int QRule_tet_1_1 (const double **wght, const Point **absc);
// Degree: 1, Points: 1
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-1-1s.html

FELIB int QRule_tet_2_4 (const double **wght, const Point **absc);
// Degree: 2, Points: 4
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-2-4bs.html

FELIB int QRule_tet_4_14 (const double **wght, const Point **absc);
// Degree: 2, Points: 14
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-4-14s.html

FELIB int QRule_tet_6_24 (const double **wght, const Point **absc);
// Degree: 4, Points: 24
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t3-6-24s.html

#endif // !__TET_QR_H
