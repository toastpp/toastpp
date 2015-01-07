// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe
// File: tri_qr
//
// Quadrature rules over the local triangle, defined by
// xi >= 0 && eta >= 0 && xi+eta <= 1
// with vertices at (0,0), (1,0), (0,1)
// All rules are from www.cs.kuleuven.ac.be/~nines/research/ecf/Formules
// Login: cubature, passwd: kepler15
// ==========================================================================

// All quadrature rules return the weights in array 'wght', the integration
// abscissae in array 'absc', and the number of points as the function
// return value

#ifndef __TRI_QR_H
#define __TRI_QR_H

#include <math.h>
#include "point.h"

FELIB int QRule_tri_1_1 (const double **wght, const Point **absc);
// Degree: 1, Points: 1
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-1-1s.html

FELIB int QRule_tri_2_3 (const double **wght, const Point **absc);
// Degree: 2, Points: 3
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-2-3as.html

FELIB int QRule_tri_3_6 (const double **wght, const Point **absc);
// Degree: 3, Points: 6
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-3-6cs.html

FELIB int QRule_tri_4_6 (const double **wght, const Point **absc);
// Degree: 4, Points: 6
// ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-4-6as.html

FELIB int QRule_tri_4_9 (const double **wght, const Point **absc);
// Degree: 4, Points: 9
// Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-4-9al.html

FELIB int QRule_tri_5_7 (const double **wght, const Point **absc);
// Degree: 5, Points: 7
// ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-5-7s.html

FELIB int QRule_tri_6_12 (const double **wght, const Point **absc);
// Degree: 6, Points: 12
// ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-6-12as.html

FELIB int QRule_tri_7_12 (const double **wght, const Point **absc);
// Degree: 7, Points: 12
// ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-7-12s.html

FELIB int QRule_tri_7_15 (const double **wght, const Point **absc);
// Degree: 7, Points: 15
// ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-7-15s.html

FELIB int QRule_tri_9_19 (const double **wght, const Point **absc);
// Degree: 9, Points: 19
// ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-9-19s.html

#endif // !__TRI_QR_H
