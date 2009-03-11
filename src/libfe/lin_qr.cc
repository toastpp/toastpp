// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe
// File: lin_qr
//
// Gaussian quadrature rules over the 1D interval [0,1]
// ==========================================================================

#include "lin_qr.h"

int QRule_lin_1 (const double **wght, const double **absc)
{
    // Region: Line
    // Degree: 1
    // Points: 1
  
    static const double w[1] = {
        1.0
    };
    static const double a[1] = {
        0.5
    };
    *wght = w;
    *absc = a;
    return 1;
}

int QRule_lin_2 (const double **wght, const double **absc)
{
    // Region: Line
    // Degree: 2
    // Points: 2

    static const double w[2] = {
        0.5, 0.5
    };
    static const double a[2] = {
        0.2113248654051871, 0.7886751345948129
    };
    *wght = w;
    *absc = a;
    return 2;
}

int QRule_lin_3 (const double **wght, const double **absc)
{
    // Region: Line
    // Degree: 3
    // Points: 3

    static const double w[3] = {
        0.277777777777777778, 0.444444444444444444, 0.277777777777777778
    };
    static const double a[3] = {
        0.1127016653792583, 0.5, 0.8872983346207416
    };
    *wght = w;
    *absc = a;
    return 3;
}

int QRule_lin_4 (const double **wght, const double **absc)
{
    // Region: Line
    // Degree: 4
    // Points: 4

    static const double w[4] = {
        0.173927422568727, 0.3260725774312729,
	0.3260725774312729, 0.173927422568727
    };
    static const double a[4] = {
        0.06943184420297372, 0.3300094782075718,
	0.6699905217924282, 0.9305681557970262
    };
    *wght = w;
    *absc = a;
    return 4;
}

int QRule_lin_5 (const double **wght, const double **absc)
{
    // Region: Line
    // Degree: 5
    // Points: 5

    static const double w[5] = {
        0.1184634425280947, 0.2393143352496831, 0.2844444444444445,
	0.2393143352496831, 0.1184634425280947
    };
    static const double a[5] = {
        0.04691007703066802, 0.2307653449471585, 0.5,
	0.7692346550528415, 0.9530899229693318
    };
    *wght = w;
    *absc = a;
    return 5;
}

int QRule_lin_6 (const double **wght, const double **absc)
{
    // Region: Line
    // Degree: 6
    // Points: 6

    static const double w[6] = {
        0.08566224618958508, 0.1803807865240693, 0.2339569672863457,
	0.2339569672863457, 0.1803807865240693, 0.08566224618958508
    };
    static const double a[6] = {
        0.03376524289842397, 0.1693953067668678, 0.3806904069584016,
	0.6193095930415985, 0.8306046932331322, 0.966234757101576
    };
    *wght = w;
    *absc = a;
    return 6;
}
