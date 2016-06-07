// ==========================================================================
// TOAST v.15                                       (c) Martin Schweiger 1999
// Library: libfe
// File: tri_qr
//
// Quadrature rules over the local triangle, defined by
// xi >= 0 && eta >= 0 && xi+eta <= 1
// with vertices at (0,0), (1,0), (0,1)
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include "felib.h"
#include "tri_qr.h"

int QRule_tri_1_1 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 1
    // Points: 1
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-1-1s.html

    static const double w[1] = {
        0.5
    };
    static const Point2D a[1] = {
        Point2D(0.33333333333333333, 0.33333333333333333)
    };
    *wght = w;
    *absc = a;
    return 1;
}

int QRule_tri_2_3 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 2
    // Points: 3
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-2-3as.html

    static const double w[3] = {
        0.16666666666666667, 0.16666666666666667, 0.16666666666666667
    };
    static const Point2D a[3] = {
        Point2D(0.5,0), Point2D(0,0.5), Point2D(0.5,0.5)
    };
    *wght = w;
    *absc = a;
    return 3;
}

int QRule_tri_3_6 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 3
    // Points: 6
    // Ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-3-6cs.html

    static const double w[6] = {
        0.016666666666666667, 0.016666666666666667, 0.016666666666666667,
	0.15, 0.15, 0.15
    };
    static const Point2D a[6] = {
        Point2D(0.5,0), Point2D(0,0.5), Point2D(0.5,0.5),
	Point2D(0.16666666666666667, 0.66666666666666667),
	Point2D(0.66666666666666667, 0.16666666666666667),
	Point2D(0.16666666666666667, 0.16666666666666667)
    };
    *wght = w;
    *absc = a;
    return 6;
}

int QRule_tri_4_6 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 4
    // Points: 6
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-4-6as.html

    static const double w[6] = {
        0.054975871827660933, 0.054975871827660933, 0.054975871827660933,
	0.11169079483900573,  0.11169079483900573,  0.11169079483900573
    };
    static const Point2D a[6] = {
        Point2D(0.091576213509770743, 0.81684757298045851),
	Point2D(0.81684757298045851, 0.091576213509770743),
	Point2D(0.091576213509770743, 0.091576213509770743),
	Point2D(0.44594849091596488, 0.1081030181680702),
	Point2D(0.1081030181680702, 0.44594849091596488),
	Point2D(0.44594849091596488, 0.44594849091596488)
    };
    *wght = w;
    *absc = a;
    return 6;
}

int QRule_tri_4_9 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 4
    // Points: 9
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-4-9al.html

    static const double w[9] = {
        0.010270067672966681, 0.010270067672966681, 0.010270067672966681,
	0.030987749434133571, 0.030987749434133571, 0.030987749434133571,
	0.12540884955956641,  0.12540884955956641,  0.12540884955956641
    };
    static const Point2D a[9] = {
        Point2D(0,0), Point2D(1,0), Point2D(0,1),
	Point2D(0.5,0), Point2D(0,0.5), Point2D(0.5,0.5),
	Point2D(0.18858048469644504, 0.62283903060711),
	Point2D(0.62283903060711, 0.18858048469644504),
	Point2D(0.18858048469644504, 0.18858048469644504)
    };
    *wght = w;
    *absc = a;
    return 9;
}

int QRule_tri_5_7 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 5
    // Points: 7
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-5-7s.html

    static const double w[7] = {
        0.1125,
	0.062969590272413576, 0.062969590272413576, 0.062969590272413576,
	0.066197076394253090, 0.066197076394253090, 0.066197076394253090
    };
    static const Point2D a[7] = {
        Point2D(0.33333333333333333, 0.33333333333333333),
	Point2D(0.10128650732345633, 0.7974269853530873),
	Point2D(0.7974269853530873, 0.10128650732345633),
	Point2D(0.10128650732345633, 0.10128650732345633),
	Point2D(0.47014206410511508, 0.0597158717897698),
	Point2D(0.0597158717897698, 0.47014206410511508),
	Point2D(0.47014206410511508, 0.47014206410511508)
    };
    *wght = w;
    *absc = a;
    return 7;
}

int QRule_tri_6_12 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 6
    // Points: 12
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-6-12as.html

    static const double w[12] = {
        0.025422453185103408, 0.025422453185103408, 0.025422453185103408,
	0.058393137863189683, 0.058393137863189683, 0.058393137863189683,
	0.041425537809186787, 0.041425537809186787, 0.041425537809186787,
	0.041425537809186787, 0.041425537809186787, 0.041425537809186787
    };
    static const Point2D a[12] = {
        Point2D(0.063089014491502228, 0.87382197101699554),
	Point2D(0.87382197101699554, 0.063089014491502228),
	Point2D(0.063089014491502228, 0.063089014491502228),
	Point2D(0.24928674517091042, 0.5014265096581792),
	Point2D(0.5014265096581792, 0.24928674517091042),
	Point2D(0.24928674517091042, 0.24928674517091042),
	Point2D(0.31035245103378440, 0.6365024991213987),
	Point2D(0.6365024991213987, 0.31035245103378440),
	Point2D(0.053145049844816947, 0.6365024991213987),
	Point2D(0.6365024991213987, 0.053145049844816947),
	Point2D(0.053145049844816947, 0.31035245103378440),
	Point2D(0.31035245103378440, 0.053145049844816947)
    };
    *wght = w;
    *absc = a;
    return 12;
}

int QRule_tri_7_12 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 7
    // Points: 12
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-7-12s.html

    static const double w[12] = {
        0.026517028157436251, 0.026517028157436251, 0.026517028157436251,
	0.043881408714446055, 0.043881408714446055, 0.043881408714446055,
	0.028775042784981585, 0.028775042784981585, 0.028775042784981585,
	0.067493187009802774, 0.067493187009802774, 0.067493187009802774
    };
    static const Point2D a[12] = {
        Point2D (0.067517867073916085, 0.87009986783168180),
	Point2D (0.062382265094402118, 0.067517867073916085),
	Point2D (0.87009986783168180, 0.062382265094402118),
	Point2D (0.32150249385198182, 0.6232720494910916),
	Point2D (0.055225456656926611, 0.32150249385198182),
	Point2D (0.6232720494910916, 0.055225456656926611),
	Point2D (0.66094919618673565, 0.3047265008681672),
	Point2D (0.034324302945097146, 0.66094919618673565),
	Point2D (0.3047265008681672, 0.034324302945097146),
	Point2D (0.27771616697639178, 0.2064414986700165),
	Point2D (0.51584233435359177, 0.27771616697639178),
	Point2D (0.2064414986700165, 0.51584233435359177)
    };
    *wght = w;
    *absc = a;
    return 12;
}

int QRule_tri_7_15 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 7
    // Points: 15
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-7-15s.html

    static const double w[15] = {
        0.026538900895116205, 0.026538900895116205, 0.026538900895116205,
	0.035426541846066783, 0.035426541846066783, 0.035426541846066783,
	0.035426541846066783, 0.035426541846066783, 0.035426541846066783,
	0.034637341039708446, 0.034637341039708446, 0.034637341039708446,
	0.034637341039708446, 0.034637341039708446, 0.034637341039708446
    };
    static const Point2D a[15] = {
        Point2D (0.064930513159164863, 0.87013897368167027),
	Point2D (0.87013897368167027, 0.064930513159164863),
	Point2D (0.064930513159164863, 0.064930513159164863),
	Point2D (0.51703993906932294, 0.2845755842491704),
	Point2D (0.2845755842491704, 0.51703993906932294),
	Point2D (0.19838447668150671, 0.2845755842491704),
	Point2D (0.2845755842491704, 0.19838447668150671),
	Point2D (0.19838447668150671, 0.51703993906932294),
	Point2D (0.51703993906932294, 0.19838447668150671),
	Point2D (0.043863471792372471, 0.3135591843849315),
	Point2D (0.3135591843849315, 0.043863471792372471),
	Point2D (0.64257734382269602, 0.3135591843849315),
	Point2D (0.3135591843849315, 0.64257734382269602),
	Point2D (0.64257734382269602, 0.043863471792372471),
	Point2D (0.043863471792372471, 0.64257734382269602)
    };
    *wght = w;
    *absc = a;
    return 15;
}

int QRule_tri_9_19 (const double **wght, const Point **absc)
{
    // Region: Triangle
    // Degree: 9
    // Points: 19
    // ref: www.cs.kuleuven.ac.be/~nines/research/ecf/Formules/t2-9-19s.html

    static const double w[19] = {
        0.048567898141399416,
	0.015667350113569535, 0.015667350113569535, 0.015667350113569535,
	0.038913770502387139, 0.038913770502387139, 0.038913770502387139,
	0.039823869463605126, 0.039823869463605126, 0.039823869463605126,
	0.012788837829349015, 0.012788837829349015, 0.012788837829349015,
	0.021641769688644688, 0.021641769688644688, 0.021641769688644688,
	0.021641769688644688, 0.021641769688644688, 0.021641769688644688
    };
    static const Point2D a[19] = {
        Point2D(0.33333333333333333, 0.33333333333333333),
	Point2D(0.48968251919873762, 0.0206349616025248),
	Point2D(0.0206349616025248, 0.48968251919873762),
	Point2D(0.48968251919873762, 0.48968251919873762),
	Point2D(0.43708959149293663, 0.1258208170141267),
	Point2D(0.1258208170141267, 0.43708959149293663),
	Point2D(0.43708959149293663, 0.43708959149293663),
	Point2D(0.18820353561903273, 0.6235929287619345),
	Point2D(0.6235929287619345, 0.18820353561903273),
	Point2D(0.18820353561903273, 0.18820353561903273),
	Point2D(0.044729513394452709, 0.91054097321109458),
	Point2D(0.91054097321109458, 0.044729513394452709),
	Point2D(0.044729513394452709, 0.044729513394452709),
	Point2D(0.036838412054736283, 0.2219629891607657),
	Point2D(0.2219629891607657, 0.036838412054736283),
	Point2D(0.74119859878449802, 0.2219629891607657),
	Point2D(0.2219629891607657, 0.74119859878449802),
	Point2D(0.74119859878449802, 0.036838412054736283),
	Point2D(0.036838412054736283, 0.74119859878449802)
    };
    *wght = w;
    *absc = a;
    return 19;
}