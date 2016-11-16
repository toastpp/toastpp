// ==========================================================================
// Module libfe
// File wdg6.cc
// Definition of class Wedge6
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"

static const double sqrt3    = sqrt(3.0);
static const double sqrt3_05 = 0.5 * sqrt(3.0);
static const double sqrt3_2  = 2.0 * sqrt(3.0);

// ***************************************************************************

Wedge6::Wedge6 (const Wedge6 &el)
: Element_Unstructured_3D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Wedge6::Copy ()
{
    return new Wedge6(*this);
}

void Wedge6::Initialise (const NodeList &nlist)
{
    Element_Unstructured_3D::Initialise (nlist);
    intff.Zero (6);
    intff = ComputeIntFF (nlist);
}

// ***************************************************************************

int Wedge6::nSideNode (int side) const
{
    dASSERT(side >= 0 && side < 5, "Invalid side number.");
    if (side <= 2) return 4;
    return 3;
}

// ***************************************************************************

int Wedge6::SideNode (int side, int node) const
{
    static const int SN[5][4] =
	{{0,1,4,3},{0,2,5,3},{1,2,5,4},{0,1,2,-1},{3,4,5,-1}};
    dASSERT(side >= 0 && side < 5, "Invalid side parameter.");
    dASSERT(node >= 0 && node < 4 && SN[side][node] >= 0,
	    "Invalid node param.");
    return SN[side][node];
}

// ***************************************************************************
// CAUTION: This works only if the top and bottom surface of the wedge are
// parallel to the xy plane, and if the other surfaces are perpendicular to
// the xy plane!

Point Wedge6::Local (const NodeList& nlist, const Point& glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid dimension of global point.");

    double x1, y1, z1, x2, y2, z2, x3, y3, x, y, z;
    static RDenseMatrix b(2,2);		// make static to avoid construction
    static RDenseMatrix inv(2,2);
    static RVector c(2), res(2);
    Point loc(3);

    x1 = nlist[Node[0]][0];  y1 = nlist[Node[0]][1];  z1 = nlist[Node[0]][2];
    x2 = nlist[Node[1]][0];  y2 = nlist[Node[1]][1];  z2 = nlist[Node[5]][2];
    x3 = nlist[Node[2]][0];  y3 = nlist[Node[2]][1];

    x = glob[0];  y = glob[1];  z = glob[2];

    b(0,0) = (2.0*x1-x2-x3)/3.0;  b(0,1) = (x3-x2)/sqrt3;
    b(1,0) = (2.0*y1-y2-y3)/3.0;  b(1,1) = (y3-y2)/sqrt3;
    c[0] = x - (x1+x2+x3)/3.0; c[1] = y - (y1+y2+y3)/3.0;

    inv = inverse (b);
    res = inv * c;
    loc[0] = res[0], loc[1] = res[1];
    loc[2] = 2.0 * (z-z1) / (z2-z1) - 1.0;
    return loc;
}

// ***************************************************************************

Point Wedge6::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 6, "Invalid value for argument 'node'.");
    Point n(3);
    switch (node) {
	case 0: n[0] = 1;    n[1] = 0;         n[2] = -1; break;
	case 1: n[0] = -0.5; n[1] = -sqrt3_05; n[2] = -1; break;
	case 2: n[0] = -0.5; n[1] =  sqrt3_05; n[2] = -1; break;
	case 3: n[0] = 1;    n[1] = 0;         n[2] =  1; break;
	case 4: n[0] = -0.5; n[1] = -sqrt3_05; n[2] =  1; break;
	case 5: n[0] = -0.5; n[1] =  sqrt3_05; n[2] =  1; break;
    }
    return n;
}

// ***************************************************************************

RVector Wedge6::DirectionCosine (int side, RDenseMatrix& jacin)
{
    dASSERT(jacin.nCols() == 3 && jacin.nRows() == 3,
	"Invalid dimension of matrix jacin.");

    double amod;
    RVector work(3);

    switch (side) {
	case 0: work[0]= 0.5, work[1]=-sqrt3_05, work[2]= 0.0; break;
	case 1: work[0]= 0.5, work[1]= sqrt3_05, work[2]= 0.0; break;
	case 2: work[0]=-1.0, work[1]= 0.0,      work[2]= 0.0; break;
	case 3: work[0]= 0.0, work[1]= 0.0,      work[2]=-1.0; break;
	case 4: work[0]= 0.0, work[1]= 0.0,      work[2]= 1.0; break;
    }
    RVector cosin = jacin * work;
    amod = sqrt (cosin[0]*cosin[0] + cosin[1]*cosin[1] + cosin[2]*cosin[2]);
    return cosin/amod;
};

const RVector &Wedge6::LNormal (int side) const
{
    static const RVector lnm0 = RVector (3, "0.5 -0.866025404 0");
    static const RVector lnm1 = RVector (3, "0.5  0.866025404 0");
    static const RVector lnm2 = RVector (3, "-1 0 0");
    static const RVector lnm3 = RVector (3, "0 0 -1");
    static const RVector lnm4 = RVector (3, "0 0 1");
    static const RVector *lnm[5] = {
      &lnm0, &lnm1, &lnm2, &lnm3, &lnm4
    };
    dASSERT(side >= 0 && side < 5, "Argument 1 index out of range");
    return *lnm[side];
}

// ***************************************************************************
// only valid for upright wedges

double Wedge6::ComputeSize (const NodeList &nlist) const
{
    double x1, y1, x2, y2, x3, y3, a, b, c, s, h;
    x1 = nlist[Node[0]][0];  y1 = nlist[Node[0]][1];
    x2 = nlist[Node[1]][0];  y2 = nlist[Node[1]][1];
    x3 = nlist[Node[2]][0];  y3 = nlist[Node[2]][1];
    a = hypot (x1-x2, y1-y2);
    b = hypot (x2-x3, y2-y3);
    c = hypot (x3-x1, y3-y1);
    s = 0.5 * (a+b+c);
    h = fabs (nlist[Node[0]][2] - nlist[Node[3]][2]);
    return sqrt (s * (s-a) * (s-b) * (s-c)) * h;
}

// ***************************************************************************

bool Wedge6::LContains (const Point& loc, bool pad) const
{
    dASSERT(loc.Dim() == 3, "Invalid dimension of local point.");
    if (pad) {
        ERROR_UNDEF;
	return false;
    } else {
        if (loc[2] < -1.0 || loc[2] > 1.0) return false;
	if (loc[0] < -0.5 || loc[0] > 1.0) return false;
	return fabs (loc[1]) <= (1.0 - loc[0]) * sqrt3_05;
    }
}

// ***************************************************************************

RVector Wedge6::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Evaluation point must be 3D");
    RVector fun(6);
    fun[0] = (1.0+2.0*loc[0]) * (1.0-loc[2]) / 6.0;
    fun[1] = (1.0-loc[0]-sqrt3*loc[1]) * (1.0-loc[2]) / 6.0;
    fun[2] = (1.0-loc[0]+sqrt3*loc[1]) * (1.0-loc[2]) / 6.0;
    fun[3] = (1.0+2.0*loc[0]) * (1.0+loc[2]) / 6.0;
    fun[4] = (1.0-loc[0]-sqrt3*loc[1]) * (1.0+loc[2]) / 6.0;
    fun[5] = (1.0-loc[0]+sqrt3*loc[1]) * (1.0+loc[2]) / 6.0;
    return fun;
}

// ***************************************************************************

RDenseMatrix Wedge6::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Evaluation point must be 3D");
    RDenseMatrix der(3,6);
    der(0,0) =  (1.0-loc[2]) / 3.0;
    der(1,0) =  0.0;
    der(2,0) = -(1.0+2.0*loc[0]) / 6.0;
    der(0,1) = -(1.0-loc[2]) / 6.0;
    der(1,1) = -(1.0-loc[2]) / sqrt3_2;
    der(2,1) = -(1.0-loc[0]-sqrt3*loc[1]) / 6.0;
    der(0,2) = -(1.0-loc[2]) / 6.0;
    der(1,2) =  (1.0-loc[2]) / sqrt3_2;
    der(2,2) = -(1.0-loc[0]+sqrt3*loc[1]) / 6.0;
    der(0,3) =  (1.0+loc[2]) / 3.0;
    der(1,3) =  0.0;
    der(2,3) =  (1.0+2.0*loc[0]) / 6.0;
    der(0,4) = -(1.0+loc[2]) / 6.0;
    der(1,4) = -(1.0+loc[2]) / sqrt3_2;
    der(2,4) =  (1.0-loc[0]-sqrt3*loc[1]) / 6.0;
    der(0,5) = -(1.0+loc[2]) / 6.0;
    der(1,5) =  (1.0+loc[2]) / sqrt3_2;
    der(2,5) =  (1.0-loc[0]+sqrt3*loc[1]) / 6.0;
    return der;
}

// ***************************************************************************

RSymMatrix Wedge6::IntFF () const {
    return intff;
}

// ***************************************************************************

double Wedge6::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6);
    return intff(i,j);
}

// ***************************************************************************

RSymMatrix Wedge6::ComputeIntFF (const NodeList &nlist) const
{
    const double *wght;
    const Point *absc;
    RSymMatrix ff(6);
    RDenseMatrix geom = Elgeom (nlist);
    int nqp = QuadRule (&wght, &absc);
    for (int iquad = 0; iquad < nqp; iquad++) {
	RVector fun = LocalShapeF (absc[iquad]);
	RDenseMatrix lder = LocalShapeD (absc[iquad]);
	RDenseMatrix jac = lder * geom;
	double djac = fabs (det (jac));
	for (int i = 0; i < 6; i++)
	    for (int j = 0; j <= i; j++)
		ff(i,j) += fun[i]*fun[j] * djac * wght[iquad];
    }
    return ff;
}

// ***************************************************************************

RSymMatrix Wedge6::ComputeIntDD (const NodeList &nlist) const
{
    const double *wght;
    const Point *absc;
    RSymMatrix dd(6);
    RDenseMatrix geom = Elgeom (nlist);
    int nqp = QuadRule (&wght, &absc);
    for (int iquad = 0; iquad < nqp; iquad++) {
	RDenseMatrix lder = LocalShapeD (absc[iquad]);
	RDenseMatrix jac = lder * geom;
	RDenseMatrix gder = inverse (jac) * lder;
	double djac = fabs (det (jac));
	dd += ATA (gder) * (djac * wght[iquad]);
    }
    return dd;
}

// ***************************************************************************

int Wedge6::QuadRule (const double **wght, const Point **absc) const
{
    static const double w = 3.0 * sqrt3_05;

    static const double local_wght[8] =
	{ 3.0/8.0*w,  3.0/8.0*w,  1.0/24.0*w, 1.0/24.0*w,
	  1.0/24.0*w, 1.0/24.0*w, 1.0/24.0*w, 1.0/24.0*w };
    static const Point3D local_absc[8] =
	{ Point3D(0,0,1),             Point3D(0,0,-1),
	  Point3D(1,0,-1),            Point3D(-0.5,sqrt3_05,-1),
	  Point3D(-0.5,-sqrt3_05,-1), Point3D(1,0,1),
	  Point3D(-0.5,sqrt3_05,1),   Point3D(-0.5,-sqrt3_05,1) };

    *wght = local_wght;
    *absc = local_absc;
    return 8;
}

// ***************************************************************************
#ifdef UNDEF

double Wedge6::UnitLength (Point& loc, RDenseMatrix& geom, int side)
{
    dASSERT(side >= 0 && side < 5, "Invalid side number.");
    dASSERT(geom.Dim(ROW) == 6 && geom.Dim(COL) == 3,
	"Invalid geometry matrix dimension.");

    double g1, g2, g3;

    RDenseMatrix lder = LocalShapeD (loc);
    RDenseMatrix jac  = lder * geom;
    switch (side) {
	case 0:
	case 1:
	    xERROR("Sides 0 and 1 not implemented yet.");
	    break;
	case 2:  // eta=const
	    g1 = jac[0][1]*jac[2][2] - jac[0][2]*jac[2][1];
	    g2 = jac[0][2]*jac[2][0] - jac[0][0]*jac[2][2];
	    g3 = jac[0][0]*jac[2][1] - jac[0][1]*jac[2][0];
	    break;
	case 3:  // bottom, zeta=const
	case 4:  // top, zeta=const
	    g1 = jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1];
	    g2 = jac[0][2]*jac[1][0] - jac[0][0]*jac[1][2];
	    g3 = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
	    break;
    }
    return sqrt (g1*g1 + g2*g2 + g3*g3);
}

// ***************************************************************************

void Wedge6::ConvLineQuadrature (Point** absc, double* labsc,
    int nqp, int side, double* coeff)
{
    ERROR_UNDEF;
}

#endif

// ***************************************************************************

int Wedge6::Intersection (const Point &p1, const Point &p2, Point *s,
    bool add_endpoints, bool boundary_only)
{
    ERROR_UNDEF;
    return 0; // dummy
}
