// ==========================================================================
// Module libfe
// File tri3.cc
// Definition of class Triangle3old
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"
#include "lin_qr.h"
#include "tri_qr.h"

using namespace std;

// local prototypes

static void Integrate_x_Cosine (double d, double dmax, double sigma, double s,
    double &res0, double &res1);
static void Integrate_u_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_u1);
static Point Triangle_SurfToLocal (int side, const Point &p);
static int Triangle_GetLocalSubsampleAbsc (const Point *&absc);
static int Triangle_GetBndSubsampleAbsc (int side, const Point *&absc);

// some global constants

static const double sqrt3    = sqrt(3.0);
static const double sqrt3_05 = 0.5 *sqrt3;

static bool subsampling_initialised = false;
static const int nsample_lin = NSUBSAMPLE; // from toastdef.h
static Point absc_bndsample[3][nsample_lin];
static const RSymMatrix sym_intff = RSymMatrix (3,
   "2 \
    1 2 \
    1 1 2") * (1.0/12.0);

static const RDenseMatrix full_intff = RDenseMatrix (3, 3,
   "2 1 1 \
    1 2 1 \
    1 1 2") * (1.0/12.0);

// index set for retrieving IntFFF(i,j,k)
static const int intfff_index[3][3][3] = {
  {{2,1,1},{1,1,0},{1,0,1}},{{1,1,0},{1,2,1},{0,1,1}},{{1,0,1},{0,1,1},{1,1,2}}
};
static const double intfff_scale[3] = {1.0/60.0, 1.0/30.0, 1.0/10.0};

Triangle3old::Triangle3old ()
{
    Node = new int[nNode()];
    intfd_0 = 0, intfd_1 = 0;
}

Triangle3old::Triangle3old (const Triangle3old& el): Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
    intfd_0 = 0, intfd_1 = 0;
}

Triangle3old::~Triangle3old ()
{
    delete []Node;
    if (intfd_0) delete []intfd_0;
    if (intfd_1) delete []intfd_1;
}

Element *Triangle3old::Copy ()
{
    return new Triangle3old(*this);
}

void Triangle3old::Initialise (const NodeList& nlist)
{
#ifndef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
#endif // otherwise class members
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];

    jac.New (2,2);
    jac(0,0) = n1x-n0x;  jac(0,1) = n1y-n0y;
    jac(1,0) = n2x-n0x;  jac(1,1) = n2y-n0y;
    djac = n0x*(n1y-n2y) + n1x*(n2y-n0y) + n2x*(n0y-n1y); // determinant

    Element_Unstructured_2D::Initialise (nlist);
    dASSERT(size > 0, "Element size not positive");

#ifdef TRI3_STORE_INTFF
    intff.New(3);
    intff = sym_intff * size;
#endif

#ifdef TRI3_STORE_INTDD
    intdd.New(3);
    intdd = ComputeIntDD (nlist);
#endif

    if (!subsampling_initialised) {
        Point loc(2);
	double bloc;
	int i, j;
        for (i = 0; i < nsample_lin; i++) {
	    bloc = (i+0.5)/(double)nsample_lin;
	    for (j = 0; j < 3; j++)
		absc_bndsample[j][i].New(2);
	    absc_bndsample[0][i][0] = bloc;
	    absc_bndsample[0][i][1] = 0.0;
	    absc_bndsample[1][i][0] = bloc;
	    absc_bndsample[1][i][1] = bloc;
	    absc_bndsample[2][i][0] = 0.0;
	    absc_bndsample[2][i][1] = bloc;
	}
	subsampling_initialised = true;
    }
}

int Triangle3old::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, "Argument 1 index out of range");
    dASSERT(node >= 0 && node < 2, "Argument 2  index out of range");
    static int SN[3][2] = {{0,1}, {1,2}, {2,0}};
    return SN[side][node];
}

double Triangle3old::SideSize (int sd, const NodeList &nlist) const
{
    return nlist[Node[SideNode(sd,0)]].Dist (nlist[Node[SideNode(sd,1)]]);
}

Point Triangle3old::Local (const NodeList &nlist, const Point& glob) const
{
    dASSERT(glob.Dim() == 2, "Argument 2 dimension must be 2");
#ifndef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
#endif
    static Point loc(2);
    loc[0] = ((glob[0]-n0x)*(n2y-n0y) - (glob[1]-n0y)*(n2x-n0x)) / djac;
    loc[1] = ((glob[1]-n0y)*(n1x-n0x) - (glob[0]-n0x)*(n1y-n0y)) / djac;
    return loc;
}

Point Triangle3old::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 3, "Argument 1 index out of range");
    Point n(2);  // note: initialised to zero
    if (node == 1)      n[0] = 1.0;
    else if (node == 2) n[1] = 1.0;
    return n;
}

Point Triangle3old::SurfToLocal (int side, const Point &p) const
{
    return Triangle_SurfToLocal (side, p);
}

RVector Triangle3old::DirectionCosine (int side, RDenseMatrix &jacin)
{
    static const double local_dc[3][2] = {
        {0,-1},{0.70710678,0.70710678},{-1,0}
    };

    dASSERT(side >= 0 && side < 3, "Argument 1 index out of range");
    dASSERT(jacin.nCols() == 2 && jacin.nRows() == 2,
	"Argument 2 invalid dimension");

    RVector cosin(2);
    cosin[0] = jacin(0,0)*local_dc[side][0] + jacin(0,1)*local_dc[side][1];
    cosin[1] = jacin(1,0)*local_dc[side][0] + jacin(1,1)*local_dc[side][1];
    return cosin / length (cosin);
};

const RVector &Triangle3old::LNormal (int side) const
{
    static const RVector lnm0 = RVector (2, "0 -1");
    static const RVector lnm1 = RVector (2, "0.7071067814 0.7071067814");
    static const RVector lnm2 = RVector (2, "-1 0");
    static const RVector *lnm[3] = {
      &lnm0, &lnm1, &lnm2
    };
    dASSERT(side >= 0 && side < 3, "Argument 1 index out of range");
    return *lnm[side];
}

double Triangle3old::ComputeSize (const NodeList &nlist) const
{
    return -0.5 * djac;   // convert local to global size
}

bool Triangle3old::LContains (const Point& loc, bool pad) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[1]+EPS >= 0.0 &&
	    loc[0]+loc[1]-EPS <= 1.0);
    } else {
	return (loc[0] >= 0.0 && loc[1] >= 0.0 && loc[0]+loc[1] <= 1.0);
    }
}

bool Triangle3old::GContains (const Point& glob, const NodeList& nlist) const
{
    dASSERT(glob.Dim() == 2, "Argument 1 invalid dimension");
    double xx = glob[0], yy = glob[1];

    // check bounding box
    if (xx < bbmin[0] || xx > bbmax[0] || yy < bbmin[1] || yy > bbmax[1])
        return false;

    double x0, x1, x2, y0, y1, y2, y0r, yyr, fac;
    const double EPS = 1e-10;

    for (int i = 0; i < 3; i++) {
        x0 = nlist[Node[i]][0],       y0 = nlist[Node[i]][1];
	x1 = nlist[Node[(i+1)%3]][0], y1 = nlist[Node[(i+1)%3]][1];
	x2 = nlist[Node[(i+2)%3]][0], y2 = nlist[Node[(i+2)%3]][1];
	if (fabs (x1-x2) < EPS) {
	    if ((x0 < x1 && xx > x1) || (x0 > x1 && xx < x1)) return false;
	} else {
	    fac = (y2-y1)/(x2-x1);
	    y0r = (x0-x1)*fac + y1;
	    yyr = (xx-x1)*fac + y1;
	    if ((y0 < y0r && yy > yyr) || (y0 > y0r && yy < yyr)) return false;
	}
    }
    return true;
    /*
    double x1, y1, x2, y2, x3, y3, tmp;
    const double EPS = 1e-10;

    x1 = nlist[Node[0]][0];  y1 = nlist[Node[0]][1];
    x2 = nlist[Node[1]][0];  y2 = nlist[Node[1]][1];
    x3 = nlist[Node[2]][0];  y3 = nlist[Node[2]][1];

    if (x1>x2) { tmp=x1; x1=x2; x2=tmp; tmp=y1; y1=y2; y2=tmp; }
    if (x1>x3) { tmp=x1; x1=x3; x3=tmp; tmp=y1; y1=y3; y3=tmp; }
    if (x2<x3) { tmp=x2; x2=x3; x3=tmp; tmp=y2; y2=y3; y3=tmp; }
    if (y3 < y1+(y2-y1)*(x3-x1)/(x2-x1)) { y1=-y1; y2=-y2; y3=-y3; yy=-yy; }
    if ((yy >= y1 + (y2-y1)*(xx-x1)/(x2-x1)) &&
	(fabs(x3-x2)<EPS || yy <= y2 + (y3-y2)*(xx-x2)/(x3-x2)) &&
	(fabs(x3-x1)<EPS || yy <= y1 + (y3-y1)*(xx-x1)/(x3-x1))) {
        return true;
    }
    return false;
    */
}

RVector Triangle3old::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static RVector fun(3);
    fun[0] = 1.0 - loc[0] - loc[1];
    fun[1] = loc[0];
    fun[2] = loc[1];
    return fun;
}

RDenseMatrix Triangle3old::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static const RDenseMatrix der (2, 3, "-1 1 0   -1 0 1");
    return der;
}

RVector Triangle3old::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{
    dASSERT(glob.Dim() == 2, "Argument 2 invalid dimension");
#ifndef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
#endif
    double scale = 1.0/djac;
    static RVector fun(3);
    fun[0] = scale * (n1x*n2y-n2x*n1y + (n1y-n2y)*glob[0] + (n2x-n1x)*glob[1]);
    fun[1] = scale * (n2x*n0y-n0x*n2y + (n2y-n0y)*glob[0] + (n0x-n2x)*glob[1]);
    fun[2] = scale * (n0x*n1y-n1x*n0y + (n0y-n1y)*glob[0] + (n1x-n0x)*glob[1]);
    //dASSERT(fun[0] >= -1e-8, Shape function 0 negative);
    //dASSERT(fun[1] >= -1e-8, Shape function 1 negative);
    //dASSERT(fun[2] >= -1e-8, Shape function 2 negative);
    return fun;
}

RDenseMatrix Triangle3old::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Argument 2 invalid dimension");
#ifndef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
#endif
    double scale = 1.0/djac;
    RDenseMatrix der (2,3);
    der(0,0) = (n1y-n2y) * scale;
    der(0,1) = (n2y-n0y) * scale;
    der(0,2) = (n0y-n1y) * scale;
    der(1,0) = (n2x-n1x) * scale;
    der(1,1) = (n0x-n2x) * scale;
    der(1,2) = (n1x-n0x) * scale;
    return der;
}

int Triangle3old::QuadRule (int order, const double **wght, const Point **absc)
    const
{
    switch (order) {
    case 1: return QRule_tri_1_1 (wght, absc);
    case 2: return QRule_tri_2_3 (wght, absc);
    case 3: return QRule_tri_3_6 (wght, absc);
    case 4: return QRule_tri_4_6 (wght, absc);
    case 6: return QRule_tri_6_12 (wght, absc);
    default: ERROR_UNDEF; return 0;
    }
}

double Triangle3old::IntF (int i) const
{
    static const double third = 1.0/3.0;
    return size*third;
}

RSymMatrix Triangle3old::IntFF () const {
#ifdef TRI3_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Triangle3old::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3);
#ifdef TRI3_STORE_INTFF
    return intff(i,j);
#else
    return full_intff(i,j) * size;
#endif
}

double Triangle3old::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3 && k >= 0 && k < 3);
    return size * intfff_scale[intfff_index[i][j][k]];
}

#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

void Triangle3old::IntFFF (double &iii, double &iij, double &ijk) const
{
    iii = size * intfff_scale[2];
    iij = size * intfff_scale[1];
    ijk = size * intfff_scale[0];
}
#endif

RSymMatrix Triangle3old::IntPFF (const RVector &P) const
{
    static RSymMatrix pff(3);
    static const double i30 = 1.0/30.0;
    double fac = size*i30;
    double p0 = P[Node[0]], p1 = P[Node[1]], p2 = P[Node[2]];
    pff(0,0) = fac * (3.0*p0 + p1 + p2);
    pff(1,0) = fac * (p0 + p1 + 0.5*p2);
    pff(2,0) = fac * (p0 + 0.5*p1 + p2);
    pff(1,1) = fac * (p0 + 3.0*p1 + p2);
    pff(2,1) = fac * (0.5*p0 + p1 + p2);
    pff(2,2) = fac * (p0 + p1 + 3.0*p2);
    return pff;
}

double Triangle3old::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3);
    if (i < j) { int tmp = i; i = j; j = tmp; }
    static const double i30 = 1.0/30.0;
    switch (i) {
    case 0:
        return size*i30 * (3.0*P[Node[0]] + P[Node[1]] + P[Node[2]]); 
    case 1:
	switch (j) {
	case 0:
	    return size*i30 * (P[Node[0]] + P[Node[1]] + 0.5*P[Node[2]]);
	case 1:
	    return size*i30 * (P[Node[0]] + 3.0*P[Node[1]] + P[Node[2]]);
	}
    case 2:
	switch (j) {
	case 0:
	    return size*i30 * (P[Node[0]] + 0.5*P[Node[1]] + P[Node[2]]);
	case 1:
	    return size*i30 * (0.5*P[Node[0]] + P[Node[1]] + P[Node[2]]);
	case 2:
	    return size*i30 * (P[Node[0]] + P[Node[1]] + 3.0*P[Node[2]]);
	}
    default:
        xERROR("Index out of range");
        return 0.0; // dummy
    }
}

RSymMatrix Triangle3old::Intdd () const
{
#ifndef TRI3_STORE_COORDS
    double b0 = 0;
    double b1 = 0;
    double b2 = 0;
    double c0 = 0;
    double c1 = 0;
    double c2 = 0;
    xERROR ("Function only available with precomputed shape parameters");
#else
    double b0 = n1y - n2y;
    double b1 = n2y - n0y;
    double b2 = n0y - n1y;
    double c0 = n2x - n1x;
    double c1 = n0x - n2x;
    double c2 = n1x - n0x;
#endif

    RSymMatrix MDD(6,6);
    MDD(0,0) = b0*b0;
    MDD(1,0) = b0*c0; MDD(1,1) = c0*c0;

    MDD(2,0) = b0*b1; MDD(2,1) = b1*c0;
    MDD(3,0) = b0*c1; MDD(3,1) = c0*c1;

    MDD(2,2) = b1*b1;
    MDD(3,2) = b1*c1; MDD(3,3) = c1*c1;

    MDD(4,0) = b0*b2; MDD(4,1) = b2*c0;
    MDD(5,0) = b0*c2; MDD(5,1) = c0*c2;

    MDD(4,2) = b1*b2; MDD(4,3) = b2*c1;
    MDD(5,2) = b1*c2; MDD(5,3) = c1*c2;

    MDD(4,4) = b2*b2;
    MDD(5,4) = b2*c2; MDD(5,5) = c2*c2;

    MDD /= 4*size;
    return MDD;
}

int Triangle3old::GetLocalSubsampleAbsc (const Point *&absc) const
{
    return Triangle_GetLocalSubsampleAbsc (absc);
}

int Triangle3old::GetBndSubsampleAbsc (int side, const Point *&absc) const
{
    dASSERT (side >= 0 && side < 3, "Argument 1 index out of range");
    return Triangle_GetBndSubsampleAbsc (side, absc);
}

void Triangle3old::ComputeIntFD (const NodeList &nlist)
{
#ifndef TRI3_STORE_COORDS
    double b0 = nlist[Node[1]][1] - nlist[Node[2]][1];
    double b1 = nlist[Node[2]][1] - nlist[Node[0]][1];
    double b2 = nlist[Node[0]][1] - nlist[Node[1]][1];
    double c0 = nlist[Node[2]][0] - nlist[Node[1]][0];
    double c1 = nlist[Node[0]][0] - nlist[Node[2]][0];
    double c2 = nlist[Node[1]][0] - nlist[Node[0]][0];
#else
    double b0 = n1y - n2y;
    double b1 = n2y - n0y;
    double b2 = n0y - n1y;
    double c0 = n2x - n1x;
    double c1 = n0x - n2x;
    double c2 = n1x - n0x;
#endif
    static const double i6 = 1.0/6.0;

    // note that Int{F_i D_j} is independent of i, so we
    // need to store only one row

    if (intfd_0) delete []intfd_0;
    if (intfd_1) delete []intfd_1;
    intfd_0 = new double[3];
    intfd_1 = new double[3];
    intfd_0[0] = b0*i6;
    intfd_1[0] = c0*i6;
    intfd_0[1] = b1*i6;
    intfd_1[1] = c1*i6;
    intfd_0[2] = b2*i6;
    intfd_1[2] = c2*i6;
}

RVector Triangle3old::IntFD (int i, int j) const
{
    dASSERT(intfd_0 && intfd_1, "FD matrix not initialised");
    dASSERT(i >= 0 && i < 3, "Parameter 1 index out of range");
    dASSERT(j >= 0 && j < 3, "Parameter 2 index out of range");

    // note that Int{F_i D_j} is independent of i

    RVector fd(2);
    fd[0] = intfd_0[j];
    fd[1] = intfd_1[j];
    return fd;
}

RSymMatrix Triangle3old::ComputeIntDD (const NodeList &nlist) const
{
    // This is the direct evaluation of:
    //    Matrix lder = LocalShapeD (Point(2));
    //    Matrix gder = inverse (jac) * lder;
    //    return ATA (gder) * size;

#ifndef TRI3_STORE_COORDS
    double b0 = nlist[Node[1]][1] - nlist[Node[2]][1];
    double b1 = nlist[Node[2]][1] - nlist[Node[0]][1];
    double b2 = nlist[Node[0]][1] - nlist[Node[1]][1];
    double c0 = nlist[Node[2]][0] - nlist[Node[1]][0];
    double c1 = nlist[Node[0]][0] - nlist[Node[2]][0];
    double c2 = nlist[Node[1]][0] - nlist[Node[0]][0];
#else
    double b0 = n1y - n2y;
    double b1 = n2y - n0y;
    double b2 = n0y - n1y;
    double c0 = n2x - n1x;
    double c1 = n0x - n2x;
    double c2 = n1x - n0x;
#endif
    static RSymMatrix dd(3);
    double scale = 0.25/size;
    dd(0,0) = scale * (b0*b0 + c0*c0);
    dd(1,0) = scale * (b1*b0 + c1*c0);
    dd(2,0) = scale * (b2*b0 + c2*c0);
    dd(1,1) = scale * (b1*b1 + c1*c1);
    dd(2,1) = scale * (b2*b1 + c2*c1);
    dd(2,2) = scale * (b2*b2 + c2*c2);
    return dd;
}

double Triangle3old::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 3, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 3, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 3: out of range");
    
    double d, dx, dy;

    int opp_nd = (sd+2) % 3;
    if (i == opp_nd || j == opp_nd) return 0.0;

    switch (sd) {
    case 0:
	dx = jac(0,0);
	dy = jac(0,1);
	break;
    case 1:
	dx = jac(1,0)-jac(0,0);
	dy = jac(1,1)-jac(0,1);
	break;
    case 2:
	dx = jac(1,0);
	dy = jac(1,1);
	break;
    }
    d  = hypot (dx, dy);
    return d / (i == j ? 3.0 : 6.0);
}

RSymMatrix Triangle3old::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(3);
    if (!bndel) return bff;  // not a boundary element -> return zero matrix

    for (int side = 0; side < nSide(); side++) {
        if (!bndside[side]) continue;  // not a boundary side
	int n0 = SideNode (side, 0);
	int n1 = SideNode (side, 1);
	double dx = nlist[Node[n0]][0] - nlist[Node[n1]][0];
	double dy = nlist[Node[n0]][1] - nlist[Node[n1]][1];
	double d  = hypot (dx, dy);
	for (int i = 0; i < 3; i++) {
	    if (i != n0 && i != n1) continue;
	    bff(i,i) += d/3.0;
	    for (int j = 0; j < i; j++)
		if (j == n0 || j == n1) bff(i,j) += d/6.0;
	}
    }
    return bff;
}

RVector Triangle3old::BndIntFX (int side, double (*func)(const Point&),
    const NodeList &nlist) const
{
    // Calculates Int [u_i(r) func(r) dr]
    // along side 'side' of the triangle
    // where func is a scalar-valued user supplied function evaluated
    // at a global coordinate r

    int p, np, n0, n1;
    double dx, dy, d, f;
    const double *wght;
    const double *absc;
    RVector bint(3);
    Point sloc(1);

    n0 = SideNode (side, 0);
    n1 = SideNode (side, 1);
    dx = nlist[Node[n0]][0] - nlist[Node[n1]][0];
    dy = nlist[Node[n0]][1] - nlist[Node[n1]][1];
    d  = hypot (dx, dy);             // length of the side

    np = QRule_lin_2 (&wght, &absc); // quadrature rule - CHECK ORDER IS OK!
    for (p = 0; p < np; p++) {
        sloc[0] = absc[p];
        Point loc   = SurfToLocal (side, sloc);    // local quadrature point
        Point glob  = Global (nlist, loc);         // global quadrature point
	RVector fun = LocalShapeF (loc);           // shape func at quad point
	f = func(glob);                   // user function at quadrature point
	             // Warning: glob is not guaranteed to sit on the surface!
	bint += fun * (d*f*wght[p]);
    }
    return bint;
}

RVector Triangle3old::BndIntFCos (int side, const Surface *surf,
    const RVector &cntcos, double sigma, double sup, const NodeList &nlist)
    const
{
    dASSERT(surf->Dimension() == 2, "Wrong surface dimension");
    Surface2D *surf2D = (Surface2D*)surf;
    double s, d0, tmp;
    bool flip;
    int n0 = Node[SideNode (side, 0)];
    int n1 = Node[SideNode (side, 1)];
    int pdim = surf2D->ParamDim();
    RVector prm0(pdim);  surf2D->Point2Param (nlist[n0], prm0);
    RVector prm1(pdim);  surf2D->Point2Param (nlist[n1], prm1);
    s  = surf2D->ChordDiff (prm1, prm0); 
    if ((flip = (s < 0))) {
        d0 = surf2D->ChordDiff (cntcos, prm1);
	s = -s;
    } else {
        d0 = surf2D->ChordDiff (cntcos, prm0);
    }
    RVector res(2);
    Integrate_x_Cosine (d0, sup, sigma, s, res[0], res[1]);
    if (flip) tmp = res[0], res[0] = res[1], res[1] = tmp;
    return res;
}

RVector Triangle3old::BndIntFCos (int side, const RVector &cntcos, double a,
    const NodeList &nlist) const
{
    // this version parametrises the surface simply as a function of
    // the geometric distance from the centre of the cosine. This
    // will cause problems for curved surfaces and large cosine support
    // radius a

    bool flip = false;
    int n0 = Node[SideNode (side, 0)];
    int n1 = Node[SideNode (side, 1)];
    double d0 = length (nlist[n0] - cntcos);
    double d1 = length (nlist[n1] - cntcos);

    // hack: if angle subtended by the nodes as seen from cosine centre > Pi/2,
    // we assume that they bracket the cosine centre and flip one of the
    // distance signs
    double cosa = dot (nlist[n0]-cntcos, nlist[n1]-cntcos) / (d0*d1);
    if (cosa < 0.0) d0 = -d0;

    RVector res(2);
    if (d0 > d1) {
	double tmp = d0; d0 = d1; d1 = tmp;
	flip = true;
    }
    Integrate_u_Cosine (0.0, a, d0, d1, res[0], res[1]);
    if (flip) {
	double tmp = res[0]; res[0] = res[1]; res[1] = tmp;
    }
    return res;
}

RVector Triangle3old::BndIntFDelta (int side, const Surface *surf,
    const RVector &pos, const NodeList &nlist) const
{
    dASSERT(surf->Dimension() == 2, "Wrong surface dimension");
    Surface2D *surf2D = (Surface2D*)surf;
    double s, d0, d1;
    int n0 = Node[SideNode (side, 0)];
    int n1 = Node[SideNode (side, 1)];
    int pdim = surf2D->ParamDim();
    RVector prm0(pdim);  surf2D->Point2Param (nlist[n0], prm0);
    RVector prm1(pdim);  surf2D->Point2Param (nlist[n1], prm1);
    d0 = surf2D->ChordDiff (pos, prm0);
    d1 = surf2D->ChordDiff (pos, prm1);
    s  = surf2D->ChordDiff (prm1, prm0);
    RVector res(2);
    if (s*d0 >= 0.0 && d0*d1 <= 0.0) { // otherwise outside support
        res[0] = -d1/s;
	res[1] =  d0/s;
    }
    return res;
}

int Triangle3old::Intersection (const Point &p1, const Point &p2,
    Point *s, bool add_endpoints, bool boundary_only)
{
    double xs, ys;
    double pmin, pmax, smin, smax;
    const double EPS = 1e-12;
    int pindx = 0;
    Point p(2);
    
    dASSERT(p1.Dim() == 2 && p2.Dim() == 2, "Points must be 2D.");
    
    // a) check whether one of the end points of the line is within the element
    if (add_endpoints) {
	if (LContains (p1)) s[pindx++]=p1;
	if (LContains (p2)) s[pindx++]=p2;
	if (pindx==2) goto Done;
    }
    
    if (boundary_only && !bndel) goto Done;
    
    // b) check whether line intersects side 0 of the triangle
    if (!boundary_only || bndside[0]) {
	if (p1[1]*p2[1] < 0.0) {
	    xs = p1[0] - p1[1]/(p2[1]-p1[1])*(p2[0]-p1[0]);
	    if (xs > -EPS && xs < 1.0+EPS) {
		p[0] = xs, p[1] = 0.0;
		s[pindx++] = p;
	    }
	}
    }
    
    // c) check whether line intersects side 1 of the triangle
    if (!boundary_only || bndside[1]) {
	if (p1[0]*p2[0] < 0.0) {
	    ys = p1[1] - p1[0]/(p2[0]-p1[0])*(p2[1]-p1[1]);
	    if (ys > -EPS && ys < 1.0+EPS) {
		p[0] = 0.0, p[1] = ys;
		s[pindx++] = p;
	    }
	}
    }
    
    // d) check whether line intersects side 2 of the triangle
    if (!boundary_only || bndside[2]) {
	if (fabs(p1[0]-p2[0]) > EPS) {
	    double a = (p1[1]-p2[1])/(p1[0]-p2[0]);
	    double c = p1[1] - a*p1[0];
	    xs = (1.0-c)/(a+1.0);
	    if (xs > -EPS && xs < 1.0+EPS) {
		p[0] = xs, p[1] = 1.0-xs;
		s[pindx++] = p;
	    }
	} else {
	    xs = p1[0];
	    if (xs > -EPS && xs < 1.0+EPS) {
		p[0] = xs, p[1] = 1.0-xs;
		s[pindx++] = p;
	    }
	}
    }
 Done:
    
    // check whether intersection points lie between line endpoints
    if (pindx) {
	double dp, a;
	int i, j, d=0;
	dp = p2[d] - p1[d];
	if (fabs(dp) < EPS) {
	    d = 1; dp = p2[d] - p1[d];
	}	    
	for (i = 0; i < pindx; i++) {
	    a = (s[i][d]-p1[d])/dp;
	    if (a < -EPS || a > 1.0+EPS) { // point outside ray
		for (j = i+1; j < pindx; j++)
		    s[j-1] = s[j];
		pindx--;
		i--;
	    }
	}
    }
    return pindx;
}

// ============================================================================
// Auxiliary functions

// res0 = Integral [(s-x)/s * sigma/2 * Cos[x*sigma-d], {x,x0,x1}] and
// res1 = Integral [x/s * sigma/2 * Cos[x*sigma-d], {x,x0,x1}]
// dmax is the radius of support of the cos (dmax >= 0) which can be used
// to limit the integration range. If not limited by dmax, then integration
// limits are x0=0, x1=s

static void Integrate_x_Cosine (double d, double dmax, double sigma, double s,
    double &res0, double &res1)
{
    if (d >= 0.0) {
        if (d-s >= dmax) { // outside support
	    res0 = res1 = 0.0;
	    return;
	}
    } else {
        if (-d >= dmax) {  // outside support
	    res0 = res1 = 0.0;
	    return;
	}
    }

    // find the integration limits
    double x0 = ::max (d-dmax, 0.0);
    double x1 = ::min (d+dmax, s);

    double arg0 = (d-x0)*sigma;
    double arg1 = (d-x1)*sigma;
    double fac  = 0.5/(s*sigma);

    double cos0 = cos(arg0), sin0 = sin(arg0);
    double cos1 = cos(arg1), sin1 = sin(arg1);

    res0 = ( cos0 - cos1 + sigma*((s-x0)*sin0 - (s-x1)*sin1))*fac;
    res1 = (-cos0 + cos1 + sigma*(x0*sin0-x1*sin1))*fac;
}


// Integration of product of shape function and cosine
// Cosine is centered at d and has radius of support a (area of positivity)
// shape function nodes are at x0 and x1 (x0 < x1)

static void Integrate_u_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_u1)
{
    if (x0 > d+a || x1 < d-a) {
	// case 1: interval completely outside support
	int_cos_u0 = 0.0;
	int_cos_u1 = 0.0;
    } else if (x0 >= d-a && x1 <= d+a) {
	// case 2: support of cosine spans full interval
	double arg1 = (2.0*a) / (Pi*Pi * (x0-x1));
	double arg2 = Pi*(d-x0)/(2.0*a);
	double arg3 = Pi*(d-x1)/(2.0*a);
	int_cos_u0 = arg1 * (-2.0*a*cos(arg2) + 2.0*a*cos(arg3) +
			     Pi*(x0-x1)*sin(arg2));
	int_cos_u1 = arg1 * ( 2.0*a*cos(arg2) - 2.0*a*cos(arg3) -
			     Pi*(x0-x1)*sin(arg3));
    } else if (x0 < d-a && x1 > d+a) {
	// case 3: support of cosine is subset of interval
	int_cos_u0 = 4.0*a*( d-x1) / (Pi * (x0-x1));
	int_cos_u1 = 4.0*a*(-d+x0) / (Pi * (x0-x1));
    } else if (x0 >= d-a) {
	// case 4: interval exceeds support of cosine to the right
	double arg1 = (2.0*a) / (Pi*Pi * (x0-x1));
	double arg2 = Pi*(d-x0)/(2.0*a);
	int_cos_u0 = arg1 * (-2.0*a*cos(arg2) +
			     Pi * (a+d-x1+(x0-x1)*sin(arg2)));
	int_cos_u1 = arg1 * (-Pi*(a+d-x0) + 2.0*a*cos(arg2));
    } else {
	// case 5: interval exceeds support of cosine to the left
	double arg1 = (2.0*a) / (Pi*Pi * (x0-x1));
	double arg2 = Pi*(d-x1)/(2.0*a);
	int_cos_u0 = arg1 * (-Pi*(a-d+x1) + 2.0*a*cos(arg2));
	int_cos_u1 = arg1 * (-2.0*a*cos(arg2) +
			     Pi * (a-d+x0-(x0-x1)*sin(arg2)));
    }
}

// =========================================================================
// Nonmember functions
// =========================================================================

static Point Triangle_SurfToLocal (int side, const Point &p)
{
    // Convert 1-D coordinate [0..1] into 2-D triangular surface coordinate
    // for the given side

    dASSERT(p.Dim() == 1, "Arg 2 wrong vector dimension");

    Point loc(2);
    switch (side) {
    case 0:  loc[0] = p[0];                     break;
    case 1:  loc[0] = 1.0 - (loc[1] = p[0]);    break;
    case 2:  loc[1] = 1.0 - p[0];               break;
    default: xERROR ("Arg 1 index out of range"); break;
    }
    return loc;
}

static int Triangle_GetLocalSubsampleAbsc (const Point *&absc)
{
    // homogeneously distributed quadrature points over triangle in
    // local coordinates

    static const int nsample = (NSUBSAMPLE * (NSUBSAMPLE+1)) / 2;
    static Point *absc_sample = 0;

    if (!absc_sample) {
	Point loc(2);
	int i, j, idx;
	absc_sample = new Point[nsample];
	for (i = idx = 0; i < NSUBSAMPLE; i++) {
	    loc[0] = (double)(i+0.5)/(double)NSUBSAMPLE;
	    for (j = 0; j < NSUBSAMPLE-i; j++) {
		loc[1] = (double)(j+0.5)/(double)NSUBSAMPLE;
		absc_sample[idx].New(2);
		absc_sample[idx] = loc;
		idx++;
	    }
	}
    }
    absc = absc_sample;
    return nsample;
}

static int Triangle_GetBndSubsampleAbsc (int side, const Point *&absc)
{
    static const int nsample = NSUBSAMPLE;
    static Point *absc_sample = 0;

    if (!absc_sample) {
	absc_sample = new Point[nsample];
	for (int i = 0; i < nsample; i++) {
	    absc_sample[i].New(1);
	    absc_sample[i][0] = (double)(i+0.5)/(double)nsample;
	}
    }
    absc = absc_sample;
    return nsample;
}
