// ==========================================================================
// Module libfe
// File tri3_qr.cc
// Definition of class Triangle3_qr
// ==========================================================================

//#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"
#include "element.h"
#include "tri3_qr.h"

// some global constants

static const double sqrt3    = sqrt(3.0);
static const double sqrt3_05 = 0.5 *sqrt3;

// index set for retrieving IntFFF(i,j,k)
static const int intfff_index[3][3][3] =
{{{2,1,1},{1,1,0},{1,0,1}},{{1,1,0},{1,2,1},{0,1,1}},{{1,0,1},{0,1,1},{1,1,2}}};

Triangle3_qr::Triangle3_qr (const Triangle3_qr& el): Element2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

void Triangle3_qr::Initialise (const NodeList &nlist)
{
#ifndef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
#endif // otherwise class members
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];

    jac.New (2,2);
    jac[0][0] = n1x-n0x;  jac[0][1] = n1y-n0y;
    jac[1][0] = n2x-n0x;  jac[1][1] = n2y-n0y;
    djac = det (jac);

    Element::Initialise (nlist);

    intfff_scale[0] = (1.0/60.0) * size;
    intfff_scale[1] = (1.0/30.0) * size;
    intfff_scale[2] = (1.0/10.0) * size;
}

int Triangle3_qr::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, Argument 1 index out of range);
    dASSERT(node >= 0 && node < 2, Argument 2  index out of range);
    static int SN[3][2] = {{0,1}, {1,2}, {2,0}};
    return SN[side][node];
}

Point Triangle3_qr::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 2, Argument 2 dimension must be 2);
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

Point Triangle3_qr::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 3, Argument 1 index out of range);
    Point n(2);  // note: initialised to zero
    if (node == 1)      n[0] = 1.0;
    else if (node == 2) n[1] = 1.0;
    return n;
}

RVector Triangle3_qr::DirectionCosine (int side, Matrix& jacin)
{
    static const double local_dc[3][2] = {
        {0,-1},{0.70710678,0.70710678},{-1,0}
    };

    dASSERT(side >= 0 && side < 3, Argument 1 index out of range);
    dASSERT(jacin.Dim(COL) == 2 && jacin.Dim(ROW) == 2,
	Argument 2 invalid dimension);

    RVector cosin(2);
    cosin[0] = jacin[0][0]*local_dc[side][0] + jacin[0][1]*local_dc[side][1];
    cosin[1] = jacin[1][0]*local_dc[side][0] + jacin[1][1]*local_dc[side][1];
    return cosin / length (cosin);
};

double Triangle3_qr::ComputeSize (const NodeList &nlist) const
{
    return 0.5 * fabs (djac);   // convert local to global size
}

bool Triangle3_qr::LContains (const Point& loc) const
{
    dASSERT(loc.Dim() == 2, Argument 1 invalid dimension);
    return (loc[0] >= 0.0 && loc[1] >= 0.0 && loc[0]+loc[1] <= 1.0);
}

bool Triangle3_qr::GContains (const Point& glob, const NodeList& nlist) const
{
    dASSERT(glob.Dim() == 2, Argument 1 invalid dimension);
    double xx = glob[0], yy = glob[1];

    // check bounding box
    if (xx < bbmin[0] || xx > bbmax[0] || yy < bbmin[1] || yy > bbmax[1])
        return FALSE;

    double x0, x1, x2, y0, y1, y2, y0r, yyr, fac;
    const double EPS = 1e-10;

    for (int i = 0; i < 3; i++) {
        x0 = nlist[Node[i]][0],       y0 = nlist[Node[i]][1];
	x1 = nlist[Node[(i+1)%3]][0], y1 = nlist[Node[(i+1)%3]][1];
	x2 = nlist[Node[(i+2)%3]][0], y2 = nlist[Node[(i+2)%3]][1];
	if (fabs (x1-x2) < EPS) {
	    if ((x0 < x1 && xx > x1) || (x0 > x1 && xx < x1)) return FALSE;
	} else {
	    fac = (y2-y1)/(x2-x1);
	    y0r = (x0-x1)*fac + y1;
	    yyr = (xx-x1)*fac + y1;
	    if ((y0 < y0r && yy > yyr) || (y0 > y0r && yy < yyr)) return FALSE;
	}
    }
    return TRUE;
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
        return TRUE;
    }
    return FALSE;
    */
}

RVector Triangle3_qr::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, Argument 1 invalid dimension);
    static RVector fun(3);
    fun[0] = 1.0 - loc[0] - loc[1];
    fun[1] = loc[0];
    fun[2] = loc[1];
    return fun;
}

Matrix Triangle3_qr::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, Argument 1 invalid dimension);
    static const Matrix der = Matrix (2, 3, "-1 1 0   -1 0 1");
    return der;
}

RVector Triangle3_qr::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{
    dASSERT(glob.Dim() == 2, Argument 2 invalid dimension);
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
    dASSERT(fun[0] >= -1e-8, Shape function 0 negative);
    dASSERT(fun[1] >= -1e-8, Shape function 1 negative);
    dASSERT(fun[2] >= -1e-8, Shape function 2 negative);
    return fun;
}

Matrix Triangle3_qr::GlobalShapeD (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 2, Argument 2 invalid dimension);
#ifndef TRI3_STORE_COORDS
    double n0x, n0y, n1x, n1y, n2x, n2y;
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
#endif
    double scale = 1.0/djac;
    static Matrix der (2,3);
    der[0][0] = (n1y-n2y) * scale;
    der[0][1] = (n2y-n0y) * scale;
    der[0][2] = (n0y-n1y) * scale;
    der[1][0] = (n2x-n1x) * scale;
    der[1][1] = (n0x-n2x) * scale;
    der[1][2] = (n1x-n0x) * scale;
    return der;
}

SymMatrix Triangle3_qr::ComputeIntFF (const NodeList&) const
{
    double val1 = size/6.0;     // or equivalently: sqrt3/8.0  * fabs(djac)
    double val2 = size/12.0;    // or equivalently: sqrt3/16.0 * fabs(djac)
    static SymMatrix ff(3);
    for (int i = 0; i < 3; i++) {
	ff[i][i] = val1;
	for (int j = 0; j < i; j++) ff[i][j] = val2;
    }
    return ff;
}

double Triangle3_qr::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3 && k >= 0 && k < 3);
    return intfff_scale[intfff_index[i][j][k]];
}


#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

void Triangle3_qr::IntFFF (double &iii, double &iij, double &ijk) const
{
    iii = intfff_scale[2];
    iij = intfff_scale[1];
    ijk = intfff_scale[0];
}
#endif

SymMatrix Triangle3_qr::IntPFF (const RVector &P) const
{
    SymMatrix pff(3);
    double fac = Size()/60.0;
    double p0 = P[Node[0]], p1 = P[Node[1]], p2 = P[Node[2]];
    pff[0][0] = fac * (p0*6.0 + p1*2.0 + p2*2.0);
    pff[1][0] = fac * (p0*2.0 + p1*2.0 + p2    );
    pff[2][0] = fac * (p0*2.0 + p1     + p2*2.0);
    pff[1][1] = fac * (p0*2.0 + p1*6.0 + p2*2.0);
    pff[2][1] = fac * (p0     + p1*2.0 + p2*2.0);
    pff[2][2] = fac * (p0*2.0 + p1*2.0 + p2*6.0);
    return pff;
}

double Triangle3_qr::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3);
    if (i < j) { int tmp = i; i = j; j = tmp; }
    switch (i) {
    case 0:
	return Size()/60.0 * (P[Node[0]]*6.0+P[Node[1]]*2.0+P[Node[2]]*2.0);
    case 1:
	switch (j) {
	case 0:
	    return Size()/60.0*(P[Node[0]]*2.0+P[Node[1]]*2.0+P[Node[2]]);
	case 1:
	    return Size()/60.0*(P[Node[0]]*2.0+P[Node[1]]*6.0+P[Node[2]]*2.0);
	}
    case 2:
	switch (j) {
	case 0:
	    return Size()/60.0*(P[Node[0]]*2.0+P[Node[1]]+P[Node[2]]*2.0);
	case 1:
	    return Size()/60.0*(P[Node[0]]+P[Node[1]]*2.0+P[Node[2]]*2.0);
	case 2:
	    return Size()/60.0*(P[Node[0]]*2.0+P[Node[1]]*2.0+P[Node[2]]*6.0);
	}
    default: // never get here, but needed to supress compiler warning
        return 0.0; // dummy
    }
}

SymMatrix Triangle3_qr::ComputeIntDD (const NodeList &nlist) const
{
    // this is the direct evaluation of:
    //    Matrix lder = LocalShapeD (Point(2));
    //    Matrix gder = inverse (jac) * lder;
    //    return ATA (gder) * size;

    SymMatrix dd(3);
    double fac = 0.25/size;
    double b0 = nlist[Node[1]][1] - nlist[Node[2]][1];
    double b1 = nlist[Node[2]][1] - nlist[Node[0]][1];
    double b2 = nlist[Node[0]][1] - nlist[Node[1]][1];
    double c0 = nlist[Node[2]][0] - nlist[Node[1]][0];
    double c1 = nlist[Node[0]][0] - nlist[Node[2]][0];
    double c2 = nlist[Node[1]][0] - nlist[Node[0]][0];
    dd[0][0] = fac * (b0*b0 + c0*c0);
    dd[1][0] = fac * (b1*b0 + c1*c0);
    dd[2][0] = fac * (b2*b0 + c2*c0);
    dd[1][1] = fac * (b1*b1 + c1*c1);
    dd[2][1] = fac * (b2*b1 + c2*c1);
    dd[2][2] = fac * (b2*b2 + c2*c2);
    return dd;
}

SymMatrix Triangle3_qr::ComputeBndIntFF (const NodeList &nlist) const
{
    SymMatrix bff(3);
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
	    bff[i][i] += d/3.0;
	    for (int j = 0; j < i; j++)
		if (j == n0 || j == n1) bff[i][j] += d/6.0;
	}
    }
    return bff;
}

#ifdef UNDEF
SymMatrix Triangle3_qr::FTF_bnd (const Matrix& geom, int side)
{
    int i, j, n0, n1, n2;
    double dx, dy, d;
    SymMatrix ftfb(3);

    n0 = SideNode (side, 0);
    n1 = SideNode (side, 1);
    for (n2 = 0; n2 == n0 || n2 == n1; n2++);

    dx = geom[n0][0] - geom[n1][0];
    dy = geom[n0][1] - geom[n1][1];
    d  = hypot (dx, dy); // length of integration side

    for (i = 0; i < 3; i++)
	for (j = 0; j <= i; j++)
	    ftfb[i][j] = (i == n2 || j == n2 ? 0.0 : (i == j ? d/3.0 : d/6.0));
    return ftfb;
}
#endif

#ifdef UNDEF
int Triangle3_qr::QuadRule (const double **wght, const Point **absc) const
{
    int i;
    static const double area = 3.0 * sqrt3 / 4.0;

    static const double local_wght[4] =
	{ 1.0/12.0*area, 1.0/12.0*area, 1.0/12.0*area, 3.0/4.0*area };
    static const Point2D local_absc[4] =
	{ Point2D(1,0), Point2D(-0.5,-sqrt3_05),
	  Point2D(-0.5,sqrt3_05), Point2D(0,0) };

    *wght = local_wght;
    *absc = local_absc;
    return 4;
}

int Triangle3_qr::BndQuadRule (int side, double **wght, Point **absc,
    double *coeff)
{
    int i, nqp = 2;
    *absc = new Point[nqp];
    *wght = new double[nqp];
    for (i = 0; i < nqp; i++) {
	(*absc)[i].New(2);
	(*wght)[i] = 1.0;
    }
    *coeff = 0.5;
    switch (side) {
    case 0:
	(*absc)[0][0] =  0.25 * (1.0 + sqrt3);
	(*absc)[0][1] =  0.25 * (1.0 - sqrt3);
	(*absc)[1][0] =  0.25 * (1.0 - sqrt3);
	(*absc)[1][1] = -0.25 * (1.0 + sqrt3);
	break;
    case 1:
	(*absc)[0][0] = -0.5;
	(*absc)[0][1] = -0.5;
	(*absc)[1][0] = -0.5;
	(*absc)[1][1] =  0.5;
	break;
    case 2:
	(*absc)[0][0] =  0.25 * (1.0 - sqrt3);
	(*absc)[0][1] =  0.25 * (1.0 + sqrt3);
	(*absc)[1][0] =  0.25 * (1.0 + sqrt3);
	(*absc)[1][1] = -0.25 * (1.0 - sqrt3);
	break;
    }
    return nqp;
}

double Triangle3_qr::UnitLength (Point&, Matrix& geom, int side)
{
    double dxdl, dydl;

    dASSERT(side >= 0 && side <= 3, Invalid side parameter.);
    dASSERT(geom.Dim(ROW) == 3 && geom.Dim(COL) == 2,
	Matrix parameter has wrong dimension.);

    switch (side) {
	case 0:
	    dxdl = geom[0][0] - geom[1][0];
	    dydl = geom[0][1] - geom[1][1];
	    break;
	case 1:
	    dxdl = geom[1][0] - geom[2][0];
	    dydl = geom[1][1] - geom[2][1];
	    break;
	case 2:
	    dxdl = geom[2][0] - geom[0][0];
	    dydl = geom[2][1] - geom[0][1];
	    break;
    }
    return hypot (dxdl, dydl);
}

inline double xi (double /*l1*/, double l2, double l3)
{ return 1.0 - 3.0/2.0 * (l2+l3); }

inline double eta (double /*l1*/, double l2, double l3)
{ return sqrt3_05 * (l3-l2); }

void Triangle3_qr::ConvLineQuadrature (Point** absc, double* labsc,
    int nqp, int side, double* coeff)
{
    int i;
    double l1, l2, l3;

    *absc = new Point[nqp];
    for (i=0; i<nqp; i++) (*absc)[i].New(2);

    dASSERT(nqp > 0, Invalid number of quadrature points.);
    dASSERT(side >= 0 && side < 3, Invalid side parameter.);

    *coeff = 0.5;
    switch (side) {
	case 0:
	    l3 = 0.0;
	    for (i = 0; i < nqp; i++) {
		l1 = 0.5 * (1.0+labsc[i]);
		l2 = 1.0 - l1;
		(*absc)[i][0] = xi (l1, l2, l3);
		(*absc)[i][1] = eta (l1, l2, l3);
	    }
	    break;
	case 1:
	    l1 = 0.0;
	    for (i = 0; i < nqp; i++) {
		l2 = 0.5 * (1.0+labsc[i]);
		l3 = 1.0 - l2;
		(*absc)[i][0] = xi (l1, l2, l3);
		(*absc)[i][1] = eta (l1, l2, l3);
	    }
	    break;
	case 2:
	    l2 = 0.0;
	    for (i = 0; i < nqp; i++) {
		l3 = 0.5 * (1.0+labsc[i]);
		l1 = 1.0 - l3;
		(*absc)[i][0] = xi (l1, l2, l3);
		(*absc)[i][1] = eta (l1, l2, l3);
	    }
	    break;
    }
}
#endif //UNDEF

int Triangle3_qr::Intersection (const Point &p1, const Point &p2,
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
