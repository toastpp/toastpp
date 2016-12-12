// ==========================================================================
// Module libfe
// File line2d2.cc
// Definition of class Line2D2
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "lin_qr.h"

// local prototypes

#ifdef DO_THE_REST
static void Integrate_x_Cosine (double d, double dmax, double sigma, double s,
    double &res0, double &res1);
static void Integrate_u_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_u1);
static Point Triangle_SurfToLocal (int side, const Point &p);
int Triangle_GetLocalSubsampleAbsc (const Point *&absc);
int Triangle_GetBndSubsampleAbsc (int side, const Point *&absc);
#endif

// some global constants
static bool subsampling_initialised = false;
static const int nsample_lin = NSUBSAMPLE; // from toastdef.h
static Point absc_bndsample[3][nsample_lin];
static const RSymMatrix sym_intff = RSymMatrix (2,
   "2 \
    1 2") * (1.0/16.0);

static const RDenseMatrix full_intff = RDenseMatrix (2, 2,
   "2 1  \
    1 2") * (1.0/6.0);

// index set for retrieving IntFFF(i,j,k)
static const int intfff_index[2][2][2] = {
  {{1,0},{0,0}},{{0,0},{0,1}}
};
static const double intfff_scale[2] = {1.0/12.0, 1.0/4.0};

Line2D2::Line2D2 ()
{
    Node = new int[nNode()];
    intfd_0 = 0, intfd_1 = 0;
}

Line2D2::Line2D2 (const Line2D2& el): Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
    intfd_0 = 0, intfd_1 = 0;
}

Line2D2::~Line2D2 ()
{
    delete []Node;
    if (intfd_0) delete []intfd_0;
    if (intfd_1) delete []intfd_1;
}

Element *Line2D2::Copy ()
{
    return new Line2D2(*this);
}

void Line2D2::Initialise (const NodeList& nlist)
{
    double x0 = nlist[Node[0]][0];
    double x1 = nlist[Node[1]][0];

    // Set up line geometry parameters (almost trivial)
    a0 = x1;  b0 = -1;
    a1 = -x0; b1 = 1;

    Element_Unstructured_2D::Initialise (nlist);
    dASSERT(size > 0, "Element size not positive");

#ifdef LINE2D2_STORE_INTFF
    intff.New(2);
    intff = sym_intff * size;
#endif

#ifdef LINE2D2_STORE_INTDD
    intdd.New(2);
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

int Line2D2::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 1, "Argument 1 index out of range");
    dASSERT(node >= 0 && node < 1, "Argument 2  index out of range");
    static int SN[2][1] = {{0}, {1}};
    return SN[side][node];
}

Point Line2D2::Local (const NodeList &nlist, const Point& glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");

    Point loc(2);
    double scale = 1.0/(a0+a1);
    loc[0] = (a1 + b1*glob[0] ) * scale;

    return loc;
}

Point Line2D2::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 2, "Argument 1 index out of range");
    Point n(2);  // note: initialised to zero
    if (node == 1)      n[0] = 1.0;
    return n;
}

//Point Line2D2::SurfToLocal (int side, const Point &p) const
//{
//    return Triangle_SurfToLocal (side, p);
//}
//
//RVector Line2D2::DirectionCosine (int side, RDenseMatrix &jacin)
//{
//    RVector cosin(2);
//    switch(side) {
//    case 0: cosin[0] = -b2, cosin[1] = -c2; break;
//    case 1: cosin[0] = -b0, cosin[1] = -c0; break;
//    case 2: cosin[0] = -b1, cosin[1] = -c1; break;
//    default: xERROR(Side index out of range);
//    }
//    return cosin/length(cosin);
//};
//
//const RVector &Line2D2::LNormal (int side) const
//{
//    static const RVector lnm0 = RVector (2, "0 -1");
//    static const RVector lnm1 = RVector (2, "0.7071067814 0.7071067814");
//    static const RVector lnm2 = RVector (2, "-1 0");
//    static const RVector *lnm[3] = {
//      &lnm0, &lnm1, &lnm2
//    };
//    dASSERT(side >= 0 && side < 3, "Argument 1 index out of range");
//    return *lnm[side];
//}

double Line2D2::ComputeSize (const NodeList &nlist) const
{
    return 0.5 * (a0+a1);
    // this should _not_ be the absolute value!
}

bool Line2D2::LContains (const Point& loc, bool pad) const
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

//bool Line2D2::GContains (const Point& glob, const NodeList& nlist) const
//{
//    dASSERT(glob.Dim() == 2, "Argument 1 invalid dimension");
//    double xx = glob[0], yy = glob[1];

    // check bounding box
//    if (xx < bbmin[0] || xx > bbmax[0] || yy < bbmin[1] || yy > bbmax[1])
//        return false;
//
//    double x0, x1, x2, y0, y1, y2, y0r, yyr, fac;
//    const double EPS = 1e-10;
//
//    for (int i = 0; i < 3; i++) {
//        x0 = nlist[Node[i]][0],       y0 = nlist[Node[i]][1];
//	x1 = nlist[Node[(i+1)%3]][0], y1 = nlist[Node[(i+1)%3]][1];
//	x2 = nlist[Node[(i+2)%3]][0], y2 = nlist[Node[(i+2)%3]][1];
//	if (fabs (x1-x2) < EPS) {
//	    if ((x0 < x1 && xx > x1) || (x0 > x1 && xx < x1)) return false;
//	} else {
//	    fac = (y2-y1)/(x2-x1);
//	    y0r = (x0-x1)*fac + y1;
//	    yyr = (xx-x1)*fac + y1;
//	    if ((y0 < y0r && yy > yyr) || (y0 > y0r && yy < yyr)) return false;
//	}
//    }
//    return true;
//}

RVector Line2D2::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static RVector fun(2);
    fun[0] = 1.0 - loc[0];
    fun[1] = loc[0];
    return fun;
}

RDenseMatrix Line2D2::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static const RDenseMatrix der (1, 2, "-1 1");
    return der;
}

RVector Line2D2::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{  // assume glob is on line. Otherwise drop perpendicular ?
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RVector fun(2);
    double scale = 1.0/(size);
    fun[0] = scale * (a0 + b0*glob[0]);
    fun[1] = scale * (a1 + b1*glob[0]);
    return fun;
}

RDenseMatrix Line2D2::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der (1,2);
    double scale = 1.0/(size);
    der(0,0) = b0*scale;
    der(0,1) = b1*scale;
    return der;
}

//int Line2D2::QuadRule (int order, const double **wght, const Point **absc)
//    const
//{
//    switch (order) {
//    case 1: return QRule_tri_1_1 (wght, absc);
//    case 2: return QRule_tri_2_3 (wght, absc);
//    case 3: return QRule_tri_3_6 (wght, absc);
//    case 4: return QRule_tri_4_6 (wght, absc);
//    case 6: return QRule_tri_6_12 (wght, absc);
//    default: xERROR(Not implemented); return 0;
//    }
//}

double Line2D2::IntF (int i) const
{
    static const double half = 1.0/2.0;
    return size*half;
}

RSymMatrix Line2D2::IntFF () const {
#ifdef LINE2D2_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Line2D2::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 2 && j >= 0 && j < 2);
#ifdef LINE2D2_STORE_INTFF
    return intff(i,j);
#else
    return full_intff(i,j) * size;
#endif
}

double Line2D2::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 2 && j >= 0 && j < 2 && k >= 0 && k < 2);
    return size * intfff_scale[intfff_index[i][j][k]];
}


RSymMatrix Line2D2::IntPFF (const RVector &P) const
{
    static RSymMatrix pff(2);
    static const double i12 = 1.0/12.0;
    double fac = size*i12;
    double p0 = P[Node[0]], p1 = P[Node[1]];
    pff(0,0) = fac * (3.0*p0 + p1);
    pff(1,0) = fac * (p0 + p1);
    pff(1,1) = fac * (p0 + 3.0*p1);
    return pff;
}

double Line2D2::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 2 && j >= 0 && j < 2);
    if (i < j) { int tmp = i; i = j; j = tmp; }
    static const double i12 = 1.0/12.0;
    switch (i) {
    case 0:
        return size*i12 * (3.0*P[Node[0]] + P[Node[1]] ); 
    case 1:
	switch (j) {
	case 0:
	    return size*i12 * (P[Node[0]] + P[Node[1]] );
	case 1:
	    return size*i12 * (P[Node[0]] + 3.0*P[Node[1]] );
	}

    default:
        xERROR("Index out of range");
        return 0.0; // dummy
    }
}

double Line2D2::IntFd (int i, int j, int k) const
{
    // Calculate Int [u_i  du_j/dx_k] dr
    // note that result is independent of i
  // NOTE only one dimension, so k index is irrelevent
    static double scale = 1.0/2.0;

    switch (j) {
    case 0:
	return ( b0) * scale;
    case 1:
	return ( b1) * scale;
    default:
	xERROR("Invalid index");
	return 0;
    }
}

double Line2D2::IntPd (const RVector &P, int j, int k) const
{
    // since the result of Int u_i du_j/dx_k is independent of
    // index i, we can pull the integral out of the sum over nodes

    double sum = 0.0;
    for (int i = 0; i < 2; i++)
	sum += P[Node[i]];
    return sum * IntFd (0, j, k);
}

double Line2D2::IntFdd (int i, int j, int k, int l, int m) const
{
    // Int [u_i du_j/dx_l du_k/dx_m] dr
    // note that result is independent of i

    double fac1, fac2;

    dASSERT (l < 2 && m < 2, "Invalid index");

    switch (j) {
    case 0: fac1 = (b0 ); break;
    case 1: fac1 = ( b1); break;
    default: xERROR ("Invalid index"); return 0;
    }
    switch (k) {
    case 0: fac2 = (b0 ); break;
    case 1: fac2 = ( b1); break;
    default: xERROR ("Invalid index"); return 0;
    }
    return fac1 * fac2 / (2.0*size);
}

double Line2D2::IntPdd (const RVector &P, int j, int k, int l, int m) const
{
    // since IntFdd doesn't depend on the first index, we can pull the
    // integral over the shape function out of the sum over nodes

    double val = 0;
    for (int i = 0; i < 2; i++) val += P[Node[i]];
    return val * IntFdd (0, j, k, l, m);
}

RSymMatrix Line2D2::Intdd () const
{
    RSymMatrix MDD(2,2);
    MDD(0,0) = b0*b0;
    MDD(1,0) = b0*b1;;
    MDD(2,2) = b1*b1;

    MDD /= size;
    return MDD;
}

//int Line2D2::GetLocalSubsampleAbsc (const Point *&absc) const
//{
//    return Triangle_GetLocalSubsampleAbsc (absc);
//}

//int Line2D2::GetBndSubsampleAbsc (int side, const Point *&absc) const
//{
//    dASSERT (side >= 0 && side < 3, "Argument 1 index out of range");
//    return Triangle_GetBndSubsampleAbsc (side, absc);
//}

void Line2D2::ComputeIntFD (const NodeList &nlist)
{
    static const double i2 = 1.0/2.0;

    // note that Int{F_i D_j} is independent of i, so we
    // need to store only one row

    if (intfd_0) delete []intfd_0;
    intfd_0 = new double[2];
    intfd_0[0] = b0*i2;
    intfd_0[1] = b1*i2;
}

RVector Line2D2::IntFD (int i, int j) const
{
    dASSERT(intfd_0 && intfd_1, "FD matrix not initialised");
    dASSERT(i >= 0 && i < 2, "Parameter 1 index out of range");
    dASSERT(j >= 0 && j < 2, "Parameter 2 index out of range");

    // note that Int{F_i D_j} is independent of i

    RVector fd(1);
    fd[0] = intfd_0[j];
    return fd;
}

RSymMatrix Line2D2::ComputeIntDD (const NodeList &nlist) const
{
    // This is the direct evaluation of:
    //    Matrix lder = LocalShapeD (Point(2));
    //    Matrix gder = inverse (jac) * lder;
    //    return ATA (gder) * size;

    static RSymMatrix dd(2);
    double scale = 1.0/size;
    dd(0,0) = scale * (b0*b0);
    dd(1,0) = scale * (b1*b0);
    dd(1,1) = scale * (b1*b1);
 
    return dd;
}
#ifdef DO_THE_REST

double Line2D2::SurfIntFF (int i, int j, int sd)
{
    double d;
    switch (sd) {
    case 0:
	if (i == 2 || j == 2) return 0.0; // u_2=0 along side 0
	else d = hypot (c2, b2);          // 0->1 node distance
	break;
    case 1:
	if (i == 0 || j == 0) return 0.0; // u_0=0 along side 1
	else d = hypot (c0, b0);          // 1->2 node distance
	break;
    case 2:
	if (i == 1 || j == 1) return 0.0; // u_1=0 along side 2
	else d = hypot (c1, b1);          // 2->0 node distance
	break;
    }
    return d / (i == j ? 3.0 : 6.0);
}

RSymMatrix Line2D2::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(3);
    if (!(bndel||interfaceel)) 
      return bff;  // not a boundary or interface element -> return zero matrix

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

RVector Line2D2::BndIntFX (int side, double (*func)(const Point&),
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

RVector Line2D2::BndIntFCos (int side, const Surface *surf,
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

RVector Line2D2::BndIntFCos (int side, const RVector &cntcos, double a,
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

RVector Line2D2::BndIntFDelta (int side, const Surface *surf,
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

int Line2D2::GlobalIntersection (const NodeList &nlist,
    const Point &p1, const Point &p2, Point **list)
{
    double bbx1 = 1e6, bbx2 = -1e6, bby1 = 1e6, bby2 = -1e6;
    double m0, m;
    int i, below, above;

    // generate bounding box
    RDenseMatrix egeom = Elgeom (nlist);
    for (i = 0; i < 3; i++) {
	if (egeom(i,0) < bbx1) bbx1 = egeom(i,0);
	if (egeom(i,0) > bbx2) bbx2 = egeom(i,0);
	if (egeom(i,1) < bby1) bby1 = egeom(i,1);
	if (egeom(i,1) > bby2) bby2 = egeom(i,1);
    }

    // test whether line is completely outside bounding box
    if (p1[0] < bbx1 && p2[0] < bbx1 ||
	p1[0] > bbx2 && p2[0] > bbx2 ||
	p1[1] < bby1 && p2[1] < bby1 ||
	p1[1] > bby2 && p2[1] > bby2) return 0;

    // test whether all points lie below or above the line
    below = 0, above = 0;
    if (p1[0] != p2[0]) {
	m0 = (p2[1] - p1[1]) / (p2[0] - p1[0]);
	for (i = 0; i < 3; i++) {
	    if (egeom(i,0) != p1[0])
		m = (egeom(i,1) - p1[1]) / (egeom(i,0) - p1[0]);
	    else if (egeom(i,1) > p1[1])
		m = 1e10;
	    else m = -1e10;
	    if (m > m0) above++;
	    else below++;
	}
    } else {
	for (i = 0; i < 3; i++) {
	    if (egeom(i,0) < p1[0]) above++;
	    else below++;
	}
    }
    if (!above || !below) return 0;

    // do the rest on local basis
    Point loc1 = Local (nlist, p1);
    Point loc2 = Local (nlist, p2);
    return Intersection (loc1, loc2, list);
}

int Line2D2::Intersection (const Point &p1, const Point &p2, Point **list)
{
    double xs, ys;
    double pmin, pmax, smin, smax;
    const double EPS = 1e-20;
    int pindx = 0;
    Point *s;

    dASSERT(p1.Dim() == 2 && p2.Dim() == 2, "Points must be 2D.");
    s = new Point[2];
    s[0].New(2), s[1].New(2);

    // a) check whether one of the end points of the line is within the element
    if (LContains (p1)) s[pindx++]=p1;
    if (LContains (p2)) s[pindx++]=p2;
    if (pindx==2) goto Done;

    // b) check whether line intersects side 1 of the triangle
    if (p1[0]-p2[0] == sqrt3*(p1[1]-p2[1])) goto DoneSide1;
    if (p1[0] == p2[0]) {
	xs = p1[0];
    } else {
	double m = (p1[1]-p2[1])/(p1[0]-p2[0]);
	xs = (-m*p1[0] + p1[1] + 1.0/sqrt3) / (1.0/sqrt3 - m);
    }
    ys = 1.0/sqrt3 * xs - 1.0/sqrt3;
    if (xs < -0.5 || xs > 1.0) goto DoneSide1;
    s[pindx][0] = xs;
    s[pindx][1] = ys;
    pindx++;
    if (pindx == 2) goto Done;
    DoneSide1:

    // c) check whether line intersects side 2 of the triangle
    if (p1[0]-p2[0] == -sqrt3 * (p1[1]-p2[1])) goto DoneSide2;
    if (p1[0] == p2[0]) {
	xs = p1[0];
    } else {
	double m = (p1[1]-p2[1])/(p1[0]-p2[0]);
	xs = (-m*p1[0] + p1[1] - 1.0/sqrt3) / (-1.0/sqrt3 - m);
    }
    ys = -1.0/sqrt3 * xs + 1.0/sqrt3;
    if (xs < -0.5 || xs > 1.0) goto DoneSide2;
    s[pindx][0] = xs;
    s[pindx][1] = ys;
    pindx++;
    if (pindx == 2) goto Done;
    DoneSide2:

    // d) check whether line intersects side 3 of the triangle
    if (p1[0] == p2[0]) goto Done;
    ys = (p1[1]-p2[1])/(p1[0]-p2[0]) * (-0.5-p1[0]) + p1[1];
    if (fabs(ys) > sqrt3_05) goto Done;
    s[pindx][0] = -0.5;
    s[pindx][1] = ys;
    pindx++;

    Done:
    if (pindx == 1) pindx=0; // if only one intersection found, forget it

    // check whether intersection points lie between line endpoints
    if (pindx == 2) {
	if (fabs(p1[0]-p2[0]) > EPS) {
	    if (p1[0] < p2[0]) pmin=p1[0], pmax=p2[0];
	    else               pmin=p2[0], pmax=p1[0];
	    if (s[0][0] < s[1][0]) smin=s[0][0], smax=s[1][0];
	    else                   smin=s[1][0], smax=s[0][0];
	    if (smax+EPS < pmin || smin-EPS > pmax) pindx=0;
	} else {
	    if (p1[1] < p2[1]) pmin=p1[1], pmax=p2[1];
	    else               pmin=p2[1], pmax=p1[1];
	    if (s[0][1] < s[1][1]) smin=s[0][1], smax=s[1][1];
	    else                   smin=s[1][1], smax=s[0][1];
	    if (smax+EPS < pmin || smin-EPS > pmax) pindx=0;
	}
    }
    xASSERT(pindx == 0 || pindx == 2, "Something went wrong...");
    if (!pindx) {
	s[0].New(0), s[1].New(0);
	delete s;
	*list = NULL;
    } else *list = s;
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
    double x0 = max (d-dmax, 0.0);
    double x1 = min (d+dmax, s);

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
	double arg1 = (2.0*a) / (M_PI*M_PI * (x0-x1));
	double arg2 = M_PI*(d-x0)/(2.0*a);
	double arg3 = M_PI*(d-x1)/(2.0*a);
	int_cos_u0 = arg1 * (-2.0*a*cos(arg2) + 2.0*a*cos(arg3) +
			     M_PI*(x0-x1)*sin(arg2));
	int_cos_u1 = arg1 * ( 2.0*a*cos(arg2) - 2.0*a*cos(arg3) -
			     M_PI*(x0-x1)*sin(arg3));
    } else if (x0 < d-a && x1 > d+a) {
	// case 3: support of cosine is subset of interval
	int_cos_u0 = 4.0*a*( d-x1) / (M_PI * (x0-x1));
	int_cos_u1 = 4.0*a*(-d+x0) / (M_PI * (x0-x1));
    } else if (x0 >= d-a) {
	// case 4: interval exceeds support of cosine to the right
	double arg1 = (2.0*a) / (M_PI*M_PI * (x0-x1));
	double arg2 = M_PI*(d-x0)/(2.0*a);
	int_cos_u0 = arg1 * (-2.0*a*cos(arg2) +
			     M_PI * (a+d-x1+(x0-x1)*sin(arg2)));
	int_cos_u1 = arg1 * (-M_PI*(a+d-x0) + 2.0*a*cos(arg2));
    } else {
	// case 5: interval exceeds support of cosine to the left
	double arg1 = (2.0*a) / (M_PI*M_PI * (x0-x1));
	double arg2 = M_PI*(d-x1)/(2.0*a);
	int_cos_u0 = arg1 * (-M_PI*(a-d+x1) + 2.0*a*cos(arg2));
	int_cos_u1 = arg1 * (-2.0*a*cos(arg2) +
			     M_PI * (a-d+x0-(x0-x1)*sin(arg2)));
    }
}
#endif




