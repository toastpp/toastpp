// ==========================================================================
// Module libfe
// File tri3.cc
// Definition of class Triangle3
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
int Triangle_GetLocalSubsampleAbsc (const Point *&absc);
int Triangle_GetBndSubsampleAbsc (int side, const Point *&absc);

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

Triangle3::Triangle3 ()
{
    Node = new int[nNode()];
    intfd_0 = 0, intfd_1 = 0;
}

Triangle3::Triangle3 (const Triangle3& el): Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
    intfd_0 = 0, intfd_1 = 0;
}

Triangle3::~Triangle3 ()
{
    delete []Node;
    if (intfd_0) delete []intfd_0;
    if (intfd_1) delete []intfd_1;
}

Element *Triangle3::Copy ()
{
    return new Triangle3(*this);
}

void Triangle3::Initialise (const NodeList& nlist)
{
    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1];

    // Set up triangle geometry parameters
    a0 = x1*y2 - x2*y1;  b0 = y1-y2;  c0 = x2-x1;
    a1 = x2*y0 - x0*y2;  b1 = y2-y0;  c1 = x0-x2;
    a2 = x0*y1 - x1*y0;  b2 = y0-y1;  c2 = x1-x0;

    Element_Unstructured_2D::Initialise (nlist);

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

int Triangle3::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, "Argument 1 index out of range");
    dASSERT(node >= 0 && node < 2, "Argument 2  index out of range");
    static int SN[3][2] = {{0,1}, {1,2}, {2,0}};
    return SN[side][node];
}

double Triangle3::SideSize (int sd, const NodeList &nlist) const
{
    return nlist[Node[SideNode(sd,0)]].Dist (nlist[Node[SideNode(sd,1)]]);
}

Point Triangle3::Local (const NodeList &nlist, const Point& glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");

    Point loc(2);
    double scale = 1.0/(a0+a1+a2);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1]) * scale;

    return loc;
}

Point Triangle3::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 3, "Argument 1 index out of range");
    Point n(2);  // note: initialised to zero
    if (node == 1)      n[0] = 1.0;
    else if (node == 2) n[1] = 1.0;
    return n;
}

Point Triangle3::SurfToLocal (int side, const Point &p) const
{
    return Triangle_SurfToLocal (side, p);
}

RVector Triangle3::DirectionCosine (int side, RDenseMatrix &jacin)
{
    RVector cosin(2);
    switch(side) {
    case 0: cosin[0] = -b2, cosin[1] = -c2; break;
    case 1: cosin[0] = -b0, cosin[1] = -c0; break;
    case 2: cosin[0] = -b1, cosin[1] = -c1; break;
    default: xERROR("Side index out of range");
    }
    return cosin/length(cosin);
};

const RVector &Triangle3::LNormal (int side) const
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

double Triangle3::ComputeSize (const NodeList &nlist) const
{
    return 0.5 * (a0+a1+a2);
    //return fabs (0.5 * (a0+a1+a2));
    // this should _not_ be the absolute value!
}

bool Triangle3::LContains (const Point& loc, bool pad) const
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

RVector Triangle3::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    RVector fun(3);
    fun[0] = 1.0 - loc[0] - loc[1];
    fun[1] = loc[0];
    fun[2] = loc[1];
    return fun;
}

void Triangle3::LocalShapeF (const Point &loc, RVector *fun) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    if (fun->Dim() != 3)
        fun->New(3);
    double *f = fun->data_buffer();
    f[0] = 1.0 - loc[0] - loc[1];
    f[1] = loc[0];
    f[2] = loc[1];
}

RDenseMatrix Triangle3::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static const RDenseMatrix der (2, 3, "-1 1 0   -1 0 1");
    return der;
}

RVector Triangle3::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RVector fun(3);
    double scale = 1.0/(2.0*size);
    fun[0] = scale * (a0 + b0*glob[0] + c0*glob[1]);
    fun[1] = scale * (a1 + b1*glob[0] + c1*glob[1]);
    fun[2] = scale * (a2 + b2*glob[0] + c2*glob[1]);
    return fun;
}

RDenseMatrix Triangle3::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der (2,3);
    double scale = 1.0/(2.0*size);
    der(0,0) = b0*scale;
    der(0,1) = b1*scale;
    der(0,2) = b2*scale;
    der(1,0) = c0*scale;
    der(1,1) = c1*scale;
    der(1,2) = c2*scale;
    return der;
}

int Triangle3::QuadRule (int order, const double **wght, const Point **absc)
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

double Triangle3::IntF (int i) const
{
    static const double third = 1.0/3.0;
    return size*third;
}

RSymMatrix Triangle3::IntFF () const {
#ifdef TRI3_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Triangle3::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3);
#ifdef TRI3_STORE_INTFF
    return intff(i,j);
#else
    return full_intff(i,j) * size;
#endif
}

double Triangle3::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 3 && j >= 0 && j < 3 && k >= 0 && k < 3);
    return size * intfff_scale[intfff_index[i][j][k]];
}

#ifdef UNDEF
  // MS 29.6.99 Removed because incompatible with quadratic shape functions

void Triangle3::IntFFF (double &iii, double &iij, double &ijk) const
{
    iii = size * intfff_scale[2];
    iij = size * intfff_scale[1];
    ijk = size * intfff_scale[0];
}
#endif

RSymMatrix Triangle3::IntPFF (const RVector &P) const
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

double Triangle3::IntPFF (int i, int j, const RVector &P) const
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

double Triangle3::IntFd (int i, int j, int k) const
{
    // Calculate Int [u_i  du_j/dx_k] dr
    // note that result is independent of i

    static double scale = 1.0/6.0;

    switch (j) {
    case 0:
	return (k ? c0 : b0) * scale;
    case 1:
	return (k ? c1 : b1) * scale;
    case 2:
	return (k ? c2 : b2) * scale;
    default:
	xERROR("Invalid index");
	return 0;
    }
}

double Triangle3::IntPd (const RVector &P, int j, int k) const
{
    // since the result of Int u_i du_j/dx_k is independent of
    // index i, we can pull the integral out of the sum over nodes

    double sum = 0.0;
    for (int i = 0; i < 3; i++)
	sum += P[Node[i]];
    return sum * IntFd (0, j, k);
}

double Triangle3::IntFdd (int i, int j, int k, int l, int m) const
{
    // Int [u_i du_j/dx_l du_k/dx_m] dr
    // note that result is independent of i

    double fac1, fac2;

    dASSERT (l < 2 && m < 2, "Invalid index");

    switch (j) {
    case 0: fac1 = (l == 0 ? b0 : c0); break;
    case 1: fac1 = (l == 0 ? b1 : c1); break;
    case 2: fac1 = (l == 0 ? b2 : c2); break;
    default: xERROR ("Invalid index"); return 0;
    }
    switch (k) {
    case 0: fac2 = (m == 0 ? b0 : c0); break;
    case 1: fac2 = (m == 0 ? b1 : c1); break;
    case 2: fac2 = (m == 0 ? b2 : c2); break;
    default: xERROR ("Invalid index"); return 0;
    }
    return fac1 * fac2 / (12.0*size);
}

double Triangle3::IntPdd (const RVector &P, int j, int k, int l, int m) const
{
    // since IntFdd doesn't depend on the first index, we can pull the
    // integral over the shape function out of the sum over nodes

    double val = 0;
    for (int i = 0; i < 3; i++) val += P[Node[i]];
    return val * IntFdd (0, j, k, l, m);
}

RSymMatrix Triangle3::Intdd () const
{
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

int Triangle3::GetLocalSubsampleAbsc (const Point *&absc) const
{
    return Triangle_GetLocalSubsampleAbsc (absc);
}

int Triangle3::GetBndSubsampleAbsc (int side, const Point *&absc) const
{
    dASSERT (side >= 0 && side < 3, "Argument 1 index out of range");
    return Triangle_GetBndSubsampleAbsc (side, absc);
}

void Triangle3::ComputeIntFD (const NodeList &nlist)
{
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

RVector Triangle3::IntFD (int i, int j) const
{
    dASSERT(i >= 0 && i < 3, "Parameter 1 index out of range");
    dASSERT(j >= 0 && j < 3, "Parameter 2 index out of range");

    // note that Int{F_i D_j} is independent of i

    static const double i6 = 1.0/6.0;

    RVector fd(2);
    fd[0] = (intfd_0 ? intfd_0[j] : j==0 ? b0*i6 : j==1 ? b1*i6 : b2*i6);
    fd[1] = (intfd_1 ? intfd_1[j] : j==0 ? c0*i6 : j==1 ? c1*i6 : c2*i6);
    return fd;
}

RSymMatrix Triangle3::ComputeIntDD (const NodeList &nlist) const
{
    // This is the direct evaluation of:
    //    Matrix lder = LocalShapeD (Point(2));
    //    Matrix gder = inverse (jac) * lder;
    //    return ATA (gder) * size;

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

double Triangle3::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 3, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 2: out of range");
    
    double d;
    switch (sd) {
    case 0:
	if (i == 2) return 0.0;
	else d = hypot (c2, b2);
	break;
    case 1:
	if (i == 0) return 0.0;
	else d = hypot (c0, b0);
	break;
    case 2:
	if (i == 1) return 0.0;
	else d = hypot (c1, b1);
	break;
    default:
	xERROR("Invalid side index");
	break;
    }
    return d * 0.5;
}

double Triangle3::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 3, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 3, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 3: out of range");
    
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
    default:
	xERROR("Invalid side index");
	break;
    }
    return d / (i == j ? 3.0 : 6.0);
}

RSymMatrix Triangle3::ComputeBndIntFF (const NodeList &nlist) const
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

RVector Triangle3::BndIntFX (int side, double (*func)(const Point&),
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

RVector Triangle3::BndIntFCos (int side, const Surface *surf,
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

RVector Triangle3::BndIntFCos (int side, const RVector &cntcos, double a,
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

RVector Triangle3::BndIntFDelta (int side, const Surface *surf,
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

int Triangle3::Intersection (const Point &p1, const Point &p2,
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

#ifdef UNDEF
void Triangle3::Subdivide (Mesh *mesh)
{
    dASSERT(sdnbhr, "Element neighbour list not initialised");
    dASSERT(subdivdata, "Element subdivision data not initialised");

    int level = subdivdata->level;

    if (level % 2 == 0) {
	// we are currently at an even (isotropic) subdivision level.
	// go to the next level by bisection


    }
}
#endif

void Triangle3::SplitSide (Mesh *mesh, int side, int newnode,
    Element *nbr1, Element *nbr2, Element *el1, Element *el2)
{
    dASSERT(subdivdata, "Subdivision data not available");
    // int level = subdivdata->level;

    if (0/*level % 2*/) {
	// we are currently at an odd subdivision level, so call for a
	// merge and resplit operation
	MergeAndResplit (mesh, side, newnode, nbr1, nbr2, el1, el2);
    } else {
	// we are currently at an even subdivision level, so call for
	// a bisection operation
	Bisect (mesh, side, newnode, nbr1, nbr2, el1, el2);
    }
}


void Triangle3::Bisect (Mesh *mesh, int side, int newnode,
    Element *nbr1, Element *nbr2, Element *el1, Element *el2)
{
    int i;

    dASSERT(subdivdata, "Subdivision data not available");
    // int level = subdivdata->level;
    //dASSERT((level % 2) == 0, "Subdivision level must be even");

    // If the side is not specified (side < 0), pick the longest side
    if (side < 0) {
	double maxlen = 0;
	for (i = 0; i < 3; i++) {
	    double len = mesh->nlist[Node[i]].Dist 
		(mesh->nlist[Node[(i+1)%3]]);
	    if (len > maxlen) maxlen = len, side = i;
	}
    }

    // Re-order the element nodes so that the bisection occurs on side 1
    int shift[3] = {1,0,2};
    int nd[3];
    Element *nbr[3];
    for (i = 0; i < 3; i++) {
	nd[i] = Node[i];
	nbr[i] = sdnbhr[i];
    }
    for (i = 0; i < 3; i++) {
	Node[(i+shift[side])%3] = nd[i];
	sdnbhr[(i+shift[side])%3] = nbr[i];
    }

    // check if the side neighbour is subdivided already
    Element *unsplit_nbr = (newnode >= 0 && nbr1 && nbr2 ? 0 : sdnbhr[1]);

#ifdef UNDEF
    // if the neighbour isn't split yet, check if we will be splitting its
    // longest side. Otherwise keep subdividing the neighbour until we do
    if (unsplit_nbr) {
	int maxside, splitside;
	do {
	    double maxlen = 0;
	    for (splitside = 0; splitside < unsplit_nbr->nSide(); splitside++)
		if (unsplit_nbr->sdnbhr[splitside] == this) break;
	    
	    for (i = 0; i < 3; i++) {
		double len = mesh->nlist[unsplit_nbr->Node[i]].Dist
		    (mesh->nlist[unsplit_nbr->Node[(i+1)%3]]);
		if (len > maxlen) maxlen = len, maxside = i;
	    }
	    if (maxside != splitside) {
		Element *e1, *e2;
		unsplit_nbr->Bisect (mesh, maxside, -1, NULL, NULL, e1, e2);
		unsplit_nbr = sdnbhr[1];
	    }
	} while (0/*maxside != splitside*/);
    }
#endif

    // if the required midpoint node doesn't exist yet, create it now
    if (newnode < 0) {
	newnode = mesh->nlist.Len();
	mesh->nlist.Append(1);
	mesh->nlist[newnode].New(2);
	for (i = 0; i < 2; i++)
	    mesh->nlist[newnode][i] =
		(mesh->nlist[Node[1]][i] + mesh->nlist[Node[2]][i]) * 0.5;
    }

    // create a new element to accommodate the sibling element
    Triangle3 *sib = new Triangle3;
    sib->Node[0] = Node[0];
    sib->Node[1] = newnode;
    sib->Node[2] = Node[2];
    mesh->elist.Append (sib);
    sib->InitNeighbourSupport();
    sib->InitSubdivisionSupport();
    sib->sdnbhr[0] = this;
    sib->sdnbhr[1] = nbr2;
    sib->sdnbhr[2] = sdnbhr[2];
    sib->subdivdata->level = subdivdata->level+1;
    sib->subdivdata->sibling = this;
    sib->subdivdata->is_sibling0 = false;

    // adjust our own node and side information to become the other sibling
    Node[2] = newnode;
    sdnbhr[1] = nbr[1];
    sdnbhr[2] = sib;
    subdivdata->level++;
    subdivdata->sibling = sib;
    subdivdata->is_sibling0 = true;

    // if neighbour subdivision data had not been provided, now subdivide
    // the neighbour
    if (unsplit_nbr) {
	Element *e1, *e2;
	for (i = 0; i < unsplit_nbr->nSide(); i++)
	    if (unsplit_nbr->sdnbhr[i] == this) break;
	unsplit_nbr->SplitSide (mesh, i, newnode, sib, this, e1, e2);
	sdnbhr[1] = e2;
	sib->sdnbhr[1] = e1;
    }
	
    // return our own subdivision data
    el1 = this;
    el2 = sib;
}

void Triangle3::MergeAndResplit (Mesh *mesh, int side, int newnode,
        Element *nbr1, Element *nbr2, Element *el1, Element *el2)
{
#ifdef UNDEF

    // This method can only be applied to triangles at odd subdivision levels:
    // - Merges the triangle with its sibling
    // - Splits the merged triangle into 4 self-similar triangles
    // In total, this operation generates two additional nodes and two
    // additional triangles
    // Recursively calls for subdivision of neighbours as required

    dASSERT(subdivdata, "Subdivision data not available");
    int level = subdivdata->level;
    dASSERT(level % 2, "Subdivision level must be odd");

    Element *sib = subdivdata->sibling; // the sibling triangle
    
    if (!subdivdata->is_sibling0) sib->MergeAndResplit (mesh, );
    // If we are not sibling0, just perform this method on the other sibling
    // so that we don't need to consider two different cases

    // We now know that we are sibling0, and that the other sibling is
    // attached at side 2

    int i;

    // Collect the vertex node indices of the merged triangle
    int nd[6];
    nd[0] = Node[0];
    nd[1] = Node[1];
    nd[2] = sib->Node[2];

    // The midpoint node for side 1 exists already
    nd[4] = Node[2];

    // We need two additional nodes to accommodate the 4-way split
    // One of those may already have been provided on the input
    int nnd = mesh->nlist.Len();
    if (side == 0) {
	if (newnode >= 0) {
	    nd[3] = newnode;
	} else {
	    mesh->nlist.Append(1);
	    mesh->nlist[nnd].New(2);
	    for (i = 0; i < 2; i++)
		mesh->nlist[nnd][i] =
		    (mesh->nlist[nd[0]][i] + mesh->nlist[nd[1]][i]) * 0.5;
	    nd[3] = nnd++;
	}
	mesh->nlist.Append(1);
	mesh->nlist[nnd].New(2);
	for (i = 0; i < 2; i++)
	    mesh->nlist[nnd][i] =
		(mesh->nlist[nd[2]][i] + mesh->nlist[nd[0]][i]) * 0.5;
	nd[5] = nnd++;
    } else if (side == 2) {
	mesh->nlist.Append(1);
	mesh->nlist[nnd].New(2);
	for (i = 0; i < 2; i++)
	    mesh->nlist[nnd][i] = 
		(mesh->nlist[nd[0]][i] + mesh->nlist[nd[1]][i]) * 0.5;
	nd[3] = nnd++;
	if (newnode >= 0) {
	    nd[5] = newnode;
	} else {
	    mesh->nlist.Append(1);
	    mesh->nlist[nnd].New(2);
	    for (i = 0; i < 2; i++)
		mesh->nlist[nnd][i] =
		    (mesh->nlist[nd[2]][i] + mesh->nlist[nd[0]][i]) * 0.5;
	    nd[5] = nnd++;
	}
    } else if (side == 1) {
	// this is a special case, because this side is not split during
	// the merge and resplit process. We will have to do another
	// subdivision after the merge and resplit
	mesh->nlist.Append(1);
	mesh->nlist[nnd].New(2);
	for (i = 0; i < 2; i++)
	    mesh->nlist[nnd][i] =
		(mesh->nlist[nd[0]][i] + mesh->nlist[nd[1]][i]) * 0.5;
	nd[3] = nnd++;
	mesh->nlist.Append(1);
	mesh->nlist[nnd].New(2);
	for (i = 0; i < 2; i++)
	    mesh->nlist[nnd][i] =
		(mesh->nlist[nd[2]][i] + mesh->nlist[nd[0]][i]) * 0.5;
	nd[5] = nnd++
    }
	
	    
    // Add two nodes to accommodate the two missing midpoints for the
    // re-split triangle
    int nnd = mesh->nlist.Len();
    mesh->nlist.Append(2);
    mesh->nlist[nnd].New(2);
    mesh->nlist[nnd+1].New(2);

    // set the vertex node indices of the merged triangle
    int nd[6];
    nd[0] = Node[0];
    nd[1] = Node[1];
    nd[2] = sib->Node[2];

    // assign the new midpoint nodes
    for (i = 0; i < 2; i++) {
	mesh->nlist[nnd][i] = 
	    (mesh->nlist[nd[0]][i] + mesh->nlist[nd[1]][i])*0.5;
	mesh->nlist[nnd+1][i] = 
	    (mesh->nlist[nd[2]][i] + mesh->nlist[nd[0]][i])*0.5;
    }

    // set the midpoint node indices
    nd[3] = nnd;
    nd[4] = Node[2];
    nd[5] = nnd+1;

    // Create two new triangles to accommodate the split
    int nel = mesh->elist.Len();
    Triangle3 *tri1, *tri2;
    tri1 = new Triangle3;
    tri1->Node[0] = nd[0];
    tri1->Node[1] = nd[3];
    tri1->Node[2] = nd[5];
    tri1->InitSubdivisionSupport();
    tri1->subdivdata->level = level+1;
    mesh->elist.Append(tri1);
    tri1->Initialise (mesh->nlist);
    tri1->sdnbhr = new Element*[3];
    tri2 = new Triangle3;
    tri2->Node[0] = nd[4];
    tri2->Node[1] = nd[5];
    tri2->Node[2] = nd[3];
    mesh->elist.Append(tri2);
    tri2->subdivdata->level = level+1;
    mesh->elist.Append(tri2);
    tri2->Initialise (mesh->nlist);
    tri2->sdnbhr = new Element*[3];

    // Re-assign the nodes for the two original siblings to form the other
    // two split triangles
    Node[0] = nd[3];
    subdivdata->level++;
    Initialise (mesh->nlist);
    sib->Node[0] = nd[5];
    sib->subdivdata->level++;
    sib->Initialise (mesh->nlist);

    // Call for subdivision of the two neighbours whose sides we have split
    // re-assign the neighbour lists
    Element *el1a, *el2a, *el1b, *el2b;
    sdnbhr[0]->Splitside(nd[0], nd[1], tri1, this, el1a, el2a, mesh);
    sib->sdnbhr[2]->Splitside(nd[2], nd[0], sib, tri1, el1b, el2b, mesh);

    tri1->sdnbhr[0] = el1a;
    tri1->sdnbhr[1] = tri2;
    tri1->sdnbhr[2] = el2b;

    tri2->sdnbhr[0] = sib;
    tri2->sdnbhr[1] = tri1;
    tri2->sdnbhr[2] = this;

    sib->sdnbhr[0] = tri2;
    // sib->sdnbhr[1] remains
    sib->sdnbhr[2] = el1b;

    sdnbhr[0] = el2a;
    //sdnbhr[1] remains
    sdnbhr[2] = tri2;
#endif
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

int Triangle_GetLocalSubsampleAbsc (const Point *&absc)
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

int Triangle_GetBndSubsampleAbsc (int side, const Point *&absc)
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
