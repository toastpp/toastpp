// ==========================================================================
// Module libfe
// File tri6.cc
// Definition of class Triangle6
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"

using namespace std;

// local prototypes

void Integrate_x2_Cosine (double d, double dmax, double sigma, double s,
    double g, double &res0, double &res1, double &res2);
void Integrate_u_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_um, double &int_cos_u1);

static bool subsampling_initialised = false;
static const int nsample_lin = NSUBSAMPLE; // from toastdef.h
static const int nsample_tot = (nsample_lin*(nsample_lin + 1)) / 2;
static Point absc_sample[nsample_tot];

static const RDenseMatrix full_intff = RDenseMatrix (6,6,
   "12 -2 -2  0 -8  0 \
    -2 12 -2  0  0 -8 \
    -2 -2 12 -8  0  0 \
     0  0 -8 64 32 32 \
    -8  0  0 32 64 32 \
     0 -8  0 32 32 64") * (2.0/720.0);

static const RSymMatrix sym_intff = RSymMatrix (6,
   "12 \
    -2 12 \
    -2 -2 12 \
     0  0 -8 64 \
    -8  0  0 32 64 \
     0 -8  0 32 32 64") * (2.0/720.0);

// below are the 6 planes of the IntFFF tensor
static const RDenseMatrix full_intf0ff = RDenseMatrix (6,6,
   " 288  -32  -32  192   64  192 \
     -32  -32   16  -64    0    0 \
     -32   16  -32    0    0  -64 \
     192  -64    0    0 -128    0 \
      64    0    0 -128 -256 -128 \
     192    0  -64    0 -128    0") * (2.0/40320.0);
static const RDenseMatrix full_intf1ff = RDenseMatrix (6,6,
   " -32  -32   16  -64    0    0 \
     -32  288  -32  192  192   64 \
      16  -32  -32    0  -64    0 \
     -64  192    0    0    0 -128 \
       0  192  -64    0    0 -128 \
       0   64    0 -128 -128 -256") * (2.0/40320.0);
static const RDenseMatrix full_intf2ff = RDenseMatrix (6,6,
   " -32   16  -32    0    0  -64 \
      16  -32  -32    0  -64    0 \
     -32  -32  288   64  192  192 \
       0    0   64 -256 -128 -128 \
       0  -64  192 -128    0    0 \
     -64    0  192 -128    0    0") * (2.0/40320.0);
static const RDenseMatrix full_intf3ff = RDenseMatrix (6,6,
   " 192  -64    0    0 -128    0 \
     -64  192    0    0    0 -128 \
       0    0   64 -256 -128 -128 \
       0    0 -256 2304  768  768 \
    -128    0 -128  768  768  512 \
       0 -128 -128  768  512  768") * (2.0/40320.0);
static const RDenseMatrix full_intf4ff = RDenseMatrix (6,6,
   "  64    0    0 -128 -256 -128 \
       0  192  -64    0    0 -128 \
       0  -64  192 -128    0    0 \
    -128    0 -128  768  768  512 \
    -256    0    0  768 2304  768 \
    -128 -128    0  512  768  768") * (2.0/40320.0);
static const RDenseMatrix full_intf5ff = RDenseMatrix (6,6,
   " 192    0  -64    0 -128    0 \
       0   64    0 -128 -128 -256 \
     -64    0  192 -128    0    0 \
       0 -128 -128  768  512  768 \
    -128 -128    0  512  768  768 \
       0 -256    0  768  768 2304") * (2.0/40320.0);

static const RDenseMatrix *full_intfff[6] = {
    &full_intf0ff,
    &full_intf1ff,
    &full_intf2ff,
    &full_intf3ff,
    &full_intf4ff,
    &full_intf5ff
};

// the same in symmetric matrix representation
static const RSymMatrix sym_intf0ff = RSymMatrix (6,
   " 288 \
     -32  -32 \
     -32   16  -32 \
     192  -64    0    0 \
      64    0    0 -128 -256 \
     192    0  -64    0 -128    0") * (2.0/40320.0);
static const RSymMatrix sym_intf1ff = RSymMatrix (6,
   " -32 \
     -32  288 \
      16  -32  -32 \
     -64  192    0    0 \
       0  192  -64    0    0 \
       0   64    0 -128 -128 -256") * (2.0/40320.0);
static const RSymMatrix sym_intf2ff = RSymMatrix (6,
   " -32 \
      16  -32 \
     -32  -32  288 \
       0    0   64 -256 \
       0  -64  192 -128    0 \
     -64    0  192 -128    0    0") * (2.0/40320.0);
static const RSymMatrix sym_intf3ff = RSymMatrix (6,
   " 192 \
     -64  192 \
       0    0   64 \
       0    0 -256 2304 \
    -128    0 -128  768  768 \
       0 -128 -128  768  512  768") * (2.0/40320.0);
static const RSymMatrix sym_intf4ff = RSymMatrix (6,
   "  64 \
       0  192 \
       0  -64  192 \
    -128    0 -128  768 \
    -256    0    0  768 2304 \
    -128 -128    0  512  768  768") * (2.0/40320.0);
static const RSymMatrix sym_intf5ff = RSymMatrix (6,
   " 192 \
       0   64 \
     -64    0  192 \
       0 -128 -128  768 \
    -128 -128    0  512  768 \
       0 -256    0  768  768 2304") * (2.0/40320.0);

// BndIntF over each side
static const RDenseMatrix bndintf = RDenseMatrix (3, 6,
   "1 1 0 4 0 0 \
    0 1 1 0 4 0 \
    1 0 1 0 0 4") * (1.0/6.0);

// BndIntFF over side 0
static const RSymMatrix sym_bndintff_sd0 = RSymMatrix (6,
   " 4 \
    -1 4 \
     0 0 0 \
     2 2 0 16 \
     0 0 0 0 0 \
     0 0 0 0 0 0") * (1.0/30.0);

// BndIntFF over side 1
static const RSymMatrix sym_bndintff_sd1 = RSymMatrix (6,
   "0 \
    0 4 \
    0 -1 4 \
    0 0 0 0 \
    0 2 2 0 16 \
    0 0 0 0 0 0") * (1.0/30.0);

// BndIntFF over side 2
static const RSymMatrix sym_bndintff_sd2 = RSymMatrix (6,
   " 4 \
     0 0 \
    -1 0 4 \
     0 0 0 0 \
     0 0 0 0 0 \
     2 0 2 0 0 16") * (1.0/30.0);


Triangle6::Triangle6 (const Triangle6 &el): Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Triangle6::Copy ()
{
    return new Triangle6(*this);
}

void Triangle6::Initialise (const NodeList &nlist)
{
    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1];

    a0 = x1*y2 - x2*y1;  b0 = y1-y2;  c0 = x2-x1;
    a1 = x2*y0 - x0*y2;  b1 = y2-y0;  c1 = x0-x2;
    a2 = x0*y1 - x1*y0;  b2 = y0-y1;  c2 = x1-x0;

    Element_Unstructured_2D::Initialise (nlist);

#ifdef TRI6_STORE_INTFF
    intff.New(6);
    intff = sym_intff * size;
#endif

    if (!subsampling_initialised) {
        Point loc(2);
	int i, j, idx;
        for (i = idx = 0; i < nsample_lin; i++) {
	    loc[0] = (double)i/(double)(nsample_lin-1);
	    for (j = 0; j < nsample_lin-i; j++) {
	        loc[1] = (double)j/(double)(nsample_lin-1);
		absc_sample[idx].New(2);
		absc_sample[idx] = loc;
		idx++;
	    }
	}
	subsampling_initialised = true;
    }
}

int Triangle6::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, "Side index out of range");
    dASSERT(node >= 0 && node < 3, "Node index out of range");
    static int SN[3][3] = {{0,1,3},{1,2,4},{2,0,5}};
    return SN[side][node];
}

Point Triangle6::Local (const NodeList &nlist, const Point& glob) const
{
    // Note: This implementation is more efficient than that of Triangle3,
    // so maybe should be used for Triangle3, but requires realignment of
    // local element

    dASSERT(glob.Dim() == 2, "Invalid point dimension");

    Point loc(2);
    double scale = 1.0/(a0+a1+a2);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1]) * scale;

    return loc;
}

Point Triangle6::NodeLocal (int node) const
{
    Point nloc(2);
    switch (node) {
        case 0:  nloc[0] = nloc[1] = 0.0; break;
        case 1:  nloc[0] = 1.0; nloc[1] = 0.0; break;
        case 2:  nloc[0] = 0.0; nloc[1] = 1.0; break;
        case 3:  nloc[0] = 0.5; nloc[1] = 0.0; break;
        case 4:  nloc[0] = nloc[1] = 0.5; break;
        case 5:  nloc[0] = 0.0; nloc[1] = 0.5; break;
        default: xERROR("Node index out of range"); break;
    }
    return nloc;
}

RVector Triangle6::DirectionCosine (int side, RDenseMatrix &/*jacin*/)
{
    RVector cosin(2);
    switch(side) {
    case 0: cosin[0] = -b2, cosin[1] = -c2; break;
    case 1: cosin[0] = -b0, cosin[1] = -c0; break;
    case 2: cosin[0] = -b1, cosin[1] = -c1; break;
    default: xERROR("Side index out of range");
    }
    return cosin/length(cosin);
}

const RVector &Triangle6::LNormal (int side) const
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

bool Triangle6::LContains (const Point& loc, bool pad) const
{
    dASSERT(loc.Dim() == 2, "Local point must be 2D.");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[1]+EPS >= 0.0 &&
	    loc[0]+loc[1]-EPS <= 1.0);
    } else {
	return (loc[0] >= 0.0 && loc[1] >= 0.0 && loc[0]+loc[1] <= 1.0);
    }
}

RVector Triangle6::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RVector fun(6);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = 4.0*L0*L1;
    fun[4] = 4.0*L1*L2;
    fun[5] = 4.0*L0*L2;
    return fun;
}

void Triangle6::LocalShapeF (const Point &loc, RVector *fun) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    if (fun->Dim() != 6)
        fun->New(6);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    double *f = fun->data_buffer();
    f[0] = L0 * (2.0*L0 - 1.0);
    f[1] = L1 * (2.0*L1 - 1.0);
    f[2] = L2 * (2.0*L2 - 1.0);
    f[3] = 4.0*L0*L1;
    f[4] = 4.0*L1*L2;
    f[5] = 4.0*L0*L2;
}

RDenseMatrix Triangle6::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der(2,6);
    double lx = loc[0], ly = loc[1];
    der(0,0) = der(1,0) = 4.0*(lx+ly) - 3.0;
    der(0,1) = 4.0*lx - 1.0;
    der(1,1) = 0.0;
    der(0,2) = 0.0;
    der(1,2) = 4.0*ly - 1.0;
    der(0,3) = 4.0*(1.0-2.0*lx-ly);
    der(1,3) = -4.0*lx;
    der(0,4) = 4.0*ly;
    der(1,4) = 4.0*lx;
    der(0,5) = -4.0*ly;
    der(1,5) = 4.0*(1.0-lx-2.0*ly);
    return der;
}

RVector Triangle6::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RVector fun(6);
    double scale = 1.0/(2.0*size);
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1]);
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = 4.0*L0*L1;
    fun[4] = 4.0*L1*L2;
    fun[5] = 4.0*L0*L2;
    return fun;
}

RDenseMatrix Triangle6::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der(2,6);
    double scale = 1.0/(2.0*size);
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1]);
    der(0,0) = b0*scale * (4.0*L0-1.0);
    der(1,0) = c0*scale * (4.0*L0-1.0);
    der(0,1) = b1*scale * (4.0*L1-1.0);
    der(1,1) = c1*scale * (4.0*L1-1.0);
    der(0,2) = b2*scale * (4.0*L2-1.0);
    der(1,2) = c2*scale * (4.0*L2-1.0);
    scale *= 4.0;
    der(0,3) = scale * (b0*L1 + b1*L0);
    der(1,3) = scale * (c0*L1 + c1*L0);
    der(0,4) = scale * (b1*L2 + b2*L1);
    der(1,4) = scale * (c1*L2 + c2*L1);
    der(0,5) = scale * (b0*L2 + b2*L0);
    der(1,5) = scale * (c0*L2 + c2*L0);
    return der;
}

double Triangle6::IntF (int i) const
{
    static const double third = 1.0/3.0;
    return (i >= 3 ? size*third : 0);
}

RSymMatrix Triangle6::IntFF () const
{
#ifdef TRI6_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Triangle6::IntFF (int i, int j) const
{
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6);
#ifdef TRI6_STORE_INTFF
    return intff(i,j);
#else
    return full_intff(i,j) * size;
#endif
}

double Triangle6::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6 && k >= 0 && k < 6);
    return full_intfff[i]->Get(j,k) * size;
}

RSymMatrix Triangle6::IntPFF (const RVector &P) const
{
    return (sym_intf0ff * P[Node[0]] +
            sym_intf1ff * P[Node[1]] +
            sym_intf2ff * P[Node[2]] +
            sym_intf3ff * P[Node[3]] +
            sym_intf4ff * P[Node[4]] +
            sym_intf5ff * P[Node[5]]) * size;
}

double Triangle6::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 6 && j >= 0 && j < 6);
    return (full_intf0ff(i,j) * P[Node[0]] +
            full_intf1ff(i,j) * P[Node[1]] +
            full_intf2ff(i,j) * P[Node[2]] +
            full_intf3ff(i,j) * P[Node[3]] +
            full_intf4ff(i,j) * P[Node[4]] +
            full_intf5ff(i,j) * P[Node[5]]) * size;
}

static const RVector intfd0d0 = RVector (6, "12 -2 -2 15 7 15")/180.0;
static const RVector intfd0d1 = RVector (6, "-3 -3 2 -1 -5 -5")/180.0;
static const RVector intfd0d2 = RVector (6, "-3 2 -3 -5 -5 -1")/180.0;
static const RDenseMatrix intfd0d3 = RDenseMatrix (6, 4,
   " 3 18  3 18 \
    -6 -5 -6 -5 \
    -1 -5 -1 -5 \
     8 24  8 24 \
    -8  4 -8  4 \
     4 24  4 24")/180.0;
static const RDenseMatrix intfd0d4 = RDenseMatrix (6, 4,
   " 3  3  3  3 \
    -1 -6 -1 -6 \
    -6 -1 -6 -1 \
     4  8  4  8 \
    -8 -8 -8 -8 \
     8  4  8  4")/180.0;
static const RDenseMatrix intfd0d5 = RDenseMatrix (6, 4,
   " 3 18 3 18 \
    -1 -5 -1 -5 \
    -6 -5 -6 -5 \
     4 24  4 24 \
    -8  4 -8  4 \
     8 24  8 24")/180.0;
static const RVector intfd1d1 = RVector (6, "-2 12 -2 15 15 7")/180.0;
static const RVector intfd1d2 = RVector (6, "2 -3 -3 -5 -1 -5")/180.0;
static const RDenseMatrix intfd1d3 = RDenseMatrix (6, 4,
   "-5 -6 -5 -6 \
    18  3 18  3 \
    -5 -1 -5 -1 \
    24  8 24  8 \
    24  4 24  4 \
     4 -8  4 -8")/180.0;
static const RDenseMatrix intfd1d4 = RDenseMatrix (6, 4,
   "-1 -5 -1 -5 \
     3 18  3 18 \
    -6 -5 -6 -5 \
     4 24  4 24 \
     8 24  8 24 \
    -8  4 -8  4")/180.0;
static const RDenseMatrix intfd1d5 = RDenseMatrix (6, 4,
   "-1 -6 -1 -6 \
     3  3  3  3 \
    -6 -1 -6 -1 \
     4  8  4  8 \
     8  4  8  4 \
    -8 -8 -8 -8")/180.0;
static const RVector intfd2d2 = RVector (6, "-2 -2 12 7 15 15")/180.0;
static const RDenseMatrix intfd2d3 = RDenseMatrix (6, 4,
   "-1 -6 -1 -6 \
    -6 -1 -6 -1 \
     3  3  3  3 \
    -8 -8 -8 -8 \
     8  4  8  4 \
     4  8  4  8")/180.0;
static const RDenseMatrix intfd2d4 = RDenseMatrix (6, 4,
   "-5 -1 -5 -1 \
    -5 -6 -5 -6 \
    18  3 18  3 \
     4 -8  4 -8 \
    24  8 24  8 \
    24  4 24  4")/180.0;
static const RDenseMatrix intfd2d5 = RDenseMatrix (6, 4,
   "-5 -6 -5 -6 \
    -5 -1 -5 -1 \
    18  3 18  3 \
     4 -8  4 -8 \
    24  4 24  4 \
    24  8 24  8")/180.0;
static const RDenseMatrix intfd3d3 = RDenseMatrix (6, 6,
   "-8  0 24 -8  0 24 \
    24  0 -8 24  0 -8 \
    -8 -8 -8 -8 -8 -8 \
    48 64 48 48 64 48 \
    48 32 16 48 32 16 \
    16 32 48 16 32 48")/180.0;
static const RDenseMatrix intfd3d4 = RDenseMatrix (6, 8,
   "-4  0 -8  0 -4  0 -8  0 \
     0 -4 24  0  0 -4 24  0 \
     0  0 -8 -4  0  0 -8 -4 \
    16 16 48 32 16 16 48 32 \
    32 16 48 16 32 16 48 16 \
    16 32 16 16 16 32 16 16")/180.0;
static const RDenseMatrix intfd3d5 = RDenseMatrix (6, 8,
   "-4  0  0 24 -4  0  0 24 \
     0 -4  0 -8  0 -4  0 -8 \
     0  0 -4 -8  0  0 -4 -8 \
    16 16 32 48 16 16 32 48 \
    32 16 16 16 32 16 16 16 \
    16 32 16 48 16 32 16 48")/180.0;
static const RDenseMatrix intfd4d4 = RDenseMatrix (6, 6,
   "-8 -8 -8 -8 -8 -8 \
    -8  0 24 -8  0 24 \
    24  0 -8 24  0 -8 \
    16 32 48 16 32 48 \
    48 64 48 48 64 48 \
    48 32 16 48 32 16")/180.0;
static const RDenseMatrix intfd4d5 = RDenseMatrix (6, 8,
   "-8 -4  0  0 -8 -4  0  0 \
    -8  0 -4  0 -8  0 -4  0 \
    24  0  0 -4 24  0  0 -4 \
    16 16 16 32 16 16 16 32 \
    48 32 16 16 48 32 16 16 \
    48 16 32 16 48 16 32 16")/180.0;
static const RDenseMatrix intfd5d5 = RDenseMatrix (6, 6,
   "-8  0 24 -8  0 24 \
    -8 -8 -8 -8 -8 -8 \
    24  0 -8 24  0 -8 \
    16 32 48 16 32 48 \
    48 32 16 48 32 16 \
    48 64 48 48 64 48")/180.0;

double Triangle6::IntFDD (int i, int j, int k) const
{
    if (j > k) { int tmp=j; j=k; k=tmp; } // symmetry
    double geom;
    switch (j) {
    case 0:
        switch (k) {
	case 0:
	    geom = (b0*b0+c0*c0)/size;
	    return intfd0d0[i] * geom;
	case 1:
	    geom = (b0*b1+c0*c1)/size;
	    return intfd0d1[i] * geom;
	case 2:
	    geom = (b0*b2+c0*c2)/size;
	    return intfd0d2[i] * geom;
	case 3:
	    return (intfd0d3(i,0)*(b0*b0) +
	            intfd0d3(i,1)*(b0*b1) +
	            intfd0d3(i,2)*(c0*c0) +
	            intfd0d3(i,3)*(c0*c1))/size;
	case 4:
	    return (intfd0d4(i,0)*(b0*b1) +
                    intfd0d4(i,1)*(b0*b2) +
		    intfd0d4(i,2)*(c0*c1) +
		    intfd0d4(i,3)*(c0*c2))/size;
	case 5:
	    return (intfd0d5(i,0)*(b0*b0) +
		    intfd0d5(i,1)*(b0*b2) +
		    intfd0d5(i,2)*(c0*c0) +
		    intfd0d5(i,3)*(c0*c2))/size;
	}
    case 1:
        switch (k) {
	case 1:
	    geom = (b1*b1+c1*c1)/size;
	    return intfd1d1[i] * geom;
	case 2:
	    geom = (b1*b2+c1*c2)/size;
	    return intfd1d2[i] * geom;
	case 3:
	    return (intfd1d3(i,0)*(b0*b1) +
		    intfd1d3(i,1)*(b1*b1) +
		    intfd1d3(i,2)*(c0*c1) +
		    intfd1d3(i,3)*(c1*c1))/size;
	case 4:
	    return (intfd1d4(i,0)*(b1*b1) +
		    intfd1d4(i,1)*(b1*b2) +
		    intfd1d4(i,2)*(c1*c1) +
		    intfd1d4(i,3)*(c1*c2))/size;
	case 5:
	    return (intfd1d5(i,0)*(b0*b1) +
		    intfd1d5(i,1)*(b1*b2) +
		    intfd1d5(i,2)*(c0*c1) +
		    intfd1d5(i,3)*(c1*c2))/size;
	}
    case 2:
        switch (k) {
	case 2:
	    geom = (b2*b2+c2*c2)/size;
	    return intfd2d2[i] * geom;
	case 3:
	    return (intfd2d3(i,0)*(b0*b2) +
		    intfd2d3(i,1)*(b1*b2) +
		    intfd2d3(i,2)*(c0*c2) +
		    intfd2d3(i,3)*(c1*c2))/size;
	case 4:
	    return (intfd2d4(i,0)*(b1*b2) +
		    intfd2d4(i,1)*(b2*b2) +
		    intfd2d4(i,2)*(c1*c2) +
		    intfd2d4(i,3)*(c2*c2))/size;
	case 5:
	    return (intfd2d5(i,0)*(b0*b2) +
		    intfd2d5(i,1)*(b2*b2) +
		    intfd2d5(i,2)*(c0*c2) +
		    intfd2d5(i,3)*(c2*c2))/size;
	}
    case 3:
        switch (k) {
	case 3:
	    return (intfd3d3(i,0)*(b0*b0) +
		    intfd3d3(i,1)*(b0*b1) +
		    intfd3d3(i,2)*(b1*b1) +
		    intfd3d3(i,3)*(c0*c0) +
		    intfd3d3(i,4)*(c0*c1) +
		    intfd3d3(i,5)*(c1*c1))/size;
	case 4:
	    return (intfd3d4(i,0)*(b0*b1) +
		    intfd3d4(i,1)*(b1*b1) +
		    intfd3d4(i,2)*(b0*b2) +
		    intfd3d4(i,3)*(b1*b2) +
		    intfd3d4(i,4)*(c0*c1) +
		    intfd3d4(i,5)*(c1*c1) +
		    intfd3d4(i,6)*(c0*c2) +
		    intfd3d4(i,7)*(c1*c2))/size;
	case 5:
	    return (intfd3d5(i,0)*(b0*b0) +
		    intfd3d5(i,1)*(b0*b1) +
		    intfd3d5(i,2)*(b0*b2) +
		    intfd3d5(i,3)*(b1*b2) +
		    intfd3d5(i,4)*(c0*c0) +
		    intfd3d5(i,5)*(c0*c1) +
		    intfd3d5(i,6)*(c0*c2) +
		    intfd3d5(i,7)*(c1*c2))/size;
	}
    case 4:
        switch (k) {
	case 4:
	    return (intfd4d4(i,0)*(b1*b1) +
		    intfd4d4(i,1)*(b1*b2) +
		    intfd4d4(i,2)*(b2*b2) +
		    intfd4d4(i,3)*(c1*c1) +
		    intfd4d4(i,4)*(c1*c2) +
		    intfd4d4(i,5)*(c2*c2))/size;
	case 5:
	    return (intfd4d5(i,0)*(b0*b1) +
		    intfd4d5(i,1)*(b0*b2) +
		    intfd4d5(i,2)*(b1*b2) +
		    intfd4d5(i,3)*(b2*b2) +
		    intfd4d5(i,4)*(c0*c1) +
		    intfd4d5(i,5)*(c0*c2) +
		    intfd4d5(i,6)*(c1*c2) +
		    intfd4d5(i,7)*(c2*c2))/size;
	}
    case 5:
        switch (k) {
	case 5:
	    return (intfd5d5(i,0)*(b0*b0) +
		    intfd5d5(i,1)*(b0*b2) +
		    intfd5d5(i,2)*(b2*b2) +
		    intfd5d5(i,3)*(c0*c0) +
		    intfd5d5(i,4)*(c0*c2) +
		    intfd5d5(i,5)*(c2*c2))/size;
	}
    }
    xERROR("Index out of range");
    return 0;
}

double Triangle6::IntPDD (int i, int j, const RVector &P) const
{
    // there should be more efficient ways to implement this

    double res = 0.0;
    for (int k = 0; k < 6; k++)
        res += P[Node[k]] * IntFDD (k,i,j);
    return res;
}

static const RDenseMatrix sd0_intf0ff = RDenseMatrix (6, 6,
   "39 -3  0 20  0  0 \
    -3 -3  0 -8  0  0 \
     0  0  0  0  0  0 \
    20 -8  0 16  0  0 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0")/420.0;
static const RDenseMatrix sd0_intf1ff = RDenseMatrix (6, 6,
   "-3 -3  0 -8  0  0 \
    -3 39  0 20  0  0 \
     0  0  0  0  0  0 \
    -8 20  0 16  0  0 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0")/420.0;
static const RDenseMatrix sd0_intf2ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd0_intf3ff = RDenseMatrix (6, 6,
   "20 -8  0 16  0  0 \
    -8 20  0 16  0  0 \
     0  0  0  0  0  0 \
    16 16  0 192 0  0 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0")/420.0;
static const RDenseMatrix sd0_intf4ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd0_intf5ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix *sd0_intfff[6] = {
    &sd0_intf0ff,
    &sd0_intf1ff,
    &sd0_intf2ff,
    &sd0_intf3ff,
    &sd0_intf4ff,
    &sd0_intf5ff
};
static const RDenseMatrix sd1_intf0ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd1_intf1ff = RDenseMatrix (6, 6,
   " 0  0  0  0  0  0 \
     0 39 -3  0 20  0 \
     0 -3 -3  0 -8  0 \
     0  0  0  0  0  0 \
     0 20 -8  0 16  0 \
     0  0  0  0  0  0")/420.0;
static const RDenseMatrix sd1_intf2ff = RDenseMatrix (6, 6,
   " 0  0  0  0  0  0 \
     0 -3 -3  0 -8  0 \
     0 -3 39  0 20  0 \
     0  0  0  0  0  0 \
     0 -8 20  0 16  0 \
     0  0  0  0  0  0")/420.0;
static const RDenseMatrix sd1_intf3ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd1_intf4ff = RDenseMatrix (6, 6,
   " 0  0  0  0  0  0 \
     0 20 -8  0 16  0 \
     0 -8 20  0 16  0 \
     0  0  0  0  0  0 \
     0 16 16  0 192 0 \
     0  0  0  0  0  0")/420.0;
static const RDenseMatrix sd1_intf5ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix *sd1_intfff[6] = {
    &sd1_intf0ff,
    &sd1_intf1ff,
    &sd1_intf2ff,
    &sd1_intf3ff,
    &sd1_intf4ff,
    &sd1_intf5ff
};
static const RDenseMatrix sd2_intf0ff = RDenseMatrix (6, 6,
   "39  0 -3  0  0 20 \
     0  0  0  0  0  0 \
    -3  0 -3  0  0 -8 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0 \
    20  0 -8  0  0 16")/420.0;
static const RDenseMatrix sd2_intf1ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd2_intf2ff = RDenseMatrix (6, 6,
   "-3  0 -3  0  0 -8 \
     0  0  0  0  0  0 \
    -3  0 39  0  0 20 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0 \
    -8  0 20  0  0 16")/420.0;
static const RDenseMatrix sd2_intf3ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd2_intf4ff = RDenseMatrix (6, 6); // = 0
static const RDenseMatrix sd2_intf5ff = RDenseMatrix (6, 6,
   "20  0 -8  0  0 16 \
     0  0  0  0  0  0 \
    -8  0 20  0  0 16 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0 \
    16  0 16  0  0 192")/420.0;
static const RDenseMatrix *sd2_intfff[6] = {
    &sd2_intf0ff,
    &sd2_intf1ff,
    &sd2_intf2ff,
    &sd2_intf3ff,
    &sd2_intf4ff,
    &sd2_intf5ff
};

double Triangle6::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    int k;
    if (bndside[0]) {
        double d0 = sqrt (b2*b2 + c2*c2); // length of side 0
        for (k = 0; k < 6; k++)
	    res += P[Node[k]] * d0 * sd0_intfff[k]->Get(i,j);
    }
    if (bndside[1]) {
        double d1 = sqrt (b0*b0 + c0*c0); // length of side 1
        for (k = 0; k < 6; k++)
	    res += P[Node[k]] * d1 * sd1_intfff[k]->Get(i,j);
    }
    if (bndside[2]) {
        double d2 = sqrt (b1*b1 + c1*c1); // length of side 2
        for (k = 0; k < 6; k++)
	    res += P[Node[k]] * d2 * sd2_intfff[k]->Get(i,j);
    }
    return res;
}

double Triangle6::IntFd (int i, int j, int k) const
{
    // Calculate Int [u_i  du_j/dx_k] dr

    static const double scale  = 1.0/30.0;
    static const double scale2 = 2.0*scale;
    static const double scale3 = 3.0*scale;
    static const double scale4 = 4.0*scale;
    static const double scale8 = 8.0*scale;

    double p0, p1, p2;
    switch (k) {
    case 0:  p0 = b0, p1 = b1, p2 = b2; break;
    case 1:  p0 = c0, p1 = c1, p2 = c2; break;
    default: xERROR ("Index out of range"); return 0.0;
    }

    switch (i) {
    case 0:
	switch (j) {
	case 0: return  p0 * scale2;
	case 1: return -p1 * scale;
	case 2: return -p2 * scale;
	case 3: return (-p0 + 2.0*p1) * scale;
	case 4: return (-p1 - p2) * scale;
	case 5: return (-p0 + 2.0*p2) * scale;
	default: xERROR ("Index out of range"); return 0.0;
	}
    case 1:
	switch (j) {
	case 0: return -p0 * scale;
	case 1: return  p1 * scale2;
	case 2: return -p2 * scale;
	case 3: return (2.0*p0 - p1) * scale;
	case 4: return (-p1 + 2.0*p2) * scale;
	case 5: return (-p0 - p2) * scale;
	default: xERROR ("Index out of range"); return 0.0;
	}
    case 2:
	switch (j) {
	case 0: return -p0 * scale;
	case 1: return -p1 * scale;
	case 2: return  p2 * scale2;
	case 3: return (-p0 - p1) * scale;
	case 4: return (2.0*p1 - p2) * scale;
	case 5: return (2.0*p0 - p2) * scale;
	default: xERROR ("Index out of range"); return 0.0;
	}
    case 3:
	switch (j) {
	case 0: return  p0 * scale3;
	case 1: return  p1 * scale3;
	case 2: return -p2 * scale;
	case 3: return (p0 + p1) * scale8;
	case 4: return (p1 + 2.0*p2) * scale4;
	case 5: return (p0 + 2.0*p2) * scale4;
	default: xERROR ("Index out of range"); return 0.0;
	}
    case 4:
	switch (j) {
	case 0: return -p0 * scale;
	case 1: return  p1 * scale3;
	case 2: return  p2 * scale3;
	case 3: return (2.0*p0 + p1) * scale4;
	case 4: return (p1 + p2) * scale8;
	case 5: return (2.0*p0 + p2) * scale4;
	default: xERROR ("Index out of range"); return 0.0;
	}
    case 5:
	switch (j) {
	case 0: return  p0 * scale3;
	case 1: return -p1 * scale;
	case 2: return  p2 * scale3;
	case 3: return (p0 + 2.0*p1) * scale4;
	case 4: return (2.0*p1 + p2) * scale4;
	case 5: return (p0 + p2) * scale8;
	default: xERROR ("Index out of range"); return 0.0;
	}
    default: xERROR ("Index out of range"); return 0.0;
    }
}

double Triangle6::IntPd (const RVector &P, int j, int k) const
{
    double sum = 0.0;
    for (int i = 0; i < 6; i++) 
	sum += P[Node[i]] * IntFd (i, j, k);
    return sum;
}

RSymMatrix Triangle6::Intdd () const
{
    RSymMatrix MDD(12,12);
    MDD(0,0) = 3.0*b0*b0;

    MDD(1,0) = 3.0*b0*c0;
    MDD(1,1) = 3.0*c0*c0;

    MDD(2,0) = -b0*b1;
    MDD(2,1) = -b1*c0;
    MDD(2,2) = 3.0*b1*b1;

    MDD(3,0) = -b0*c1;
    MDD(3,1) = -c0*c1;
    MDD(3,2) = 3.0*b1*c1;
    MDD(3,3) = 3.0*c1*c1;

    MDD(4,0) = -b0*b2;
    MDD(4,1) = -b2*c0;
    MDD(4,2) = -b1*b2;
    MDD(4,3) = -b2*c1;
    MDD(4,4) = 3.0*b2*b2;

    MDD(5,0) = -b0*c2;
    MDD(5,1) = -c0*c2;
    MDD(5,2) = -b1*c2;
    MDD(5,3) = -c1*c2;
    MDD(5,4) = 3.0*b2*c2;
    MDD(5,5) = 3.0*c2*c2;

    MDD(6,0) = 4.0*b0*b1;
    MDD(6,1) = 4.0*b1*c0;
    MDD(6,2) = 4.0*b0*b1;
    MDD(6,3) = 4.0*b0*c1;
    MDD(6,4) = 0.0;
    MDD(6,5) = 0.0;
    MDD(6,6) = 8.0*(b0*b0 + b0*b1 + b1*b1);

    MDD(7,0) = 4.0*b0*c1;
    MDD(7,1) = 4.0*c0*c1;
    MDD(7,2) = 4.0*b1*c0;
    MDD(7,3) = 4.0*c0*c1;
    MDD(7,4) = 0.0;
    MDD(7,5) = 0.0;
    MDD(7,6) = 4.0*(b0*(2.0*c0+c1) + b1*(c0+2.0*c1));
    MDD(7,7) = 8.0*(c0*c0 + c0*c1 + c1*c1);

    MDD(8,0) = 0.0;
    MDD(8,1) = 0.0;
    MDD(8,2) = 4.0*b1*b2;
    MDD(8,3) = 4.0*b2*c1;
    MDD(8,4) = 4.0*b1*b2;
    MDD(8,5) = 4.0*b1*c2;
    MDD(8,6) = 4.0*(b1*(b1+b2) + b0*(b1+2.0*b2));
    MDD(8,7) = 4.0*(b1*(c0+c1) + b2*(2.0*c0+c1));
    MDD(8,8) = 8.0*(b1*b1 + b1*b2 + b2*b2);

    MDD(9,0) = 0.0;
    MDD(9,1) = 0.0;
    MDD(9,2) = 4.0*b1*c2;
    MDD(9,3) = 4.0*c1*c2;
    MDD(9,4) = 4.0*b2*c1;
    MDD(9,5) = 4.0*c1*c2;
    MDD(9,6) = 4.0*(b1*(c1+c2) + b0*(c1+2.0*c2));
    MDD(9,7) = 4.0*(c1*(c1+c2) + c0*(c1+2.0*c2));
    MDD(9,8) = 4.0*(b1*(2.0*c1+c2) + b2*(c1+2.0*c2));
    MDD(9,9) = 8.0*(c1*c1 + c1*c2 + c2*c2);

    MDD(10,0) = 4.0*b0*b2;
    MDD(10,1) = 4.0*b2*c0;
    MDD(10,2) = 0.0;
    MDD(10,3) = 0.0;
    MDD(10,4) = 4.0*b0*b2;
    MDD(10,5) = 4.0*b0*c2;
    MDD(10,6) = 4.0*(b0*b0 + 2.0*b1*b2 + b0*(b1+b2));
    MDD(10,7) = 4.0*(b0*(c0+c1) + b2*(c0+2.0*c1));
    MDD(10,8) = 4.0*(b2*(b1+b2) + b0*(2.0*b1+b2));
    MDD(10,9) = 4.0*(b2*(c1+c2) + b0*(2.0*c1+c2));
    MDD(10,10) = 8.0*(b0*b0 + b0*b2 + b2*b2);

    MDD(11,0) = 4.0*b0*c2;
    MDD(11,1) = 4.0*c0*c2;
    MDD(11,2) = 0.0;
    MDD(11,3) = 0.0;
    MDD(11,4) = 4.0*b2*c0;
    MDD(11,5) = 4.0*c0*c2;
    MDD(11,6) = 4.0*(b0*(c0+c2) + b1*(c0+2.0*c2));
    MDD(11,7) = 4.0*(c0*c0 + 2.0*c1*c2 + c0*(c1+c2));
    MDD(11,8) = 4.0*(b2*(c0+c2) + b1*(2.0*c0+c2));
    MDD(11,9) = 4.0*(c2*(c1+c2) + c0*(2.0*c1+c2));
    MDD(11,10) = 4.0*(b0*(2.0*c0+c2) + b2*(c0+2.0*c2));
    MDD(11,11) = 8.0*(c0*c0 + c0*c2 + c2*c2);

    MDD /= 12.0*size;
    return MDD;
}

int Triangle6::GetLocalSubsampleAbsc (const Point *&absc) const
{
    absc = absc_sample;
    return nsample_tot;
}

double Triangle6::ComputeSize (const NodeList &nlist) const
{
    return 0.5 * (a0+a1+a2);
}

RVector Triangle6::IntFD (int i, int j) const
{
    dASSERT(i >= 0 && i < 6, "Parameter 1 index out of range");
    dASSERT(j >= 0 && j < 6, "Parameter 2 index out of range");
    RVector fd(2);
    fd[0] = intfd_0(i,j);
    fd[1] = intfd_1(i,j);
    return fd;
}

void Triangle6::ComputeIntFD (const NodeList &nlist)
{
    static const double i30 = 1.0/30.0;

    intfd_0.New(6,6);
    intfd_1.New(6,6);
    intfd_0(0,0) =  2.0*b0;         intfd_1(0,0) =  2.0*c0;
    intfd_0(0,1) = -b1;             intfd_1(0,1) = -c1;
    intfd_0(0,2) = -b2;             intfd_1(0,2) = -c2;
    intfd_0(0,3) = -b0+2.0*b1;      intfd_1(0,3) = -c0+2.0*c1;
    intfd_0(0,4) = -b1-b2;          intfd_1(0,4) = -c1-c2;
    intfd_0(0,5) = -b0+2.0*b2;      intfd_1(0,5) = -c0+2.0*c2;
    intfd_0(1,0) = -b0;             intfd_1(1,0) = -c0;
    intfd_0(1,1) =  2.0*b1;         intfd_1(1,1) =  2.0*c1;
    intfd_0(1,2) = -b2;             intfd_1(1,2) = -c2;
    intfd_0(1,3) =  2.0*b0-b1;      intfd_1(1,3) =  2.0*c0-c1;
    intfd_0(1,4) = -b1+2.0*b2;      intfd_1(1,4) = -c1+2.0*c2;
    intfd_0(1,5) = -b0-b2;          intfd_1(1,5) = -c0-c2;
    intfd_0(2,0) = -b0;             intfd_1(2,0) = -c0;
    intfd_0(2,1) = -b1;             intfd_1(2,1) = -c1;
    intfd_0(2,2) =  2.0*b2;         intfd_1(2,2) =  2.0*c2;
    intfd_0(2,3) = -b0-b1;          intfd_1(2,3) = -c0-c1;
    intfd_0(2,4) =  2.0*b1-b2;      intfd_1(2,4) =  2.0*c1-c2;
    intfd_0(2,5) =  2.0*b0-b2;      intfd_1(2,5) =  2.0*c0-c2;
    intfd_0(3,0) =  3.0*b0;         intfd_1(3,0) =  3.0*c0;
    intfd_0(3,1) =  3.0*b1;         intfd_1(3,1) =  3.0*c1;
    intfd_0(3,2) = -b2;             intfd_1(3,2) = -c2;
    intfd_0(3,3) =  8.0*(b0+b1);    intfd_1(3,3) =  8.0*(c0+c1);
    intfd_0(3,4) =  4.0*b1+8.0*b2;  intfd_1(3,4) =  4.0*c1+8.0*c2;
    intfd_0(3,5) =  4.0*b0+8.0*b2;  intfd_1(3,5) =  4.0*c0+8.0*c2;
    intfd_0(4,0) = -b0;             intfd_1(4,0) = -c0;
    intfd_0(4,1) =  3.0*b1;         intfd_1(4,1) =  3.0*c1;
    intfd_0(4,2) =  3.0*b2;         intfd_1(4,2) =  3.0*c2;
    intfd_0(4,3) =  8.0*b0+4.0*b1;  intfd_1(4,3) =  8.0*c0+4.0*c1;
    intfd_0(4,4) =  8.0*(b1+b2);    intfd_1(4,4) =  8.0*(c1+c2);
    intfd_0(4,5) =  8.0*b0+4.0*b2;  intfd_1(4,5) =  8.0*c0+4.0*c2;
    intfd_0(5,0) =  3.0*b0;         intfd_1(5,0) =  3.0*c0;
    intfd_0(5,1) = -b1;             intfd_1(5,1) = -c1;
    intfd_0(5,2) =  3.0*b2;         intfd_1(5,2) =  3.0*c2;
    intfd_0(5,3) =  4.0*b0+8.0*b1;  intfd_1(5,3) =  4.0*c0+8.0*c1;
    intfd_0(5,4) =  8.0*b1+4.0*b2;  intfd_1(5,4) =  8.0*c1+4.0*c2;
    intfd_0(5,5) =  8.0*(b0+b2);    intfd_1(5,5) =  8.0*(c0+c2);

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
	    intfd_0(i,j) *= i30;
	    intfd_1(i,j) *= i30;
	}
}

RSymMatrix Triangle6::ComputeIntDD (const NodeList &nlist) const
{
    static RSymMatrix dd(6);
    double bc00, bc01, bc11, bc02, bc12, bc22;
    double scale = 1.0/(12.0*size), scale4, scale8;

    dd(0,0) =  3.0*scale * (bc00 = b0*b0 + c0*c0);
    dd(1,0) =     -scale * (bc01 = b0*b1 + c0*c1);
    dd(1,1) =  3.0*scale * (bc11 = b1*b1 + c1*c1);
    dd(2,0) =     -scale * (bc02 = b0*b2 + c0*c2);
    dd(2,1) =     -scale * (bc12 = b1*b2 + c1*c2);
    dd(2,2) =  3.0*scale * (bc22 = b2*b2 + c2*c2);
    dd(3,0) =
    dd(3,1) =  (scale4 = 4.0*scale) * bc01;
    // dd(3,2) =  0.0;
    dd(3,3) =  (scale8 = 8.0*scale) * (bc00 + bc01 + bc11);
    // dd(4,0) =  0.0;
    dd(4,1) =
    dd(4,2) =  scale4 * bc12;
    dd(4,3) =  scale4 * (bc01 + 2.0*bc02 + bc11 + bc12);
    dd(4,4) =  scale8 * (bc11 + bc12 + bc22);
    dd(5,0) =
    dd(5,2) =  scale4 * bc02;
    // dd(5,1) =  0.0;
    dd(5,3) =  scale4 * (bc00 + bc02 + bc01 + 2.0*bc12);
    dd(5,4) =  scale4 * (2.0*bc01 + bc12 + bc02 + bc22);
    dd(5,5) =  scale8 * (bc00 + bc02 + bc22);

    return dd;
}

double Triangle6::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 6, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 2: out of range");

    double f = bndintf(sd,i);
    if (!f) return f;
    switch (sd) {
    case 0: return f * hypot (c2, b2);
    case 1: return f * hypot (c0, b0);
    case 2: return f * hypot (c1, b1);
    }
    return 0.0;
}

double Triangle6::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 6, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 6, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 3: out of range");
    
    switch (sd) {
    case 0:
	return hypot (c2, b2) * sym_bndintff_sd0(i,j);
    case 1:
	return hypot (c0, b0) * sym_bndintff_sd1(i,j);
    case 2:
	return hypot (c1, b1) * sym_bndintff_sd2(i,j);
    default:
	xERROR("Invalid side index");
	return 0.0;
    }
}

RSymMatrix Triangle6::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(6);
    if (!bndel) return bff;  // not a boundary element -> return zero matrix

    if (bndside[0]) {
        double d0 = nlist[Node[0]].Dist (nlist[Node[1]]);
	bff += sym_bndintff_sd0 * d0;
    }
    if (bndside[1]) {
        double d1 = nlist[Node[1]].Dist (nlist[Node[2]]);
	bff += sym_bndintff_sd1 * d1;
    }
    if (bndside[2]) {
        double d2 = nlist[Node[2]].Dist (nlist[Node[0]]);
	bff += sym_bndintff_sd2 * d2;
    }
    return bff;
}

RVector Triangle6::BndIntFDelta (int side, const Surface *surf,
    const RVector &pos, const NodeList &nlist) const
{
    double d0, d1, s, g;
    dASSERT(surf->Dimension() == 2, "Wrong surface dimension");
    Surface2D *surf2D = (Surface2D*)surf;
    int n0 = Node[SideNode (side, 0)];
    int n1 = Node[SideNode (side, 1)];
    int n2 = Node[SideNode (side, 2)]; // centre node
    int pdim = surf2D->ParamDim();
    RVector prm0(pdim);  surf2D->Point2Param (nlist[n0], prm0);
    RVector prm1(pdim);  surf2D->Point2Param (nlist[n1], prm1);
    RVector prm2(pdim);  surf2D->Point2Param (nlist[n2], prm2);
    RVector res(3);
    d0 = surf2D->ChordDiff (pos, prm0);
    d1 = surf2D->ChordDiff (pos, prm1);
    s  = surf2D->ChordDiff (prm1, prm0);
    g  = surf2D->ChordDiff (prm2, prm0);
    if (s*d0 >= 0.0 && d0*d1 <= 0.0) { // otherwise outside support
        res[0] = ( d0*d0 - g*s*d0 + 1.0)/(g*s);
	res[1] = (-d0*d0 + g*d0)/((g-s)*s);
	res[2] = ( d0*d0 - s*d0)/((g-s)*g);
    }
    return res;
}

RVector Triangle6::BndIntFCos (int side, const Surface *surf,
    const RVector &cntcos, double sigma, double sup, const NodeList &nlist)
    const
{
    dASSERT(surf->Dimension() == 2, "Wrong surface dimension");
    Surface2D *surf2D = (Surface2D*)surf;
    double d0, s, g, tmp;
    bool flip;
    int n0 = Node[SideNode (side, 0)];
    int n1 = Node[SideNode (side, 1)];
    int n2 = Node[SideNode (side, 2)];
    int pdim = surf2D->ParamDim();
    RVector prm0(pdim);  surf2D->Point2Param (nlist[n0], prm0);
    RVector prm1(pdim);  surf2D->Point2Param (nlist[n1], prm1);
    RVector prm2(pdim);  surf2D->Point2Param (nlist[n2], prm2);
    s  = surf2D->ChordDiff (prm1, prm0);
    if ((flip = (s < 0))) {
        g = surf2D->ChordDiff (prm2, prm1);
	d0 = surf2D->ChordDiff (cntcos, prm1);
	s = -s;
    } else {
        g  = surf2D->ChordDiff (prm2, prm0);
	d0 = surf2D->ChordDiff (cntcos, prm0);
    }
    RVector res(3);
    Integrate_x2_Cosine (d0, sup, sigma, s, g, res[0], res[1], res[2]);
    if (flip) tmp = res[0], res[0] = res[1], res[1] = tmp;
    return res;
}

RVector Triangle6::BndIntFCos (int side, const RVector &cntcos, double a,
    const NodeList &nlist) const
{
    // this version parametrises the surface simply as a function of
    // the geometric distance from the centre of the cosine. This
    // will cause problems for curved surfaces and large cosine support
    // radius a
    
    bool flip = false;
    int n0 = Node[SideNode (side,0)];  // end-node 0
    int n1 = Node[SideNode (side,1)];  // end-node 1
    double d0 = length (nlist[n0] - cntcos);
    double d1 = length (nlist[n1] - cntcos);

    // hack: if angle subtended by the nodes as seen from cosine centre > Pi/2,
    // we assume that they bracket the cosine centre and flip one of the
    // distance signs
    double cosa = dot (nlist[n0]-cntcos, nlist[n1]-cntcos) / (d0*d1);
    if (cosa < 0.0) d0 = -d0;

    RVector res(3);
    if (d0 > d1) {
	double tmp = d0; d0 = d1; d1 = tmp;
	flip = true;
    }
    Integrate_u_Cosine (0.0, a, d0, d1, res[0], res[1], res[2]);
    if (flip) {
	double tmp = res[0]; res[0] = res[2]; res[2] = tmp;
    }
    return res;
}

void Integrate_x2_Cosine (double d, double dmax, double sigma, double s,
    double g, double &res0, double &res1, double &res2)
{
    double arg0, arg1, sin0, sin1, cos0, cos1;

    if (d >= 0.0) {
        if (d-s >= dmax) {
	    res0 = res1 = res2 = 0.0;
	    return;
	}
    } else {
        if (-d >= dmax) {
	    res0 = res1 = res2 = 0.0;
	    return;
	}
    }

    // find the integration limits
    double x0 = ::max (d-dmax, 0.0);
    double x1 = ::min (d+dmax, s);

    arg0 = sigma * (d-x0);
    arg1 = sigma * (d-x1);
    sin0 = sin(arg0), cos0 = cos(arg0);
    sin1 = sin(arg1), cos1 = cos(arg1);

    res0 = 1.0/(g*s*sigma*sigma*sigma) *
      (sigma*(g+s-x0-x0)*cos0 - sigma*(g+s-x1-x1)*cos1 -
       2.0*sin0 + g*s*sigma*sigma*sin0 - g*sigma*sigma*x0*sin0 -
       s*sigma*sigma*x0*sin0 + sigma*sigma*x0*x0*sin0 +
       2.0*sin1 - g*s*sigma*sigma*sin1 + g*sigma*sigma*x1*sin1 +
       s*sigma*sigma*x1*sin1 - sigma*sigma*x1*x1*sin1);

    res1 = 1.0/((g-s)*s*sigma*sigma*sigma) *
      (-sigma*(g-x0-x0)*cos0 + sigma*(g-x1-x1)*cos1 +
       2.0*sin0 + g*sigma*sigma*x0*sin0 -
       sigma*sigma*x0*x0*sin0 - 2.0*sin1 -
       g*sigma*sigma*x1*sin1 + sigma*sigma*x1*x1*sin1);

    res2 = 1.0/(g*(g-s)*sigma*sigma*sigma) *
      (sigma*(s-x0-x0)*cos0 - sigma*(s-x1-x1)*cos1 -
       2.0*sin0 - s*sigma*sigma*x0*sin1 +
       sigma*sigma*x0*x0*sin0 + 2.0*sin1 +
       s*sigma*sigma*sin1 - sigma*sigma*x1*x1*sin1);
}

// Integration of product of shape function and cosine
// Cosine is centered at d and has radius of support a (area of positivity)
// shape function nodes are at x0 and x1 (x0 < x1)

void Integrate_u_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_um, double &int_cos_u1)
{
    if (x0 > d+a || x1 < d-a) {
	// case 1: interval completely outside support
	int_cos_u0 = 0.0;
	int_cos_um = 0.0;
	int_cos_u1 = 0.0;
    } else if (x0 >= d-a && x1 <= d+a) {
	// case 2: support of cosine spans full interval
	double arg1 = 2.0*a / (Pi*Pi*Pi * (x0-x1)*(x0-x1));
	double arg2 = Pi * (d-x0) / (2.0*a);
	double arg3 = Pi * (d-x1) / (2.0*a);
	int_cos_u0 = arg1 * (-6.0*a*Pi*(x0-x1)*cos(arg2) -
			     2.0*a*Pi*(x0-x1)*cos(arg3) -
			     16.0*a*a*sin(arg2) +
			     Pi*Pi*x0*x0*sin(arg2) -
			     2.0*Pi*Pi*x0*x1*sin(arg2) +
			     Pi*Pi*x1*x1*sin(arg2) +
			     16.0*a*a*sin(arg3));
	int_cos_um = arg1 * (16.0*a*cos(Pi*(2.0*d-x0-x1)/(4.0*a)) *
			     (Pi*(x0-x1)*cos(Pi*(x0-x1)/(4.0*a)) -
			      4.0*a*sin(Pi*(x0-x1)/(4.0*a))));
	int_cos_u1 =-arg1 * (2.0*a*Pi*(x0-x1)*cos(arg2) +
			     6.0*a*Pi*(x0-x1)*cos(arg3) +
			     16.0*a*a*sin(arg2) -
			     16.0*a*a*sin(arg3) +
			     Pi*Pi*x0*x0*sin(arg3) -
			     2.0*Pi*Pi*x0*x1*sin(arg3) +
			     Pi*Pi*x1*x1*sin(arg3));
    } else if (x0 < d-a && x1 > d+a) {
	// case 3: support of cosine is subset of interval
	double arg1 = 4.0*a / (Pi*Pi*Pi * (x1-x0)*(x1-x0));
	int_cos_u0 = arg1 * (-16.0*a*a +
			     2.0*a*a*Pi*Pi +
			     2.0*d*d*Pi*Pi -
			     d*Pi*Pi*x0 -
			     3.0*d*Pi*Pi*x1 +
			     Pi*Pi*x0*x1 +
			     Pi*Pi*x1*x1);
	int_cos_um = arg1*4.0 * (8.0*a*a -
			     a*a*Pi*Pi -
			     d*d*Pi*Pi +
			     d*Pi*Pi*x0 +
			     d*Pi*Pi*x1 -
			     Pi*Pi*x0*x1);
	int_cos_u1 = arg1 * (-16.0*a*a +
			     2.0*a*a*Pi*Pi +
			     2.0*d*d*Pi*Pi -
			     3.0*d*Pi*Pi*x0 +
			     Pi*Pi*x0*x0 -
			     d*Pi*Pi*x1 +
			     Pi*Pi*x0*x1);
    } else if (x0 >= d-a) {
	// case 4: interval exceeds support of cosine to the right
	double arg1 = 2.0*a / (Pi*Pi*Pi * (x0-x1)*(x0-x1));
	double arg2 = Pi * (d-x0) / (2.0*a);
	int_cos_u0 = arg1 * (-16.0*a*a +
			     2.0*a*a*Pi*Pi +
			     4.0*a*d*Pi*Pi +
			     2.0*d*d*Pi*Pi -
			     a*Pi*Pi*x0 -
			     d*Pi*Pi*x0 -
			     3.0*a*Pi*Pi*x1 -
			     3.0*d*Pi*Pi*x1 +
			     Pi*Pi*x0*x1 +
			     Pi*Pi*x1*x1 -
			     6.0*a*Pi*(x0-x1)*cos(arg2) +
			     (-16.0*a*a+Pi*Pi*(x0-x1)*(x0-x1))*sin(arg2));
	int_cos_um = arg1*4.0 * (8.0*a*a -
			     a*a*Pi*Pi -
			     2.0*a*d*Pi*Pi -
			     d*d*Pi*Pi +
			     a*Pi*Pi*x0 +
			     d*Pi*Pi*x0 +
			     a*Pi*Pi*x1 +
			     d*Pi*Pi*x1 -
			     Pi*Pi*x0*x1 +
			     2.0*a*Pi*(x0-x1)*cos(arg2) +
			     8.0*a*a*sin(arg2));
	int_cos_u1 = -arg1 * (16.0*a*a -
			     2.0*a*a*Pi*Pi -
			     4.0*a*d*Pi*Pi -
			     2.0*d*d*Pi*Pi +
			     3.0*a*Pi*Pi*x0 +
			     3.0*d*Pi*Pi*x0 -
			     Pi*Pi*x0*x0 +
			     a*Pi*Pi*x1 +
			     d*Pi*Pi*x1 -
			     Pi*Pi*x0*x1 +
			     2.0*a*Pi*(x0-x1)*cos(arg2) +
			     16.0*a*a*sin(arg2));
    } else {
	// case 5: interval exceeds support of cosine to the left
	double arg1 = 2.0*a / (Pi*Pi*Pi * (x0-x1)*(x0-x1));
	double arg2 = Pi * (d-x1) / (2.0*a);
	int_cos_u0 = arg1 * (-16.0*a*a +
			     2.0*a*a*Pi*Pi -
			     4.0*a*d*Pi*Pi +
			     2.0*d*d*Pi*Pi +
			     a*Pi*Pi*x0 -
			     d*Pi*Pi*x0 +
			     3.0*a*Pi*Pi*x1 -
			     3.0*d*Pi*Pi*x1 +
			     Pi*Pi*x0*x1 +
			     Pi*Pi*x1*x1 -
			     2.0*a*Pi*(x0-x1)*cos(arg2) +
			     16.0*a*a*sin(arg2));
	int_cos_um = -arg1*4.0 * (-8.0*a*a +
			     a*a*Pi*Pi -
			     2.0*a*d*Pi*Pi +
			     d*d*Pi*Pi +
			     a*Pi*Pi*x0 -
			     d*Pi*Pi*x0 +
			     a*Pi*Pi*x1 -
			     d*Pi*Pi*x1 +
			     Pi*Pi*x0*x1 -
			     2.0*a*Pi*(x0-x1)*cos(arg2) +
			     8.0*a*a*sin(arg2));
	int_cos_u1 = arg1 * (-16.0*a*a +
			     2.0*a*a*Pi*Pi -
			     4.0*a*d*Pi*Pi +
			     2.0*d*d*Pi*Pi +
			     3.0*a*Pi*Pi*x0 -
			     3.0*d*Pi*Pi*x0 +
			     Pi*Pi*x0*x0 +
			     a*Pi*Pi*x1 -
			     d*Pi*Pi*x1 +
			     Pi*Pi*x0*x1 -
			     6.0*a*Pi*(x0-x1)*cos(arg2) +
			     (16.0*a*a - Pi*Pi*(x0-x1)*(x0-x1))*sin(arg2));
    }
}

int Triangle6::Intersection (const Point &p1, const Point &p2,
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
