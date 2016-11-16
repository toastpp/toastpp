// ==========================================================================
// Module libfe
// File tet10.cc
// Definition of class Tetrahedron10
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "tri_qr.h"

// some global constants

static bool subsampling_initialised = false;
static const int nsample_lin = NSUBSAMPLE; // from toastdef.h
static const int nsample_tot = (((nsample_lin+3)*nsample_lin+2)*nsample_lin)/6;
static Point absc_sample[nsample_tot];

static const RSymMatrix sym_intff = RSymMatrix (10,
   " 6 \
     1  6 \
     1  1  6 \
     1  1  1  6 \
    -4 -4 -6 -6 32 \
    -4 -6 -4 -6 16 32 \
    -4 -6 -6 -4 16 16 32 \
    -6 -4 -4 -6 16 16  8 32 \
    -6 -4 -6 -4 16  8 16 16 32 \
    -6 -6 -4 -4  8 16 16 16 16 32") * (1.0/420.0);

static const RDenseMatrix full_intff = RDenseMatrix (10, 10,
   " 6  1  1  1 -4 -4 -4 -6 -6 -6 \
     1  6  1  1 -4 -6 -6 -4 -4 -6 \
     1  1  6  1 -6 -4 -6 -4 -6 -4 \
     1  1  1  6 -6 -6 -4 -6 -4 -4 \
    -4 -4 -6 -6 32 16 16 16 16  8 \
    -4 -6 -4 -6 16 32 16 16  8 16 \
    -4 -6 -6 -4 16 16 32  8 16 16 \
    -6 -4 -4 -6 16 16  8 32 16 16 \
    -6 -4 -6 -4 16  8 16 16 32 16 \
    -6 -6 -4 -4  8 16 16 16 16 32") * (1.0/420.0);

// below are the 10 planes of the IntFFF tensor
static const RDenseMatrix full_intf0ff = RDenseMatrix (10,10,
   "18 -6 -6 -6  24  24  24  12  12  12 \
    -6 -6 -1 -1   0   6   6   6   6   8 \
    -6 -1 -6 -1   6   0   6   6   8   6 \
    -6 -1 -1 -6   6   6   0   8   6   6 \
    24  0  6  6 -24 -12 -12 -24 -24 -12 \
    24  6  0  6 -12 -24 -12 -24 -12 -24 \
    24  6  6  0 -12 -12 -24 -12 -24 -24 \
    12  6  6  8 -24 -24 -12 -40 -20 -20 \
    12  6  8  6 -24 -12 -24 -20 -40 -20 \
    12  8  6  6 -12 -24 -24 -20 -20 -40") / 7560.0;
static const RDenseMatrix full_intf1ff = RDenseMatrix (10,10,
   "-6 -6 -1 -1   0   6   6   6   6   8 \
    -6 18 -6 -6  24  12  12  24  24  12 \
    -1 -6 -6 -1   6   6   8   0   6   6 \
    -1 -6 -1 -6   6   8   6   6   0   6 \
     0 24  6  6 -24 -24 -24 -12 -12 -12 \
     6 12  6  8 -24 -40 -20 -24 -12 -20 \
     6 12  8  6 -24 -20 -40 -12 -24 -20 \
     6 24  0  6 -12 -24 -12 -24 -12 -24 \
     6 24  6  0 -12 -12 -24 -12 -24 -24 \
     8 12  6  6 -12 -20 -20 -24 -24 -40") / 7560.0;
static const RDenseMatrix full_intf2ff = RDenseMatrix (10,10,
   "-6 -1 -6 -1  6  0  6  6  8  6 \
    -1 -6 -6 -1  6  6  8  0  6  6 \
    -6 -6 18 -6 12 24 12 24 12 24 \
    -1 -1 -6 -6  8  6  6  6  6  0 \
     6  6 12  8 -40 -24 -20 -24 -20 -12 \
     0  6 24  6 -24 -24 -24 -12 -12 -12 \
     6  8 12  6 -20 -24 -40 -12 -20 -24 \
     6  0 24  6 -24 -12 -12 -24 -24 -12 \
     8  6 12  6 -20 -12 -20 -24 -40 -24 \
     6  6 24  0 -12 -12 -24 -12 -24 -24") / 7560.0;
static const RDenseMatrix full_intf3ff = RDenseMatrix (10,10,
   "-6 -1 -1 -6  6  6  0  8  6  6 \
    -1 -6 -1 -6  6  8  6  6  0  6 \
    -1 -1 -6 -6  8  6  6  6  6  0 \
    -6 -6 -6 18 12 12 24 12 24 24 \
     6  6  8 12 -40 -20 -24 -20 -24 -12 \
     6  8  6 12 -20 -40 -24 -20 -12 -24 \
     0  6  6 24 -24 -24 -24 -12 -12 -12 \
     8  6  6 12 -20 -20 -12 -40 -24 -24 \
     6  0  6 24 -24 -12 -12 -24 -24 -12 \
     6  6  0 24 -12 -24 -12 -24 -12 -24") / 7560.0;
static const RDenseMatrix full_intf4ff = RDenseMatrix (10,10,
   "24  0  6  6 -24 -12 -12 -24 -24 -12 \
     0 24  6  6 -24 -24 -24 -12 -12 -12 \
     6  6 12  8 -40 -24 -20 -24 -20 -12 \
     6  6  8 12 -40 -20 -24 -20 -24 -12 \
    -24 -24 -40 -40 288 96 96 96 96 32 \
    -12 -24 -24 -20 96 96 48 64 32 32 \
    -12 -24 -20 -24 96 48 96 32 64 32 \
    -24 -12 -24 -20 96 64 32 96 48 32 \
    -24 -12 -20 -24 96 32 64 48 96 32 \
    -12 -12 -12 -12 32 32 32 32 32 32") / 7560.0;
static const RDenseMatrix full_intf5ff = RDenseMatrix (10,10,
   "24  6  0  6 -12 -24 -12 -24 -12 -24 \
     6 12  6  8 -24 -40 -20 -24 -12 -20 \
     0  6 24  6 -24 -24 -24 -12 -12 -12 \
     6  8  6 12 -20 -40 -24 -20 -12 -24 \
    -12 -24 -24 -20 96 96 48 64 32 32 \
    -24 -40 -24 -40 96 288 96 96 32 96 \
    -12 -20 -24 -24 48 96 96 32 32 64 \
    -24 -24 -12 -20 64 96 32 96 32 48 \
    -12 -12 -12 -12 32 32 32 32 32 32 \
    -24 -20 -12 -24 32 96 64 48 32 96") / 7560.0;
static const RDenseMatrix full_intf6ff = RDenseMatrix (10,10,
   "24  6  6  0 -12 -12 -24 -12 -24 -24 \
     6 12  8  6 -24 -20 -40 -12 -24 -20 \
     6  8 12  6 -20 -24 -40 -12 -20 -24 \
     0  6  6 24 -24 -24 -24 -12 -12 -12 \
    -12 -24 -20 -24 96 48 96 32 64 32 \
    -12 -20 -24 -24 48 96 96 32 32 64 \
    -24 -40 -40 -24 96 96 288 32 96 96 \
    -12 -12 -12 -12 32 32 32 32 32 32 \
    -24 -24 -20 -12 64 32 96 32 96 48 \
    -24 -20 -24 -12 32 64 96 32 48 96") / 7560.0;
static const RDenseMatrix full_intf7ff = RDenseMatrix (10,10,
   "12  6  6  8 -24 -24 -12 -40 -20 -20 \
     6 24  0  6 -12 -24 -12 -24 -12 -24 \
     6  0 24  6 -24 -12 -12 -24 -24 -12 \
     8  6  6 12 -20 -20 -12 -40 -24 -24 \
    -24 -12 -24 -20 96 64 32 96 48 32 \
    -24 -24 -12 -20 64 96 32 96 32 48 \
    -12 -12 -12 -12 32 32 32 32 32 32 \
    -40 -24 -24 -40 96 96 32 288 96 96 \
    -20 -12 -24 -24 48 32 32 96 96 64 \
    -20 -24 -12 -24 32 48 32 96 64 96") / 7560.0;
static const RDenseMatrix full_intf8ff = RDenseMatrix (10,10,
   "12  6  8  6 -24 -12 -24 -20 -40 -20 \
     6 24  6  0 -12 -12 -24 -12 -24 -24 \
     8  6 12  6 -20 -12 -20 -24 -40 -24 \
     6  0  6 24 -24 -12 -12 -24 -24 -12 \
    -24 -12 -20 -24 96 32 64 48 96 32 \
    -12 -12 -12 -12 32 32 32 32 32 32 \
    -24 -24 -20 -12 64 32 96 32 96 48 \
    -20 -12 -24 -24 48 32 32 96 96 64 \
    -40 -24 -40 -24 96 32 96 96 288 96 \
    -20 -24 -24 -12 32 32 48 64 96 96") / 7560.0;
static const RDenseMatrix full_intf9ff = RDenseMatrix (10,10,
   "12  8  6  6 -12 -24 -24 -20 -20 -40 \
     8 12  6  6 -12 -20 -20 -24 -24 -40 \
     6  6 24  0 -12 -12 -24 -12 -24 -24 \
     6  6  0 24 -12 -24 -12 -24 -12 -24 \
    -12 -12 -12 -12 32 32 32 32 32 32 \
    -24 -20 -12 -24 32 96 64 48 32 96 \
    -24 -20 -24 -12 32 64 96 32 48 96 \
    -20 -24 -12 -24 32 48 32 96 64 96 \
    -20 -24 -24 -12 32 32 48 64 96 96 \
    -40 -40 -24 -24 32 96 96 96 96 288") / 7560.0;

static const RDenseMatrix *full_intfff[10] = {
    &full_intf0ff,
    &full_intf1ff,
    &full_intf2ff,
    &full_intf3ff,
    &full_intf4ff,
    &full_intf5ff,
    &full_intf6ff,
    &full_intf7ff,
    &full_intf8ff,
    &full_intf9ff
};

// the same in symmetric matrix representation
static const RSymMatrix sym_intf0ff = RSymMatrix (10,
   "18 \
    -6 -6 \
    -6 -1 -6 \
    -6 -1 -1 -6 \
    24  0  6  6 -24 \
    24  6  0  6 -12 -24 \
    24  6  6  0 -12 -12 -24 \
    12  6  6  8 -24 -24 -12 -40 \
    12  6  8  6 -24 -12 -24 -20 -40 \
    12  8  6  6 -12 -24 -24 -20 -20 -40") * (1.0/7560.0);
static const RSymMatrix sym_intf1ff = RSymMatrix (10,
   "-6 \
    -6 18 \
    -1 -6 -6 \
    -1 -6 -1 -6 \
     0 24  6  6 -24 \
     6 12  6  8 -24 -40 \
     6 12  8  6 -24 -20 -40 \
     6 24  0  6 -12 -24 -12 -24 \
     6 24  6  0 -12 -12 -24 -12 -24 \
     8 12  6  6 -12 -20 -20 -24 -24 -40") * (1.0/7560.0);
static const RSymMatrix sym_intf2ff = RSymMatrix (10,
   "-6 \
    -1 -6 \
    -6 -6 18 \
    -1 -1 -6 -6 \
     6  6 12  8 -40 \
     0  6 24  6 -24 -24 \
     6  8 12  6 -20 -24 -40 \
     6  0 24  6 -24 -12 -12 -24 \
     8  6 12  6 -20 -12 -20 -24 -40 \
     6  6 24  0 -12 -12 -24 -12 -24 -24") * (1.0/7560.0);
static const RSymMatrix sym_intf3ff = RSymMatrix (10,
   "-6 \
    -1 -6 \
    -1 -1 -6 \
    -6 -6 -6 18 \
     6  6  8 12 -40 \
     6  8  6 12 -20 -40 \
     0  6  6 24 -24 -24 -24 \
     8  6  6 12 -20 -20 -12 -40 \
     6  0  6 24 -24 -12 -12 -24 -24 \
     6  6  0 24 -12 -24 -12 -24 -12 -24") * (1.0/7560.0);
static const RSymMatrix sym_intf4ff = RSymMatrix (10,
   "24 \
     0 24 \
     6  6 12 \
     6  6  8 12 \
    -24 -24 -40 -40 288 \
    -12 -24 -24 -20 96 96 \
    -12 -24 -20 -24 96 48 96 \
    -24 -12 -24 -20 96 64 32 96 \
    -24 -12 -20 -24 96 32 64 48 96 \
    -12 -12 -12 -12 32 32 32 32 32 32") * (1.0/7560.0);
static const RSymMatrix sym_intf5ff = RSymMatrix (10,
   "24 \
     6 12 \
     0  6 24 \
     6  8  6 12 \
    -12 -24 -24 -20 96 \
    -24 -40 -24 -40 96 288 \
    -12 -20 -24 -24 48 96 96 \
    -24 -24 -12 -20 64 96 32 96 \
    -12 -12 -12 -12 32 32 32 32 32 \
    -24 -20 -12 -24 32 96 64 48 32 96") * (1.0/7560.0);
static const RSymMatrix sym_intf6ff = RSymMatrix (10,
   "24 \
     6 12 \
     6  8 12 \
     0  6  6 24 \
    -12 -24 -20 -24 96 \
    -12 -20 -24 -24 48 96 \
    -24 -40 -40 -24 96 96 288 \
    -12 -12 -12 -12 32 32 32 32 \
    -24 -24 -20 -12 64 32 96 32 96 \
    -24 -20 -24 -12 32 64 96 32 48 96") * (1.0/7560.0);
static const RSymMatrix sym_intf7ff = RSymMatrix (10,
   "12 \
     6 24 \
     6  0 24 \
     8  6  6 12 \
    -24 -12 -24 -20 96 \
    -24 -24 -12 -20 64 96 \
    -12 -12 -12 -12 32 32 32 \
    -40 -24 -24 -40 96 96 32 288 \
    -20 -12 -24 -24 48 32 32 96 96 \
    -20 -24 -12 -24 32 48 32 96 64 96") * (1.0/7560.0);
static const RSymMatrix sym_intf8ff = RSymMatrix (10,
   "12 \
     6 24 \
     8  6 12 \
     6  0  6 24 \
    -24 -12 -20 -24 96 \
    -12 -12 -12 -12 32 32 \
    -24 -24 -20 -12 64 32 96 \
    -20 -12 -24 -24 48 32 32 96 \
    -40 -24 -40 -24 96 32 96 96 288 \
    -20 -24 -24 -12 32 32 48 64 96 96") * (1.0/7560.0);
static const RSymMatrix sym_intf9ff = RSymMatrix (10,
   "12 \
     8 12 \
     6  6 24 \
     6  6  0 24 \
    -12 -12 -12 -12 32 \
    -24 -20 -12 -24 32 96 \
    -24 -20 -24 -12 32 64 96 \
    -20 -24 -12 -24 32 48 32 96 \
    -20 -24 -24 -12 32 32 48 64 96 \
    -40 -40 -24 -24 32 96 96 96 96 288") * (1.0/7560.0);

Tetrahedron10::Tetrahedron10 (const Tetrahedron10 &el)
: Element_Unstructured_3D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Tetrahedron10::Copy ()
{
    return new Tetrahedron10(*this);
}

void Tetrahedron10::Initialise (const NodeList &nlist)
{
    extern double TriangleArea (const Point &p1,
				const Point &p2, const Point &p3);

    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1],
           z0 = nlist[Node[0]][2];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1],
           z1 = nlist[Node[1]][2];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1],
           z2 = nlist[Node[2]][2];
    double x3 = nlist[Node[3]][0], y3 = nlist[Node[3]][1],
           z3 = nlist[Node[3]][2];

    a0 = x1*(y2*z3-y3*z2) - x2*(y1*z3-y3*z1) + x3*(y1*z2-y2*z1);
    a1 = x0*(y3*z2-y2*z3) - x2*(y3*z0-y0*z3) + x3*(y2*z0-y0*z2);
    a2 = x0*(y1*z3-y3*z1) - x1*(y0*z3-y3*z0) + x3*(y0*z1-y1*z0);
    a3 = x0*(y2*z1-y1*z2) - x1*(y2*z0-y0*z2) + x2*(y1*z0-y0*z1);
    b0 = y1*(z3-z2) + y2*(z1-z3) + y3*(z2-z1);
    b1 = y0*(z2-z3) + y2*(z3-z0) + y3*(z0-z2);
    b2 = y0*(z3-z1) + y1*(z0-z3) + y3*(z1-z0);
    b3 = y0*(z1-z2) + y1*(z2-z0) + y2*(z0-z1);
    c0 = x1*(z2-z3) + x2*(z3-z1) + x3*(z1-z2);
    c1 = x0*(z3-z2) + x2*(z0-z3) + x3*(z2-z0);
    c2 = x0*(z1-z3) + x1*(z3-z0) + x3*(z0-z1);
    c3 = x0*(z2-z1) + x1*(z0-z2) + x2*(z1-z0);
    d0 = x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1);
    d1 = x0*(y2-y3) + x2*(y3-y0) + x3*(y0-y2);
    d2 = x0*(y3-y1) + x1*(y0-y3) + x3*(y1-y0);
    d3 = x0*(y1-y2) + x1*(y2-y0) + x2*(y0-y1);

    // calculate surface triangle areas
    for (int sd = 0; sd < 4; sd++)
        side_size[sd] = TriangleArea (nlist[Node[SideNode(sd,0)]],
			 	      nlist[Node[SideNode(sd,1)]],
				      nlist[Node[SideNode(sd,2)]]);

    Element_Unstructured_3D::Initialise (nlist);

#ifdef TET10_STORE_INTFF
    intff.New(10);
    intff = sym_intff * size;
#endif

    if (!subsampling_initialised) {
        int i, j, k, idx;
        Point loc(3);
	for (k = idx = 0; k < nsample_lin; k++) {
	    loc[2] = (double)k/(double)(nsample_lin-1);
	    for (j = 0; j < nsample_lin-k; j++) {
	        loc[1] = (double)j/(double)(nsample_lin-1);
		for (i = 0; i < nsample_lin-k-j; i++) {
		    loc[0] = (double)i/(double)(nsample_lin-1);
		    absc_sample[idx].New (3);
		    absc_sample[idx] = loc;
		    idx++;
		}
	    }
	}
	subsampling_initialised = true;
    }
}

int Tetrahedron10::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 4, "Side index out of range");
    dASSERT(node >= 0 && node < 6, "Node index out of range");
    static int SN[4][6] = {{0,1,2,4,7,5},{0,3,1,6,8,4},
			   {0,2,3,5,9,6},{1,3,2,8,9,7}};
    return SN[side][node];
}

double Tetrahedron10::SideSize (int sd, const NodeList &nlist) const
{
    return side_size[sd];
}

Point Tetrahedron10::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");

    Point loc(3);
    double scale = 1.0/(a0+a1+a2+a3);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1] + d1*glob[2]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1] + d2*glob[2]) * scale;
    loc[2] = (a3 + b3*glob[0] + c3*glob[1] + d3*glob[2]) * scale;
    return loc;
}

Point Tetrahedron10::NodeLocal (int node) const
{
    Point nloc(3);
    switch (node) {
    case 0: nloc[0] = nloc[1] = nloc[2] = 0.0; break;
    case 1: nloc[0] = 1.0; nloc[1] = nloc[2] = 0.0; break;
    case 2: nloc[0] = nloc[2] = 0.0; nloc[1] = 1.0; break;
    case 3: nloc[0] = nloc[1] = 0.0; nloc[2] = 1.0; break;
    case 4: nloc[0] = 0.5; nloc[1] = nloc[2] = 0.0; break;
    case 5: nloc[0] = nloc[2] = 0.0; nloc[1] = 0.5; break;
    case 6: nloc[0] = nloc[1] = 0.0; nloc[2] = 0.5; break;
    case 7: nloc[0] = nloc[1] = 0.5; nloc[2] = 0.0; break;
    case 8: nloc[0] = nloc[2] = 0.5; nloc[1] = 0.0; break;
    case 9: nloc[0] = 0.0; nloc[1] = nloc[2] = 0.5; break;
    default: xERROR("Node index out of range");
    }
    return nloc;
}

Point Tetrahedron10::SurfToLocal (int side, const Point &p) const
{
    dASSERT(p.Dim() == 2, "Arg 2 wrong vector dimension");

    Point loc(3);
    switch (side) {
    case 0: loc[0] = p[0]; loc[1] = p[1]; loc[2] = 0.0;  break;
    case 1: loc[0] = p[0]; loc[1] = 0.0;  loc[2] = p[1]; break;
    case 2: loc[0] = 0.0;  loc[1] = p[0]; loc[2] = p[1]; break;
    case 3: loc[0] = p[0]; loc[1] = p[1];
            loc[2] = 1.0-p[0]-p[1];                      break;
    default: xERROR("Arg 1 index out of range");
    }
    return loc;
}

RVector Tetrahedron10::DirectionCosine (int side, RDenseMatrix &jacin)
{
    RVector cosin(3);
    switch (side) {
    case 0: cosin[0] = -b3, cosin[1] = -c3, cosin[2] = -d3; break;
    case 1: cosin[0] = -b2, cosin[1] = -c2, cosin[2] = -d2; break;
    case 2: cosin[0] = -b1, cosin[1] = -c1, cosin[2] = -d1; break;
    case 3: cosin[0] = -b0, cosin[1] = -c0, cosin[2] = -d0; break;
    default: xERROR("Side index out of range");
    }
    return cosin/length(cosin);
}

const RVector &Tetrahedron10::LNormal (int side) const
{
    static const RVector lnm0 = RVector (3, "0 0 -1");
    static const RVector lnm1 = RVector (3, "0 -1 0");
    static const RVector lnm2 = RVector (3, "-1 0 0");
    static const RVector lnm3 = RVector (3,
         "0.577350269 0.577350269 0.577350269");
    static const RVector *lnm[4] = {
      &lnm0, &lnm1, &lnm2, &lnm3
    };
    dASSERT(side >= 0 && side < 4, "Argument 1 index out of range");
    return *lnm[side];
}

bool Tetrahedron10::LContains (const Point &loc, bool pad) const
{
    dASSERT(loc.Dim() == 3, "Local point must be 3D");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[1]+EPS >= 0.0 &&
	    loc[2]+EPS >= 0.0 && loc[0]+loc[1]+loc[2]-EPS <= 1.0);
    } else {
	return (loc[0] >= 0.0 && loc[1] >= 0.0 &&
	    loc[2] >= 0.0 && loc[0]+loc[1]+loc[2] <= 1.0);
    }
}

RVector Tetrahedron10::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Invalid point dimension");
    RVector fun(10);
    double L0 = 1.0-loc[0]-loc[1]-loc[2];
    double L1 = loc[0];
    double L2 = loc[1];
    double L3 = loc[2];
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = L3 * (2.0*L3 - 1.0);
    fun[4] = 4.0*L0*L1;
    fun[5] = 4.0*L0*L2;
    fun[6] = 4.0*L0*L3;
    fun[7] = 4.0*L1*L2;
    fun[8] = 4.0*L1*L3;
    fun[9] = 4.0*L2*L3;
    return fun;
}

RDenseMatrix Tetrahedron10::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Invalid point dimension");
    RDenseMatrix der(3,10);
    double lx = loc[0], ly = loc[1], lz = loc[2];

    // elements not set are zero
    der(0,0) = der(1,0) = der(2,0) = 4.0*(lx+ly+lz)-3.0;
    der(0,1) = 4.0*lx - 1.0;
    der(1,2) = 4.0*ly - 1.0;
    der(2,3) = 4.0*lz - 1.0;
    der(0,4) = 4.0*(1.0-2.0*lx-ly-lz);
    der(1,4) = der(2,4) = -4.0*lx;
    der(0,5) = der(2,5) = -4.0*ly;
    der(1,5) = 4.0*(1.0-lx-2.0*ly-lz);
    der(0,6) = der(1,6) = -4.0*lz;
    der(2,6) = 4.0*(1.0-lx-ly-2.0*lz);
    der(0,7) = der(2,9) = 4.0*ly;
    der(1,7) = der(2,8) = 4.0*lx;
    der(0,8) = der(1,9) = 4.0*lz;
    return der;
}

RVector Tetrahedron10::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");
    RVector fun(10);
    double scale = 1.0/(6.0*size);
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1] + d0*glob[2]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1] + d1*glob[2]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1] + d2*glob[2]);
    double L3 = scale * (a3 + b3*glob[0] + c3*glob[1] + d3*glob[2]);
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = L3 * (2.0*L3 - 1.0);
    fun[4] = 4.0*L0*L1;
    fun[5] = 4.0*L0*L2;
    fun[6] = 4.0*L0*L3;
    fun[7] = 4.0*L1*L2;
    fun[8] = 4.0*L1*L3;
    fun[9] = 4.0*L2*L3;
    return fun;
}

RDenseMatrix Tetrahedron10::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");
    RDenseMatrix der(3,10);
    double scale = 1.0/(6.0*size);
    double scale4 = 4.0*scale;
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1] + d0*glob[2]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1] + d1*glob[2]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1] + d2*glob[2]);
    double L3 = scale * (a3 + b3*glob[0] + c3*glob[1] + d3*glob[2]);

    der(0,0) = b0*scale * (4.0*L0-1.0);
    der(1,0) = c0*scale * (4.0*L0-1.0);
    der(2,0) = d0*scale * (4.0*L0-1.0);
    der(0,1) = b1*scale * (4.0*L1-1.0);
    der(1,1) = c1*scale * (4.0*L1-1.0);
    der(2,1) = d1*scale * (4.0*L1-1.0);
    der(0,2) = b2*scale * (4.0*L2-1.0);
    der(1,2) = c2*scale * (4.0*L2-1.0);
    der(2,2) = d2*scale * (4.0*L2-1.0);
    der(0,3) = b3*scale * (4.0*L3-1.0);
    der(1,3) = c3*scale * (4.0*L3-1.0);
    der(2,3) = d3*scale * (4.0*L3-1.0);
    der(0,4) = scale4 * (b0*L1 + b1*L0);
    der(1,4) = scale4 * (c0*L1 + c1*L0);
    der(2,4) = scale4 * (d0*L1 + d1*L0);
    der(0,5) = scale4 * (b0*L2 + b2*L0);
    der(1,5) = scale4 * (c0*L2 + c2*L0);
    der(2,5) = scale4 * (d0*L2 + d2*L0);
    der(0,6) = scale4 * (b0*L3 + b3*L0);
    der(1,6) = scale4 * (c0*L3 + c3*L0);
    der(2,6) = scale4 * (d0*L3 + d3*L0);
    der(0,7) = scale4 * (b1*L2 + b2*L1);
    der(1,7) = scale4 * (c1*L2 + c2*L1);
    der(2,7) = scale4 * (d1*L2 + d2*L1);
    der(0,8) = scale4 * (b1*L3 + b3*L1);
    der(1,8) = scale4 * (c1*L3 + c3*L1);
    der(2,8) = scale4 * (d1*L3 + d3*L1);
    der(0,9) = scale4 * (b2*L3 + b3*L2);
    der(1,9) = scale4 * (c2*L3 + c3*L2);
    der(2,9) = scale4 * (d2*L3 + d3*L2);
    return der;
}

double Tetrahedron10::IntF (int i) const
{
    return 0.25*size;
}

RSymMatrix Tetrahedron10::IntFF () const {
#ifdef TET10_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Tetrahedron10::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
#ifdef TET10_STORE_INTFF
    return intff(i,j);
#else
    return full_intff(i,j) * size;
#endif
}

double Tetrahedron10::IntFFF (int i, int j, int k) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10 && k >= 0 && k < 10);
    return full_intfff[i]->Get(j,k) * size;
}

RSymMatrix Tetrahedron10::IntPFF (const RVector &P) const
{
    return (sym_intf0ff * P[Node[0]] +
	    sym_intf1ff * P[Node[1]] +
	    sym_intf2ff * P[Node[2]] +
	    sym_intf3ff * P[Node[3]] +
	    sym_intf4ff * P[Node[4]] +
	    sym_intf5ff * P[Node[5]] +
	    sym_intf6ff * P[Node[6]] +
	    sym_intf7ff * P[Node[7]] +
	    sym_intf8ff * P[Node[8]] +
	    sym_intf9ff * P[Node[9]]) * size;
}

double Tetrahedron10::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
    return (full_intf0ff(i,j) * P[Node[0]] +
	    full_intf1ff(i,j) * P[Node[1]] +
	    full_intf2ff(i,j) * P[Node[2]] +
	    full_intf3ff(i,j) * P[Node[3]] +
	    full_intf4ff(i,j) * P[Node[4]] +
	    full_intf5ff(i,j) * P[Node[5]] +
	    full_intf6ff(i,j) * P[Node[6]] +
	    full_intf7ff(i,j) * P[Node[7]] +
	    full_intf8ff(i,j) * P[Node[8]] +
	    full_intf9ff(i,j) * P[Node[9]]) * size;
}

static const RVector intfd0d0 = RVector (10,
   "27 -13 -13 -13 52 52 52 36 36 36")/15120.0;
static const RVector intfd0d1 = RVector (10,
   "-9 -9 11 11 -12 -20 -20 -20 -20 4")/15120.0;
static const RVector intfd0d2 = RVector (10,
   "-9 11 -9 11 -20 -12 -20 -20 4 -20")/15120.0;
static const RVector intfd0d3 = RVector (10,
   "-9 11 11 -9 -20 -20 -12 4 -20 -20")/15120.0;
static const RDenseMatrix intfd0d4 = RDenseMatrix (10, 6,
   "12 48 12 48 12 48 \
    -16 -20 -16 -20 -16 -20 \
      4 -20   4 -20   4 -20 \
      4 -20   4 -20   4 -20 \
     16  80  16  80  16  80 \
      8  80   8  80   8  80 \
      8  80   8  80   8  80 \
    -48   8 -48   8 -48   8 \
    -48   8 -48   8 -48   8 \
    -24   8 -24   8 -24   8")/15120.0;
static const RDenseMatrix intfd0d5 = RDenseMatrix (10, 6,
   " 12  48  12  48  12  48 \
      4 -20   4 -20   4 -20 \
    -16 -20 -16 -20 -16 -20 \
      4 -20   4 -20   4 -20 \
      8  80   8  80   8  80 \
     16  80  16  80  16  80 \
      8  80   8  80   8  80 \
    -48   8 -48   8 -48   8 \
    -24   8 -24   8 -24   8 \
    -48   8 -48   8 -48   8")/15120.0;
static const RDenseMatrix intfd0d6 = RDenseMatrix (10, 6,
   " 12  48  12  48  12  48 \
      4 -20   4 -20   4 -20 \
      4 -20   4 -20   4 -20 \
    -16 -20 -16 -20 -16 -20 \
      8  80   8  80   8  80 \
      8  80   8  80   8  80 \
     16  80  16  80  16  80 \
    -24   8 -24   8 -24   8 \
    -48   8 -48   8 -48   8 \
    -48   8 -48   8 -48   8")/15120.0;
static const RDenseMatrix intfd0d7 = RDenseMatrix (10, 6,
   " 12  12  12  12  12  12 \
      4 -16   4 -16   4 -16 \
    -16   4 -16   4 -16   4 \
      4   4   4   4   4   4 \
      8  16   8  16   8  16 \
     16   8  16   8  16   8 \
      8   8   8   8   8   8 \
    -48 -48 -48 -48 -48 -48 \
    -24 -48 -24 -48 -24 -48 \
    -48 -24 -48 -24 -48 -24")/15120.0;
static const RDenseMatrix intfd0d8 = RDenseMatrix (10, 6,
   " 12  12  12  12  12  12 \
      4 -16   4 -16   4 -16 \
      4   4   4   4   4   4 \
    -16   4 -16   4 -16   4 \
      8  16   8  16   8  16 \
      8   8   8   8   8   8 \
     16   8  16   8  16   8 \
    -24 -48 -24 -48 -24 -48 \
    -48 -48 -48 -48 -48 -48 \
    -48 -24 -48 -24 -48 -24")/15120.0;
static const RDenseMatrix intfd0d9 = RDenseMatrix (10, 6,
   " 12  12  12  12  12  12 \
      4   4   4   4   4   4 \
      4 -16   4 -16   4 -16 \
    -16   4 -16   4 -16   4 \
      8   8   8   8   8   8 \
      8  16   8  16   8  16 \
     16   8  16   8  16   8 \
    -24 -48 -24 -48 -24 -48 \
    -48 -24 -48 -24 -48 -24 \
    -48 -48 -48 -48 -48 -48")/15120.0;
static const RVector intfd1d1 = RVector (10,
   "-13 27 -13 -13 52 36 36 52 52 36")/15120.0;
static const RVector intfd1d2 = RVector (10,
   "11 -9 -9 11 -20 -20 4 -12 -20 -20")/15120.0;
static const RVector intfd1d3 = RVector (10,
   "11 -9 11 -9 -20 4 -20 -20 -12 -20")/15120.0;
static const RDenseMatrix intfd1d4 = RDenseMatrix (10, 6,
   "-20 -16 -20 -16 -20 -16 \
     48  12  48  12  48  12 \
    -20   4 -20   4 -20   4 \
    -20   4 -20   4 -20   4 \
     80  16  80  16  80  16 \
      8 -48   8 -48   8 -48 \
      8 -48   8 -48   8 -48 \
     80   8  80   8  80   8 \
     80   8  80   8  80   8 \
      8 -24   8 -24   8 -24")/15120.0;
static const RDenseMatrix intfd1d5 = RDenseMatrix (10, 6,
   "  4 -16   4 -16   4 -16 \
     12  12  12  12  12  12 \
    -16   4 -16   4 -16   4 \
      4   4   4   4   4   4 \
      8  16   8  16   8  16 \
    -48 -48 -48 -48 -48 -48 \
    -24 -48 -24 -48 -24 -48 \
     16   8  16   8  16   8 \
      8   8   8   8   8   8 \
    -48 -24 -48 -24 -48 -24")/15120.0;
static const RDenseMatrix intfd1d6 = RDenseMatrix (10, 6,
   "  4 -16   4 -16   4 -16 \
     12  12  12  12  12  12 \
      4   4   4   4   4   4 \
    -16   4 -16   4 -16   4 \
      8  16   8  16   8  16 \
    -24 -48 -24 -48 -24 -48 \
    -48 -48 -48 -48 -48 -48 \
      8   8   8   8   8   8 \
     16   8  16   8  16   8 \
    -48 -24 -48 -24 -48 -24")/15120.0;
static const RDenseMatrix intfd1d7 = RDenseMatrix (10, 6,
   "  4 -20   4 -20   4 -20 \
     12  48  12  48  12  48 \
    -16 -20 -16 -20 -16 -20 \
      4 -20   4 -20   4 -20 \
      8  80   8  80   8  80 \
    -48   8 -48   8 -48   8 \
    -24   8 -24   8 -24   8 \
     16  80  16  80  16  80 \
      8  80   8  80   8  80 \
    -48   8 -48   8 -48   8")/15120.0;
static const RDenseMatrix intfd1d8 = RDenseMatrix (10, 6,
   "  4 -20   4 -20   4 -20 \
     12  48  12  48  12  48 \
      4 -20   4 -20   4 -20 \
    -16 -20 -16 -20 -16 -20 \
      8  80   8  80   8  80 \
    -24   8 -24   8 -24   8 \
    -48   8 -48   8 -48   8 \
      8  80   8  80   8  80 \
     16  80  16  80  16  80 \
    -48   8 -48   8 -48   8")/15120.0;
static const RDenseMatrix intfd1d9 = RDenseMatrix (10, 6,
   "  4   4   4   4   4   4 \
     12  12  12  12  12  12 \
      4 -16   4 -16   4 -16 \
    -16   4 -16   4 -16   4 \
      8   8   8   8   8   8 \
    -24 -48 -24 -48 -24 -48 \
    -48 -24 -48 -24 -48 -24 \
      8  16   8  16   8  16 \
     16   8  16   8  16   8 \
    -48 -48 -48 -48 -48 -48")/15120.0;
static const RVector intfd2d2 = RVector (10,
   "-13 -13 27 -13 36 52 36 52 36 52")/15120.0;
static const RVector intfd2d3 = RVector (10,
   "11 11 -9 -9 4 -20 -20 -20 -20 -12")/15120.0;
static const RDenseMatrix intfd2d4 = RDenseMatrix (10, 6,
   "  4 -16   4 -16   4 -16 \
    -16   4 -16   4 -16   4 \
     12  12  12  12  12  12 \
      4   4   4   4   4   4 \
    -48 -48 -48 -48 -48 -48 \
      8  16   8  16   8  16 \
    -24 -48 -24 -48 -24 -48 \
     16   8  16   8  16   8 \
    -48 -24 -48 -24 -48 -24 \
      8   8   8   8   8   8")/15120.0;
static const RDenseMatrix intfd2d5 = RDenseMatrix (10, 6,
   "-20 -16 -20 -16 -20 -16 \
    -20   4 -20   4 -20   4 \
     48  12  48  12  48  12 \
    -20   4 -20   4 -20   4 \
      8 -48   8 -48   8 -48 \
     80  16  80  16  80  16 \
      8 -48   8 -48   8 -48 \
     80   8  80   8  80   8 \
      8 -24   8 -24   8 -24 \
     80   8  80   8  80   8")/15120.0;
static const RDenseMatrix intfd2d6 = RDenseMatrix (10, 6,
   "  4 -16   4 -16   4 -16 \
      4   4   4   4   4   4 \
     12  12  12  12  12  12 \
    -16   4 -16   4 -16   4 \
    -24 -48 -24 -48 -24 -48 \
      8  16   8  16   8  16 \
    -48 -48 -48 -48 -48 -48 \
      8   8   8   8   8   8 \
    -48 -24 -48 -24 -48 -24 \
     16   8  16   8  16   8")/15120.0;
static const RDenseMatrix intfd2d7 = RDenseMatrix (10, 6,
   "-20   4 -20   4 -20   4 \
    -20 -16 -20 -16 -20 -16 \
     48  12  48  12  48  12 \
    -20   4 -20   4 -20   4 \
      8 -48   8 -48   8 -48 \
     80   8  80   8  80   8 \
      8 -24   8 -24   8 -24 \
     80  16  80  16  80  16 \
      8 -48   8 -48   8 -48 \
     80   8  80   8  80   8")/15120.0;
static const RDenseMatrix intfd2d8 = RDenseMatrix (10, 6,
   "  4   4   4   4   4   4 \
      4 -16   4 -16   4 -16 \
     12  12  12  12  12  12 \
    -16   4 -16   4 -16   4 \
    -24 -48 -24 -48 -24 -48 \
      8   8   8   8   8   8 \
    -48 -24 -48 -24 -48 -24 \
      8  16   8  16   8  16 \
    -48 -48 -48 -48 -48 -48 \
     16   8  16   8  16   8")/15120.0;
static const RDenseMatrix intfd2d9 = RDenseMatrix (10, 6,
   "  4 -20   4 -20   4 -20 \
      4 -20   4 -20   4 -20 \
     12  48  12  48  12  48 \
    -16 -20 -16 -20 -16 -20 \
    -24   8 -24   8 -24   8 \
      8  80   8  80   8  80 \
    -48   8 -48   8 -48   8 \
      8  80   8  80   8  80 \
    -48   8 -48   8 -48   8 \
     16  80  16  80  16  80")/15120.0;
static const RVector intfd3d3 = RVector (10,
   "-13 -13 -13 27 36 36 52 36 52 52")/15120.0;
static const RDenseMatrix intfd3d4 = RDenseMatrix (10, 6,
   "  4 -16   4 -16   4 -16 \
    -16   4 -16   4 -16   4 \
      4   4   4   4   4   4 \
     12  12  12  12  12  12 \
    -48 -48 -48 -48 -48 -48 \
    -24 -48 -24 -48 -24 -48 \
      8  16   8  16   8  16 \
    -48 -24 -48 -24 -48 -24 \
     16   8  16   8  16   8 \
      8   8   8   8   8   8")/15120.0;
static const RDenseMatrix intfd3d5 = RDenseMatrix (10, 6,
   "  4 -16   4 -16   4 -16 \
      4   4   4   4   4   4 \
    -16   4 -16   4 -16   4 \
     12  12  12  12  12  12 \
    -24 -48 -24 -48 -24 -48 \
    -48 -48 -48 -48 -48 -48 \
      8  16   8  16   8  16 \
    -48 -24 -48 -24 -48 -24 \
      8   8   8   8   8   8 \
     16   8  16   8  16   8")/15120.0;
static const RDenseMatrix intfd3d6 = RDenseMatrix (10, 6,
   "-20 -16 -20 -16 -20 -16 \
    -20   4 -20   4 -20   4 \
    -20   4 -20   4 -20   4 \
     48  12  48  12  48  12 \
      8 -48   8 -48   8 -48 \
      8 -48   8 -48   8 -48 \
     80  16  80  16  80  16 \
      8 -24   8 -24   8 -24 \
     80   8  80   8  80   8 \
     80   8  80   8  80   8")/15120.0;
static const RDenseMatrix intfd3d7 = RDenseMatrix (10, 6,
   "  4   4   4   4   4   4 \
      4 -16   4 -16   4 -16 \
    -16   4 -16   4 -16   4 \
     12  12  12  12  12  12 \
    -24 -48 -24 -48 -24 -48 \
    -48 -24 -48 -24 -48 -24 \
      8   8   8   8   8   8 \
    -48 -48 -48 -48 -48 -48 \
      8  16   8  16   8  16 \
     16   8  16   8  16   8")/15120.0;
static const RDenseMatrix intfd3d8 = RDenseMatrix (10, 6,
   "-20   4 -20   4 -20   4 \
    -20 -16 -20 -16 -20 -16 \
    -20   4 -20   4 -20   4 \
     48  12  48  12  48  12 \
      8 -48   8 -48   8 -48 \
      8 -24   8 -24   8 -24 \
     80   8  80   8  80   8 \
      8 -48   8 -48   8 -48 \
     80  16  80  16  80  16 \
     80   8  80   8  80   8")/15120.0;
static const RDenseMatrix intfd3d9 = RDenseMatrix (10, 6,
   "-20   4 -20   4 -20   4 \
    -20   4 -20   4 -20   4 \
    -20 -16 -20 -16 -20 -16 \
     48  12  48  12  48  12 \
      8 -24   8 -24   8 -24 \
      8 -48   8 -48   8 -48 \
     80   8  80   8  80   8 \
      8 -48   8 -48   8 -48 \
     80   8  80   8  80   8 \
     80  16  80  16  80  16")/15120.0;
static const RDenseMatrix intfd4d4 = RDenseMatrix (10, 9,
   "-48 -32  48 -48 -32  48 -48 -32  48 \
     48 -32 -48  48 -32 -48  48 -32 -48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
    192 256 192 192 256 192 192 256 192 \
     64 128 192  64 128 192  64 128 192 \
     64 128 192  64 128 192  64 128 192 \
    192 128  64 192 128  64 192 128  64 \
    192 128  64 192 128  64 192 128  64 \
     64  64  64  64  64  64  64  64  64")/15120.0;
static const RDenseMatrix intfd4d5 = RDenseMatrix (10, 12,
   "-24 -16 -16  48 -24 -16 -16  48 -24 -16 -16  48 \
    -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 \
    -16 -16 -24 -48 -16 -16 -24 -48 -16 -16 -24 -48 \
    -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 \
     64  64 128 192  64  64 128 192  64  64 128 192 \
     64 128  64 192  64 128  64 192  64 128  64 192 \
     32  64  64 192  32  64  64 192  32  64  64 192 \
    128  64  64  64 128  64  64  64 128  64  64  64 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
     64  64  32  64  64  64  32  64  64  64  32  64")/15120.0;
static const RDenseMatrix intfd4d6 = RDenseMatrix (10, 12,
   "-24 -16 -16  48 -24 -16 -16  48 -24 -16 -16  48 \
    -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 \
    -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 \
    -16 -16 -24 -48 -16 -16 -24 -48 -16 -16 -24 -48 \
     64  64 128 192  64  64 128 192  64  64 128 192 \
     32  64  64 192  32  64  64 192  32  64  64 192 \
     64 128  64 192  64 128  64 192  64 128  64 192 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
    128  64  64  64 128  64  64  64 128  64  64  64 \
     64  64  32  64  64  64  32  64  64  64  32  64")/15120.0;
static const RDenseMatrix intfd4d7 = RDenseMatrix (10, 12,
   "-24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 \
    -16 -24  48 -16 -16 -24  48 -16 -16 -24  48 -16 \
    -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 -48 -24 \
    -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 \
     64  64 192 128  64  64 192 128  64  64 192 128 \
     64 128  64  64  64 128  64  64  64 128  64  64 \
     32  64  64  64  32  64  64  64  32  64  64  64 \
    128  64 192  64 128  64 192  64 128  64 192  64 \
     64  32 192  64  64  32 192  64  64  32 192  64 \
     64  64  64  32  64  64  64  32  64  64  64  32")/15120.0;
static const RDenseMatrix intfd4d8 = RDenseMatrix (10, 12,
   "-24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 \
    -16 -24  48 -16 -16 -24  48 -16 -16 -24  48 -16 \
    -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 \
    -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 -48 -24 \
     64  64 192 128  64  64 192 128  64  64 192 128 \
     32  64  64  64  32  64  64  64  32  64  64  64 \
     64 128  64  64  64 128  64  64  64 128  64  64 \
     64  32 192  64  64  32 192  64  64  32 192  64 \
    128  64 192  64 128  64 192  64 128  64 192  64 \
     64  64  64  32  64  64  64  32  64  64  64  32")/15120.0;
static const RDenseMatrix intfd4d9 = RDenseMatrix (10, 12,
   "-24 -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 \
    -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 \
    -24 -24 -16 -16 -24 -24 -16 -16 -24 -24 -16 -16 \
    -16 -16 -24 -24 -16 -16 -24 -24 -16 -16 -24 -24 \
     64  64  64  64  64  64  64  64  64  64  64  64 \
     32  64  64 128  32  64  64 128  32  64  64 128 \
     64 128  32  64  64 128  32  64  64 128  32  64 \
     64  32 128  64  64  32 128  64  64  32 128  64 \
    128  64  64  32 128  64  64  32 128  64  64  32 \
     64  64  64  64  64  64  64  64  64  64  64  64")/15120.0;
static const RDenseMatrix intfd5d5 = RDenseMatrix (10, 9,
   "-48 -32  48 -48 -32  48 -48 -32  48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
     48 -32 -48  48 -32 -48  48 -32 -48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
     64 128 192  64 128 192  64 128 192 \
    192 256 192 192 256 192 192 256 192 \
     64 128 192  64 128 192  64 128 192 \
    192 128  64 192 128  64 192 128  64 \
     64  64  64  64  64  64  64  64  64 \
    192 128  64 192 128  64 192 128  64")/15120.0;
static const RDenseMatrix intfd5d6 = RDenseMatrix (10, 12,
   "-24 -16 -16  48 -24 -16 -16  48 -24 -16 -16  48 \
    -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 \
    -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 \
    -16 -16 -24 -48 -16 -16 -24 -48 -16 -16 -24 -48 \
     32  64  64 192  32  64  64 192  32  64  64 192 \
     64  64 128 192  64  64 128 192  64  64 128 192 \
     64 128  64 192  64 128  64 192  64 128  64 192 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
     64  64  32  64  64  64  32  64  64  64  32  64 \
    128  64  64  64 128  64  64  64 128  64  64  64")/15120.0;
static const RDenseMatrix intfd5d7 = RDenseMatrix (10, 12,
   "-48 -24 -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 \
    -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 \
     48 -16 -16 -24  48 -16 -16 -24  48 -16 -16 -24 \
    -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 \
     64  64  64 128  64  64  64 128  64  64  64 128 \
    192  64 128  64 192  64 128  64 192  64 128  64 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
    192 128  64  64 192 128  64  64 192 128  64  64 \
     64  64  32  64  64  64  32  64  64  64  32  64 \
    192  64  64  32 192  64  64  32 192  64  64  32")/15120.0;
static const RDenseMatrix intfd5d8 = RDenseMatrix (10, 12,
   "-24 -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 \
    -24 -24 -16 -16 -24 -24 -16 -16 -24 -24 -16 -16 \
    -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 \
    -16 -16 -24 -24 -16 -16 -24 -24 -16 -16 -24 -24 \
     32  64  64 128  32  64  64 128  32  64  64 128 \
     64  64  64  64  64  64  64  64  64  64  64  64 \
     64 128  32  64  64 128  32  64  64 128  32  64 \
     64  32 128  64  64  32 128  64  64  32 128  64 \
     64  64  64  64  64  64  64  64  64  64  64  64 \
    128  64  64  32 128  64  64  32 128  64  64  32")/15120.0;
static const RDenseMatrix intfd5d9 = RDenseMatrix (10, 12,
   "-24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 \
    -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 \
    -16 -24  48 -16 -16 -24  48 -16 -16 -24  48 -16 \
    -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 -48 -24 \
     32  64  64  64  32  64  64  64  32  64  64  64 \
     64  64 192 128  64  64 192 128  64  64 192 128 \
     64 128  64  64  64 128  64  64  64 128  64  64 \
     64  32 192  64  64  32 192  64  64  32 192  64 \
     64  64  64  32  64  64  64  32  64  64  64  32 \
    128  64 192  64 128  64 192  64 128  64 192  64")/15120.0;
static const RDenseMatrix intfd6d6 = RDenseMatrix (10, 9,
   "-48 -32  48 -48 -32  48 -48 -32  48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
     48 -32 -48  48 -32 -48  48 -32 -48 \
     64 128 192  64 128 192  64 128 192 \
     64 128 192  64 128 192  64 128 192 \
    192 256 192 192 256 192 192 256 192 \
     64  64  64  64  64  64  64  64  64 \
    192 128  64 192 128  64 192 128  64 \
    192 128  64 192 128  64 192 128  64")/15120.0;
static const RDenseMatrix intfd6d7 = RDenseMatrix (10, 12,
   "-24 -24 -16 -16 -24 -24 -16 -16 -24 -24 -16 -16 \
    -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 \
    -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 -16 -24 \
    -16 -16 -24 -24 -16 -16 -24 -24 -16 -16 -24 -24 \
     32  64  64 128  32  64  64 128  32  64  64 128 \
     64  32 128  64  64  32 128  64  64  32 128  64 \
     64  64  64  64  64  64  64  64  64  64  64  64 \
     64  64  64  64  64  64  64  64  64  64  64  64 \
     64 128  32  64  64 128  32  64  64 128  32  64 \
    128  64  64  32 128  64  64  32 128  64  64  32")/15120.0;
static const RDenseMatrix intfd6d8 = RDenseMatrix (10, 12,
   "-48 -24 -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 \
    -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 \
    -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 \
     48 -16 -16 -24  48 -16 -16 -24  48 -16 -16 -24 \
     64  64  64 128  64  64  64 128  64  64  64 128 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
    192  64 128  64 192  64 128  64 192  64 128  64 \
     64  64  32  64  64  64  32  64  64  64  32  64 \
    192 128  64  64 192 128  64  64 192 128  64  64 \
    192  64  64  32 192  64  64  32 192  64  64  32")/15120.0;
static const RDenseMatrix intfd6d9 = RDenseMatrix (10, 12,
   "-48 -24 -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 \
    -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 \
    -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 \
     48 -16 -16 -24  48 -16 -16 -24  48 -16 -16 -24 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
     64  64  64 128  64  64  64 128  64  64  64 128 \
    192  64 128  64 192  64 128  64 192  64 128  64 \
     64  64  32  64  64  64  32  64  64  64  32  64 \
    192  64  64  32 192  64  64  32 192  64  64  32 \
    192 128  64  64 192 128  64  64 192 128  64  64")/15120.0;
static const RDenseMatrix intfd7d7 = RDenseMatrix (10, 9,
   "-48 -48 -48 -48 -48 -48 -48 -48 -48 \
    -48 -32  48 -48 -32  48 -48 -32  48 \
     48 -32 -48  48 -32 -48  48 -32 -48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
     64 128 192  64 128 192  64 128 192 \
    192 128  64 192 128  64 192 128  64 \
     64  64  64  64  64  64  64  64  64 \
    192 256 192 192 256 192 192 256 192 \
     64 128 192  64 128 192  64 128 192 \
    192 128  64 192 128  64 192 128  64")/15120.0;
static const RDenseMatrix intfd7d8 = RDenseMatrix (10, 12,
   "-24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 \
    -24 -16 -16  48 -24 -16 -16  48 -24 -16 -16  48 \
    -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 \
    -16 -16 -24 -48 -16 -16 -24 -48 -16 -16 -24 -48 \
     32  64  64 192  32  64  64 192  32  64  64 192 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
     64  64  32  64  64  64  32  64  64  64  32  64 \
     64  64 128 192  64  64 128 192  64  64 128 192 \
     64 128  64 192  64 128  64 192  64 128  64 192 \
    128  64  64  64 128  64  64  64 128  64  64  64")/15120.0;
static const RDenseMatrix intfd7d9 = RDenseMatrix (10, 12,
   "-24 -24 -48 -24 -24 -24 -48 -24 -24 -24 -48 -24 \
    -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 \
    -16 -24  48 -16 -16 -24  48 -16 -16 -24  48 -16 \
    -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 -48 -24 \
     32  64  64  64  32  64  64  64  32  64  64  64 \
     64  32 192  64  64  32 192  64  64  32 192  64 \
     64  64  64  32  64  64  64  32  64  64  64  32 \
     64  64 192 128  64  64 192 128  64  64 192 128 \
     64 128  64  64  64 128  64  64  64 128  64  64 \
    128  64 192  64 128  64 192  64 128  64 192  64")/15120.0;
static const RDenseMatrix intfd8d8 = RDenseMatrix (10, 9,
   "-48 -48 -48 -48 -48 -48 -48 -48 -48 \
    -48 -32  48 -48 -32  48 -48 -32  48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
     48 -32 -48  48 -32 -48  48 -32 -48 \
     64 128 192  64 128 192  64 128 192 \
     64  64  64  64  64  64  64  64  64 \
    192 128  64 192 128  64 192 128  64 \
     64 128 192  64 128 192  64 128 192 \
    192 256 192 192 256 192 192 256 192 \
    192 128  64 192 128  64 192 128  64")/15120.0;
static const RDenseMatrix intfd8d9 = RDenseMatrix (10, 12,
   "-48 -24 -24 -24 -48 -24 -24 -24 -48 -24 -24 -24 \
    -48 -24 -16 -16 -48 -24 -16 -16 -48 -24 -16 -16 \
    -48 -16 -24 -16 -48 -16 -24 -16 -48 -16 -24 -16 \
     48 -16 -16 -24  48 -16 -16 -24  48 -16 -16 -24 \
     64  32  64  64  64  32  64  64  64  32  64  64 \
     64  64  32  64  64  64  32  64  64  64  32  64 \
    192  64  64  32 192  64  64  32 192  64  64  32 \
     64  64  64 128  64  64  64 128  64  64  64 128 \
    192  64 128  64 192  64 128  64 192  64 128  64 \
    192 128  64  64 192 128  64  64 192 128  64  64")/15120.0;
static const RDenseMatrix intfd9d9 = RDenseMatrix (10, 9,
   "-48 -48 -48 -48 -48 -48 -48 -48 -48 \
    -48 -48 -48 -48 -48 -48 -48 -48 -48 \
    -48 -32  48 -48 -32  48 -48 -32  48 \
     48 -32 -48  48 -32 -48  48 -32 -48 \
     64  64  64  64  64  64  64  64  64 \
     64 128 192  64 128 192  64 128 192 \
    192 128  64 192 128  64 192 128  64 \
     64 128 192  64 128 192  64 128 192 \
    192 128  64 192 128  64 192 128  64 \
    192 256 192 192 256 192 192 256 192")/15120.0;

bool debug_flag = false;
double Tetrahedron10::IntFDD (int i, int j, int k) const 
{
    if (j > k) { int tmp=j; j=k; k=tmp; } // symmetry
    double geom;
    switch (j) {
    case 0:
        switch (k) {
	case 0:
	    geom = (b0*b0+c0*c0+d0*d0)/size;
	    return intfd0d0[i] * geom;
	case 1:
	    geom = (b0*b1+c0*c1+d0*d1)/size;
	    return intfd0d1[i] * geom;
	case 2:
	    geom = (b0*b2+c0*c2+d0*d2)/size;
	    return intfd0d2[i] * geom;
	case 3:
	    geom = (b0*b3+c0*c3+d0*d3)/size;
	    return intfd0d3[i] * geom;
	case 4:
	    return (intfd0d4(i,0)*(b0*b0) +
                    intfd0d4(i,1)*(b0*b1) +
		    intfd0d4(i,2)*(c0*c0) +
		    intfd0d4(i,3)*(c0*c1) +
		    intfd0d4(i,4)*(d0*d0) +
		    intfd0d4(i,5)*(d0*d1))/size;
	case 5:
	    return (intfd0d5(i,0)*(b0*b0) +
		    intfd0d5(i,1)*(b0*b2) +
		    intfd0d5(i,2)*(c0*c0) +
		    intfd0d5(i,3)*(c0*c2) +
		    intfd0d5(i,4)*(d0*d0) +
		    intfd0d5(i,5)*(d0*d2))/size;
	case 6:
	    return (intfd0d6(i,0)*(b0*b0) +
		    intfd0d6(i,1)*(b0*b3) +
		    intfd0d6(i,2)*(c0*c0) +
		    intfd0d6(i,3)*(c0*c3) +
		    intfd0d6(i,4)*(d0*d0) +
		    intfd0d6(i,5)*(d0*d3))/size;
	case 7:
	    return (intfd0d7(i,0)*(b0*b1) +
		    intfd0d7(i,1)*(b0*b2) +
		    intfd0d7(i,2)*(c0*c1) +
		    intfd0d7(i,3)*(c0*c2) +
		    intfd0d7(i,4)*(d0*d1) +
		    intfd0d7(i,5)*(d0*d2))/size;
	case 8:
	    return (intfd0d8(i,0)*(b0*b1) +
		    intfd0d8(i,1)*(b0*b3) +
		    intfd0d8(i,2)*(c0*c1) +
		    intfd0d8(i,3)*(c0*c3) +
		    intfd0d8(i,4)*(d0*d1) +
		    intfd0d8(i,5)*(d0*d3))/size;
	case 9:
	    return (intfd0d9(i,0)*(b0*b2) +
		    intfd0d9(i,1)*(b0*b3) +
		    intfd0d9(i,2)*(c0*c2) +
		    intfd0d9(i,3)*(c0*c3) +
		    intfd0d9(i,4)*(d0*d2) +
		    intfd0d9(i,5)*(d0*d3))/size;
	}
    case 1:
        switch (k) {
	case 1:
	    geom = (b1*b1+c1*c1+d1*d1)/size;
	    return intfd1d1[i] * geom;
	case 2:
	    geom = (b1*b2+c1*c2+d1*d2)/size;
	    return intfd1d2[i] * geom;
	case 3:
	    geom = (b1*b3+c1*c3+d1*d3)/size;
	    return intfd1d3[i] * geom;
	case 4:
	    return (intfd1d4(i,0)*(b0*b1) +
		    intfd1d4(i,1)*(b1*b1) +
		    intfd1d4(i,2)*(c0*c1) +
		    intfd1d4(i,3)*(c1*c1) +
		    intfd1d4(i,4)*(d0*d1) +
		    intfd1d4(i,5)*(d1*d1))/size;
	case 5:
	    return (intfd1d5(i,0)*(b0*b1) +
		    intfd1d5(i,1)*(b1*b2) +
		    intfd1d5(i,2)*(c0*c1) +
		    intfd1d5(i,3)*(c1*c2) +
		    intfd1d5(i,4)*(d0*d1) +
		    intfd1d5(i,5)*(d1*d2))/size;
	case 6:
	    return (intfd1d6(i,0)*(b0*b1) +
		    intfd1d6(i,1)*(b1*b3) +
		    intfd1d6(i,2)*(c0*c1) +
		    intfd1d6(i,3)*(c1*c3) +
		    intfd1d6(i,4)*(d0*d1) +
		    intfd1d6(i,5)*(d1*d3))/size;
	case 7:
	    return (intfd1d7(i,0)*(b1*b1) +
		    intfd1d7(i,1)*(b1*b2) +
		    intfd1d7(i,2)*(c1*c1) +
		    intfd1d7(i,3)*(c1*c2) +
		    intfd1d7(i,4)*(d1*d1) +
		    intfd1d7(i,5)*(d1*d2))/size;
	case 8:
	    return (intfd1d8(i,0)*(b1*b1) +
		    intfd1d8(i,1)*(b1*b3) +
		    intfd1d8(i,2)*(c1*c1) +
		    intfd1d8(i,3)*(c1*c3) +
		    intfd1d8(i,4)*(d1*d1) +
		    intfd1d8(i,5)*(d1*d3))/size;
	case 9:
	    return (intfd1d9(i,0)*(b1*b2) +
		    intfd1d9(i,1)*(b1*b3) +
		    intfd1d9(i,2)*(c1*c2) +
		    intfd1d9(i,3)*(c1*c3) +
		    intfd1d9(i,4)*(d1*d2) +
		    intfd1d9(i,5)*(d1*d3))/size;
	}
    case 2:
        switch (k) {
	case 2:
	    geom = (b2*b2+c2*c2+d2*d2)/size;
	    return intfd2d2[i] * geom;
	case 3:
	    geom = (b2*b3+c2*c3+d2*d3)/size;
	    return intfd2d3[i] * geom;
	case 4:
	    return (intfd2d4(i,0)*(b0*b2) +
		    intfd2d4(i,1)*(b1*b2) +
		    intfd2d4(i,2)*(c0*c2) +
		    intfd2d4(i,3)*(c1*c2) +
		    intfd2d4(i,4)*(d0*d2) +
		    intfd2d4(i,5)*(d1*d2))/size;
	case 5:
	    return (intfd2d5(i,0)*(b0*b2) +
		    intfd2d5(i,1)*(b2*b2) +
		    intfd2d5(i,2)*(c0*c2) +
		    intfd2d5(i,3)*(c2*c2) +
		    intfd2d5(i,4)*(d0*d2) +
		    intfd2d5(i,5)*(d2*d2))/size;
	case 6:
	    return (intfd2d6(i,0)*(b0*b2) +
		    intfd2d6(i,1)*(b2*b3) +
		    intfd2d6(i,2)*(c0*c2) +
		    intfd2d6(i,3)*(c2*c3) +
		    intfd2d6(i,4)*(d0*d2) +
		    intfd2d6(i,5)*(d2*d3))/size;
	case 7:
	    return (intfd2d7(i,0)*(b1*b2) +
		    intfd2d7(i,1)*(b2*b2) +
		    intfd2d7(i,2)*(c1*c2) +
		    intfd2d7(i,3)*(c2*c2) +
		    intfd2d7(i,4)*(d1*d2) +
		    intfd2d7(i,5)*(d2*d2))/size;
	case 8:
	    return (intfd2d8(i,0)*(b1*b2) +
		    intfd2d8(i,1)*(b2*b3) +
		    intfd2d8(i,2)*(c1*c2) +
		    intfd2d8(i,3)*(c2*c3) +
		    intfd2d8(i,4)*(d1*d2) +
		    intfd2d8(i,5)*(d2*d3))/size;
	case 9:
	    return (intfd2d9(i,0)*(b2*b2) +
		    intfd2d9(i,1)*(b2*b3) +
		    intfd2d9(i,2)*(c2*c2) +
		    intfd2d9(i,3)*(c2*c3) +
		    intfd2d9(i,4)*(d2*d2) +
		    intfd2d9(i,5)*(d2*d3))/size;
	}
    case 3:
        switch (k) {
	case 3:
	    geom = (b3*b3+c3*c3+d3*d3)/size;
	    return intfd3d3[i] * geom;
	case 4:
	    return (intfd3d4(i,0)*(b0*b3) +
		    intfd3d4(i,1)*(b1*b3) +
		    intfd3d4(i,2)*(c0*c3) +
		    intfd3d4(i,3)*(c1*c3) +
		    intfd3d4(i,4)*(d0*d3) +
		    intfd3d4(i,5)*(d1*d3))/size;
	case 5:
	    return (intfd3d5(i,0)*(b0*b3) +
		    intfd3d5(i,1)*(b2*b3) +
		    intfd3d5(i,2)*(c0*c3) +
		    intfd3d5(i,3)*(c2*c3) +
		    intfd3d5(i,4)*(d0*d3) +
		    intfd3d5(i,5)*(d2*d3))/size;
	case 6:
	    return (intfd3d6(i,0)*(b0*b3) +
		    intfd3d6(i,1)*(b3*b3) +
		    intfd3d6(i,2)*(c0*c3) +
		    intfd3d6(i,3)*(c3*c3) +
		    intfd3d6(i,4)*(d0*d3) +
		    intfd3d6(i,5)*(d3*d3))/size;
	case 7:
	    return (intfd3d7(i,0)*(b1*b3) +
		    intfd3d7(i,1)*(b2*b3) +
		    intfd3d7(i,2)*(c1*c3) +
		    intfd3d7(i,3)*(c2*c3) +
		    intfd3d7(i,4)*(d1*d3) +
		    intfd3d7(i,5)*(d2*d3))/size;
	case 8:
	    return (intfd3d8(i,0)*(b1*b3) +
		    intfd3d8(i,1)*(b3*b3) +
		    intfd3d8(i,2)*(c1*c3) +
		    intfd3d8(i,3)*(c3*c3) +
		    intfd3d8(i,4)*(d1*d3) +
		    intfd3d8(i,5)*(d3*d3))/size;
	case 9:
	    return (intfd3d9(i,0)*(b2*b3) +
		    intfd3d9(i,1)*(b3*b3) +
		    intfd3d9(i,2)*(c2*c3) +
		    intfd3d9(i,3)*(c3*c3) +
		    intfd3d9(i,4)*(d2*d3) +
		    intfd3d9(i,5)*(d3*d3))/size;
	}
    case 4:
        switch (k) {
	case 4:
	    return (intfd4d4(i,0)*(b0*b0) +
		    intfd4d4(i,1)*(b0*b1) +
		    intfd4d4(i,2)*(b1*b1) +
		    intfd4d4(i,3)*(c0*c0) +
		    intfd4d4(i,4)*(c0*c1) +
		    intfd4d4(i,5)*(c1*c1) +
		    intfd4d4(i,6)*(d0*d0) +
		    intfd4d4(i,7)*(d0*d1) +
		    intfd4d4(i,8)*(d1*d1))/size;
	case 5:
	    return (intfd4d5(i,0)*(b0*b0) +
		    intfd4d5(i,1)*(b0*b1) +
		    intfd4d5(i,2)*(b0*b2) +
		    intfd4d5(i,3)*(b1*b2) +
		    intfd4d5(i,4)*(c0*c0) +
		    intfd4d5(i,5)*(c0*c1) +
		    intfd4d5(i,6)*(c0*c2) +
		    intfd4d5(i,7)*(c1*c2) +
		    intfd4d5(i,8)*(d0*d0) +
		    intfd4d5(i,9)*(d0*d1) +
		    intfd4d5(i,10)*(d0*d2) +
		    intfd4d5(i,11)*(d1*d2))/size;
	case 6:
	    return (intfd4d6(i,0)*(b0*b0) +
		    intfd4d6(i,1)*(b0*b1) +
		    intfd4d6(i,2)*(b0*b3) +
		    intfd4d6(i,3)*(b1*b3) +
		    intfd4d6(i,4)*(c0*c0) +
		    intfd4d6(i,5)*(c0*c1) +
		    intfd4d6(i,6)*(c0*c3) +
		    intfd4d6(i,7)*(c1*c3) +
		    intfd4d6(i,8)*(d0*d0) +
		    intfd4d6(i,9)*(d0*d1) +
		    intfd4d6(i,10)*(d0*d3) +
		    intfd4d6(i,11)*(d1*d3))/size;
	case 7:
	    return (intfd4d7(i,0)*(b0*b1) +
		    intfd4d7(i,1)*(b1*b1) +
		    intfd4d7(i,2)*(b0*b2) +
		    intfd4d7(i,3)*(b1*b2) +
		    intfd4d7(i,4)*(c0*c1) +
		    intfd4d7(i,5)*(c1*c1) +
		    intfd4d7(i,6)*(c0*c2) +
		    intfd4d7(i,7)*(c1*c2) +
		    intfd4d7(i,8)*(d0*d1) +
		    intfd4d7(i,9)*(d1*d1) +
		    intfd4d7(i,10)*(d0*d2) +
		    intfd4d7(i,11)*(d1*d2))/size;
	case 8:
	    return (intfd4d8(i,0)*(b0*b1) +
		    intfd4d8(i,1)*(b1*b1) +
		    intfd4d8(i,2)*(b0*b3) +
		    intfd4d8(i,3)*(b1*b3) +
		    intfd4d8(i,4)*(c0*c1) +
		    intfd4d8(i,5)*(c1*c1) +
		    intfd4d8(i,6)*(c0*c3) +
		    intfd4d8(i,7)*(c1*c3) +
		    intfd4d8(i,8)*(d0*d1) +
		    intfd4d8(i,9)*(d1*d1) +
		    intfd4d8(i,10)*(d0*d3) +
		    intfd4d8(i,11)*(d1*d3))/size;
	case 9:
	    return (intfd4d9(i,0)*(b0*b2) +
		    intfd4d9(i,1)*(b1*b2) +
		    intfd4d9(i,2)*(b0*b3) +
		    intfd4d9(i,3)*(b1*b3) +
		    intfd4d9(i,4)*(c0*c2) +
		    intfd4d9(i,5)*(c1*c2) +
		    intfd4d9(i,6)*(c0*c3) +
		    intfd4d9(i,7)*(c1*c3) +
		    intfd4d9(i,8)*(d0*d2) +
		    intfd4d9(i,9)*(d1*d2) +
		    intfd4d9(i,10)*(d0*d3) +
		    intfd4d9(i,11)*(d1*d3))/size;
	}
    case 5:
        switch (k) {
	case 5:
	    return (intfd5d5(i,0)*(b0*b0) +
		    intfd5d5(i,1)*(b0*b2) +
		    intfd5d5(i,2)*(b2*b2) +
		    intfd5d5(i,3)*(c0*c0) +
		    intfd5d5(i,4)*(c0*c2) +
		    intfd5d5(i,5)*(c2*c2) +
		    intfd5d5(i,6)*(d0*d0) +
		    intfd5d5(i,7)*(d0*d2) +
		    intfd5d5(i,8)*(d2*d2))/size;
	case 6:
	    return (intfd5d6(i,0)*(b0*b0) +
		    intfd5d6(i,1)*(b0*b2) +
		    intfd5d6(i,2)*(b0*b3) +
		    intfd5d6(i,3)*(b2*b3) +
		    intfd5d6(i,4)*(c0*c0) +
		    intfd5d6(i,5)*(c0*c2) +
		    intfd5d6(i,6)*(c0*c3) +
		    intfd5d6(i,7)*(c2*c3) +
		    intfd5d6(i,8)*(d0*d0) +
		    intfd5d6(i,9)*(d0*d2) +
		    intfd5d6(i,10)*(d0*d3) +
		    intfd5d6(i,11)*(d2*d3))/size;
	case 7:
	    return (intfd5d7(i,0)*(b0*b1) +
		    intfd5d7(i,1)*(b0*b2) +
		    intfd5d7(i,2)*(b1*b2) +
		    intfd5d7(i,3)*(b2*b2) +
		    intfd5d7(i,4)*(c0*c1) +
		    intfd5d7(i,5)*(c0*c2) +
		    intfd5d7(i,6)*(c1*c2) +
		    intfd5d7(i,7)*(c2*c2) +
		    intfd5d7(i,8)*(d0*d1) +
		    intfd5d7(i,9)*(d0*d2) +
		    intfd5d7(i,10)*(d1*d2) +
		    intfd5d7(i,11)*(d2*d2))/size;
	case 8:
	    return (intfd5d8(i,0)*(b0*b1) +
		    intfd5d8(i,1)*(b1*b2) +
		    intfd5d8(i,2)*(b0*b3) +
		    intfd5d8(i,3)*(b2*b3) +
		    intfd5d8(i,4)*(c0*c1) +
		    intfd5d8(i,5)*(c1*c2) +
		    intfd5d8(i,6)*(c0*c3) +
		    intfd5d8(i,7)*(c2*c3) +
		    intfd5d8(i,8)*(d0*d1) +
		    intfd5d8(i,9)*(d1*d2) +
		    intfd5d8(i,10)*(d0*d3) +
		    intfd5d8(i,11)*(d2*d3))/size;
	case 9:
	    return (intfd5d9(i,0)*(b0*b2) +
		    intfd5d9(i,1)*(b2*b2) +
		    intfd5d9(i,2)*(b0*b3) +
		    intfd5d9(i,3)*(b2*b3) +
		    intfd5d9(i,4)*(c0*c2) +
		    intfd5d9(i,5)*(c2*c2) +
		    intfd5d9(i,6)*(c0*c3) +
		    intfd5d9(i,7)*(c2*c3) +
		    intfd5d9(i,8)*(d0*d2) +
		    intfd5d9(i,9)*(d2*d2) +
		    intfd5d9(i,10)*(d0*d3) +
		    intfd5d9(i,11)*(d2*d3))/size;
	}
    case 6:
        switch (k) {
	case 6:
	    return (intfd6d6(i,0)*(b0*b0) +
		    intfd6d6(i,1)*(b0*b3) +
		    intfd6d6(i,2)*(b3*b3) +
		    intfd6d6(i,3)*(c0*c0) +
		    intfd6d6(i,4)*(c0*c3) +
		    intfd6d6(i,5)*(c3*c3) +
		    intfd6d6(i,6)*(d0*d0) +
		    intfd6d6(i,7)*(d0*d3) +
		    intfd6d6(i,8)*(d3*d3))/size;
	case 7:
	    return (intfd6d7(i,0)*(b0*b1) +
		    intfd6d7(i,1)*(b0*b2) +
		    intfd6d7(i,2)*(b1*b3) +
		    intfd6d7(i,3)*(b2*b3) +
		    intfd6d7(i,4)*(c0*c1) +
		    intfd6d7(i,5)*(c0*c2) +
		    intfd6d7(i,6)*(c1*c3) +
		    intfd6d7(i,7)*(c2*c3) +
		    intfd6d7(i,8)*(d0*d1) +
		    intfd6d7(i,9)*(d0*d2) +
		    intfd6d7(i,10)*(d1*d3) +
		    intfd6d7(i,11)*(d2*d3))/size;
	case 8:
	    return (intfd6d8(i,0)*(b0*b1) +
		    intfd6d8(i,1)*(b0*b3) +
		    intfd6d8(i,2)*(b1*b3) +
		    intfd6d8(i,3)*(b3*b3) +
		    intfd6d8(i,4)*(c0*c1) +
		    intfd6d8(i,5)*(c0*c3) +
		    intfd6d8(i,6)*(c1*c3) +
		    intfd6d8(i,7)*(c3*c3) +
		    intfd6d8(i,8)*(d0*d1) +
		    intfd6d8(i,9)*(d0*d3) +
		    intfd6d8(i,10)*(d1*d3) +
		    intfd6d8(i,11)*(d3*d3))/size;
	case 9:
	    return (intfd6d9(i,0)*(b0*b2) +
		    intfd6d9(i,1)*(b0*b3) +
		    intfd6d9(i,2)*(b2*b3) +
		    intfd6d9(i,3)*(b3*b3) +
		    intfd6d9(i,4)*(c0*c2) +
		    intfd6d9(i,5)*(c0*c3) +
		    intfd6d9(i,6)*(c2*c3) +
		    intfd6d9(i,7)*(c3*c3) +
		    intfd6d9(i,8)*(d0*d2) +
		    intfd6d9(i,9)*(d0*d3) +
		    intfd6d9(i,10)*(d2*d3) +
		    intfd6d9(i,11)*(d3*d3))/size;
	}
    case 7:
        switch (k) {
	case 7:
	    return (intfd7d7(i,0)*(b1*b1) +
		    intfd7d7(i,1)*(b1*b2) +
		    intfd7d7(i,2)*(b2*b2) +
		    intfd7d7(i,3)*(c1*c1) +
		    intfd7d7(i,4)*(c1*c2) +
		    intfd7d7(i,5)*(c2*c2) +
		    intfd7d7(i,6)*(d1*d1) +
		    intfd7d7(i,7)*(d1*d2) +
		    intfd7d7(i,8)*(d2*d2))/size;
	case 8:
	    return (intfd7d8(i,0)*(b1*b1) +
		    intfd7d8(i,1)*(b1*b2) +
		    intfd7d8(i,2)*(b1*b3) +
		    intfd7d8(i,3)*(b2*b3) +
		    intfd7d8(i,4)*(c1*c1) +
		    intfd7d8(i,5)*(c1*c2) +
		    intfd7d8(i,6)*(c1*c3) +
		    intfd7d8(i,7)*(c2*c3) +
		    intfd7d8(i,8)*(d1*d1) +
		    intfd7d8(i,9)*(d1*d2) +
		    intfd7d8(i,10)*(d1*d3) +
		    intfd7d8(i,11)*(d2*d3))/size;
	case 9:
	    return (intfd7d9(i,0)*(b1*b2) +
		    intfd7d9(i,1)*(b2*b2) +
		    intfd7d9(i,2)*(b1*b3) +
		    intfd7d9(i,3)*(b2*b3) +
		    intfd7d9(i,4)*(c1*c2) +
		    intfd7d9(i,5)*(c2*c2) +
		    intfd7d9(i,6)*(c1*c3) +
		    intfd7d9(i,7)*(c2*c3) +
		    intfd7d9(i,8)*(d1*d2) +
		    intfd7d9(i,9)*(d2*d2) +
		    intfd7d9(i,10)*(d1*d3) +
		    intfd7d9(i,11)*(d2*d3))/size;
	}
    case 8:
        switch (k) {
	case 8:
	    return (intfd8d8(i,0)*(b1*b1) +
		    intfd8d8(i,1)*(b1*b3) +
		    intfd8d8(i,2)*(b3*b3) +
		    intfd8d8(i,3)*(c1*c1) +
		    intfd8d8(i,4)*(c1*c3) +
		    intfd8d8(i,5)*(c3*c3) +
		    intfd8d8(i,6)*(d1*d1) +
		    intfd8d8(i,7)*(d1*d3) +
		    intfd8d8(i,8)*(d3*d3))/size;
	case 9:
	    return (intfd8d9(i,0)*(b1*b2) +
		    intfd8d9(i,1)*(b1*b3) +
		    intfd8d9(i,2)*(b2*b3) +
		    intfd8d9(i,3)*(b3*b3) +
		    intfd8d9(i,4)*(c1*c2) +
		    intfd8d9(i,5)*(c1*c3) +
		    intfd8d9(i,6)*(c2*c3) +
		    intfd8d9(i,7)*(c3*c3) +
		    intfd8d9(i,8)*(d1*d2) +
		    intfd8d9(i,9)*(d1*d3) +
		    intfd8d9(i,10)*(d2*d3) +
		    intfd8d9(i,11)*(d3*d3))/size;
	}
    case 9:
        switch (k) {
	case 9:
	    return (intfd9d9(i,0)*(b2*b2) +
		    intfd9d9(i,1)*(b2*b3) +
		    intfd9d9(i,2)*(b3*b3) +
		    intfd9d9(i,3)*(c2*c2) +
		    intfd9d9(i,4)*(c2*c3) +
		    intfd9d9(i,5)*(c3*c3) +
		    intfd9d9(i,6)*(d2*d2) +
		    intfd9d9(i,7)*(d2*d3) +
		    intfd9d9(i,8)*(d3*d3))/size;
	}
    }
    xERROR("Index out of range");
    return 0;
}

double Tetrahedron10::IntPDD (int i, int j, const RVector &P) const
{
    double res = 0.0;
    for (int k = 0; k < 10; k++)
        res += P[Node[k]] * IntFDD (k,i,j);
    return res;
}

// boundary integrals of shape functions on local element
static const RDenseMatrix bndintf = RDenseMatrix (4, 10,
    "0 0 0 0 1 1 0 1 0 0\
     0 0 0 0 1 0 1 0 1 0\
     0 0 0 0 0 1 1 0 0 1\
     0 0 0 0 0 0 0 1 1 1") * (2.0/6.0);

static const RSymMatrix sym_bndintff_sd0 = RSymMatrix (10,
   " 6 \
    -1  6 \
    -1 -1  6 \
     0  0  0  0 \
     0  0 -4  0 32 \
     0 -4  0  0 16 32 \
     0  0  0  0  0  0  0 \
    -4  0  0  0 16 16  0 32 \
     0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/360.0);
// see tet4.cc for reason for 2 in numerator
static const RSymMatrix sym_bndintff_sd1 = RSymMatrix (10,
   " 6 \
    -1  6 \
     0  0  0 \
    -1 -1  0  6 \
     0  0  0 -4 32 \
     0  0  0  0  0  0 \
     0 -4  0  0 16  0 32 \
     0  0  0  0  0  0  0  0 \
    -4  0  0  0 16  0 16  0 32 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/360.0);
static const RSymMatrix sym_bndintff_sd2 = RSymMatrix (10,
   " 6 \
     0  0 \
    -1  0  6 \
    -1  0 -1  6 \
     0  0  0  0  0 \
     0  0  0 -4  0 32 \
     0  0 -4  0  0 16 32 \
     0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0 \
    -4  0  0  0  0 16 16  0  0 32") * (2.0/360.0);
static const RSymMatrix sym_bndintff_sd3 = RSymMatrix (10,
   " 0 \
     0  6 \
     0 -1  6 \
     0 -1 -1  6 \
     0  0  0  0  0 \
     0  0  0  0  0  0 \
     0  0  0  0  0  0  0 \
     0  0  0 -4  0  0  0 32 \
     0  0 -4  0  0  0  0 16 32 \
     0 -4  0  0  0  0  0 16 16 32") * (2.0/360.0);
static const RSymMatrix *sym_bndintff[4] = {
   &sym_bndintff_sd0,
   &sym_bndintff_sd1,
   &sym_bndintff_sd2,
   &sym_bndintff_sd3
};

// boundary integrals of products of 3 shape functions on local element
static const RDenseMatrix sd0_intf0ff = RDenseMatrix (10, 10,
   "18 -2 -2  0 12 12  0  4  0  0 \
    -2 -2  1  0 -4  0  0  0  0  0 \
    -2  1 -2  0  0 -4  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    12 -4  0  0  0  0  0 -8  0  0 \
    12  0 -4  0  0  0  0 -8  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     4  0  0  0 -8 -8  0 -16 0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd0_intf1ff = RDenseMatrix (10, 10,
   "-2 -2  1  0 -4  0  0  0  0  0 \
    -2 18 -2  0 12  4  0 12  0  0 \
     1 -2 -2  0  0  0  0 -4  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -4 12  0  0  0 -8  0  0  0  0 \
     0  4  0  0 -8 -16 0 -8  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 12 -4  0  0 -8  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd0_intf2ff = RDenseMatrix (10, 10,
   "-2  1 -2  0  0 -4  0  0  0  0 \
     1 -2 -2  0  0  0  0 -4  0  0 \
    -2 -2 18  0  4 12  0 12  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  4  0 -16 -8 0 -8  0  0 \
    -4  0 12  0 -8  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -4 12  0 -8  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd0_intf3ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd0_intf4ff = RDenseMatrix (10, 10,
   "12 -4  0  0  0  0  0 -8  0  0 \
    -4 12  0  0  0 -8  0  0  0  0 \
     0  0  4  0 -16 -8 0 -8  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0 -16 0 144 48 0 48  0  0 \
     0 -8 -8  0 48 48  0 32  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8  0 -8  0 48 32  0 48  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd0_intf5ff = RDenseMatrix (10, 10,
   "12  0 -4  0  0  0  0 -8  0  0 \
     0  4  0  0 -8 -16 0 -8  0  0 \
    -4  0 12  0 -8  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -8 -8  0 48 48  0 32  0  0 \
     0 -16 0  0 48 144 0 48  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8 -8  0  0 32 48  0 48  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd0_intf6ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd0_intf7ff = RDenseMatrix (10, 10,
   " 4  0  0  0 -8 -8  0 -16 0  0 \
     0 12 -4  0  0 -8  0  0  0  0 \
     0 -4 12  0 -8  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8  0 -8  0 48 32  0 48  0  0 \
    -8 -8  0  0 32 48  0 48  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -16 0  0  0 48 48  0 144 0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd0_intf8ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd0_intf9ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix *sd0_intfff[10] = {
    &sd0_intf0ff,
    &sd0_intf1ff,
    &sd0_intf2ff,
    &sd0_intf3ff,
    &sd0_intf4ff,
    &sd0_intf5ff,
    &sd0_intf6ff,
    &sd0_intf7ff,
    &sd0_intf8ff,
    &sd0_intf9ff
};
static const RDenseMatrix sd1_intf0ff = RDenseMatrix (10, 10,
   "18 -2  0 -2 12  0 12  0  4  0 \
    -2 -2  0  1 -4  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -2  1  0 -2  0  0 -4  0  0  0 \
    12 -4  0  0  0  0  0  0 -8  0 \
     0  0  0  0  0  0  0  0  0  0 \
    12  0  0 -4  0  0  0  0 -8  0 \
     0  0  0  0  0  0  0  0  0  0 \
     4  0  0  0 -8  0 -8  0 -16 0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd1_intf1ff = RDenseMatrix (10, 10,
   "-2 -2  0  1 -4  0  0  0  0  0 \
    -2 18  0 -2 12  0  4  0 12  0 \
     0  0  0  0  0  0  0  0  0  0 \
     1 -2  0 -2  0  0  0  0 -4  0 \
    -4 12  0  0  0  0 -8  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  4  0  0 -8  0 -16 0 -8  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 12  0 -4  0  0 -8  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd1_intf2ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd1_intf3ff = RDenseMatrix (10, 10,
   "-2  1  0 -2  0  0 -4  0  0  0 \
     1 -2  0 -2  0  0  0  0 -4  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -2 -2  0 18  4  0 12  0 12  0 \
     0  0  0  4 -16 0 -8  0 -8  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -4  0  0 12 -8  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -4  0 12 -8  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd1_intf4ff = RDenseMatrix (10, 10,
   "12 -4  0  0  0  0  0  0 -8  0 \
    -4 12  0  0  0  0 -8  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  4 -16 0 -8  0 -8  0 \
     0  0  0 -16 144 0 48 0 48  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -8  0 -8 48  0 48  0 32  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8  0  0 -8 48  0 32  0 48  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd1_intf5ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd1_intf6ff = RDenseMatrix (10, 10,
   "12  0  0 -4  0  0  0  0 -8  0 \
     0  4  0  0 -8  0 -16 0 -8  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -4  0  0 12 -8  0  0  0  0  0 \
     0 -8  0 -8 48  0 48  0 32  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -16 0  0 48  0 144 0 48  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8 -8  0  0 32  0 48  0 48  0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd1_intf7ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd1_intf8ff = RDenseMatrix (10, 10,
   " 4  0  0  0 -8  0 -8  0 -16 0 \
     0 12  0 -4  0  0 -8  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -4  0 12 -8  0  0  0  0  0 \
    -8  0  0 -8 48  0 32  0 48  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8 -8  0  0 32  0 48  0 48  0 \
     0  0  0  0  0  0  0  0  0  0 \
   -16  0  0  0 48  0 48  0 144 0 \
     0  0  0  0  0  0  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd1_intf9ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix *sd1_intfff[10] = {
    &sd1_intf0ff,
    &sd1_intf1ff,
    &sd1_intf2ff,
    &sd1_intf3ff,
    &sd1_intf4ff,
    &sd1_intf5ff,
    &sd1_intf6ff,
    &sd1_intf7ff,
    &sd1_intf8ff,
    &sd1_intf9ff
};
static const RDenseMatrix sd2_intf0ff = RDenseMatrix (10, 10,
   "18  0 -2 -2  0 12 12  0  0  4 \
     0  0  0  0  0  0  0  0  0  0 \
    -2  0 -2  1  0 -4  0  0  0  0 \
    -2  0  1 -2  0  0 -4  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    12  0 -4  0  0  0  0  0  0 -8 \
    12  0  0 -4  0  0  0  0  0 -8 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     4  0  0  0  0 -8 -8  0  0 -16") * (2.0/2520.0);
static const RDenseMatrix sd2_intf1ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd2_intf2ff = RDenseMatrix (10, 10,
   "-2  0 -2  1  0 -4  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -2  0 18 -2  0 12  4  0  0 12 \
     1  0 -2 -2  0  0  0  0  0 -4 \
     0  0  0  0  0  0  0  0  0  0 \
    -4  0 12  0  0  0 -8  0  0  0 \
     0  0  4  0  0 -8 -16 0  0 -8 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0 12 -4  0  0 -8  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd2_intf3ff = RDenseMatrix (10, 10,
   "-2  0  1 -2  0  0 -4  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     1  0 -2 -2  0  0  0  0  0 -4 \
    -2  0 -2 18  0  4 12  0  0 12 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  4  0 -16 -8 0  0 -8 \
    -4  0  0 12  0 -8  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0 -4 12  0 -8  0  0  0  0") * (2.0/2520.0);
static const RDenseMatrix sd2_intf4ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd2_intf5ff = RDenseMatrix (10, 10,
   "12  0 -4  0  0  0  0  0  0 -8 \
     0  0  0  0  0  0  0  0  0  0 \
    -4  0 12  0  0  0 -8  0  0  0 \
     0  0  0  4  0 -16 -8 0  0 -8 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0 -16 0 144 48 0  0 48 \
     0  0 -8 -8  0 48 48  0  0 32 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8  0  0 -8  0 48 32  0  0 48") * (2.0/2520.0);
static const RDenseMatrix sd2_intf6ff = RDenseMatrix (10, 10,
   "12  0  0 -4  0  0  0  0  0 -8 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  4  0  0 -8 -16 0  0 -8 \
    -4  0  0 12  0 -8  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0 -8 -8  0 48 48  0  0 32 \
     0  0 -16 0  0 48 144 0  0 48 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8  0 -8  0  0 32 48  0  0 48") * (2.0/2520.0);
static const RDenseMatrix sd2_intf7ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd2_intf8ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd2_intf9ff = RDenseMatrix (10, 10,
   " 4  0  0  0  0 -8 -8  0  0 -16 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0 12 -4  0  0 -8  0  0  0 \
     0  0 -4 12  0 -8  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
    -8  0  0 -8  0 48 32  0  0 48 \
    -8  0 -8  0  0 32 48  0  0 48 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
   -16  0  0  0  0 48 48  0  0 144") * (2.0/2520.0);
static const RDenseMatrix *sd2_intfff[10] = {
    &sd2_intf0ff,
    &sd2_intf1ff,
    &sd2_intf2ff,
    &sd2_intf3ff,
    &sd2_intf4ff,
    &sd2_intf5ff,
    &sd2_intf6ff,
    &sd2_intf7ff,
    &sd2_intf8ff,
    &sd2_intf9ff
};
static const RDenseMatrix sd3_intf0ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd3_intf1ff = RDenseMatrix (10, 10,
   " 0  0  0  0  0  0  0  0  0  0 \
     0 18 -2 -2  0  0  0 12 12  4 \
     0 -2 -2  1  0  0  0 -4  0  0 \
     0 -2  1 -2  0  0  0  0 -4  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 12 -4  0  0  0  0  0  0 -8 \
     0 12  0 -4  0  0  0  0  0 -8 \
     0  4  0  0  0  0  0 -8 -8 -16") * (2.0/2520.0);
static const RDenseMatrix sd3_intf2ff = RDenseMatrix (10, 10,
   " 0  0  0  0  0  0  0  0  0  0 \
     0 -2 -2  1  0  0  0 -4  0  0 \
     0 -2 18 -2  0  0  0 12  4 12 \
     0  1 -2 -2  0  0  0  0  0 -4 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -4 12  0  0  0  0  0 -8  0 \
     0  0  4  0  0  0  0 -8 -16 -8 \
     0  0 12 -4  0  0  0  0 -8  0") * (2.0/2520.0);
static const RDenseMatrix sd3_intf3ff = RDenseMatrix (10, 10,
   " 0  0  0  0  0  0  0  0  0  0 \
     0 -2  1 -2  0  0  0  0 -4  0 \
     0  1 -2 -2  0  0  0  0  0 -4 \
     0 -2 -2 18  0  0  0  4 12 12 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  4  0  0  0 -16 -8 -8 \
     0 -4  0 12  0  0  0 -8  0  0 \
     0  0 -4 12  0  0  0 -8  0  0") * (2.0/2520.0);
static const RDenseMatrix sd3_intf4ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd3_intf5ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd3_intf6ff = RDenseMatrix (10, 10); // = 0
static const RDenseMatrix sd3_intf7ff = RDenseMatrix (10, 10,
   " 0  0  0  0  0  0  0  0  0  0 \
     0 12 -4  0  0  0  0  0  0 -8 \
     0 -4 12  0  0  0  0  0 -8  0 \
     0  0  0  4  0  0  0 -16 -8 -8 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0 -16 0  0  0 144 48 48 \
     0  0 -8 -8  0  0  0 48 48 32 \
     0 -8  0 -8  0  0  0 48 32 48") * (2.0/2520.0);
static const RDenseMatrix sd3_intf8ff = RDenseMatrix (10, 10,
   " 0  0  0  0  0  0  0  0  0  0 \
     0 12  0 -4  0  0  0  0  0 -8 \
     0  0  4  0  0  0  0 -8 -16 -8 \
     0 -4  0 12  0  0  0 -8  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0 -8 -8  0  0  0 48 48 32 \
     0  0 -16 0  0  0  0 48 144 48 \
     0 -8 -8  0  0  0  0 32 48 48") * (2.0/2520.0);
static const RDenseMatrix sd3_intf9ff = RDenseMatrix (10, 10,
   " 0  0  0  0  0  0  0  0  0  0 \
     0  4  0  0  0  0  0 -8 -8 -16 \
     0  0 12 -4  0  0  0  0 -8  0 \
     0  0 -4 12  0  0  0 -8  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0  0  0  0  0  0  0  0  0  0 \
     0 -8  0 -8  0  0  0 48 32 48 \
     0 -8 -8  0  0  0  0 32 48 48 \
     0 -16 0  0  0  0  0 48 48 144") * (2.0/2520.0);
static const RDenseMatrix *sd3_intfff[10] = {
    &sd3_intf0ff,
    &sd3_intf1ff,
    &sd3_intf2ff,
    &sd3_intf3ff,
    &sd3_intf4ff,
    &sd3_intf5ff,
    &sd3_intf6ff,
    &sd3_intf7ff,
    &sd3_intf8ff,
    &sd3_intf9ff
};

static const RDenseMatrix **bndintfff[4] = {
    sd0_intfff,
    sd1_intfff,
    sd2_intfff,
    sd3_intfff
};

double Tetrahedron10::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 10, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 4, "Argument 2: out of range");
    return bndintf(sd,i) * side_size[sd];
}

double Tetrahedron10::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 10, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 10, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 4, "Argument 3: out of range");
    return sym_bndintff[sd]->Get(i,j) * side_size[sd];
}

double Tetrahedron10::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    for (int sd = 0; sd < 4; sd++) {
        if (bndside[sd]) {
	    double sres = 0.0;
	    for (int k = 0; k < 10; k++)
	        sres += P[Node[k]] * bndintfff[sd][k]->Get(i,j);
	    res += sres * side_size[sd];
	}
    }
    return res;
}

int Tetrahedron10::GetLocalSubsampleAbsc (const Point *&absc) const
{
    absc = absc_sample;
    return nsample_tot;
}

int Tetrahedron10::GetBndSubsampleAbsc (int side, const Point *&absc) const
{
    extern int Triangle_GetLocalSubsampleAbsc (const Point *&);
    return Triangle_GetLocalSubsampleAbsc (absc);
}

inline double Power (double x, int y)
{
    if (y == 2) return x*x;
    else        return pow(x,(double)y);
}

RDenseMatrix Tetrahedron10::ElasticityStiffnessMatrix (double E, double nu)
    const
{
    int i, j;

    double mu = E/(2.0*(1.0+nu));                 // shear modulus
    double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu)); // Lame modulus

    double r = lambda + 2.0*mu;
    double s = lambda;
    double t = mu;
    
    RDenseMatrix BTDB(30,30);
    double b0b0 = b0*b0, c0c0 = c0*c0, d0d0 = d0*d0;
    double b0c0 = b0*c0, b0d0 = b0*d0, c0d0 = c0*d0;
    double b0b1 = b0*b1, b0c1 = b0*c1, b0d1 = b0*d1;
    double c0c1 = c0*c1, c0d1 = c0*d1, d0d1 = d0*d1;

    BTDB( 0, 0) = 3*b0b0*r + (3*c0c0 + 3*d0d0)*t;
    BTDB( 0, 1) = 3*b0c0*s + 3*b0c0*t;
    BTDB( 0, 2) = 3*b0d0*s + 3*b0d0*t;
    BTDB( 0, 3) = -(b0b1*r) + (-c0c1 - d0d1)*t;
    BTDB( 0, 4) = -(b0c1*s) - b1*c0*t;
    BTDB( 0, 5) = -(b0d1*s) - b1*d0*t;
    BTDB( 0, 6) = -(b0*b2*r) + (-(c0*c2) - d0*d2)*t;
    BTDB( 0, 7) = -(b0*c2*s) - b2*c0*t;
    BTDB( 0, 8) = -(b0*d2*s) - b2*d0*t;
    BTDB( 0, 9) = -(b0*b3*r) + (-(c0*c3) - d0*d3)*t;
    BTDB( 0,10) = -(b0*c3*s) - b3*c0*t;
    BTDB( 0,11) = -(b0*d3*s) - b3*d0*t;
    BTDB( 0,12) = (-b0b0 + 3*b0b1)*r + (-c0c0 + 3*c0c1 - d0d0 + 3*d0d1)*t;
    BTDB( 0,13) = (-b0c0 + 3*b0c1)*s + (-b0c0 + 3*b1*c0)*t;
    BTDB( 0,14) = (-b0d0 + 3*b0d1)*s + (-b0d0 + 3*b1*d0)*t;
    BTDB( 0,15) = (-b0b0 + 3*b0*b2)*r + (-c0c0 + 3*c0*c2 - d0d0 + 3*d0*d2)*t;
    BTDB( 0,16) = (-b0c0 + 3*b0*c2)*s + (-b0c0 + 3*b2*c0)*t;
    BTDB( 0,17) = (-b0d0 + 3*b0*d2)*s + (-b0d0 + 3*b2*d0)*t;
    BTDB( 0,18) = (-b0b0 + 3*b0*b3)*r + (-c0c0 + 3*c0*c3 - d0d0 + 3*d0*d3)*t;
    BTDB( 0,19) = (-b0c0 + 3*b0*c3)*s + (-b0c0 + 3*b3*c0)*t;
    BTDB( 0,20) = (-b0d0 + 3*b0*d3)*s + (-b0d0 + 3*b3*d0)*t;
    BTDB( 0,21) = (-b0b1 - b0*b2)*r + (-c0c1 - c0*c2 - d0d1 - d0*d2)*t;
    BTDB( 0,22) = (-b0c1 - b0*c2)*s + (-(b1*c0) - b2*c0)*t;
    BTDB( 0,23) = (-b0d1 - b0*d2)*s + (-(b1*d0) - b2*d0)*t;
    BTDB( 0,24) = (-b0b1 - b0*b3)*r + (-c0c1 - c0*c3 - d0d1 - d0*d3)*t;
    BTDB( 0,25) = (-b0c1 - b0*c3)*s + (-(b1*c0) - b3*c0)*t;
    BTDB( 0,26) = (-b0d1 - b0*d3)*s + (-(b1*d0) - b3*d0)*t;
    BTDB( 0,27) = (-(b0*b2) - b0*b3)*r + (-(c0*c2) - c0*c3 - d0*d2 - d0*d3)*t;
    BTDB( 0,28) = (-(b0*c2) - b0*c3)*s + (-(b2*c0) - b3*c0)*t;
    BTDB( 0,29) = (-(b0*d2) - b0*d3)*s + (-(b2*d0) - b3*d0)*t;

    BTDB( 1, 1) = 3*c0c0*r + (3*b0b0 + 3*d0d0)*t;
    BTDB( 1, 2) = 3*c0d0*s + 3*c0d0*t;
    BTDB( 1, 3) = -(b1*c0*s) - b0c1*t;
    BTDB( 1, 4) = -(c0c1*r) + (-b0b1 - d0d1)*t;
    BTDB( 1, 5) = -(c0d1*s) - c1*d0*t;
    BTDB( 1, 6) = -(b2*c0*s) - b0*c2*t;
    BTDB( 1, 7) = -(c0*c2*r) + (-(b0*b2) - d0*d2)*t;
    BTDB( 1, 8) = -(c0*d2*s) - c2*d0*t;
    BTDB( 1, 9) = -(b3*c0*s) - b0*c3*t;
    BTDB( 1,10) = -(c0*c3*r) + (-(b0*b3) - d0*d3)*t;
    BTDB( 1,11) = -(c0*d3*s) - c3*d0*t;
    BTDB( 1,12) = (-b0c0 + 3*b1*c0)*s + (-b0c0 + 3*b0c1)*t;
    BTDB( 1,13) = (-c0c0 + 3*c0c1)*r + (-b0b0 + 3*b0b1 - d0d0 + 3*d0d1)*t;
    BTDB( 1,14) = (-c0d0 + 3*c0d1)*s + (-c0d0 + 3*c1*d0)*t;
    BTDB( 1,15) = (-b0c0 + 3*b2*c0)*s + (-b0c0 + 3*b0*c2)*t;
    BTDB( 1,16) = (-c0c0 + 3*c0*c2)*r + (-b0b0 + 3*b0*b2 - d0d0 + 3*d0*d2)*t;
    BTDB( 1,17) = (-c0d0 + 3*c0*d2)*s + (-c0d0 + 3*c2*d0)*t;
    BTDB( 1,18) = (-b0c0 + 3*b3*c0)*s + (-b0c0 + 3*b0*c3)*t;
    BTDB( 1,19) = (-c0c0 + 3*c0*c3)*r + (-b0b0 + 3*b0*b3 - d0d0 + 3*d0*d3)*t;
    BTDB( 1,20) = (-c0d0 + 3*c0*d3)*s + (-c0d0 + 3*c3*d0)*t;
    BTDB( 1,21) = (-(b1*c0) - b2*c0)*s + (-b0c1 - b0*c2)*t;
    BTDB( 1,22) = (-c0c1 - c0*c2)*r + (-b0b1 - b0*b2 - d0d1 - d0*d2)*t;
    BTDB( 1,23) = (-c0d1 - c0*d2)*s + (-(c1*d0) - c2*d0)*t;
    BTDB( 1,24) = (-(b1*c0) - b3*c0)*s + (-b0c1 - b0*c3)*t;
    BTDB( 1,25) = (-c0c1 - c0*c3)*r + (-b0b1 - b0*b3 - d0d1 - d0*d3)*t;
    BTDB( 1,26) = (-c0d1 - c0*d3)*s + (-(c1*d0) - c3*d0)*t;
    BTDB( 1,27) = (-(b2*c0) - b3*c0)*s + (-(b0*c2) - b0*c3)*t;
    BTDB( 1,28) = (-(c0*c2) - c0*c3)*r + (-(b0*b2) - b0*b3 - d0*d2 - d0*d3)*t;
    BTDB( 1,29) = (-(c0*d2) - c0*d3)*s + (-(c2*d0) - c3*d0)*t;

    BTDB( 2, 2) = 3*d0d0*r + (3*b0b0 + 3*c0c0)*t;
    BTDB( 2, 3) = -(b1*d0*s) - b0d1*t;
    BTDB( 2, 4) = -(c1*d0*s) - c0d1*t;
    BTDB( 2, 5) = -(d0d1*r) + (-b0b1 - c0c1)*t;
    BTDB( 2, 6) = -(b2*d0*s) - b0*d2*t;
    BTDB( 2, 7) = -(c2*d0*s) - c0*d2*t;
    BTDB( 2, 8) = -(d0*d2*r) + (-(b0*b2) - c0*c2)*t;
    BTDB( 2, 9) = -(b3*d0*s) - b0*d3*t;
    BTDB( 2,10) = -(c3*d0*s) - c0*d3*t;
    BTDB( 2,11) = -(d0*d3*r) + (-(b0*b3) - c0*c3)*t;
    BTDB( 2,12) = (-b0d0 + 3*b1*d0)*s + (-b0d0 + 3*b0d1)*t;
    BTDB( 2,13) = (-c0d0 + 3*c1*d0)*s + (-c0d0 + 3*c0d1)*t;
    BTDB( 2,14) = (-d0d0 + 3*d0d1)*r + (-b0b0 + 3*b0b1 - c0c0 + 3*c0c1)*t;
    BTDB( 2,15) = (-b0d0 + 3*b2*d0)*s + (-b0d0 + 3*b0*d2)*t;
    BTDB( 2,16) = (-c0d0 + 3*c2*d0)*s + (-c0d0 + 3*c0*d2)*t;
    BTDB( 2,17) = (-d0d0 + 3*d0*d2)*r + (-b0b0 + 3*b0*b2 - c0c0 + 3*c0*c2)*t;
    BTDB( 2,18) = (-b0d0 + 3*b3*d0)*s + (-b0d0 + 3*b0*d3)*t;
    BTDB( 2,19) = (-c0d0 + 3*c3*d0)*s + (-c0d0 + 3*c0*d3)*t;
    BTDB( 2,20) = (-d0d0 + 3*d0*d3)*r + (-b0b0 + 3*b0*b3 - c0c0 + 3*c0*c3)*t;
    BTDB( 2,21) = (-(b1*d0) - b2*d0)*s + (-b0d1 - b0*d2)*t;
    BTDB( 2,22) = (-(c1*d0) - c2*d0)*s + (-c0d1 - c0*d2)*t;
    BTDB( 2,23) = (-d0d1 - d0*d2)*r + (-b0b1 - b0*b2 - c0c1 - c0*c2)*t;
    BTDB( 2,24) = (-(b1*d0) - b3*d0)*s + (-b0d1 - b0*d3)*t;
    BTDB( 2,25) = (-(c1*d0) - c3*d0)*s + (-c0d1 - c0*d3)*t;
    BTDB( 2,26) = (-d0d1 - d0*d3)*r + (-b0b1 - b0*b3 - c0c1 - c0*c3)*t;
    BTDB( 2,27) = (-(b2*d0) - b3*d0)*s + (-(b0*d2) - b0*d3)*t;
    BTDB( 2,28) = (-(c2*d0) - c3*d0)*s + (-(c0*d2) - c0*d3)*t;
    BTDB( 2,29) = (-(d0*d2) - d0*d3)*r + (-(b0*b2) - b0*b3 - c0*c2 - c0*c3)*t;

    BTDB( 3, 3) = 3*b1*b1*r + (3*c1*c1 + 3*d1*d1)*t;
    BTDB( 3, 4) = 3*b1*c1*s + 3*b1*c1*t;
    BTDB( 3, 5) = 3*b1*d1*s + 3*b1*d1*t;
    BTDB( 3, 6) = -(b1*b2*r) + (-(c1*c2) - d1*d2)*t;
    BTDB( 3, 7) = -(b1*c2*s) - b2*c1*t;
    BTDB( 3, 8) = -(b1*d2*s) - b2*d1*t;
    BTDB( 3, 9) = -(b1*b3*r) + (-(c1*c3) - d1*d3)*t;
    BTDB( 3,10) = -(b1*c3*s) - b3*c1*t;
    BTDB( 3,11) = -(b1*d3*s) - b3*d1*t;
    BTDB( 3,12) = (3*b0b1 - b1*b1)*r + (3*c0c1 - c1*c1 + 3*d0d1 - d1*d1)*t;
    BTDB( 3,13) = (3*b1*c0 - b1*c1)*s + (3*b0c1 - b1*c1)*t;
    BTDB( 3,14) = (3*b1*d0 - b1*d1)*s + (3*b0d1 - b1*d1)*t;
    BTDB( 3,15) = (-b0b1 - b1*b2)*r + (-c0c1 - c1*c2 - d0d1 - d1*d2)*t;
    BTDB( 3,16) = (-(b1*c0) - b1*c2)*s + (-b0c1 - b2*c1)*t;
    BTDB( 3,17) = (-(b1*d0) - b1*d2)*s + (-b0d1 - b2*d1)*t;
    BTDB( 3,18) = (-b0b1 - b1*b3)*r + (-c0c1 - c1*c3 - d0d1 - d1*d3)*t;
    BTDB( 3,19) = (-(b1*c0) - b1*c3)*s + (-b0c1 - b3*c1)*t;
    BTDB( 3,20) = (-(b1*d0) - b1*d3)*s + (-b0d1 - b3*d1)*t;
    BTDB( 3,21) = (-b1*b1 + 3*b1*b2)*r + 
      (-c1*c1 + 3*c1*c2 - d1*d1 + 3*d1*d2)*t;
    BTDB( 3,22) = (-(b1*c1) + 3*b1*c2)*s + (-(b1*c1) + 3*b2*c1)*t;
    BTDB( 3,23) = (-(b1*d1) + 3*b1*d2)*s + (-(b1*d1) + 3*b2*d1)*t;
    BTDB( 3,24) = (-b1*b1 + 3*b1*b3)*r + 
      (-c1*c1 + 3*c1*c3 - d1*d1 + 3*d1*d3)*t;
    BTDB( 3,25) = (-(b1*c1) + 3*b1*c3)*s + (-(b1*c1) + 3*b3*c1)*t;
    BTDB( 3,26) = (-(b1*d1) + 3*b1*d3)*s + (-(b1*d1) + 3*b3*d1)*t;
    BTDB( 3,27) = (-(b1*b2) - b1*b3)*r + 
      (-(c1*c2) - c1*c3 - d1*d2 - d1*d3)*t;
    BTDB( 3,28) = (-(b1*c2) - b1*c3)*s + (-(b2*c1) - b3*c1)*t;
    BTDB( 3,29) = (-(b1*d2) - b1*d3)*s + (-(b2*d1) - b3*d1)*t;
    
    BTDB( 4, 4) = 3*c1*c1*r + (3*b1*b1 + 3*d1*d1)*t;
    BTDB( 4, 5) = 3*c1*d1*s + 3*c1*d1*t;
    BTDB( 4, 6) = -(b2*c1*s) - b1*c2*t;
    BTDB( 4, 7) = -(c1*c2*r) + (-(b1*b2) - d1*d2)*t;
    BTDB( 4, 8) = -(c1*d2*s) - c2*d1*t;
    BTDB( 4, 9) = -(b3*c1*s) - b1*c3*t;
    BTDB( 4,10) = -(c1*c3*r) + (-(b1*b3) - d1*d3)*t;
    BTDB( 4,11) = -(c1*d3*s) - c3*d1*t;
    BTDB( 4,12) = (3*b0c1 - b1*c1)*s + (3*b1*c0 - b1*c1)*t;
    BTDB( 4,13) = (3*c0c1 - c1*c1)*r + 
      (3*b0b1 - b1*b1 + 3*d0d1 - d1*d1)*t;
    BTDB( 4,14) = (3*c1*d0 - c1*d1)*s + (3*c0d1 - c1*d1)*t;
    BTDB( 4,15) = (-b0c1 - b2*c1)*s + (-(b1*c0) - b1*c2)*t;
    BTDB( 4,16) = (-c0c1 - c1*c2)*r + 
      (-b0b1 - b1*b2 - d0d1 - d1*d2)*t;
    BTDB( 4,17) = (-(c1*d0) - c1*d2)*s + (-c0d1 - c2*d1)*t;
    BTDB( 4,18) = (-b0c1 - b3*c1)*s + (-(b1*c0) - b1*c3)*t;
    BTDB( 4,19) = (-c0c1 - c1*c3)*r + 
      (-b0b1 - b1*b3 - d0d1 - d1*d3)*t;
    BTDB( 4,20) = (-(c1*d0) - c1*d3)*s + (-c0d1 - c3*d1)*t;
    BTDB( 4,21) = (-(b1*c1) + 3*b2*c1)*s + (-(b1*c1) + 3*b1*c2)*t;
    BTDB( 4,22) = (-c1*c1 + 3*c1*c2)*r + 
      (-b1*b1 + 3*b1*b2 - d1*d1 + 3*d1*d2)*t;
    BTDB( 4,23) = (-(c1*d1) + 3*c1*d2)*s + (-(c1*d1) + 3*c2*d1)*t;
    BTDB( 4,24) = (-(b1*c1) + 3*b3*c1)*s + (-(b1*c1) + 3*b1*c3)*t;
    BTDB( 4,25) = (-c1*c1 + 3*c1*c3)*r + 
      (-b1*b1 + 3*b1*b3 - d1*d1 + 3*d1*d3)*t;
    BTDB( 4,26) = (-(c1*d1) + 3*c1*d3)*s + (-(c1*d1) + 3*c3*d1)*t;
    BTDB( 4,27) = (-(b2*c1) - b3*c1)*s + (-(b1*c2) - b1*c3)*t;
    BTDB( 4,28) = (-(c1*c2) - c1*c3)*r + 
      (-(b1*b2) - b1*b3 - d1*d2 - d1*d3)*t;
    BTDB( 4,29) = (-(c1*d2) - c1*d3)*s + (-(c2*d1) - c3*d1)*t;
    
    BTDB( 5, 5) = 3*d1*d1*r + (3*b1*b1 + 3*c1*c1)*t;
    BTDB( 5, 6) = -(b2*d1*s) - b1*d2*t;
    BTDB( 5, 7) = -(c2*d1*s) - c1*d2*t;
    BTDB( 5, 8) = -(d1*d2*r) + (-(b1*b2) - c1*c2)*t;
    BTDB( 5, 9) = -(b3*d1*s) - b1*d3*t;
    BTDB( 5,10) = -(c3*d1*s) - c1*d3*t;
    BTDB( 5,11) = -(d1*d3*r) + (-(b1*b3) - c1*c3)*t;
    BTDB( 5,12) = (3*b0d1 - b1*d1)*s + (3*b1*d0 - b1*d1)*t;
    BTDB( 5,13) = (3*c0d1 - c1*d1)*s + (3*c1*d0 - c1*d1)*t;
    BTDB( 5,14) = (3*d0d1 - d1*d1)*r + 
      (3*b0b1 - b1*b1 + 3*c0c1 - c1*c1)*t;
    BTDB( 5,15) = (-b0d1 - b2*d1)*s + (-(b1*d0) - b1*d2)*t;
    BTDB( 5,16) = (-c0d1 - c2*d1)*s + (-(c1*d0) - c1*d2)*t;
    BTDB( 5,17) = (-d0d1 - d1*d2)*r + 
      (-b0b1 - b1*b2 - c0c1 - c1*c2)*t;
    BTDB( 5,18) = (-b0d1 - b3*d1)*s + (-(b1*d0) - b1*d3)*t;
    BTDB( 5,19) = (-c0d1 - c3*d1)*s + (-(c1*d0) - c1*d3)*t;
    BTDB( 5,20) = (-d0d1 - d1*d3)*r + 
      (-b0b1 - b1*b3 - c0c1 - c1*c3)*t;
    BTDB( 5,21) = (-(b1*d1) + 3*b2*d1)*s + (-(b1*d1) + 3*b1*d2)*t;
    BTDB( 5,22) = (-(c1*d1) + 3*c2*d1)*s + (-(c1*d1) + 3*c1*d2)*t;
    BTDB( 5,23) = (-d1*d1 + 3*d1*d2)*r + 
      (-b1*b1 + 3*b1*b2 - c1*c1 + 3*c1*c2)*t;
    BTDB( 5,24) = (-(b1*d1) + 3*b3*d1)*s + (-(b1*d1) + 3*b1*d3)*t;
    BTDB( 5,25) = (-(c1*d1) + 3*c3*d1)*s + (-(c1*d1) + 3*c1*d3)*t;
    BTDB( 5,26) = (-d1*d1 + 3*d1*d3)*r + 
      (-b1*b1 + 3*b1*b3 - c1*c1 + 3*c1*c3)*t;
    BTDB( 5,27) = (-(b2*d1) - b3*d1)*s + (-(b1*d2) - b1*d3)*t;
    BTDB( 5,28) = (-(c2*d1) - c3*d1)*s + (-(c1*d2) - c1*d3)*t;
    BTDB( 5,29) = (-(d1*d2) - d1*d3)*r + 
      (-(b1*b2) - b1*b3 - c1*c2 - c1*c3)*t;
    
    BTDB( 6, 6) = 3*b2*b2*r + (3*c2*c2 + 3*d2*d2)*t;
    BTDB( 6, 7) = 3*b2*c2*s + 3*b2*c2*t;
    BTDB( 6, 8) = 3*b2*d2*s + 3*b2*d2*t;
    BTDB( 6, 9) = -(b2*b3*r) + (-(c2*c3) - d2*d3)*t;
    BTDB( 6,10) = -(b2*c3*s) - b3*c2*t;
    BTDB( 6,11) = -(b2*d3*s) - b3*d2*t;
    BTDB( 6,12) = (-(b0*b2) - b1*b2)*r + 
      (-(c0*c2) - c1*c2 - d0*d2 - d1*d2)*t;
    BTDB( 6,13) = (-(b2*c0) - b2*c1)*s + (-(b0*c2) - b1*c2)*t;
    BTDB( 6,14) = (-(b2*d0) - b2*d1)*s + (-(b0*d2) - b1*d2)*t;
    BTDB( 6,15) = (3*b0*b2 - b2*b2)*r + 
      (3*c0*c2 - c2*c2 + 3*d0*d2 - d2*d2)*t;
    BTDB( 6,16) = (3*b2*c0 - b2*c2)*s + (3*b0*c2 - b2*c2)*t;
    BTDB( 6,17) = (3*b2*d0 - b2*d2)*s + (3*b0*d2 - b2*d2)*t;
    BTDB( 6,18) = (-(b0*b2) - b2*b3)*r + 
      (-(c0*c2) - c2*c3 - d0*d2 - d2*d3)*t;
    BTDB( 6,19) = (-(b2*c0) - b2*c3)*s + (-(b0*c2) - b3*c2)*t;
    BTDB( 6,20) = (-(b2*d0) - b2*d3)*s + (-(b0*d2) - b3*d2)*t;
    BTDB( 6,21) = (3*b1*b2 - b2*b2)*r + 
      (3*c1*c2 - c2*c2 + 3*d1*d2 - d2*d2)*t;
    BTDB( 6,22) = (3*b2*c1 - b2*c2)*s + (3*b1*c2 - b2*c2)*t;
    BTDB( 6,23) = (3*b2*d1 - b2*d2)*s + (3*b1*d2 - b2*d2)*t;
    BTDB( 6,24) = (-(b1*b2) - b2*b3)*r + 
      (-(c1*c2) - c2*c3 - d1*d2 - d2*d3)*t;
    BTDB( 6,25) = (-(b2*c1) - b2*c3)*s + (-(b1*c2) - b3*c2)*t;
    BTDB( 6,26) = (-(b2*d1) - b2*d3)*s + (-(b1*d2) - b3*d2)*t;
    BTDB( 6,27) = (-b2*b2 + 3*b2*b3)*r + 
      (-c2*c2 + 3*c2*c3 - d2*d2 + 3*d2*d3)*t;
    BTDB( 6,28) = (-(b2*c2) + 3*b2*c3)*s + (-(b2*c2) + 3*b3*c2)*t;
    BTDB( 6,29) = (-(b2*d2) + 3*b2*d3)*s + (-(b2*d2) + 3*b3*d2)*t;
    
    BTDB( 7, 7) = 3*c2*c2*r + (3*b2*b2 + 3*d2*d2)*t;
    BTDB( 7, 8) = 3*c2*d2*s + 3*c2*d2*t;
    BTDB( 7, 9) = -(b3*c2*s) - b2*c3*t;
    BTDB( 7,10) = -(c2*c3*r) + (-(b2*b3) - d2*d3)*t;
    BTDB( 7,11) = -(c2*d3*s) - c3*d2*t;
    BTDB( 7,12) = (-(b0*c2) - b1*c2)*s + (-(b2*c0) - b2*c1)*t;
    BTDB( 7,13) = (-(c0*c2) - c1*c2)*r + 
      (-(b0*b2) - b1*b2 - d0*d2 - d1*d2)*t;
    BTDB( 7,14) = (-(c2*d0) - c2*d1)*s + (-(c0*d2) - c1*d2)*t;
    BTDB( 7,15) = (3*b0*c2 - b2*c2)*s + (3*b2*c0 - b2*c2)*t;
    BTDB( 7,16) = (3*c0*c2 - c2*c2)*r + 
      (3*b0*b2 - b2*b2 + 3*d0*d2 - d2*d2)*t;
    BTDB( 7,17) = (3*c2*d0 - c2*d2)*s + (3*c0*d2 - c2*d2)*t;
    BTDB( 7,18) = (-(b0*c2) - b3*c2)*s + (-(b2*c0) - b2*c3)*t;
    BTDB( 7,19) = (-(c0*c2) - c2*c3)*r + 
      (-(b0*b2) - b2*b3 - d0*d2 - d2*d3)*t;
    BTDB( 7,20) = (-(c2*d0) - c2*d3)*s + (-(c0*d2) - c3*d2)*t;
    BTDB( 7,21) = (3*b1*c2 - b2*c2)*s + (3*b2*c1 - b2*c2)*t;
    BTDB( 7,22) = (3*c1*c2 - c2*c2)*r + 
      (3*b1*b2 - b2*b2 + 3*d1*d2 - d2*d2)*t;
    BTDB( 7,23) = (3*c2*d1 - c2*d2)*s + (3*c1*d2 - c2*d2)*t;
    BTDB( 7,24) = (-(b1*c2) - b3*c2)*s + (-(b2*c1) - b2*c3)*t;
    BTDB( 7,25) = (-(c1*c2) - c2*c3)*r + 
      (-(b1*b2) - b2*b3 - d1*d2 - d2*d3)*t;
    BTDB( 7,26) = (-(c2*d1) - c2*d3)*s + (-(c1*d2) - c3*d2)*t;
    BTDB( 7,27) = (-(b2*c2) + 3*b3*c2)*s + (-(b2*c2) + 3*b2*c3)*t;
    BTDB( 7,28) = (-c2*c2 + 3*c2*c3)*r + 
      (-b2*b2 + 3*b2*b3 - d2*d2 + 3*d2*d3)*t;
    BTDB( 7,29) = (-(c2*d2) + 3*c2*d3)*s + (-(c2*d2) + 3*c3*d2)*t;
    
    BTDB( 8, 8) = 3*d2*d2*r + (3*b2*b2 + 3*c2*c2)*t;
    BTDB( 8, 9) = -(b3*d2*s) - b2*d3*t;
    BTDB( 8,10) = -(c3*d2*s) - c2*d3*t;
    BTDB( 8,11) = -(d2*d3*r) + (-(b2*b3) - c2*c3)*t;
    BTDB( 8,12) = (-(b0*d2) - b1*d2)*s + (-(b2*d0) - b2*d1)*t;
    BTDB( 8,13) = (-(c0*d2) - c1*d2)*s + (-(c2*d0) - c2*d1)*t;
    BTDB( 8,14) = (-(d0*d2) - d1*d2)*r + 
      (-(b0*b2) - b1*b2 - c0*c2 - c1*c2)*t;
    BTDB( 8,15) = (3*b0*d2 - b2*d2)*s + (3*b2*d0 - b2*d2)*t;
    BTDB( 8,16) = (3*c0*d2 - c2*d2)*s + (3*c2*d0 - c2*d2)*t;
    BTDB( 8,17) = (3*d0*d2 - d2*d2)*r + 
      (3*b0*b2 - b2*b2 + 3*c0*c2 - c2*c2)*t;
    BTDB( 8,18) = (-(b0*d2) - b3*d2)*s + (-(b2*d0) - b2*d3)*t;
    BTDB( 8,19) = (-(c0*d2) - c3*d2)*s + (-(c2*d0) - c2*d3)*t;
    BTDB( 8,20) = (-(d0*d2) - d2*d3)*r + 
      (-(b0*b2) - b2*b3 - c0*c2 - c2*c3)*t;
    BTDB( 8,21) = (3*b1*d2 - b2*d2)*s + (3*b2*d1 - b2*d2)*t;
    BTDB( 8,22) = (3*c1*d2 - c2*d2)*s + (3*c2*d1 - c2*d2)*t;
    BTDB( 8,23) = (3*d1*d2 - d2*d2)*r + 
      (3*b1*b2 - b2*b2 + 3*c1*c2 - c2*c2)*t;
    BTDB( 8,24) = (-(b1*d2) - b3*d2)*s + (-(b2*d1) - b2*d3)*t;
    BTDB( 8,25) = (-(c1*d2) - c3*d2)*s + (-(c2*d1) - c2*d3)*t;
    BTDB( 8,26) = (-(d1*d2) - d2*d3)*r + 
      (-(b1*b2) - b2*b3 - c1*c2 - c2*c3)*t;
    BTDB( 8,27) = (-(b2*d2) + 3*b3*d2)*s + (-(b2*d2) + 3*b2*d3)*t;
    BTDB( 8,28) = (-(c2*d2) + 3*c3*d2)*s + (-(c2*d2) + 3*c2*d3)*t;
    BTDB( 8,29) = (-d2*d2 + 3*d2*d3)*r + 
      (-b2*b2 + 3*b2*b3 - c2*c2 + 3*c2*c3)*t;
    
    BTDB( 9, 9) = 3*b3*b3*r + (3*c3*c3 + 3*d3*d3)*t;
    BTDB( 9,10) = 3*b3*c3*s + 3*b3*c3*t;
    BTDB( 9,11) = 3*b3*d3*s + 3*b3*d3*t;
    BTDB( 9,12) = (-(b0*b3) - b1*b3)*r + 
      (-(c0*c3) - c1*c3 - d0*d3 - d1*d3)*t;
    BTDB( 9,13) = (-(b3*c0) - b3*c1)*s + (-(b0*c3) - b1*c3)*t;
    BTDB( 9,14) = (-(b3*d0) - b3*d1)*s + (-(b0*d3) - b1*d3)*t;
    BTDB( 9,15) = (-(b0*b3) - b2*b3)*r + 
      (-(c0*c3) - c2*c3 - d0*d3 - d2*d3)*t;
    BTDB( 9,16) = (-(b3*c0) - b3*c2)*s + (-(b0*c3) - b2*c3)*t;
    BTDB( 9,17) = (-(b3*d0) - b3*d2)*s + (-(b0*d3) - b2*d3)*t;
    BTDB( 9,18) = (3*b0*b3 - b3*b3)*r + 
      (3*c0*c3 - c3*c3 + 3*d0*d3 - d3*d3)*t;
    BTDB( 9,19) = (3*b3*c0 - b3*c3)*s + (3*b0*c3 - b3*c3)*t;
    BTDB( 9,20) = (3*b3*d0 - b3*d3)*s + (3*b0*d3 - b3*d3)*t;
    BTDB( 9,21) = (-(b1*b3) - b2*b3)*r + 
      (-(c1*c3) - c2*c3 - d1*d3 - d2*d3)*t;
    BTDB( 9,22) = (-(b3*c1) - b3*c2)*s + (-(b1*c3) - b2*c3)*t;
    BTDB( 9,23) = (-(b3*d1) - b3*d2)*s + (-(b1*d3) - b2*d3)*t;
    BTDB( 9,24) = (3*b1*b3 - b3*b3)*r + 
      (3*c1*c3 - c3*c3 + 3*d1*d3 - d3*d3)*t;
    BTDB( 9,25) = (3*b3*c1 - b3*c3)*s + (3*b1*c3 - b3*c3)*t;
    BTDB( 9,26) = (3*b3*d1 - b3*d3)*s + (3*b1*d3 - b3*d3)*t;
    BTDB( 9,27) = (3*b2*b3 - b3*b3)*r + 
      (3*c2*c3 - c3*c3 + 3*d2*d3 - d3*d3)*t;
    BTDB( 9,28) = (3*b3*c2 - b3*c3)*s + (3*b2*c3 - b3*c3)*t;
    BTDB( 9,29) = (3*b3*d2 - b3*d3)*s + (3*b2*d3 - b3*d3)*t;
    
    BTDB(10,10) = 3*c3*c3*r + (3*b3*b3 + 3*d3*d3)*t;
    BTDB(10,11) = 3*c3*d3*s + 3*c3*d3*t;
    BTDB(10,12) = (-(b0*c3) - b1*c3)*s + (-(b3*c0) - b3*c1)*t;
    BTDB(10,13) = (-(c0*c3) - c1*c3)*r + 
      (-(b0*b3) - b1*b3 - d0*d3 - d1*d3)*t;
    BTDB(10,14) = (-(c3*d0) - c3*d1)*s + (-(c0*d3) - c1*d3)*t;
    BTDB(10,15) = (-(b0*c3) - b2*c3)*s + (-(b3*c0) - b3*c2)*t;
    BTDB(10,16) = (-(c0*c3) - c2*c3)*r + 
      (-(b0*b3) - b2*b3 - d0*d3 - d2*d3)*t;
    BTDB(10,17) = (-(c3*d0) - c3*d2)*s + (-(c0*d3) - c2*d3)*t;
    BTDB(10,18) = (3*b0*c3 - b3*c3)*s + (3*b3*c0 - b3*c3)*t;
    BTDB(10,19) = (3*c0*c3 - c3*c3)*r + 
      (3*b0*b3 - b3*b3 + 3*d0*d3 - d3*d3)*t;
    BTDB(10,20) = (3*c3*d0 - c3*d3)*s + (3*c0*d3 - c3*d3)*t;
    BTDB(10,21) = (-(b1*c3) - b2*c3)*s + (-(b3*c1) - b3*c2)*t;
    BTDB(10,22) = (-(c1*c3) - c2*c3)*r + 
      (-(b1*b3) - b2*b3 - d1*d3 - d2*d3)*t;
    BTDB(10,23) = (-(c3*d1) - c3*d2)*s + (-(c1*d3) - c2*d3)*t;
    BTDB(10,24) = (3*b1*c3 - b3*c3)*s + (3*b3*c1 - b3*c3)*t;
    BTDB(10,25) = (3*c1*c3 - c3*c3)*r + 
      (3*b1*b3 - b3*b3 + 3*d1*d3 - d3*d3)*t;
    BTDB(10,26) = (3*c3*d1 - c3*d3)*s + (3*c1*d3 - c3*d3)*t;
    BTDB(10,27) = (3*b2*c3 - b3*c3)*s + (3*b3*c2 - b3*c3)*t;
    BTDB(10,28) = (3*c2*c3 - c3*c3)*r + 
      (3*b2*b3 - b3*b3 + 3*d2*d3 - d3*d3)*t;
    BTDB(10,29) = (3*c3*d2 - c3*d3)*s + (3*c2*d3 - c3*d3)*t;
    
    BTDB(11,11) = 3*d3*d3*r + (3*b3*b3 + 3*c3*c3)*t;
    BTDB(11,12) = (-(b0*d3) - b1*d3)*s + (-(b3*d0) - b3*d1)*t;
    BTDB(11,13) = (-(c0*d3) - c1*d3)*s + (-(c3*d0) - c3*d1)*t;
    BTDB(11,14) = (-(d0*d3) - d1*d3)*r + 
      (-(b0*b3) - b1*b3 - c0*c3 - c1*c3)*t;
    BTDB(11,15) = (-(b0*d3) - b2*d3)*s + (-(b3*d0) - b3*d2)*t;
    BTDB(11,16) = (-(c0*d3) - c2*d3)*s + (-(c3*d0) - c3*d2)*t;
    BTDB(11,17) = (-(d0*d3) - d2*d3)*r + 
      (-(b0*b3) - b2*b3 - c0*c3 - c2*c3)*t;
    BTDB(11,18) = (3*b0*d3 - b3*d3)*s + (3*b3*d0 - b3*d3)*t;
    BTDB(11,19) = (3*c0*d3 - c3*d3)*s + (3*c3*d0 - c3*d3)*t;
    BTDB(11,20) = (3*d0*d3 - d3*d3)*r + 
      (3*b0*b3 - b3*b3 + 3*c0*c3 - c3*c3)*t;
    BTDB(11,21) = (-(b1*d3) - b2*d3)*s + (-(b3*d1) - b3*d2)*t;
    BTDB(11,22) = (-(c1*d3) - c2*d3)*s + (-(c3*d1) - c3*d2)*t;
    BTDB(11,23) = (-(d1*d3) - d2*d3)*r + 
      (-(b1*b3) - b2*b3 - c1*c3 - c2*c3)*t;
    BTDB(11,24) = (3*b1*d3 - b3*d3)*s + (3*b3*d1 - b3*d3)*t;
    BTDB(11,25) = (3*c1*d3 - c3*d3)*s + (3*c3*d1 - c3*d3)*t;
    BTDB(11,26) = (3*d1*d3 - d3*d3)*r + 
      (3*b1*b3 - b3*b3 + 3*c1*c3 - c3*c3)*t;
    BTDB(11,27) = (3*b2*d3 - b3*d3)*s + (3*b3*d2 - b3*d3)*t;
    BTDB(11,28) = (3*c2*d3 - c3*d3)*s + (3*c3*d2 - c3*d3)*t;
    BTDB(11,29) = (3*d2*d3 - d3*d3)*r + 
      (3*b2*b3 - b3*b3 + 3*c2*c3 - c3*c3)*t;
    
    BTDB(12,12) = (8*b0b0 + 8*b0b1 + 8*b1*b1)*r + 
      (8*c0c0 + 8*c0c1 + 8*c1*c1 + 
       8*d0d0 + 8*d0d1 + 8*d1*d1)*t;
    BTDB(12,13) = (8*b0c0 + 4*b1*c0 + 4*b0c1 + 8*b1*c1)*s + 
      (8*b0c0 + 4*b1*c0 + 4*b0c1 + 8*b1*c1)*t;
    BTDB(12,14) = (8*b0d0 + 4*b1*d0 + 4*b0d1 + 8*b1*d1)*s + 
      (8*b0d0 + 4*b1*d0 + 4*b0d1 + 8*b1*d1)*t;
    BTDB(12,15) = (4*b0b0 + 4*b0b1 + 4*b0*b2 + 8*b1*b2)*r + 
      (4*c0c0 + 4*c0c1 + 4*c0*c2 + 8*c1*c2 + 
       4*d0d0 + 4*d0d1 + 4*d0*d2 + 8*d1*d2)*t;
    BTDB(12,16) = (4*b0c0 + 4*b1*c0 + 4*b0*c2 + 8*b1*c2)*s + 
      (4*b0c0 + 4*b2*c0 + 4*b0c1 + 8*b2*c1)*t;
    BTDB(12,17) = (4*b0d0 + 4*b1*d0 + 4*b0*d2 + 8*b1*d2)*s + 
      (4*b0d0 + 4*b2*d0 + 4*b0d1 + 8*b2*d1)*t;
    BTDB(12,18) = (4*b0b0 + 4*b0b1 + 4*b0*b3 + 8*b1*b3)*r + 
      (4*c0c0 + 4*c0c1 + 4*c0*c3 + 8*c1*c3 + 
       4*d0d0 + 4*d0d1 + 4*d0*d3 + 8*d1*d3)*t;
    BTDB(12,19) = (4*b0c0 + 4*b1*c0 + 4*b0*c3 + 8*b1*c3)*s + 
      (4*b0c0 + 4*b3*c0 + 4*b0c1 + 8*b3*c1)*t;
    BTDB(12,20) = (4*b0d0 + 4*b1*d0 + 4*b0*d3 + 8*b1*d3)*s + 
      (4*b0d0 + 4*b3*d0 + 4*b0d1 + 8*b3*d1)*t;
    BTDB(12,21) = (4*b0b1 + 4*b1*b1 + 8*b0*b2 + 4*b1*b2)*r + 
      (4*c0c1 + 4*c1*c1 + 8*c0*c2 + 4*c1*c2 + 
       4*d0d1 + 4*d1*d1 + 8*d0*d2 + 4*d1*d2)*t;
    BTDB(12,22) = (4*b0c1 + 4*b1*c1 + 8*b0*c2 + 4*b1*c2)*s + 
      (4*b1*c0 + 8*b2*c0 + 4*b1*c1 + 4*b2*c1)*t;
    BTDB(12,23) = (4*b0d1 + 4*b1*d1 + 8*b0*d2 + 4*b1*d2)*s + 
      (4*b1*d0 + 8*b2*d0 + 4*b1*d1 + 4*b2*d1)*t;
    BTDB(12,24) = (4*b0b1 + 4*b1*b1 + 8*b0*b3 + 4*b1*b3)*r + 
      (4*c0c1 + 4*c1*c1 + 8*c0*c3 + 4*c1*c3 + 
       4*d0d1 + 4*d1*d1 + 8*d0*d3 + 4*d1*d3)*t;
    BTDB(12,25) = (4*b0c1 + 4*b1*c1 + 8*b0*c3 + 4*b1*c3)*s + 
      (4*b1*c0 + 8*b3*c0 + 4*b1*c1 + 4*b3*c1)*t;
    BTDB(12,26) = (4*b0d1 + 4*b1*d1 + 8*b0*d3 + 4*b1*d3)*s + 
      (4*b1*d0 + 8*b3*d0 + 4*b1*d1 + 4*b3*d1)*t;
    BTDB(12,27) = (4*b0*b2 + 4*b1*b2 + 4*b0*b3 + 4*b1*b3)*r + 
      (4*c0*c2 + 4*c1*c2 + 4*c0*c3 + 4*c1*c3 + 4*d0*d2 + 
       4*d1*d2 + 4*d0*d3 + 4*d1*d3)*t;
    BTDB(12,28) = (4*b0*c2 + 4*b1*c2 + 4*b0*c3 + 4*b1*c3)*s + 
      (4*b2*c0 + 4*b3*c0 + 4*b2*c1 + 4*b3*c1)*t;
    BTDB(12,29) = (4*b0*d2 + 4*b1*d2 + 4*b0*d3 + 4*b1*d3)*s + 
      (4*b2*d0 + 4*b3*d0 + 4*b2*d1 + 4*b3*d1)*t;
    
    BTDB(13,13) = (8*c0c0 + 8*c0c1 + 8*c1*c1)*r + 
      (8*b0b0 + 8*b0b1 + 8*b1*b1 + 
       8*d0d0 + 8*d0d1 + 8*d1*d1)*t;
    BTDB(13,14) = (8*c0d0 + 4*c1*d0 + 4*c0d1 + 8*c1*d1)*s + 
      (8*c0d0 + 4*c1*d0 + 4*c0d1 + 8*c1*d1)*t;
    BTDB(13,15) = (4*b0c0 + 4*b2*c0 + 4*b0c1 + 8*b2*c1)*s + 
      (4*b0c0 + 4*b1*c0 + 4*b0*c2 + 8*b1*c2)*t;
    BTDB(13,16) = (4*c0c0 + 4*c0c1 + 4*c0*c2 + 8*c1*c2)*r + 
      (4*b0b0 + 4*b0b1 + 4*b0*b2 + 8*b1*b2 + 
       4*d0d0 + 4*d0d1 + 4*d0*d2 + 8*d1*d2)*t;
    BTDB(13,17) = (4*c0d0 + 4*c1*d0 + 4*c0*d2 + 8*c1*d2)*s + 
      (4*c0d0 + 4*c2*d0 + 4*c0d1 + 8*c2*d1)*t;
    BTDB(13,18) = (4*b0c0 + 4*b3*c0 + 4*b0c1 + 8*b3*c1)*s + 
      (4*b0c0 + 4*b1*c0 + 4*b0*c3 + 8*b1*c3)*t;
    BTDB(13,19) = (4*c0c0 + 4*c0c1 + 4*c0*c3 + 8*c1*c3)*r + 
      (4*b0b0 + 4*b0b1 + 4*b0*b3 + 8*b1*b3 + 
       4*d0d0 + 4*d0d1 + 4*d0*d3 + 8*d1*d3)*t;
    BTDB(13,20) = (4*c0d0 + 4*c1*d0 + 4*c0*d3 + 8*c1*d3)*s + 
      (4*c0d0 + 4*c3*d0 + 4*c0d1 + 8*c3*d1)*t;
    BTDB(13,21) = (4*b1*c0 + 8*b2*c0 + 4*b1*c1 + 4*b2*c1)*s + 
      (4*b0c1 + 4*b1*c1 + 8*b0*c2 + 4*b1*c2)*t;
    BTDB(13,22) = (4*c0c1 + 4*c1*c1 + 8*c0*c2 + 4*c1*c2)*r + 
      (4*b0b1 + 4*b1*b1 + 8*b0*b2 + 4*b1*b2 + 
       4*d0d1 + 4*d1*d1 + 8*d0*d2 + 4*d1*d2)*t;
    BTDB(13,23) = (4*c0d1 + 4*c1*d1 + 8*c0*d2 + 4*c1*d2)*s + 
      (4*c1*d0 + 8*c2*d0 + 4*c1*d1 + 4*c2*d1)*t;
    BTDB(13,24) = (4*b1*c0 + 8*b3*c0 + 4*b1*c1 + 4*b3*c1)*s + 
      (4*b0c1 + 4*b1*c1 + 8*b0*c3 + 4*b1*c3)*t;
    BTDB(13,25) = (4*c0c1 + 4*c1*c1 + 8*c0*c3 + 4*c1*c3)*r + 
      (4*b0b1 + 4*b1*b1 + 8*b0*b3 + 4*b1*b3 + 
       4*d0d1 + 4*d1*d1 + 8*d0*d3 + 4*d1*d3)*t;
    BTDB(13,26) = (4*c0d1 + 4*c1*d1 + 8*c0*d3 + 4*c1*d3)*s + 
      (4*c1*d0 + 8*c3*d0 + 4*c1*d1 + 4*c3*d1)*t;
    BTDB(13,27) = (4*b2*c0 + 4*b3*c0 + 4*b2*c1 + 4*b3*c1)*s + 
      (4*b0*c2 + 4*b1*c2 + 4*b0*c3 + 4*b1*c3)*t;
    BTDB(13,28) = (4*c0*c2 + 4*c1*c2 + 4*c0*c3 + 4*c1*c3)*r + 
      (4*b0*b2 + 4*b1*b2 + 4*b0*b3 + 4*b1*b3 + 4*d0*d2 + 
       4*d1*d2 + 4*d0*d3 + 4*d1*d3)*t;
    BTDB(13,29) = (4*c0*d2 + 4*c1*d2 + 4*c0*d3 + 4*c1*d3)*s + 
      (4*c2*d0 + 4*c3*d0 + 4*c2*d1 + 4*c3*d1)*t;
    
    BTDB(14,14) = (8*d0d0 + 8*d0d1 + 8*d1*d1)*r + 
      (8*b0b0 + 8*b0b1 + 8*b1*b1 + 
       8*c0c0 + 8*c0c1 + 8*c1*c1)*t;
    BTDB(14,15) = (4*b0d0 + 4*b2*d0 + 4*b0d1 + 8*b2*d1)*s + 
      (4*b0d0 + 4*b1*d0 + 4*b0*d2 + 8*b1*d2)*t;
    BTDB(14,16) = (4*c0d0 + 4*c2*d0 + 4*c0d1 + 8*c2*d1)*s + 
      (4*c0d0 + 4*c1*d0 + 4*c0*d2 + 8*c1*d2)*t;
    BTDB(14,17) = (4*d0d0 + 4*d0d1 + 4*d0*d2 + 8*d1*d2)*r + 
      (4*b0b0 + 4*b0b1 + 4*b0*b2 + 8*b1*b2 + 
       4*c0c0 + 4*c0c1 + 4*c0*c2 + 8*c1*c2)*t;
    BTDB(14,18) = (4*b0d0 + 4*b3*d0 + 4*b0d1 + 8*b3*d1)*s + 
      (4*b0d0 + 4*b1*d0 + 4*b0*d3 + 8*b1*d3)*t;
    BTDB(14,19) = (4*c0d0 + 4*c3*d0 + 4*c0d1 + 8*c3*d1)*s + 
      (4*c0d0 + 4*c1*d0 + 4*c0*d3 + 8*c1*d3)*t;
    BTDB(14,20) = (4*d0d0 + 4*d0d1 + 4*d0*d3 + 8*d1*d3)*r + 
      (4*b0b0 + 4*b0b1 + 4*b0*b3 + 8*b1*b3 + 
       4*c0c0 + 4*c0c1 + 4*c0*c3 + 8*c1*c3)*t;
    BTDB(14,21) = (4*b1*d0 + 8*b2*d0 + 4*b1*d1 + 4*b2*d1)*s + 
      (4*b0d1 + 4*b1*d1 + 8*b0*d2 + 4*b1*d2)*t;
    BTDB(14,22) = (4*c1*d0 + 8*c2*d0 + 4*c1*d1 + 4*c2*d1)*s + 
      (4*c0d1 + 4*c1*d1 + 8*c0*d2 + 4*c1*d2)*t;
    BTDB(14,23) = (4*d0d1 + 4*d1*d1 + 8*d0*d2 + 4*d1*d2)*r + 
      (4*b0b1 + 4*b1*b1 + 8*b0*b2 + 4*b1*b2 + 
       4*c0c1 + 4*c1*c1 + 8*c0*c2 + 4*c1*c2)*t;
    BTDB(14,24) = (4*b1*d0 + 8*b3*d0 + 4*b1*d1 + 4*b3*d1)*s + 
      (4*b0d1 + 4*b1*d1 + 8*b0*d3 + 4*b1*d3)*t;
    BTDB(14,25) = (4*c1*d0 + 8*c3*d0 + 4*c1*d1 + 4*c3*d1)*s + 
      (4*c0d1 + 4*c1*d1 + 8*c0*d3 + 4*c1*d3)*t;
    BTDB(14,26) = (4*d0d1 + 4*d1*d1 + 8*d0*d3 + 4*d1*d3)*r + 
      (4*b0b1 + 4*b1*b1 + 8*b0*b3 + 4*b1*b3 + 
       4*c0c1 + 4*c1*c1 + 8*c0*c3 + 4*c1*c3)*t;
    BTDB(14,27) = (4*b2*d0 + 4*b3*d0 + 4*b2*d1 + 4*b3*d1)*s + 
      (4*b0*d2 + 4*b1*d2 + 4*b0*d3 + 4*b1*d3)*t;
    BTDB(14,28) = (4*c2*d0 + 4*c3*d0 + 4*c2*d1 + 4*c3*d1)*s + 
      (4*c0*d2 + 4*c1*d2 + 4*c0*d3 + 4*c1*d3)*t;
    BTDB(14,29) = (4*d0*d2 + 4*d1*d2 + 4*d0*d3 + 4*d1*d3)*r + 
      (4*b0*b2 + 4*b1*b2 + 4*b0*b3 + 4*b1*b3 + 4*c0*c2 + 
       4*c1*c2 + 4*c0*c3 + 4*c1*c3)*t;
    
    BTDB(15,15) = (8*b0b0 + 8*b0*b2 + 8*b2*b2)*r + 
      (8*c0c0 + 8*c0*c2 + 8*c2*c2 + 
       8*d0d0 + 8*d0*d2 + 8*d2*d2)*t;
    BTDB(15,16) = (8*b0c0 + 4*b2*c0 + 4*b0*c2 + 8*b2*c2)*s + 
      (8*b0c0 + 4*b2*c0 + 4*b0*c2 + 8*b2*c2)*t;
    BTDB(15,17) = (8*b0d0 + 4*b2*d0 + 4*b0*d2 + 8*b2*d2)*s + 
      (8*b0d0 + 4*b2*d0 + 4*b0*d2 + 8*b2*d2)*t;
    BTDB(15,18) = (4*b0b0 + 4*b0*b2 + 4*b0*b3 + 8*b2*b3)*r + 
      (4*c0c0 + 4*c0*c2 + 4*c0*c3 + 8*c2*c3 + 
       4*d0d0 + 4*d0*d2 + 4*d0*d3 + 8*d2*d3)*t;
    BTDB(15,19) = (4*b0c0 + 4*b2*c0 + 4*b0*c3 + 8*b2*c3)*s + 
      (4*b0c0 + 4*b3*c0 + 4*b0*c2 + 8*b3*c2)*t;
    BTDB(15,20) = (4*b0d0 + 4*b2*d0 + 4*b0*d3 + 8*b2*d3)*s + 
      (4*b0d0 + 4*b3*d0 + 4*b0*d2 + 8*b3*d2)*t;
    BTDB(15,21) = (8*b0b1 + 4*b0*b2 + 4*b1*b2 + 4*b2*b2)*r + 
      (8*c0c1 + 4*c0*c2 + 4*c1*c2 + 4*c2*c2 + 
       8*d0d1 + 4*d0*d2 + 4*d1*d2 + 4*d2*d2)*t;
    BTDB(15,22) = (8*b0c1 + 4*b2*c1 + 4*b0*c2 + 4*b2*c2)*s + 
      (8*b1*c0 + 4*b2*c0 + 4*b1*c2 + 4*b2*c2)*t;
    BTDB(15,23) = (8*b0d1 + 4*b2*d1 + 4*b0*d2 + 4*b2*d2)*s + 
      (8*b1*d0 + 4*b2*d0 + 4*b1*d2 + 4*b2*d2)*t;
    BTDB(15,24) = (4*b0b1 + 4*b1*b2 + 4*b0*b3 + 4*b2*b3)*r + 
      (4*c0c1 + 4*c1*c2 + 4*c0*c3 + 4*c2*c3 + 4*d0d1 + 
       4*d1*d2 + 4*d0*d3 + 4*d2*d3)*t;
    BTDB(15,25) = (4*b0c1 + 4*b2*c1 + 4*b0*c3 + 4*b2*c3)*s + 
      (4*b1*c0 + 4*b3*c0 + 4*b1*c2 + 4*b3*c2)*t;
    BTDB(15,26) = (4*b0d1 + 4*b2*d1 + 4*b0*d3 + 4*b2*d3)*s + 
      (4*b1*d0 + 4*b3*d0 + 4*b1*d2 + 4*b3*d2)*t;
    BTDB(15,27) = (4*b0*b2 + 4*b2*b2 + 8*b0*b3 + 4*b2*b3)*r + 
      (4*c0*c2 + 4*c2*c2 + 8*c0*c3 + 4*c2*c3 + 
       4*d0*d2 + 4*d2*d2 + 8*d0*d3 + 4*d2*d3)*t;
    BTDB(15,28) = (4*b0*c2 + 4*b2*c2 + 8*b0*c3 + 4*b2*c3)*s + 
      (4*b2*c0 + 8*b3*c0 + 4*b2*c2 + 4*b3*c2)*t;
    BTDB(15,29) = (4*b0*d2 + 4*b2*d2 + 8*b0*d3 + 4*b2*d3)*s + 
      (4*b2*d0 + 8*b3*d0 + 4*b2*d2 + 4*b3*d2)*t;
    
    BTDB(16,16) = (8*c0c0 + 8*c0*c2 + 8*c2*c2)*r + 
      (8*b0b0 + 8*b0*b2 + 8*b2*b2 + 
       8*d0d0 + 8*d0*d2 + 8*d2*d2)*t;
    BTDB(16,17) = (8*c0d0 + 4*c2*d0 + 4*c0*d2 + 8*c2*d2)*s + 
      (8*c0d0 + 4*c2*d0 + 4*c0*d2 + 8*c2*d2)*t;
    BTDB(16,18) = (4*b0c0 + 4*b3*c0 + 4*b0*c2 + 8*b3*c2)*s + 
      (4*b0c0 + 4*b2*c0 + 4*b0*c3 + 8*b2*c3)*t;
    BTDB(16,19) = (4*c0c0 + 4*c0*c2 + 4*c0*c3 + 8*c2*c3)*r + 
      (4*b0b0 + 4*b0*b2 + 4*b0*b3 + 8*b2*b3 + 
       4*d0d0 + 4*d0*d2 + 4*d0*d3 + 8*d2*d3)*t;
    BTDB(16,20) = (4*c0d0 + 4*c2*d0 + 4*c0*d3 + 8*c2*d3)*s + 
      (4*c0d0 + 4*c3*d0 + 4*c0*d2 + 8*c3*d2)*t;
    BTDB(16,21) = (8*b1*c0 + 4*b2*c0 + 4*b1*c2 + 4*b2*c2)*s + 
      (8*b0c1 + 4*b2*c1 + 4*b0*c2 + 4*b2*c2)*t;
    BTDB(16,22) = (8*c0c1 + 4*c0*c2 + 4*c1*c2 + 4*c2*c2)*r + 
      (8*b0b1 + 4*b0*b2 + 4*b1*b2 + 4*b2*b2 + 
       8*d0d1 + 4*d0*d2 + 4*d1*d2 + 4*d2*d2)*t;
    BTDB(16,23) = (8*c0d1 + 4*c2*d1 + 4*c0*d2 + 4*c2*d2)*s + 
      (8*c1*d0 + 4*c2*d0 + 4*c1*d2 + 4*c2*d2)*t;
    BTDB(16,24) = (4*b1*c0 + 4*b3*c0 + 4*b1*c2 + 4*b3*c2)*s + 
      (4*b0c1 + 4*b2*c1 + 4*b0*c3 + 4*b2*c3)*t;
    BTDB(16,25) = (4*c0c1 + 4*c1*c2 + 4*c0*c3 + 4*c2*c3)*r + 
      (4*b0b1 + 4*b1*b2 + 4*b0*b3 + 4*b2*b3 + 4*d0d1 + 
       4*d1*d2 + 4*d0*d3 + 4*d2*d3)*t;
    BTDB(16,26) = (4*c0d1 + 4*c2*d1 + 4*c0*d3 + 4*c2*d3)*s + 
      (4*c1*d0 + 4*c3*d0 + 4*c1*d2 + 4*c3*d2)*t;
    BTDB(16,27) = (4*b2*c0 + 8*b3*c0 + 4*b2*c2 + 4*b3*c2)*s + 
      (4*b0*c2 + 4*b2*c2 + 8*b0*c3 + 4*b2*c3)*t;
    BTDB(16,28) = (4*c0*c2 + 4*c2*c2 + 8*c0*c3 + 4*c2*c3)*r + 
      (4*b0*b2 + 4*b2*b2 + 8*b0*b3 + 4*b2*b3 + 
       4*d0*d2 + 4*d2*d2 + 8*d0*d3 + 4*d2*d3)*t;
    BTDB(16,29) = (4*c0*d2 + 4*c2*d2 + 8*c0*d3 + 4*c2*d3)*s + 
      (4*c2*d0 + 8*c3*d0 + 4*c2*d2 + 4*c3*d2)*t;
    
    BTDB(17,17) = (8*d0d0 + 8*d0*d2 + 8*d2*d2)*r + 
      (8*b0b0 + 8*b0*b2 + 8*b2*b2 + 
       8*c0c0 + 8*c0*c2 + 8*c2*c2)*t;
    BTDB(17,18) = (4*b0d0 + 4*b3*d0 + 4*b0*d2 + 8*b3*d2)*s + 
      (4*b0d0 + 4*b2*d0 + 4*b0*d3 + 8*b2*d3)*t;
    BTDB(17,19) = (4*c0d0 + 4*c3*d0 + 4*c0*d2 + 8*c3*d2)*s + 
      (4*c0d0 + 4*c2*d0 + 4*c0*d3 + 8*c2*d3)*t;
    BTDB(17,20) = (4*d0d0 + 4*d0*d2 + 4*d0*d3 + 8*d2*d3)*r + 
      (4*b0b0 + 4*b0*b2 + 4*b0*b3 + 8*b2*b3 + 
       4*c0c0 + 4*c0*c2 + 4*c0*c3 + 8*c2*c3)*t;
    BTDB(17,21) = (8*b1*d0 + 4*b2*d0 + 4*b1*d2 + 4*b2*d2)*s + 
      (8*b0d1 + 4*b2*d1 + 4*b0*d2 + 4*b2*d2)*t;
    BTDB(17,22) = (8*c1*d0 + 4*c2*d0 + 4*c1*d2 + 4*c2*d2)*s + 
      (8*c0d1 + 4*c2*d1 + 4*c0*d2 + 4*c2*d2)*t;
    BTDB(17,23) = (8*d0d1 + 4*d0*d2 + 4*d1*d2 + 4*d2*d2)*r + 
      (8*b0b1 + 4*b0*b2 + 4*b1*b2 + 4*b2*b2 + 
       8*c0c1 + 4*c0*c2 + 4*c1*c2 + 4*c2*c2)*t;
    BTDB(17,24) = (4*b1*d0 + 4*b3*d0 + 4*b1*d2 + 4*b3*d2)*s + 
      (4*b0d1 + 4*b2*d1 + 4*b0*d3 + 4*b2*d3)*t;
    BTDB(17,25) = (4*c1*d0 + 4*c3*d0 + 4*c1*d2 + 4*c3*d2)*s + 
      (4*c0d1 + 4*c2*d1 + 4*c0*d3 + 4*c2*d3)*t;
    BTDB(17,26) = (4*d0d1 + 4*d1*d2 + 4*d0*d3 + 4*d2*d3)*r + 
      (4*b0b1 + 4*b1*b2 + 4*b0*b3 + 4*b2*b3 + 4*c0c1 + 
       4*c1*c2 + 4*c0*c3 + 4*c2*c3)*t;
    BTDB(17,27) = (4*b2*d0 + 8*b3*d0 + 4*b2*d2 + 4*b3*d2)*s + 
      (4*b0*d2 + 4*b2*d2 + 8*b0*d3 + 4*b2*d3)*t;
    BTDB(17,28) = (4*c2*d0 + 8*c3*d0 + 4*c2*d2 + 4*c3*d2)*s + 
      (4*c0*d2 + 4*c2*d2 + 8*c0*d3 + 4*c2*d3)*t;
    BTDB(17,29) = (4*d0*d2 + 4*d2*d2 + 8*d0*d3 + 4*d2*d3)*r + 
      (4*b0*b2 + 4*b2*b2 + 8*b0*b3 + 4*b2*b3 + 
       4*c0*c2 + 4*c2*c2 + 8*c0*c3 + 4*c2*c3)*t;
    
    BTDB(18,18) = (8*b0b0 + 8*b0*b3 + 8*b3*b3)*r + 
      (8*c0c0 + 8*c0*c3 + 8*c3*c3 + 
       8*d0d0 + 8*d0*d3 + 8*d3*d3)*t;
    BTDB(18,19) = (8*b0c0 + 4*b3*c0 + 4*b0*c3 + 8*b3*c3)*s + 
      (8*b0c0 + 4*b3*c0 + 4*b0*c3 + 8*b3*c3)*t;
    BTDB(18,20) = (8*b0d0 + 4*b3*d0 + 4*b0*d3 + 8*b3*d3)*s + 
      (8*b0d0 + 4*b3*d0 + 4*b0*d3 + 8*b3*d3)*t;
    BTDB(18,21) = (4*b0b1 + 4*b0*b2 + 4*b1*b3 + 4*b2*b3)*r + 
      (4*c0c1 + 4*c0*c2 + 4*c1*c3 + 4*c2*c3 + 4*d0d1 + 
       4*d0*d2 + 4*d1*d3 + 4*d2*d3)*t;
    BTDB(18,22) = (4*b0c1 + 4*b3*c1 + 4*b0*c2 + 4*b3*c2)*s + 
      (4*b1*c0 + 4*b2*c0 + 4*b1*c3 + 4*b2*c3)*t;
    BTDB(18,23) = (4*b0d1 + 4*b3*d1 + 4*b0*d2 + 4*b3*d2)*s + 
      (4*b1*d0 + 4*b2*d0 + 4*b1*d3 + 4*b2*d3)*t;
    BTDB(18,24) = (8*b0b1 + 4*b0*b3 + 4*b1*b3 + 4*b3*b3)*r + 
      (8*c0c1 + 4*c0*c3 + 4*c1*c3 + 4*c3*c3 + 
       8*d0d1 + 4*d0*d3 + 4*d1*d3 + 4*d3*d3)*t;
    BTDB(18,25) = (8*b0c1 + 4*b3*c1 + 4*b0*c3 + 4*b3*c3)*s + 
      (8*b1*c0 + 4*b3*c0 + 4*b1*c3 + 4*b3*c3)*t;
    BTDB(18,26) = (8*b0d1 + 4*b3*d1 + 4*b0*d3 + 4*b3*d3)*s + 
      (8*b1*d0 + 4*b3*d0 + 4*b1*d3 + 4*b3*d3)*t;
    BTDB(18,27) = (8*b0*b2 + 4*b0*b3 + 4*b2*b3 + 4*b3*b3)*r + 
      (8*c0*c2 + 4*c0*c3 + 4*c2*c3 + 4*c3*c3 + 
       8*d0*d2 + 4*d0*d3 + 4*d2*d3 + 4*d3*d3)*t;
    BTDB(18,28) = (8*b0*c2 + 4*b3*c2 + 4*b0*c3 + 4*b3*c3)*s + 
      (8*b2*c0 + 4*b3*c0 + 4*b2*c3 + 4*b3*c3)*t;
    BTDB(18,29) = (8*b0*d2 + 4*b3*d2 + 4*b0*d3 + 4*b3*d3)*s + 
      (8*b2*d0 + 4*b3*d0 + 4*b2*d3 + 4*b3*d3)*t;
    
    BTDB(19,19) = (8*c0c0 + 8*c0*c3 + 8*c3*c3)*r + 
      (8*b0b0 + 8*b0*b3 + 8*b3*b3 + 
       8*d0d0 + 8*d0*d3 + 8*d3*d3)*t;
    BTDB(19,20) = (8*c0d0 + 4*c3*d0 + 4*c0*d3 + 8*c3*d3)*s + 
      (8*c0d0 + 4*c3*d0 + 4*c0*d3 + 8*c3*d3)*t;
    BTDB(19,21) = (4*b1*c0 + 4*b2*c0 + 4*b1*c3 + 4*b2*c3)*s + 
      (4*b0c1 + 4*b3*c1 + 4*b0*c2 + 4*b3*c2)*t;
    BTDB(19,22) = (4*c0c1 + 4*c0*c2 + 4*c1*c3 + 4*c2*c3)*r + 
      (4*b0b1 + 4*b0*b2 + 4*b1*b3 + 4*b2*b3 + 4*d0d1 + 
       4*d0*d2 + 4*d1*d3 + 4*d2*d3)*t;
    BTDB(19,23) = (4*c0d1 + 4*c3*d1 + 4*c0*d2 + 4*c3*d2)*s + 
      (4*c1*d0 + 4*c2*d0 + 4*c1*d3 + 4*c2*d3)*t;
    BTDB(19,24) = (8*b1*c0 + 4*b3*c0 + 4*b1*c3 + 4*b3*c3)*s + 
      (8*b0c1 + 4*b3*c1 + 4*b0*c3 + 4*b3*c3)*t;
    BTDB(19,25) = (8*c0c1 + 4*c0*c3 + 4*c1*c3 + 4*c3*c3)*r + 
      (8*b0b1 + 4*b0*b3 + 4*b1*b3 + 4*b3*b3 + 
       8*d0d1 + 4*d0*d3 + 4*d1*d3 + 4*d3*d3)*t;
    BTDB(19,26) = (8*c0d1 + 4*c3*d1 + 4*c0*d3 + 4*c3*d3)*s + 
      (8*c1*d0 + 4*c3*d0 + 4*c1*d3 + 4*c3*d3)*t;
    BTDB(19,27) = (8*b2*c0 + 4*b3*c0 + 4*b2*c3 + 4*b3*c3)*s + 
      (8*b0*c2 + 4*b3*c2 + 4*b0*c3 + 4*b3*c3)*t;
    BTDB(19,28) = (8*c0*c2 + 4*c0*c3 + 4*c2*c3 + 4*c3*c3)*r + 
      (8*b0*b2 + 4*b0*b3 + 4*b2*b3 + 4*b3*b3 + 
       8*d0*d2 + 4*d0*d3 + 4*d2*d3 + 4*d3*d3)*t;
    BTDB(19,29) = (8*c0*d2 + 4*c3*d2 + 4*c0*d3 + 4*c3*d3)*s + 
      (8*c2*d0 + 4*c3*d0 + 4*c2*d3 + 4*c3*d3)*t;
    
    BTDB(20,20) = (8*d0d0 + 8*d0*d3 + 8*d3*d3)*r + 
      (8*b0b0 + 8*b0*b3 + 8*b3*b3 + 
       8*c0c0 + 8*c0*c3 + 8*c3*c3)*t;
    BTDB(20,21) = (4*b1*d0 + 4*b2*d0 + 4*b1*d3 + 4*b2*d3)*s + 
      (4*b0d1 + 4*b3*d1 + 4*b0*d2 + 4*b3*d2)*t;
    BTDB(20,22) = (4*c1*d0 + 4*c2*d0 + 4*c1*d3 + 4*c2*d3)*s + 
      (4*c0d1 + 4*c3*d1 + 4*c0*d2 + 4*c3*d2)*t;
    BTDB(20,23) = (4*d0d1 + 4*d0*d2 + 4*d1*d3 + 4*d2*d3)*r + 
      (4*b0b1 + 4*b0*b2 + 4*b1*b3 + 4*b2*b3 + 4*c0c1 + 
       4*c0*c2 + 4*c1*c3 + 4*c2*c3)*t;
    BTDB(20,24) = (8*b1*d0 + 4*b3*d0 + 4*b1*d3 + 4*b3*d3)*s + 
      (8*b0d1 + 4*b3*d1 + 4*b0*d3 + 4*b3*d3)*t;
    BTDB(20,25) = (8*c1*d0 + 4*c3*d0 + 4*c1*d3 + 4*c3*d3)*s + 
      (8*c0d1 + 4*c3*d1 + 4*c0*d3 + 4*c3*d3)*t;
    BTDB(20,26) = (8*d0d1 + 4*d0*d3 + 4*d1*d3 + 4*d3*d3)*r + 
      (8*b0b1 + 4*b0*b3 + 4*b1*b3 + 4*b3*b3 + 
       8*c0c1 + 4*c0*c3 + 4*c1*c3 + 4*c3*c3)*t;
    BTDB(20,27) = (8*b2*d0 + 4*b3*d0 + 4*b2*d3 + 4*b3*d3)*s + 
      (8*b0*d2 + 4*b3*d2 + 4*b0*d3 + 4*b3*d3)*t;
    BTDB(20,28) = (8*c2*d0 + 4*c3*d0 + 4*c2*d3 + 4*c3*d3)*s + 
      (8*c0*d2 + 4*c3*d2 + 4*c0*d3 + 4*c3*d3)*t;
    BTDB(20,29) = (8*d0*d2 + 4*d0*d3 + 4*d2*d3 + 4*d3*d3)*r + 
      (8*b0*b2 + 4*b0*b3 + 4*b2*b3 + 4*b3*b3 + 
       8*c0*c2 + 4*c0*c3 + 4*c2*c3 + 4*c3*c3)*t;
    
    BTDB(21,21) = (8*b1*b1 + 8*b1*b2 + 8*b2*b2)*r + 
      (8*c1*c1 + 8*c1*c2 + 8*c2*c2 + 
       8*d1*d1 + 8*d1*d2 + 8*d2*d2)*t;
    BTDB(21,22) = (8*b1*c1 + 4*b2*c1 + 4*b1*c2 + 8*b2*c2)*s + 
      (8*b1*c1 + 4*b2*c1 + 4*b1*c2 + 8*b2*c2)*t;
    BTDB(21,23) = (8*b1*d1 + 4*b2*d1 + 4*b1*d2 + 8*b2*d2)*s + 
      (8*b1*d1 + 4*b2*d1 + 4*b1*d2 + 8*b2*d2)*t;
    BTDB(21,24) = (4*b1*b1 + 4*b1*b2 + 4*b1*b3 + 8*b2*b3)*r + 
      (4*c1*c1 + 4*c1*c2 + 4*c1*c3 + 8*c2*c3 + 
       4*d1*d1 + 4*d1*d2 + 4*d1*d3 + 8*d2*d3)*t;
    BTDB(21,25) = (4*b1*c1 + 4*b2*c1 + 4*b1*c3 + 8*b2*c3)*s + 
      (4*b1*c1 + 4*b3*c1 + 4*b1*c2 + 8*b3*c2)*t;
    BTDB(21,26) = (4*b1*d1 + 4*b2*d1 + 4*b1*d3 + 8*b2*d3)*s + 
      (4*b1*d1 + 4*b3*d1 + 4*b1*d2 + 8*b3*d2)*t;
    BTDB(21,27) = (4*b1*b2 + 4*b2*b2 + 8*b1*b3 + 4*b2*b3)*r + 
      (4*c1*c2 + 4*c2*c2 + 8*c1*c3 + 4*c2*c3 + 
       4*d1*d2 + 4*d2*d2 + 8*d1*d3 + 4*d2*d3)*t;
    BTDB(21,28) = (4*b1*c2 + 4*b2*c2 + 8*b1*c3 + 4*b2*c3)*s + 
      (4*b2*c1 + 8*b3*c1 + 4*b2*c2 + 4*b3*c2)*t;
    BTDB(21,29) = (4*b1*d2 + 4*b2*d2 + 8*b1*d3 + 4*b2*d3)*s + 
      (4*b2*d1 + 8*b3*d1 + 4*b2*d2 + 4*b3*d2)*t;
    
    BTDB(22,22) = (8*c1*c1 + 8*c1*c2 + 8*c2*c2)*r + 
      (8*b1*b1 + 8*b1*b2 + 8*b2*b2 + 
       8*d1*d1 + 8*d1*d2 + 8*d2*d2)*t;
    BTDB(22,23) = (8*c1*d1 + 4*c2*d1 + 4*c1*d2 + 8*c2*d2)*s + 
      (8*c1*d1 + 4*c2*d1 + 4*c1*d2 + 8*c2*d2)*t;
    BTDB(22,24) = (4*b1*c1 + 4*b3*c1 + 4*b1*c2 + 8*b3*c2)*s + 
      (4*b1*c1 + 4*b2*c1 + 4*b1*c3 + 8*b2*c3)*t;
    BTDB(22,25) = (4*c1*c1 + 4*c1*c2 + 4*c1*c3 + 8*c2*c3)*r + 
      (4*b1*b1 + 4*b1*b2 + 4*b1*b3 + 8*b2*b3 + 
       4*d1*d1 + 4*d1*d2 + 4*d1*d3 + 8*d2*d3)*t;
    BTDB(22,26) = (4*c1*d1 + 4*c2*d1 + 4*c1*d3 + 8*c2*d3)*s + 
      (4*c1*d1 + 4*c3*d1 + 4*c1*d2 + 8*c3*d2)*t;
    BTDB(22,27) = (4*b2*c1 + 8*b3*c1 + 4*b2*c2 + 4*b3*c2)*s + 
      (4*b1*c2 + 4*b2*c2 + 8*b1*c3 + 4*b2*c3)*t;
    BTDB(22,28) = (4*c1*c2 + 4*c2*c2 + 8*c1*c3 + 4*c2*c3)*r + 
      (4*b1*b2 + 4*b2*b2 + 8*b1*b3 + 4*b2*b3 + 
       4*d1*d2 + 4*d2*d2 + 8*d1*d3 + 4*d2*d3)*t;
    BTDB(22,29) = (4*c1*d2 + 4*c2*d2 + 8*c1*d3 + 4*c2*d3)*s + 
      (4*c2*d1 + 8*c3*d1 + 4*c2*d2 + 4*c3*d2)*t;
    
    BTDB(23,23) = (8*d1*d1 + 8*d1*d2 + 8*d2*d2)*r + 
      (8*b1*b1 + 8*b1*b2 + 8*b2*b2 + 
       8*c1*c1 + 8*c1*c2 + 8*c2*c2)*t;
    BTDB(23,24) = (4*b1*d1 + 4*b3*d1 + 4*b1*d2 + 8*b3*d2)*s + 
      (4*b1*d1 + 4*b2*d1 + 4*b1*d3 + 8*b2*d3)*t;
    BTDB(23,25) = (4*c1*d1 + 4*c3*d1 + 4*c1*d2 + 8*c3*d2)*s + 
      (4*c1*d1 + 4*c2*d1 + 4*c1*d3 + 8*c2*d3)*t;
    BTDB(23,26) = (4*d1*d1 + 4*d1*d2 + 4*d1*d3 + 8*d2*d3)*r + 
      (4*b1*b1 + 4*b1*b2 + 4*b1*b3 + 8*b2*b3 + 
       4*c1*c1 + 4*c1*c2 + 4*c1*c3 + 8*c2*c3)*t;
    BTDB(23,27) = (4*b2*d1 + 8*b3*d1 + 4*b2*d2 + 4*b3*d2)*s + 
      (4*b1*d2 + 4*b2*d2 + 8*b1*d3 + 4*b2*d3)*t;
    BTDB(23,28) = (4*c2*d1 + 8*c3*d1 + 4*c2*d2 + 4*c3*d2)*s + 
      (4*c1*d2 + 4*c2*d2 + 8*c1*d3 + 4*c2*d3)*t;
    BTDB(23,29) = (4*d1*d2 + 4*d2*d2 + 8*d1*d3 + 4*d2*d3)*r + 
      (4*b1*b2 + 4*b2*b2 + 8*b1*b3 + 4*b2*b3 + 
       4*c1*c2 + 4*c2*c2 + 8*c1*c3 + 4*c2*c3)*t;
    
    BTDB(24,24) = (8*b1*b1 + 8*b1*b3 + 8*b3*b3)*r + 
      (8*c1*c1 + 8*c1*c3 + 8*c3*c3 + 
       8*d1*d1 + 8*d1*d3 + 8*d3*d3)*t;
    BTDB(24,25) = (8*b1*c1 + 4*b3*c1 + 4*b1*c3 + 8*b3*c3)*s + 
      (8*b1*c1 + 4*b3*c1 + 4*b1*c3 + 8*b3*c3)*t;
    BTDB(24,26) = (8*b1*d1 + 4*b3*d1 + 4*b1*d3 + 8*b3*d3)*s + 
      (8*b1*d1 + 4*b3*d1 + 4*b1*d3 + 8*b3*d3)*t;
    BTDB(24,27) = (8*b1*b2 + 4*b1*b3 + 4*b2*b3 + 4*b3*b3)*r + 
      (8*c1*c2 + 4*c1*c3 + 4*c2*c3 + 4*c3*c3 + 
       8*d1*d2 + 4*d1*d3 + 4*d2*d3 + 4*d3*d3)*t;
    BTDB(24,28) = (8*b1*c2 + 4*b3*c2 + 4*b1*c3 + 4*b3*c3)*s + 
      (8*b2*c1 + 4*b3*c1 + 4*b2*c3 + 4*b3*c3)*t;
    BTDB(24,29) = (8*b1*d2 + 4*b3*d2 + 4*b1*d3 + 4*b3*d3)*s + 
      (8*b2*d1 + 4*b3*d1 + 4*b2*d3 + 4*b3*d3)*t;
    
    BTDB(25,25) = (8*c1*c1 + 8*c1*c3 + 8*c3*c3)*r + 
      (8*b1*b1 + 8*b1*b3 + 8*b3*b3 + 
       8*d1*d1 + 8*d1*d3 + 8*d3*d3)*t;
    BTDB(25,26) = (8*c1*d1 + 4*c3*d1 + 4*c1*d3 + 8*c3*d3)*s + 
      (8*c1*d1 + 4*c3*d1 + 4*c1*d3 + 8*c3*d3)*t;
    BTDB(25,27) = (8*b2*c1 + 4*b3*c1 + 4*b2*c3 + 4*b3*c3)*s + 
      (8*b1*c2 + 4*b3*c2 + 4*b1*c3 + 4*b3*c3)*t;
    BTDB(25,28) = (8*c1*c2 + 4*c1*c3 + 4*c2*c3 + 4*c3*c3)*r + 
      (8*b1*b2 + 4*b1*b3 + 4*b2*b3 + 4*b3*b3 + 
       8*d1*d2 + 4*d1*d3 + 4*d2*d3 + 4*d3*d3)*t;
    BTDB(25,29) = (8*c1*d2 + 4*c3*d2 + 4*c1*d3 + 4*c3*d3)*s + 
      (8*c2*d1 + 4*c3*d1 + 4*c2*d3 + 4*c3*d3)*t;
    
    BTDB(26,26) = (8*d1*d1 + 8*d1*d3 + 8*d3*d3)*r + 
      (8*b1*b1 + 8*b1*b3 + 8*b3*b3 + 
       8*c1*c1 + 8*c1*c3 + 8*c3*c3)*t;
    BTDB(26,27) = (8*b2*d1 + 4*b3*d1 + 4*b2*d3 + 4*b3*d3)*s + 
      (8*b1*d2 + 4*b3*d2 + 4*b1*d3 + 4*b3*d3)*t;
    BTDB(26,28) = (8*c2*d1 + 4*c3*d1 + 4*c2*d3 + 4*c3*d3)*s + 
      (8*c1*d2 + 4*c3*d2 + 4*c1*d3 + 4*c3*d3)*t;
    BTDB(26,29) = (8*d1*d2 + 4*d1*d3 + 4*d2*d3 + 4*d3*d3)*r + 
      (8*b1*b2 + 4*b1*b3 + 4*b2*b3 + 4*b3*b3 + 
       8*c1*c2 + 4*c1*c3 + 4*c2*c3 + 4*c3*c3)*t;
    
    BTDB(27,27) = (8*b2*b2 + 8*b2*b3 + 8*b3*b3)*r + 
      (8*c2*c2 + 8*c2*c3 + 8*c3*c3 + 
       8*d2*d2 + 8*d2*d3 + 8*d3*d3)*t;
    BTDB(27,28) = (8*b2*c2 + 4*b3*c2 + 4*b2*c3 + 8*b3*c3)*s + 
      (8*b2*c2 + 4*b3*c2 + 4*b2*c3 + 8*b3*c3)*t;
    BTDB(27,29) = (8*b2*d2 + 4*b3*d2 + 4*b2*d3 + 8*b3*d3)*s + 
      (8*b2*d2 + 4*b3*d2 + 4*b2*d3 + 8*b3*d3)*t;
    
    BTDB(28,28) = (8*c2*c2 + 8*c2*c3 + 8*c3*c3)*r + 
      (8*b2*b2 + 8*b2*b3 + 8*b3*b3 + 
       8*d2*d2 + 8*d2*d3 + 8*d3*d3)*t;
    BTDB(28,29) = (8*c2*d2 + 4*c3*d2 + 4*c2*d3 + 8*c3*d3)*s + 
      (8*c2*d2 + 4*c3*d2 + 4*c2*d3 + 8*c3*d3)*t;
    
    BTDB(29,29) = (8*d2*d2 + 8*d2*d3 + 8*d3*d3)*r + 
      (8*b2*b2 + 8*b2*b3 + 8*b3*b3 + 
       8*c2*c2 + 8*c2*c3 + 8*c3*c3)*t;
    
    double scale = 1.0/(180.0*size);
    for (i = 0; i < 30; i++)
        for (j = i; j < 30; j++)
	    BTDB(i,j) *= scale;

    for (i = 1; i < 30; i++)
        for (j = 0; j < i; j++)
	    BTDB(i,j) = BTDB(j,i);

    return BTDB;
}

RVector Tetrahedron10::InitialStrainVector (double E, double nu,
    const RVector &e)
{
    double mu = E/(2.0*(1.0+nu));                 // shear modulus
    double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu)); // Lame modulus

    double r = lambda + 2.0*mu;
    double s = lambda;
    double t = mu;
    double e0=e[0], e1=e[1], e2=e[2], e3=e[3], e4=e[4], e5=e[5];
    
    RVector res(30);
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
    res[4] = 0;
    res[5] = 0;
    res[6] = 0;
    res[7] = 0;
    res[8] = 0;
    res[9] = 0;
    res[10] = 0;
    res[11] = 0;
    res[12] = (b0+b1)*e0*r + ((b0+b1)*e1 + (b0+b1)*e2)*s + ((c0+c1)*e3 + (d0+d1)*e5)*t;
    res[13] = (c0+c1)*e1*r + ((c0+c1)*e0 + (c0+c1)*e2)*s + ((b0+b1)*e3 + (d0+d1)*e4)*t;
    res[14] = (d0+d1)*e2*r + ((d0+d1)*e0 + (d0+d1)*e1)*s + ((c0+c1)*e4 + (b0+b1)*e5)*t;
    res[15] = (b0+b2)*e0*r + ((b0+b2)*e1 + (b0+b2)*e2)*s + ((c0+c2)*e3 + (d0+d2)*e5)*t;
    res[16] = (c0+c2)*e1*r + ((c0+c2)*e0 + (c0+c2)*e2)*s + ((b0+b2)*e3 + (d0+d2)*e4)*t;
    res[17] = (d0+d2)*e2*r + ((d0+d2)*e0 + (d0+d2)*e1)*s + ((c0+c2)*e4 + (b0+b2)*e5)*t;
    res[18] = (b0+b3)*e0*r + ((b0+b3)*e1 + (b0+b3)*e2)*s + ((c0+c3)*e3 + (d0+d3)*e5)*t;
    res[19] = (c0+c3)*e1*r + ((c0+c3)*e0 + (c0+c3)*e2)*s + ((b0+b3)*e3 + (d0+d3)*e4)*t;
    res[20] = (d0+d3)*e2*r + ((d0+d3)*e0 + (d0+d3)*e1)*s + ((c0+c3)*e4 + (b0+b3)*e5)*t;
    res[21] = (b1+b2)*e0*r + ((b1+b2)*e1 + (b1+b2)*e2)*s + ((c1+c2)*e3 + (d1+d2)*e5)*t;
    res[22] = (c1+c2)*e1*r + ((c1+c2)*e0 + (c1+c2)*e2)*s + ((b1+b2)*e3 + (d1+d2)*e4)*t;
    res[23] = (d1+d2)*e2*r + ((d1+d2)*e0 + (d1+d2)*e1)*s + ((c1+c2)*e4 + (b1+b2)*e5)*t;
    res[24] = (b1+b3)*e0*r + ((b1+b3)*e1 + (b1+b3)*e2)*s + ((c1+c3)*e3 + (d1+d3)*e5)*t;
    res[25] = (c1+c3)*e1*r + ((c1+c3)*e0 + (c1+c3)*e2)*s + ((b1+b3)*e3 + (d1+d3)*e4)*t;
    res[26] = (d1+d3)*e2*r + ((d1+d3)*e0 + (d1+d3)*e1)*s + ((c1+c3)*e4 + (b1+b3)*e5)*t;
    res[27] = (b2+b3)*e0*r + ((b2+b3)*e1 + (b2+b3)*e2)*s + ((c2+c3)*e3 + (d2+d3)*e5)*t;
    res[28] = (c2+c3)*e1*r + ((c2+c3)*e0 + (c2+c3)*e2)*s + ((b2+b3)*e3 + (d2+d3)*e4)*t;
    res[29] = (d2+d3)*e2*r + ((d2+d3)*e0 + (d2+d3)*e1)*s + ((c2+c3)*e4 + (b2+b3)*e5)*t;

    return res/6.0;
}

RVector Tetrahedron10::ThermalExpansionVector (double E, double nu,
    double alpha, double dT)
{
    RVector e0(6);
    e0[0] = e0[1] = e0[2] = alpha*dT;
    return InitialStrainVector (E, nu, e0);
}

double Tetrahedron10::ComputeSize (const NodeList &nlist) const
{
    return (1.0/6.0) * (a0+a1+a2+a3);
}

RSymMatrix Tetrahedron10::ComputeIntDD (const NodeList &nlist) const
{
    double b0b0, b0b1, b1b1, b0b2, b1b2, b0b3, b1b3, b2b2, b2b3, b3b3;
    double c0c0, c0c1, c1c1, c0c2, c1c2, c0c3, c1c3, c2c2, c2c3, c3c3;
    double d0d0, d0d1, d1d1, d0d2, d1d2, d0d3, d1d3, d2d2, d2d3, d3d3;

    RSymMatrix dd(10);
    dd(0,0) = 3.0*((b0b0=b0*b0) + (c0c0=c0*c0) + (d0d0=d0*d0));
    dd(1,0) = -(b0b1=b0*b1) - (c0c1=c0*c1) - (d0d1=d0*d1);
    dd(1,1) = 3.0*((b1b1=b1*b1) + (c1c1=c1*c1) + (d1d1=d1*d1));
    dd(2,0) = -(b0b2=b0*b2) - (c0c2=c0*c2) - (d0d2=d0*d2);
    dd(2,1) = -(b1b2=b1*b2) - (c1c2=c1*c2) - (d1d2=d1*d2);
    dd(2,2) = 3.0*((b2b2=b2*b2) + (c2c2=c2*c2) + (d2d2=d2*d2));
    dd(3,0) = -(b0b3=b0*b3) - (c0c3=c0*c3) - (d0d3=d0*d3);
    dd(3,1) = -(b1b3=b1*b3) - (c1c3=c1*c3) - (d1d3=d1*d3);
    dd(3,2) = -(b2b3=b2*b3) - (c2c3=c2*c3) - (d2d3=d2*d3);
    dd(3,3) = 3.0*((b3b3=b3*b3) + (c3c3=c3*c3) + (d3d3=d3*d3));
    dd(4,0) = -b0b0 - c0c0 - d0d0 + 3.0*(b0b1 + c0c1 + d0d1);
    dd(4,1) = -b1b1 - c1c1 - d1d1 + 3.0*(b0b1 + c0c1 + d0d1);
    dd(4,2) = -b0b2 - b1b2 - c0c2 - c1c2 - d0d2 - d1d2;
    dd(4,3) = -b0b3 - b1b3 - c0c3 - c1c3 - d0d3 - d1d3;
    dd(4,4) = 8.0*(b0b0 + b0b1 + b1b1 + c0c0 + c0c1 + c1c1 + d0d0 + d0d1 +
                    d1d1);
    dd(5,0) = -b0b0 - c0c0 - d0d0 + 3.0*(b0b2 + c0c2 + d0d2);
    dd(5,1) = -b0b1 - b1b2 - c0c1 - c1c2 - d0d1 - d1d2;
    dd(5,2) = -b2b2 - c2c2 - d2d2 + 3.0*(b0b2 + c0c2 + d0d2);
    dd(5,3) = -b0b3 - b2b3 - c0c3 - c2c3 - d0d3 - d2d3;
    dd(5,4) = 4.0*(b0b0 + b0b1 + b0b2 + c0c0 + c0c1 + c0c2 + d0d0 + d0d1 +
                    d0d2 + 2.0*(b1b2 + c1c2 + d1d2));
    dd(5,5) = 8.0*(b0b0 + b0b2 + b2b2 + c0c0 + c0c2 + c2c2 + d0d0 + d0d2 +
                    d2d2);
    dd(6,0) = -b0b0 - c0c0 - d0d0 + 3.0*(b0b3 + c0c3 + d0d3);
    dd(6,1) = -b0b1 - b1b3 - c0c1 - c1c3 - d0d1 - d1d3;
    dd(6,2) = -b0b2 - b2b3 - c0c2 - c2c3 - d0d2 - d2d3;
    dd(6,3) = -b3b3 - c3c3 - d3d3 + 3.0*(b0b3 + c0c3 + d0d3);
    dd(6,4) = 4.0*(b0b0 + b0b1 + b0b3 + c0c0 + c0c1 + c0c3 + d0d0 + d0d1 +
                    d0d3 + 2.0*(b1b3 + c1c3 + d1d3));
    dd(6,5) = 4.0*(b0b0 + b0b2 + b0b3 + c0c0 + c0c2 + c0c3 + d0d0 + d0d2 +
                    d0d3 + 2.0*(b2b3 + c2c3 + d2d3));
    dd(6,6) = 8.0*(b0b0 + b0b3 + b3b3 + c0c0 + c0c3 + c3c3 + d0d0 + d0d3 +
                    d3d3);
    dd(7,0) = -b0b1 - b0b2 - c0c1 - c0c2 - d0d1 - d0d2;
    dd(7,1) = -b1b1 - c1c1 - d1d1 + 3.0*(b1b2 + c1c2 + d1d2);
    dd(7,2) = -b2b2 - c2c2 - d2d2 + 3.0*(b1b2 + c1c2 + d1d2);
    dd(7,3) = -b1b3 - b2b3 - c1c3 - c2c3 - d1d3 - d2d3;
    dd(7,4) = 4.0*(b1b1 + b1b2 + b0b1 + c0c1 + c1c1 + c1c2 + d0d1 + d1d1 +
                    d1d2 + 2.0*(b0b2 + c0c2 + d0d2));
    dd(7,5) = 4.0*(b1b2 + b2b2 + b0b2 + c0c2 + c1c2 + c2c2 + d0d2 + d1d2 +
                    d2d2 + 2.0*(b0b1 + c0c1 + d0d1));
    dd(7,6) = 4.0*(b0b1 + b0b2 + b1b3 + b2b3 + c0c1 + c0c2 + c1c3 + c2c3 +
                    d0d1 + d0d2 + d1d3 + d2d3);
    dd(7,7) = 8.0*(b1b1 + b1b2 + b2b2 + c1c1 + c1c2 + c2c2 + d1d1 + d1d2 +
                    d2d2);
    dd(8,0) = -b0b1 - b0b3 - c0c1 - c0c3 - d0d1 - d0d3;
    dd(8,1) = -b1b1 - c1c1 - d1d1 + 3.0*(b1b3 + c1c3 + d1d3);
    dd(8,2) = -b1b2 - b2b3 - c1c2 - c2c3 - d1d2 - d2d3;
    dd(8,3) = -b3b3 - c3c3 - d3d3 + 3.0*(b1b3 + c1c3 + d1d3);
    dd(8,4) = 4.0*(b1b1 + b1b3 + b0b1 + c0c1 + c1c1 + c1c3 + d0d1 + d1d1 +
                    d1d3 + 2.0*(b0b3 + c0c3 + d0d3));
    dd(8,5) = 4.0*(b1b2 + b2b3 + b0b1 + b0b3 + c0c1 + c1c2 + c0c3 + c2c3 +
                    d0d1 + d1d2 + d0d3 + d2d3);
    dd(8,6) = 4.0*(b1b3 + b3b3 + b0b3 + c0c3 + c1c3 + c3c3 + d0d3 + d3d3 +
                    d1d3 + 2.0*(b0b1 + c0c1 + d0d1));
    dd(8,7) = 4.0*(b1b1 + b1b2 + b1b3 + c1c1 + c1c2 + c1c3 + d1d1 + d1d2 +
                    d1d3 + 2.0*(b2b3 + c2c3 + d2d3));
    dd(8,8) = 8.0*(b1b1 + b1b3 + b3b3 + c1c1 + c1c3 + c3c3 + d1d1 + d1d3 +
                    d3d3);
    dd(9,0) = -b0b2 - b0b3 - c0c2 - c0c3 - d0d2 - d0d3;
    dd(9,1) = -b1b2 - b1b3 - c1c2 - c1c3 - d1d2 - d1d3;
    dd(9,2) = -b2b2 - c2c2 - d2d2 + 3.0*(b2b3 + c2c3 + d2d3);
    dd(9,3) = -b3b3 - c3c3 - d3d3 + 3.0*(b2b3 + c2c3 + d2d3);
    dd(9,4) = 4.0*(b0b2 + b0b3 + b1b2 + b1b3 + c0c2 + c1c2 + c0c3 + c1c3 +
                    d0d2 + d1d2 + d0d3 + d1d3);
    dd(9,5) = 4.0*(b2b2 + b2b3 + b0b2 + c0c2 + c2c2 + c2c3 + d0d2 + d2d2 +
                    d2d3 + 2.0*(b0b3 + c0c3 + d0d3));
    dd(9,6) = 4.0*(b2b3 + b3b3 + b0b3 + c0c3 + c2c3 + c3c3 + d0d3 + d2d3 +
                    d3d3 + 2.0*(b0b2 + c0c2 + d0d2));
    dd(9,7) = 4.0*(b2b2 + b2b3 + b1b2 + c1c2 + c2c2 + c2c3 + d1d2 + d2d2 +
                    d2d3 + 2.0*(b1b3 + c1c3 + d1d3));
    dd(9,8) = 4.0*(b2b3 + b3b3 + b1b3 + c1c3 + c2c3 + c3c3 + d1d3 + d2d3 +
                    d3d3 + 2.0*(b1b2 + c1c2 + d1d2));
    dd(9,9) = 8.0*(b2b2 + b2b3 + b3b3 + c2c2 + c2c3 + c3c3 + d2d2 + d2d3 +
                    d3d3);

    return dd * (1.0/(180.0*size));
}

RSymMatrix Tetrahedron10::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(10);
    if (!bndel) return bff;

    for (int sd = 0; sd < 4; sd++) {
        if (bndside[sd]) {
	    bff += *sym_bndintff[sd] * side_size[sd];
	}
    }
    return bff;
}

RVector Tetrahedron10::BndIntFX (int side, double (*func)(const Point&),
    const NodeList &nlist) const
{
    int p, np;
    double f;
    const double *wght;
    const Point  *absc;
    RVector bint(4);

    np = QRule_tri_4_6 (&wght, &absc); // quadrature rule
    for (p = 0; p < np; p++) {
        Point loc  = SurfToLocal (side, absc[p]); // local quadrature point
	Point glob = Global (nlist, loc);         // global quadrature point
	RVector fun = LocalShapeF (loc);          // shape func at quad point
	f = func(glob);                  // user function at quadrature point
	bint += fun * (f*side_size[side]*wght[p]);
    }
    return bint;
}

RVector Tetrahedron10::BndIntFCos (int side, const RVector &cntcos, double a,
    const NodeList &nlist) const
{
    // Note: this is copied directly from Tetrahedron4.
    // Needs to be checked.
    int j, p;
    double d, f;
    RVector sum(10);
    double ssize = SideSize (side, nlist);

    static bool needsetup = true;
    static int pp;
    static Point *loc;

    if (needsetup) {
	const Point *ab;
	pp = GetBndSubsampleAbsc (side, ab);
	loc = new Point[pp];
	for (p = 0; p < pp; p++) {
	    loc[p] = SurfToLocal (side, ab[p]);
	}
	needsetup = false;
    }

    for (p = 0; p < pp; p++) {
	//Point loc  = SurfToLocal (side, ab[p]);
	Point glob = Global (nlist, loc[p]);
	d = glob.Dist (cntcos);
	if (d < a) { // sample is withing support of cosine
	    RVector F  = LocalShapeF (loc[p]);
	    f = cos ((d/a)*(0.5*Pi));
	    for (j = 0; j < 4; j++)
		sum[j] += f*F[j];
	}
    }
    sum *= ssize;
    return sum;
}

int Tetrahedron10::Intersection (const Point &p1, const Point &p2,
    Point *s, bool add_endpoints, bool boundary_only)
{
    int i, n = 0;
    double a, rx, ry, rz;
    double sx = p1[0], sy = p1[1], sz = p1[2];
    double dx = p2[0]-p1[0], dy = p2[1]-p1[1], dz = p2[2]-p1[2];
    Point p(3);
    
    // intersection with plane z=0
    if ((!boundary_only || bndside[0]) && dz) {
	a = -sz/dz;
	rx = sx + a*dx;
	ry = sy + a*dy;
	if (rx >= 0 && ry >= 0 && rx+ry <= 1) {
	    p[0] = rx;
	    p[1] = ry;
	    p[2] = 0.0;
	    s[n++] = p;
	}
    }

    // intersection with plane y=0
    if ((!boundary_only || bndside[1]) && dy) {
	a = -sy/dy;
	rx = sx + a*dx;
	rz = sz + a*dz;
	if (rx >= 0 && rz >= 0 && rx+rz <= 1) {
	    p[0] = rx;
	    p[1] = 0.0;
	    p[2] = rz;
	    s[n++] = p;
	}
    }

    // intersection with plane x=0
    if ((!boundary_only || bndside[2]) && dx) {
	a = -sx/dx;
	ry = sy + a*dy;
	rz = sz + a*dz;
	if (ry >= 0 && rz >= 0 && ry+rz <= 1) {
	    p[0] = 0.0;
	    p[1] = ry;
	    p[2] = rz;
	    s[n++] = p;
	}
    }

    // intersection with plane 1-x-y-z=0
    if (!boundary_only || bndside[3]) {
	a = (1-sx-sy-sz)/(dx+dy+dz);
	rx = sx + a*dx;
	ry = sy + a*dy;
	rz = sz + a*dz;
	if (rx >= 0 && ry >= 0 && rx+ry <= 1) {
	    p[0] = rx;
	    p[1] = ry;
	    p[2] = rz;
	    s[n++] = p;
	}
    }
    return n;
}
