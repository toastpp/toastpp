// ========================================================================= //
// TOAST v.15                                      (c) Martin Schweiger 1999 //
// Library: libfe     File: tri10.cc                                         //
//                                                                           //
// 10-noded isoparametric triangle for cubic shape functions and straight    //
// boundaries.                                                               //
// ========================================================================= //

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "lin_qr.h"
#include "tri_qr.h"

using namespace std;

static bool subsampling_initialised = false;
static const int nsample_lin = NSUBSAMPLE; // from toastdef.h
static const int nsample_tot = (nsample_lin*(nsample_lin + 1)) / 2;
static Point absc_sample[nsample_tot];

static const RSymMatrix sym_intff = RSymMatrix (10,
    " 76 \
      11  76 \
      11  11  76 \
      18   0  27  540 \
       0  18  27 -189  540 \
      27  18   0 -135  270  540 \
      27   0  18  -54 -135 -189  540 \
       0  27  18 -135  -54 -135  270  540 \
      18  27   0  270 -135  -54 -135 -189 540 \
      36  36  36  162  162  162  162  162 162 1944") * (1.0/6720.0);

static const RDenseMatrix full_intf0ff = RDenseMatrix (10, 10,
    "3720 280 280 2340 -900 60 60 -900 2340 1080 \
     280 280 80 219 219 96 -21 -21 96 -18 \
     280 80 280 96 -21 -21 96 219 219 -18 \
     2340 219 96 3240 -2430 -702 27 -810 1620 -1620 \
     -900 219 -21 -2430 1944 1161 -459 324 -810 972 \
     60 96 -21 -702 1161 2214 -864 -459 27 1458 \
     60 -21 96 27 -459 -864 2214 1161 -702 1458 \
     -900 -21 219 -810 324 -459 1161 1944 -2430 972 \
     2340 96 219 1620 -810 27 -702 -2430 3240 -1620 \
     1080 -18 -18 -1620 972 1458 1458 972 -1620 1296") * (1.0/739200.0);

static const RDenseMatrix full_intf1ff = RDenseMatrix (10, 10,
    "280 280 80 219 219 96 -21 -21 96 -18 \
     280 3720 280 -900 2340 2340 -900 60 60 1080 \
     80 280 280 -21 96 219 219 96 -21 -18 \
     219 -900 -21 1944 -2430 -810 324 -459 1161 972 \
     219 2340 96 -2430 3240 1620 -810 27 -702 -1620 \
     96 2340 219 -810 1620 3240 -2430 -702 27 -1620 \
     -21 -900 219 324 -810 -2430 1944 1161 -459 972 \
     -21 60 96 -459 27 -702 1161 2214 -864 1458 \
     96 60 -21 1161 -702 27 -459 -864 2214 1458 \
     -18 1080 -18 972 -1620 -1620 972 1458 1458 1296") * (1.0/739200.0);

static const RDenseMatrix full_intf2ff = RDenseMatrix (10, 10,
    "280 80 280 96 -21 -21 96 219 219 -18 \
     80 280 280 -21 96 219 219 96 -21 -18 \
     280 280 3720 60 60 -900 2340 2340 -900 1080 \
     96 -21 60 2214 -864 -459 27 -702 1161 1458 \
     -21 96 60 -864 2214 1161 -702 27 -459 1458 \
     -21 219 -900 -459 1161 1944 -2430 -810 324 972 \
     96 219 2340 27 -702 -2430 3240 1620 -810 -1620 \
     219 96 2340 -702 27 -810 1620 3240 -2430 -1620 \
     219 -21 -900 1161 -459 324 -810 -2430 1944 972 \
     -18 -18 1080 1458 1458 972 -1620 -1620 972 1296") * (1.0/739200.0);

static const RDenseMatrix full_intf3ff = RDenseMatrix (10, 10,
    "2340 219 96 3240 -2430 -702 27 -810 1620 -1620 \
     219 -900 -21 1944 -2430 -810 324 -459 1161 972 \
     96 -21 60 2214 -864 -459 27 -702 1161 1458 \
     3240 1944 2214 32805 -2187 1701 -3159 -7047 10935 18954 \
     -2430 -2430 -864 -2187 -2187 -3888 2187 2187 -3888 -7290 \
     -702 -810 -459 1701 -3888 -7047 2187 972 -486 -6318 \
     27 324 27 -3159 2187 2187 -3159 -486 -486 -3402 \
     -810 -459 -702 -7047 2187 972 -486 1701 -3888 -6318 \
     1620 1161 1161 10935 -3888 -486 -486 -3888 10935 12636 \
     -1620 972 1458 18954 -7290 -6318 -3402 -6318 12636 8748") *(1.0/739200.0);

static const RDenseMatrix full_intf4ff = RDenseMatrix (10, 10,
    "-900 219 -21 -2430 1944 1161 -459 324 -810 972 \
     219 2340 96 -2430 3240 1620 -810 27 -702 -1620 \
     -21 96 60 -864 2214 1161 -702 27 -459 1458 \
     -2430 -2430 -864 -2187 -2187 -3888 2187 2187 -3888 -7290 \
     1944 3240 2214 -2187 32805 10935 -7047 -3159 1701 18954 \
     1161 1620 1161 -3888 10935 10935 -3888 -486 -486 12636 \
     -459 -810 -702 2187 -7047 -3888 1701 -486 972 -6318 \
     324 27 27 2187 -3159 -486 -486 -3159 2187 -3402 \
     -810 -702 -459 -3888 1701 -486 972 2187 -7047 -6318 \
     972 -1620 1458 -7290 18954 12636 -6318 -3402 -6318 8748") *(1.0/739200.0);

static const RDenseMatrix full_intf5ff = RDenseMatrix (10, 10,
    "60 96 -21 -702 1161 2214 -864 -459 27 1458 \
     96 2340 219 -810 1620 3240 -2430 -702 27 -1620 \
     -21 219 -900 -459 1161 1944 -2430 -810 324 972 \
     -702 -810 -459 1701 -3888 -7047 2187 972 -486 -6318 \
     1161 1620 1161 -3888 10935 10935 -3888 -486 -486 12636 \
     2214 3240 1944 -7047 10935 32805 -2187 1701 -3159 18954 \
     -864 -2430 -2430 2187 -3888 -2187 -2187 -3888 2187 -7290 \
     -459 -702 -810 972 -486 1701 -3888 -7047 2187 -6318 \
     27 27 324 -486 -486 -3159 2187 2187 -3159 -3402 \
     1458 -1620 972 -6318 12636 18954 -7290 -6318 -3402 8748") *(1.0/739200.0);

static const RDenseMatrix full_intf6ff = RDenseMatrix (10, 10,
    "60 -21 96 27 -459 -864 2214 1161 -702 1458 \
     -21 -900 219 324 -810 -2430 1944 1161 -459 972 \
     96 219 2340 27 -702 -2430 3240 1620 -810 -1620 \
     27 324 27 -3159 2187 2187 -3159 -486 -486 -3402 \
     -459 -810 -702 2187 -7047 -3888 1701 -486 972 -6318 \
     -864 -2430 -2430 2187 -3888 -2187 -2187 -3888 2187 -7290 \
     2214 1944 3240 -3159 1701 -2187 32805 10935 -7047 18954 \
     1161 1161 1620 -486 -486 -3888 10935 10935 -3888 12636 \
     -702 -459 -810 -486 972 2187 -7047 -3888 1701 -6318 \
     1458 972 -1620 -3402 -6318 -7290 18954 12636 -6318 8748") *(1.0/739200.0);

static const RDenseMatrix full_intf7ff = RDenseMatrix (10, 10,
    "-900 -21 219 -810 324 -459 1161 1944 -2430 972 \
     -21 60 96 -459 27 -702 1161 2214 -864 1458 \
     219 96 2340 -702 27 -810 1620 3240 -2430 -1620 \
     -810 -459 -702 -7047 2187 972 -486 1701 -3888 -6318 \
     324 27 27 2187 -3159 -486 -486 -3159 2187 -3402 \
     -459 -702 -810 972 -486 1701 -3888 -7047 2187 -6318 \
     1161 1161 1620 -486 -486 -3888 10935 10935 -3888 12636 \
     1944 2214 3240 1701 -3159 -7047 10935 32805 -2187 18954 \
     -2430 -864 -2430 -3888 2187 2187 -3888 -2187 -2187 -7290 \
     972 1458 -1620 -6318 -3402 -6318 12636 18954 -7290 8748") *(1.0/739200.0);

static const RDenseMatrix full_intf8ff = RDenseMatrix (10, 10,
    "2340 96 219 1620 -810 27 -702 -2430 3240 -1620 \
     96 60 -21 1161 -702 27 -459 -864 2214 1458 \
     219 -21 -900 1161 -459 324 -810 -2430 1944 972 \
     1620 1161 1161 10935 -3888 -486 -486 -3888 10935 12636 \
     -810 -702 -459 -3888 1701 -486 972 2187 -7047 -6318 \
     27 27 324 -486 -486 -3159 2187 2187 -3159 -3402 \
     -702 -459 -810 -486 972 2187 -7047 -3888 1701 -6318 \
     -2430 -864 -2430 -3888 2187 2187 -3888 -2187 -2187 -7290 \
     3240 2214 1944 10935 -7047 -3159 1701 -2187 32805 18954 \
     -1620 1458 972 12636 -6318 -3402 -6318 -7290 18954 8748") *(1.0/739200.0);

static const RDenseMatrix full_intf9ff = RDenseMatrix (10, 10,
    "1080 -18 -18 -1620 972 1458 1458 972 -1620 1296 \
     -18 1080 -18 972 -1620 -1620 972 1458 1458 1296 \
     -18 -18 1080 1458 1458 972 -1620 -1620 972 1296 \
     -1620 972 1458 18954 -7290 -6318 -3402 -6318 12636 8748 \
     972 -1620 1458 -7290 18954 12636 -6318 -3402 -6318 8748 \
     1458 -1620 972 -6318 12636 18954 -7290 -6318 -3402 8748 \
     1458 972 -1620 -3402 -6318 -7290 18954 12636 -6318 8748 \
     972 1458 -1620 -6318 -3402 -6318 12636 18954 -7290 8748 \
     -1620 1458 972 12636 -6318 -3402 -6318 -7290 18954 8748 \
     1296 1296 1296 8748 8748 8748 8748 8748 8748 157464") * (1.0/739200.0);

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

static const RSymMatrix sym_intf0ff = RSymMatrix (10,
    "3720 \
     280 280 \
     280 80 280 \
     2340 219 96 3240 \
     -900 219 -21 -2430 1944 \
     60 96 -21 -702 1161 2214 \
     60 -21 96 27 -459 -864 2214 \
     -900 -21 219 -810 324 -459 1161 1944 \
     2340 96 219 1620 -810 27 -702 -2430 3240 \
     1080 -18 -18 -1620 972 1458 1458 972 -1620 1296") * (1.0/739200.0);

static const RSymMatrix sym_intf1ff = RSymMatrix (10,
    "280 \
     280 3720 \
     80 280 280 \
     219 -900 -21 1944 \
     219 2340 96 -2430 3240 \
     96 2340 219 -810 1620 3240 \
     -21 -900 219 324 -810 -2430 1944 \
     -21 60 96 -459 27 -702 1161 2214 \
     96 60 -21 1161 -702 27 -459 -864 2214 \
     -18 1080 -18 972 -1620 -1620 972 1458 1458 1296") * (1.0/739200.0);

static const RSymMatrix sym_intf2ff = RSymMatrix (10,
    "280 \
     80 280 \
     280 280 3720 \
     96 -21 60 2214 \
     -21 96 60 -864 2214 \
     -21 219 -900 -459 1161 1944 \
     96 219 2340 27 -702 -2430 3240 \
     219 96 2340 -702 27 -810 1620 3240 \
     219 -21 -900 1161 -459 324 -810 -2430 1944 \
     -18 -18 1080 1458 1458 972 -1620 -1620 972 1296") * (1.0/739200.0);

static const RSymMatrix sym_intf3ff = RSymMatrix (10,
    "2340 \
     219 -900 \
     96 -21 60 \
     3240 1944 2214 32805 \
     -2430 -2430 -864 -2187 -2187 \
     -702 -810 -459 1701 -3888 -7047 \
     27 324 27 -3159 2187 2187 -3159 \
     -810 -459 -702 -7047 2187 972 -486 1701 \
     1620 1161 1161 10935 -3888 -486 -486 -3888 10935 \
     -1620 972 1458 18954 -7290 -6318 -3402 -6318 12636 8748") *(1.0/739200.0);

static const RSymMatrix sym_intf4ff = RSymMatrix (10,
    "-900 \
     219 2340 \
     -21 96 60 \
     -2430 -2430 -864 -2187 \
     1944 3240 2214 -2187 32805 \
     1161 1620 1161 -3888 10935 10935 \
     -459 -810 -702 2187 -7047 -3888 1701 \
     324 27 27 2187 -3159 -486 -486 -3159 \
     -810 -702 -459 -3888 1701 -486 972 2187 -7047 \
     972 -1620 1458 -7290 18954 12636 -6318 -3402 -6318 8748") *(1.0/739200.0);

static const RSymMatrix sym_intf5ff = RSymMatrix (10,
    "60 \
     96 2340 \
     -21 219 -900 \
     -702 -810 -459 1701 \
     1161 1620 1161 -3888 10935 \
     2214 3240 1944 -7047 10935 32805 \
     -864 -2430 -2430 2187 -3888 -2187 -2187 \
     -459 -702 -810 972 -486 1701 -3888 -7047 \
     27 27 324 -486 -486 -3159 2187 2187 -3159 \
     1458 -1620 972 -6318 12636 18954 -7290 -6318 -3402 8748") *(1.0/739200.0);

static const RSymMatrix sym_intf6ff = RSymMatrix (10,
    "60 \
     -21 -900 \
     96 219 2340 \
     27 324 27 -3159 \
     -459 -810 -702 2187 -7047 \
     -864 -2430 -2430 2187 -3888 -2187 \
     2214 1944 3240 -3159 1701 -2187 32805 \
     1161 1161 1620 -486 -486 -3888 10935 10935 \
     -702 -459 -810 -486 972 2187 -7047 -3888 1701 \
     1458 972 -1620 -3402 -6318 -7290 18954 12636 -6318 8748") *(1.0/739200.0);

static const RSymMatrix sym_intf7ff = RSymMatrix (10,
    "-900 \
     -21 60 \
     219 96 2340 \
     -810 -459 -702 -7047 \
     324 27 27 2187 -3159 \
     -459 -702 -810 972 -486 1701 \
     1161 1161 1620 -486 -486 -3888 10935 \
     1944 2214 3240 1701 -3159 -7047 10935 32805 \
     -2430 -864 -2430 -3888 2187 2187 -3888 -2187 -2187 \
     972 1458 -1620 -6318 -3402 -6318 12636 18954 -7290 8748") *(1.0/739200.0);

static const RSymMatrix sym_intf8ff = RSymMatrix (10,
    "2340 \
     96 60 \
     219 -21 -900 \
     1620 1161 1161 10935 \
     -810 -702 -459 -3888 1701 \
     27 27 324 -486 -486 -3159 \
     -702 -459 -810 -486 972 2187 -7047 \
     -2430 -864 -2430 -3888 2187 2187 -3888 -2187 \
     3240 2214 1944 10935 -7047 -3159 1701 -2187 32805 \
     -1620 1458 972 12636 -6318 -3402 -6318 -7290 18954 8748") *(1.0/739200.0);

static const RSymMatrix sym_intf9ff = RSymMatrix (10,
    "1080 \
     -18 1080 \
     -18 -18 1080 \
     -1620 972 1458 18954 \
     972 -1620 1458 -7290 18954 \
     1458 -1620 972 -6318 12636 18954 \
     1458 972 -1620 -3402 -6318 -7290 18954 \
     972 1458 -1620 -6318 -3402 -6318 12636 18954 \
     -1620 1458 972 12636 -6318 -3402 -6318 -7290 18954 \
     1296 1296 1296 8748 8748 8748 8748 8748 8748 157464") * (1.0/739200.0);

static const RDenseMatrix bndintf = RDenseMatrix (3, 10,
   "1 1 0 3 3 0 0 0 0 0 \
    0 1 1 0 0 3 3 0 0 0 \
    1 0 1 0 0 0 0 3 3 0") * (1.0/8.0);

// BndIntFF over side 0
static const RSymMatrix sym_bndintff_sd0 = RSymMatrix (10,
   "128 \
     19 128 \
      0   0   0 \
     99 -36   0 648 \
    -36  99   0 -81 648 \
      0   0   0   0   0   0 \
      0   0   0   0   0   0   0 \
      0   0   0   0   0   0   0   0 \
      0   0   0   0   0   0   0   0   0 \
      0   0   0   0   0   0   0   0   0   0") * (1.0/1680.0);

// BndIntFF over side 1
static const RSymMatrix sym_bndintff_sd1 = RSymMatrix (10,
   "  0 \
      0 128 \
      0  19 128 \
      0   0   0   0 \
      0   0   0   0   0 \
      0  99 -36   0   0 648 \
      0 -36  99   0   0 -81 648 \
      0   0   0   0   0   0   0   0 \
      0   0   0   0   0   0   0   0   0 \
      0   0   0   0   0   0   0   0   0   0") * (1.0/1680.0);

// BndIntFF over side 2
static const RSymMatrix sym_bndintff_sd2 = RSymMatrix (10,
   "128 \
      0   0 \
     19   0 128 \
      0   0   0   0 \
      0   0   0   0   0 \
      0   0   0   0   0   0 \
      0   0   0   0   0   0   0 \
    -36   0  99   0   0   0   0 648 \
     99   0 -36   0   0   0   0 -81 648 \
      0   0   0   0   0   0   0   0   0   0") * (1.0/1680.0);

Triangle10::Triangle10 (const Triangle10 &el): Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Triangle10::Copy ()
{
    return new Triangle10(*this);
}

void Triangle10::Initialise (const NodeList &nlist)
{
#ifdef TRI10_STORE_COORDS
    n0x = nlist[Node[0]][0], n0y = nlist[Node[0]][1];
    n1x = nlist[Node[1]][0], n1y = nlist[Node[1]][1];
    n2x = nlist[Node[2]][0], n2y = nlist[Node[2]][1];
    n3x = nlist[Node[3]][0], n3y = nlist[Node[3]][1];
    n4x = nlist[Node[4]][0], n4y = nlist[Node[4]][1];
    n5x = nlist[Node[5]][0], n5y = nlist[Node[5]][1];
    n6x = nlist[Node[6]][0], n6y = nlist[Node[6]][1];
    n7x = nlist[Node[7]][0], n7y = nlist[Node[7]][1];
    n8x = nlist[Node[8]][0], n8y = nlist[Node[8]][1];
    n9x = nlist[Node[9]][0], n9y = nlist[Node[9]][1];
#endif

    double x0 = nlist[Node[0]][0], y0 = nlist[Node[0]][1];
    double x1 = nlist[Node[1]][0], y1 = nlist[Node[1]][1];
    double x2 = nlist[Node[2]][0], y2 = nlist[Node[2]][1];

    a0 = x1*y2 - x2*y1;  b0 = y1-y2;  c0 = x2-x1;
    a1 = x2*y0 - x0*y2;  b1 = y2-y0;  c1 = x0-x2;
    a2 = x0*y1 - x1*y0;  b2 = y0-y1;  c2 = x1-x0;

    Element_Unstructured_2D::Initialise (nlist);

#ifdef TRI10_STORE_INTFF
    intff.New(10);
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

int Triangle10::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 3, "Side index out of range");
    dASSERT(node >= 0 && node < 4, "Node index out of range");
    static int SN[3][4] = {{0,1,3,4},{1,2,5,6},{2,0,7,8}};
    return SN[side][node];
}

Point Triangle10::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");

    Point loc(2);
    double scale = 1.0/(a0+a1+a2);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1]) * scale;

    return loc;
}

Point Triangle10::NodeLocal (int node) const
{
    dASSERT(node >= 0 && node < 10, "Node index out of range");
    static Point2D nloc[10] = {
        Point2D(0,0), Point2D(1,0), Point2D(0,1),
	Point2D(1.0/3.0,0),       Point2D(2.0/3.0,0),
	Point2D(2.0/3.0,1.0/3.0), Point2D(1.0/3.0,2.0/3.0),
	Point2D(0,2.0/3.0),       Point2D(0,1.0/3.0),
	Point2D(1.0/3.0,1.0/3.0)
    };
    return nloc[node];
}

RVector Triangle10::DirectionCosine (int side, RDenseMatrix &/*jacin*/)
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

const RVector &Triangle10::LNormal (int side) const
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

bool Triangle10::LContains (const Point& loc, bool pad) const
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

RVector Triangle10::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RVector fun(10);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    double f0 = L0*(3.0*L0-1.0);
    double f1 = L1*(3.0*L1-1.0);
    double f2 = L2*(3.0*L2-1.0);
    fun[0] = 0.5 * f0 * (3.0*L0-2.0);
    fun[1] = 0.5 * f1 * (3.0*L1-2.0);
    fun[2] = 0.5 * f2 * (3.0*L2-2.0);
    fun[3] = 4.5 * L1 * f0;
    fun[4] = 4.5 * L0 * f1;
    fun[5] = 4.5 * L2 * f1;
    fun[6] = 4.5 * L1 * f2;
    fun[7] = 4.5 * L0 * f2;
    fun[8] = 4.5 * L2 * f0;
    fun[9] = 27.0 * L0 * L1 * L2;
    return fun;
}

void Triangle10::LocalShapeF (const Point &loc, RVector *fun) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    if (fun->Dim() != 10)
        fun->New(10);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    double f0 = L0*(3.0*L0-1.0);
    double f1 = L1*(3.0*L1-1.0);
    double f2 = L2*(3.0*L2-1.0);
    double *f = fun->data_buffer();
    f[0] = 0.5 * f0 * (3.0*L0-2.0);
    f[1] = 0.5 * f1 * (3.0*L1-2.0);
    f[2] = 0.5 * f2 * (3.0*L2-2.0);
    f[3] = 4.5 * L1 * f0;
    f[4] = 4.5 * L0 * f1;
    f[5] = 4.5 * L2 * f1;
    f[6] = 4.5 * L1 * f2;
    f[7] = 4.5 * L0 * f2;
    f[8] = 4.5 * L2 * f0;
    f[9] = 27.0 * L0 * L1 * L2;
}

RDenseMatrix Triangle10::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der (2,10);
    double xi = loc[0], eta = loc[1];
    double s = 1.0-xi-eta;
    double f1 = 3.0*s-1.0,   g1 = 3.0*s-2.0,    h1 = 13.5*xi*eta;
    double f2 = 3.0*xi-1.0,  g2 = 3.0*xi-2.0,   h2 = 13.5*s*xi;
    double f3 = 3.0*eta-1.0, g3 = 3.0*eta-2.0,  h3 = 13.5*s*eta;

    der(0,0) = der(1,0) = -0.5*g1*f1 - 1.5*s*(g1+f1);
    der(0,1) = 0.5*f2*g2 + 1.5*xi*(f2+g2);
    // der(1,1) = 0.0;
    // der(0,2) = 0.0;
    der(1,2) = 0.5*g3*f3 + 1.5*eta*(g3+f3);

    f1 *= 4.5, f2 *= 4.5, f3 *= 4.5;

    der(0,3) =  f1*(s-xi) - h2;
    der(1,3) = -f1*xi - h2;
    der(0,4) =  f2*(s-xi) + h2;
    der(1,4) = -f2*xi;
    der(0,5) =  f2*eta + h1;
    der(1,5) =  f2*xi;
    der(0,6) =  f3*eta;
    der(1,6) =  f3*xi + h1;
    der(0,7) = -f3*eta;
    der(1,7) =  f3*(s-eta) + h3;
    der(0,8) = -f1*eta - h3;
    der(1,8) =  f1*(s-eta) - h3;
    der(0,9) =  27.0*eta*(s-xi);
    der(1,9) =  27.0*xi*(s-eta);
    return der;
}

RVector Triangle10::GlobalShapeF (const NodeList& nlist, const Point& glob)
    const
{
    dASSERT(glob.Dim() == 2, "Invalid point dimension");
    RVector fun(10);
    double scale = 1.0/(2.0*size);
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1]);
    double f0 = 3.0*L0-1.0;
    double f1 = 3.0*L1-1.0;
    double f2 = 3.0*L2-1.0;
    fun[0] = 0.5*L0*f0*(3.0*L0-2.0);
    fun[1] = 0.5*L1*f1*(3.0*L1-2.0);
    fun[2] = 0.5*L2*f2*(3.0*L2-2.0);
    f0 *= 4.5, f1 *= 4.5, f2 *= 4.5;
    fun[3] = L0*L1*f0;
    fun[4] = L0*L1*f1;
    fun[5] = L1*L2*f1;
    fun[6] = L1*L2*f2;
    fun[7] = L2*L0*f2;
    fun[8] = L2*L0*f0;
    fun[9] = 27.0*L0*L1*L2;
    return fun;
}

double Triangle10::IntF (int i) const
{
    static const double fac = 1.0/120.0;
    return size*fac * (i<3 ? 4.0 : i<9 ? 9.0 : 54.0);
}

RSymMatrix Triangle10::IntFF () const
{
#ifdef TRI10_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Triangle10::IntFF (int i, int j) const
{
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10);
#ifdef TRI10_STORE_INTFF
    return intff(i,j);
#else
    return sym_intff(i,j) * size;
#endif
}

double Triangle10::IntFFF (int i, int j, int k) const {
    RANGE_CHECK(i >= 0 && i < 10 && j >= 0 && j < 10 && k >= 0 && k < 10);
    return full_intfff[i]->Get(j,k) * size;
}

RSymMatrix Triangle10::IntPFF (const RVector &P) const
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

double Triangle10::IntPFF (int i, int j, const RVector &P) const
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

#define a (1.0/53760.0)
static const double intfd0d0[10] =
{3544*a,484*a,484*a,3762*a,-1422*a,198*a,198*a,-1422*a,3762*a,1836*a};

static const double intfd0d1[10] =
{556*a,556*a,160*a,-45*a,-45*a,765*a,-369*a,-369*a,765*a,378*a};

static const double intfd0d2[10] =
{556*a,160*a,556*a,765*a,-369*a,-369*a,765*a,-45*a,-45*a,378*a};

static const double intfd0d3a[10] =
{1062*a,-90*a,117*a,4212*a,-2511*a,-1701*a,-324*a,-1053*a,2106*a,-810*a};
static const double intfd0d3b[10] =
{6444*a,765*a,765*a,8100*a,-3402*a,-405*a,-405*a,-3402*a,8100*a,2592*a};

static const double intfd0d4a[10] =
{-126*a,936*a,-45*a,-729*a,-162*a,1782*a,-243*a,324*a,-405*a,-324*a};
static const double intfd0d4b[10] =
{-2520*a,9*a,-369*a,-810*a,-648*a,-81*a,405*a,1620*a,-3402*a,-2268*a};

static const double intfd0d5a[10] =
{-126*a,171*a,-414*a,-324*a,-81*a,1782*a,891*a,729*a,-810*a,-810*a};
static const double intfd0d5b[10] =
{-126*a,936*a,-45*a,-729*a,-162*a,1782*a,-243*a,324*a,-405*a,-324*a};

static const double intfd0d6a[10] =
{-126*a,-45*a,936*a,-405*a,324*a,-243*a,1782*a,-162*a,-729*a,-324*a};
static const double intfd0d6b[10] =
{-126*a,-414*a,171*a,-810*a,729*a,891*a,1782*a,-81*a,-324*a,-810*a};

static const double intfd0d7a[10] =
{-126*a,-45*a,936*a,-405*a,324*a,-243*a,1782*a,-162*a,-729*a,-324*a};
static const double intfd0d7b[10] =
{-2520*a,-369*a,9*a,-3402*a,1620*a,405*a,-81*a,-648*a,-810*a,-2268*a};

static const double intfd0d8a[10] =
{1062*a,117*a,-90*a,2106*a,-1053*a,-324*a,-1701*a,-2511*a,4212*a,-810*a};
static const double intfd0d8b[10] =
{6444*a,765*a,765*a,8100*a,-3402*a,-405*a,-405*a,-3402*a,8100*a,2592*a};

static const double intfd0d9a[10] =
{216*a,270*a,270*a,162*a,-324*a,2268*a,2268*a,-324*a,162*a,-1944*a};
static const double intfd0d9b[10] =
{1404*a,216*a,594*a,2592*a,-1296*a,162*a,-324*a,-3564*a,5184*a,-1944*a};
static const double intfd0d9c[10] =
{1404*a,594*a,216*a,5184*a,-3564*a,-324*a,162*a,-1296*a,2592*a,-1944*a};

static const double intfd1d1[10] = 
{484*a,3544*a,484*a,-1422*a,3762*a,3762*a,-1422*a,198*a,198*a,1836*a};

static const double intfd1d2[10] =
{160*a,556*a,556*a,-369*a,765*a,-45*a,-45*a,765*a,-369*a,378*a};

static const double intfd1d3a[10] =
{9*a,-2520*a,-369*a,-648*a,-810*a,-3402*a,1620*a,405*a,-81*a,-2268*a};
static const double intfd1d3b[10] =
{936*a,-126*a,-45*a,-162*a,-729*a,-405*a,324*a,-243*a,1782*a,-324*a};

static const double intfd1d4a[10] =
{765*a,6444*a,765*a,-3402*a,8100*a,8100*a,-3402*a,-405*a,-405*a,2592*a};
static const double intfd1d4b[10] =
{-90*a,1062*a,117*a,-2511*a,4212*a,2106*a,-1053*a,-324*a,-1701*a,-810*a};

static const double intfd1d5a[10] =
{117*a,1062*a,-90*a,-1053*a,2106*a,4212*a,-2511*a,-1701*a,-324*a,-810*a};
static const double intfd1d5b[10] =
{765*a,6444*a,765*a,-3402*a,8100*a,8100*a,-3402*a,-405*a,-405*a,2592*a};

static const double intfd1d6a[10] =
{-45*a,-126*a,936*a,324*a,-405*a,-729*a,-162*a,1782*a,-243*a,-324*a};
static const double intfd1d6b[10] =
{-369*a,-2520*a,9*a,1620*a,-3402*a,-810*a,-648*a,-81*a,405*a,-2268*a};

static const double intfd1d7a[10] =
{-45*a,-126*a,936*a,324*a,-405*a,-729*a,-162*a,1782*a,-243*a,-324*a};
static const double intfd1d7b[10] =
{-414*a,-126*a,171*a,729*a,-810*a,-324*a,-81*a,1782*a,891*a,-810*a};

static const double intfd1d8a[10] =
{171*a,-126*a,-414*a,-81*a,-324*a,-810*a,729*a,891*a,1782*a,-810*a};
static const double intfd1d8b[10] =
{936*a,-126*a,-45*a,-162*a,-729*a,-405*a,324*a,-243*a,1782*a,-324*a};

static const double intfd1d9a[10] =
{216*a,1404*a,594*a,-1296*a,2592*a,5184*a,-3564*a,-324*a,162*a,-1944*a};
static const double intfd1d9b[10] =
{270*a,216*a,270*a,-324*a,162*a,162*a,-324*a,2268*a,2268*a,-1944*a};
static const double intfd1d9c[10] =
{594*a,1404*a,216*a,-3564*a,5184*a,2592*a,-1296*a,162*a,-324*a,-1944*a};

static const double intfd2d2[10] =
{484*a,484*a,3544*a,198*a,198*a,-1422*a,3762*a,3762*a,-1422*a,1836*a};

static const double intfd2d3a[10] =
{171*a,-414*a,-126*a,1782*a,891*a,729*a,-810*a,-324*a,-81*a,-810*a};
static const double intfd2d3b[10] =
{936*a,-45*a,-126*a,1782*a,-243*a,324*a,-405*a,-729*a,-162*a,-324*a};

static const double intfd2d4a[10] =
{-45*a,936*a,-126*a,-243*a,1782*a,-162*a,-729*a,-405*a,324*a,-324*a};
static const double intfd2d4b[10] =
{-414*a,171*a,-126*a,891*a,1782*a,-81*a,-324*a,-810*a,729*a,-810*a};
static const double intfd2d5a[10] =

{-369*a,9*a,-2520*a,405*a,-81*a,-648*a,-810*a,-3402*a,1620*a,-2268*a};
static const double intfd2d5b[10] =
{-45*a,936*a,-126*a,-243*a,1782*a,-162*a,-729*a,-405*a,324*a,-324*a};

static const double intfd2d6a[10] =
{765*a,765*a,6444*a,-405*a,-405*a,-3402*a,8100*a,8100*a,-3402*a,2592*a};
static const double intfd2d6b[10] =
{117*a,-90*a,1062*a,-324*a,-1701*a,-2511*a,4212*a,2106*a,-1053*a,-810*a};

static const double intfd2d7a[10] =
{765*a,765*a,6444*a,-405*a,-405*a,-3402*a,8100*a,8100*a,-3402*a,2592*a};
static const double intfd2d7b[10] =
{-90*a,117*a,1062*a,-1701*a,-324*a,-1053*a,2106*a,4212*a,-2511*a,-810*a};

static const double intfd2d8a[10] =
{9*a,-369*a,-2520*a,-81*a,405*a,1620*a,-3402*a,-810*a,-648*a,-2268*a};
static const double intfd2d8b[10] =
{936*a,-45*a,-126*a,1782*a,-243*a,324*a,-405*a,-729*a,-162*a,-324*a};

static const double intfd2d9a[10] =
{216*a,594*a,1404*a,162*a,-324*a,-3564*a,5184*a,2592*a,-1296*a,-1944*a};
static const double intfd2d9b[10] =
{594*a,216*a,1404*a,-324*a,162*a,-1296*a,2592*a,5184*a,-3564*a,-1944*a};
static const double intfd2d9c[10] =
{270*a,270*a,216*a,2268*a,2268*a,-324*a,162*a,162*a,-324*a,-1944*a};

static const double intfd3d3a[10] =
{0,1944*a,1620*a,14580*a,5832*a,5832*a,-2916*a,-3888*a,4860*a,17496*a};
static const double intfd3d3b[10] =
{3564*a,1620*a,1782*a,25272*a,-5346*a,-486*a,-1944*a,-6318*a,12636*a,14580*a};
static const double intfd3d3c[10] =
{11664*a,1782*a,1782*a,18954*a,-7290*a,-972*a,-972*a,-7290*a,18954*a,8748*a};

static const double intfd3d4a[10] =
{-567*a,-4860*a,-243*a,1458*a,1458*a,-7290*a,1458*a,-243*a,-243*a,0};
static const double intfd3d4b[10] =
{-1782*a,-1782*a,810*a,7776*a,7776*a,-486*a,-1944*a,-1944*a,-486*a,10206*a};
static const double intfd3d4c[10] =
{-4860*a,-567*a,-243*a,1458*a,1458*a,-243*a,-243*a,1458*a,-7290*a,0};

static const double intfd3d5a[10] =
{-486*a,-972*a,324*a,1944*a,486*a,-6318*a,-1944*a,-972*a,486*a,2916*a};
static const double intfd3d5b[10] =
{-486*a,-81*a,162*a,486*a,-243*a,-972*a,-243*a,243*a,-1944*a,-1458*a};
static const double intfd3d5c[10] =
{-567*a,-4860*a,-243*a,1458*a,1458*a,-7290*a,1458*a,-243*a,-243*a,0};
static const double intfd3d5d[10] =
{-486*a,-486*a,-162*a,-486*a,-486*a,-972*a,486*a,486*a,-972*a,-1458*a};

static const double intfd3d6a[10] =
{-81*a,162*a,-486*a,-972*a,-243*a,243*a,-1944*a,486*a,-243*a,-1458*a};
static const double intfd3d6b[10] =
{-486*a,-162*a,-486*a,-972*a,486*a,486*a,-972*a,-486*a,-486*a,-1458*a};
static const double intfd3d6c[10] =
{-162*a,1944*a,-162*a,-2916*a,-2916*a,-2916*a,-2916*a,486*a,486*a,0};
static const double intfd3d6d[10] =
{-486*a,162*a,-81*a,-1944*a,243*a,-243*a,-972*a,-243*a,486*a,-1458*a};

static const double intfd3d7a[10] =
{-81*a,162*a,-486*a,-972*a,-243*a,243*a,-1944*a,486*a,-243*a,-1458*a};
static const double intfd3d7b[10] =
{-486*a,-162*a,-486*a,-972*a,486*a,486*a,-972*a,-486*a,-486*a,-1458*a};
static const double intfd3d7c[10] =
{-972*a,324*a,-486*a,-6318*a,-1944*a,-972*a,486*a,1944*a,486*a,2916*a};
static const double intfd3d7d[10] =
{-4860*a,-243*a,-567*a,-7290*a,1458*a,-243*a,-243*a,1458*a,1458*a,0};

static const double intfd3d8a[10] =
{0,648*a,648*a,4860*a,-972*a,972*a,972*a,-972*a,4860*a,11664*a};
static const double intfd3d8b[10] =
{1782*a,891*a,810*a,6318*a,-3159*a,-972*a,-243*a,-2673*a,12636*a,7290*a};
static const double intfd3d8c[10] =
{1782*a,810*a,891*a,12636*a,-2673*a,-243*a,-972*a,-3159*a,6318*a,7290*a};
static const double intfd3d8d[10] =
{11664*a,1782*a,1782*a,18954*a,-7290*a,-972*a,-972*a,-7290*a,18954*a,8748*a};

static const double intfd3d9a[10] =
{-810*a,-972*a,324*a,4374*a,0,-8748*a,-4374*a,-1458*a,2916*a,8748*a};
static const double intfd3d9b[10] =
{-648*a,810*a,810*a,10206*a,-2916*a,-2916*a,-2916*a,-2916*a,10206*a,17496*a};
static const double intfd3d9c[10] =
{1944*a,1134*a,810*a,8748*a,-4374*a,-1458*a,-1458*a,-4374*a,17496*a,8748*a};
static const double intfd3d9d[10] =
{-648*a,-972*a,2106*a,21870*a,8748*a,0,-4374*a,-5832*a,7290*a,26244*a};
static const double intfd3d9e[10] =
{1944*a,810*a,1134*a,17496*a,-4374*a,-1458*a,-1458*a,-4374*a,8748*a,8748*a};

static const double intfd4d4a[10] =
{1782*a,11664*a,1782*a,-7290*a,18954*a,18954*a,-7290*a,-972*a,-972*a,8748*a};
static const double intfd4d4b[10] =
{1620*a,3564*a,1782*a,-5346*a,25272*a,12636*a,-6318*a,-1944*a,-486*a,14580*a};
static const double intfd4d4c[10] =
{1944*a,0,1620*a,5832*a,14580*a,4860*a,-3888*a,-2916*a,5832*a,17496*a};

static const double intfd4d5a[10] = // b0b1
{891*a,1782*a,810*a,-3159*a,6318*a,12636*a,-2673*a,-243*a,-972*a,7290*a};
static const double intfd4d5b[10] = // b1b1
{648*a,0,648*a,-972*a,4860*a,4860*a,-972*a,972*a,972*a,11664*a};
static const double intfd4d5c[10] = // b0b2
{1782*a,11664*a,1782*a,-7290*a,18954*a,18954*a,-7290*a,-972*a,-972*a,8748*a};
static const double intfd4d5d[10] = // b1b2
{810*a,1782*a,891*a,-2673*a,12636*a,6318*a,-3159*a,-972*a,-243*a,7290*a};

static const double intfd4d6a[10] = //b0b1
{-162*a,-486*a,-486*a,486*a,-972*a,-486*a,-486*a,-972*a,486*a,-1458*a};
static const double intfd4d6b[10] = //b1b1
{162*a,-81*a,-486*a,-243*a,-972*a,-243*a,486*a,-1944*a,243*a,-1458*a};
static const double intfd4d6c[10] = //b0b2
{-243*a,-4860*a,-567*a,1458*a,-7290*a,1458*a,1458*a,-243*a,-243*a,0};
static const double intfd4d6d[10] = //b1b2
{324*a,-972*a,-486*a,-1944*a,-6318*a,486*a,1944*a,486*a,-972*a,2916*a};

static const double intfd4d7a[10] = //b0b0
{-162*a,-486*a,-486*a,486*a,-972*a,-486*a,-486*a,-972*a,486*a,-1458*a};
static const double intfd4d7b[10] = //b0b1
{162*a,-81*a,-486*a,-243*a,-972*a,-243*a,486*a,-1944*a,243*a,-1458*a};
static const double intfd4d7c[10] = //b0b2
{162*a,-486*a,-81*a,243*a,-1944*a,486*a,-243*a,-972*a,-243*a,-1458*a};
static const double intfd4d7d[10] = //b1b2
{1944*a,-162*a,-162*a,-2916*a,-2916*a,486*a,486*a,-2916*a,-2916*a,0};

static const double intfd4d8a[10] = //b0b0
{-81*a,-486*a,162*a,-243*a,486*a,-1944*a,243*a,-243*a,-972*a,-1458*a};
static const double intfd4d8b[10] = //b0b1
{-972*a,-486*a,324*a,486*a,1944*a,486*a,-972*a,-1944*a,-6318*a,2916*a};
static const double intfd4d8c[10] = //b0b2
{-486*a,-486*a,-162*a,-486*a,-486*a,-972*a,486*a,486*a,-972*a,-1458*a};
static const double intfd4d8d[10] = //b1b2
{-4860*a,-567*a,-243*a,1458*a,1458*a,-243*a,-243*a,1458*a,-7290*a,0};

static const double intfd4d9a[10] = //b0b0
{1134*a,1944*a,810*a,-4374*a,8748*a,17496*a,-4374*a,-1458*a,-1458*a,8748*a};
static const double intfd4d9b[10] = //b0b1
{810*a,-648*a,810*a,-2916*a,10206*a,10206*a,-2916*a,-2916*a,-2916*a,17496*a};
static const double intfd4d9c[10] = //b1b1
{-972*a,-810*a,324*a,0,4374*a,2916*a,-1458*a,-4374*a,-8748*a,8748*a};
static const double intfd4d9d[10] = //b0b2
{810*a,1944*a,1134*a,-4374*a,17496*a,8748*a,-4374*a,-1458*a,-1458*a,8748*a};
static const double intfd4d9e[10] = //b1b2
{-972*a,-648*a,2106*a,8748*a,21870*a,7290*a,-5832*a,-4374*a,0,26244*a};

static const double intfd5d5a[10] = //b1b1
{1620*a,0,1944*a,-3888*a,4860*a,14580*a,5832*a,5832*a,-2916*a,17496*a};
static const double intfd5d5b[10] = //b1b2
{1782*a,3564*a,1620*a,-6318*a,12636*a,25272*a,-5346*a,-486*a,-1944*a,14580*a};
static const double intfd5d5c[10] = //b2b2
{1782*a,11664*a,1782*a,-7290*a,18954*a,18954*a,-7290*a,-972*a,-972*a,8748*a};

static const double intfd5d6a[10] = //b1b1
{-243*a,-567*a,-4860*a,-243*a,-243*a,1458*a,1458*a,-7290*a,1458*a,0};
static const double intfd5d6b[10] = //b1b2
{810*a,-1782*a,-1782*a,-1944*a,-486*a,7776*a,7776*a,-486*a,-1944*a,10206*a};
static const double intfd5d6c[10] = //b2b2
{-243*a,-4860*a,-567*a,1458*a,-7290*a,1458*a,1458*a,-243*a,-243*a,0};

static const double intfd5d7a[10] = //b0b1
{-243*a,-567*a,-4860*a,-243*a,-243*a,1458*a,1458*a,-7290*a,1458*a,0};
static const double intfd5d7b[10] = //b0b2
{-162*a,-486*a,-486*a,486*a,-972*a,-486*a,-486*a,-972*a,486*a,-1458*a};
static const double intfd5d7c[10] = //b1b2
{324*a,-486*a,-972*a,-972*a,486*a,1944*a,486*a,-6318*a,-1944*a,2916*a};
static const double intfd5d7d[10] = //b2b2
{162*a,-486*a,-81*a,243*a,-1944*a,486*a,-243*a,-972*a,-243*a,-1458*a};

static const double intfd5d8a[10] = //b0b1
{-162*a,-162*a,1944*a,486*a,486*a,-2916*a,-2916*a,-2916*a,-2916*a,0};
static const double intfd5d8b[10] = //b0b2
{-81*a,-486*a,162*a,-243*a,486*a,-1944*a,243*a,-243*a,-972*a,-1458*a};
static const double intfd5d8c[10] = //b1b2
{-486*a,-81*a,162*a,486*a,-243*a,-972*a,-243*a,243*a,-1944*a,-1458*a};
static const double intfd5d8d[10] = //b2b2
{-486*a,-486*a,-162*a,-486*a,-486*a,-972*a,486*a,486*a,-972*a,-1458*a};

static const double intfd5d9a[10] = //b0b1
{2106*a,-648*a,-972*a,-5832*a,7290*a,21870*a,8748*a,0,-4374*a,26244*a};
static const double intfd5d9b[10] = //b1b1
{324*a,-810*a,-972*a,-1458*a,2916*a,4374*a,0,-8748*a,-4374*a,8748*a};
static const double intfd5d9c[10] = //b0b2
{1134*a,1944*a,810*a,-4374*a,8748*a,17496*a,-4374*a,-1458*a,-1458*a,8748*a};
static const double intfd5d9d[10] = //b1b2
{810*a,-648*a,810*a,-2916*a,10206*a,10206*a,-2916*a,-2916*a,-2916*a,17496*a};
static const double intfd5d9e[10] = //b2b2
{810*a,1944*a,1134*a,-4374*a,17496*a,8748*a,-4374*a,-1458*a,-1458*a,8748*a};

static const double intfd6d6a[10] = //b1b1
{1782*a,1782*a,11664*a,-972*a,-972*a,-7290*a,18954*a,18954*a,-7290*a,8748*a};
static const double intfd6d6b[10] = //b1b2
{1782*a,1620*a,3564*a,-1944*a,-486*a,-5346*a,25272*a,12636*a,-6318*a,14580*a};
static const double intfd6d6c[10] = //b2b2
{1620*a,1944*a,0,-2916*a,5832*a,5832*a,14580*a,4860*a,-3888*a,17496*a};

static const double intfd6d7a[10] = //b0b1
{1782*a,1782*a,11664*a,-972*a,-972*a,-7290*a,18954*a,18954*a,-7290*a,8748*a};
static const double intfd6d7b[10] = //b0b2
{891*a,810*a,1782*a,-972*a,-243*a,-2673*a,12636*a,6318*a,-3159*a,7290*a};
static const double intfd6d7c[10] = //b1b2
{810*a,891*a,1782*a,-243*a,-972*a,-3159*a,6318*a,12636*a,-2673*a,7290*a};
static const double intfd6d7d[10] = //b2b2
{648*a,648*a,0,972*a,972*a,-972*a,4860*a,4860*a,-972*a,11664*a};

static const double intfd6d8a[10] = //b0b1
{-567*a,-243*a,-4860*a,-243*a,-243*a,1458*a,-7290*a,1458*a,1458*a,0};
static const double intfd6d8b[10] = //b0b2
{-486*a,324*a,-972*a,486*a,-972*a,-1944*a,-6318*a,486*a,1944*a,2916*a};
static const double intfd6d8c[10] = //b1b2
{-486*a,-162*a,-486*a,-972*a,486*a,486*a,-972*a,-486*a,-486*a,-1458*a};
static const double intfd6d8d[10] = //b2b2
{-486*a,162*a,-81*a,-1944*a,243*a,-243*a,-972*a,-243*a,486*a,-1458*a};

static const double intfd6d9a[10] = //b0b1
{1134*a,810*a,1944*a,-1458*a,-1458*a,-4374*a,17496*a,8748*a,-4374*a,8748*a};
static const double intfd6d9b[10] = //b1b1
{810*a,1134*a,1944*a,-1458*a,-1458*a,-4374*a,8748*a,17496*a,-4374*a,8748*a};
static const double intfd6d9c[10] = //b0b2
{2106*a,-972*a,-648*a,-4374*a,0,8748*a,21870*a,7290*a,-5832*a,26244*a};
static const double intfd6d9d[10] = //b1b2
{810*a,810*a,-648*a,-2916*a,-2916*a,-2916*a,10206*a,10206*a,-2916*a,17496*a};
static const double intfd6d9e[10] = //b2b2
{324*a,-972*a,-810*a,-4374*a,-8748*a,0,4374*a,2916*a,-1458*a,8748*a};

static const double intfd7d7a[10] = //b0b0
{1782*a,1782*a,11664*a,-972*a,-972*a,-7290*a,18954*a,18954*a,-7290*a,8748*a};
static const double intfd7d7b[10] = //b0b2
{1620*a,1782*a,3564*a,-486*a,-1944*a,-6318*a,12636*a,25272*a,-5346*a,14580*a};
static const double intfd7d7c[10] = //b2b2
{1944*a,1620*a,0,5832*a,-2916*a,-3888*a,4860*a,14580*a,5832*a,17496*a};

static const double intfd7d8a[10] = //b0b0
{-567*a,-243*a,-4860*a,-243*a,-243*a,1458*a,-7290*a,1458*a,1458*a,0};
static const double intfd7d8b[10] = //b0b2
{-1782*a,810*a,-1782*a,-486*a,-1944*a,-1944*a,-486*a,7776*a,7776*a,10206*a};
static const double intfd7d8c[10] = //b2b2
{-4860*a,-243*a,-567*a,-7290*a,1458*a,-243*a,-243*a,1458*a,1458*a,0};

static const double intfd7d9a[10] = //b0b0
{1134*a,810*a,1944*a,-1458*a,-1458*a,-4374*a,17496*a,8748*a,-4374*a,8748*a};
static const double intfd7d9b[10] = //b0b1
{810*a,1134*a,1944*a,-1458*a,-1458*a,-4374*a,8748*a,17496*a,-4374*a,8748*a};
static const double intfd7d9c[10] = //b0b2
{810*a,810*a,-648*a,-2916*a,-2916*a,-2916*a,10206*a,10206*a,-2916*a,17496*a};
static const double intfd7d9d[10] = //b1b2
{-972*a,2106*a,-648*a,0,-4374*a,-5832*a,7290*a,21870*a,8748*a,26244*a};
static const double intfd7d9e[10] = //b2b2
{-972*a,324*a,-810*a,-8748*a,-4374*a,-1458*a,2916*a,4374*a,0,8748*a};

static const double intfd8d8a[10] = //b0b0
{0,1620*a,1944*a,4860*a,-3888*a,-2916*a,5832*a,5832*a,14580*a,17496*a};
static const double intfd8d8b[10] = //b0b2
{3564*a,1782*a,1620*a,12636*a,-6318*a,-1944*a,-486*a,-5346*a,25272*a,14580*a};
static const double intfd8d8c[10] = //b2b2
{11664*a,1782*a,1782*a,18954*a,-7290*a,-972*a,-972*a,-7290*a,18954*a,8748*a};

static const double intfd8d9a[10] = //b0b0
{-810*a,324*a,-972*a,2916*a,-1458*a,-4374*a,-8748*a,0,4374*a,8748*a};
static const double intfd8d9b[10] = //b0b1
{-648*a,2106*a,-972*a,7290*a,-5832*a,-4374*a,0,8748*a,21870*a,26244*a};
static const double intfd8d9c[10] = //b0b2
{-648*a,810*a,810*a,10206*a,-2916*a,-2916*a,-2916*a,-2916*a,10206*a,17496*a};
static const double intfd8d9d[10] = //b1b2
{1944*a,1134*a,810*a,8748*a,-4374*a,-1458*a,-1458*a,-4374*a,17496*a,8748*a};
static const double intfd8d9e[10] = //b2b2
{1944*a,810*a,1134*a,17496*a,-4374*a,-1458*a,-1458*a,-4374*a,8748*a,8748*a};

static const double intfd9d9a[10] = //b0b0
{3888*a,0,0,-8748*a,8748*a,26244*a,26244*a,8748*a,-8748*a,52488*a};
static const double intfd9d9b[10] = //b0b1
{1944*a,1944*a,0,0,0,0,17496*a,17496*a,0,69984*a};
static const double intfd9d9c[10] = //b1b1
{0,3888*a,0,8748*a,-8748*a,-8748*a,8748*a,26244*a,26244*a,52488*a};
static const double intfd9d9d[10] = //b0b2
{1944*a,0,1944*a,0,17496*a,17496*a,0,0,0,69984*a};
static const double intfd9d9e[10] = //b1b2
{0,1944*a,1944*a,17496*a,0,0,0,0,17496*a,69984*a};
static const double intfd9d9f[10] = //b2b2
{0,0,3888*a,26244*a,26244*a,8748*a,-8748*a,-8748*a,8748*a,52488*a};
#undef a

double Triangle10::IntFDD (int i, int j, int k) const
{
    if (j > k) { int tmp=j; j=k; k=tmp; } // symmetry
    switch (j) {
    case 0:
        switch (k) {
	case 0:
	    return intfd0d0[i] * (b0*b0+c0*c0)/size;
	case 1:
	    return intfd0d1[i] * (b0*b1+c0*c1)/size;
	case 2:
	    return intfd0d2[i] * (b0*b2+c0*c2)/size;
	case 3:
	    return (intfd0d3a[i] * (b0*b0+c0*c0) +
		    intfd0d3b[i] * (b0*b1+c0*c1))/size;
	case 4:
	    return (intfd0d4a[i] * (b0*b0+c0*c0) +
		    intfd0d4b[i] * (b0*b1+c0*c1))/size;
	case 5:
	    return (intfd0d5a[i] * (b0*b1+c0*c1) +
		    intfd0d5b[i] * (b0*b2+c0*c2))/size;
	case 6:
	    return (intfd0d6a[i] * (b0*b1+c0*c1) +
		    intfd0d6b[i] * (b0*b2+c0*c2))/size;
	case 7:
	    return (intfd0d7a[i] * (b0*b0+c0*c0) +
		    intfd0d7b[i] * (b0*b2+c0*c2))/size;
	case 8:
	    return (intfd0d8a[i] * (b0*b0+c0*c0) +
		    intfd0d8b[i] * (b0*b2+c0*c2))/size;
	case 9:
	    return (intfd0d9a[i] * (b0*b0+c0*c0) +
		    intfd0d9b[i] * (b0*b1+c0*c1) +
		    intfd0d9c[i] * (b0*b2+c0*c2))/size;
	}
    case 1:
        switch (k) {
	case 1:
	    return intfd1d1[i] * (b1*b1+c1*c1)/size;
	case 2:
	    return intfd1d2[i] * (b1*b2+c1*c2)/size;
	case 3:
	    return (intfd1d3a[i] * (b0*b1+c0*c1) +
		    intfd1d3b[i] * (b1*b1+c1*c1))/size;
	case 4:
	    return (intfd1d4a[i] * (b0*b1+c0*c1) +
		    intfd1d4b[i] * (b1*b1+c1*c1))/size;
	case 5:
	    return (intfd1d5a[i] * (b1*b1+c1*c1) +
		    intfd1d5b[i] * (b1*b2+c1*c2))/size;
	case 6:
	    return (intfd1d6a[i] * (b1*b1+c1*c1) +
		    intfd1d6b[i] * (b1*b2+c1*c2))/size;
	case 7:
	    return (intfd1d7a[i] * (b0*b1+c0*c1) +
		    intfd1d7b[i] * (b1*b2+c1*c2))/size;
	case 8:
	    return (intfd1d8a[i] * (b0*b1+c0*c1) +
		    intfd1d8b[i] * (b1*b2+c1*c2))/size;
	case 9:
	    return (intfd1d9a[i] * (b0*b1+c0*c1) +
		    intfd1d9b[i] * (b1*b1+c1*c1) +
		    intfd1d9c[i] * (b1*b2+c1*c2))/size;
	}
    case 2:
        switch (k) {
	case 2:
	    return intfd2d2[i] * (b2*b2+c2*c2)/size;
	case 3:
	    return (intfd2d3a[i] * (b0*b2+c0*c2) +
		    intfd2d3b[i] * (b1*b2+c1*c2))/size;
	case 4:
	    return (intfd2d4a[i] * (b0*b2+c0*c2) +
		    intfd2d4b[i] * (b1*b2+c1*c2))/size;
	case 5:
	    return (intfd2d5a[i] * (b1*b2+c1*c2) +
		    intfd2d5b[i] * (b2*b2+c2*c2))/size;
	case 6:
	    return (intfd2d6a[i] * (b1*b2+c1*c2) +
		    intfd2d6b[i] * (b2*b2+c2*c2))/size;
	case 7:
	    return (intfd2d7a[i] * (b0*b2+c0*c2) +
		    intfd2d7b[i] * (b2*b2+c2*c2))/size;
	case 8:
	    return (intfd2d8a[i] * (b0*b2+c0*c2) +
		    intfd2d8b[i] * (b2*b2+c2*c2))/size;
	case 9:
	    return (intfd2d9a[i] * (b0*b2+c0*c2) +
		    intfd2d9b[i] * (b1*b2+c1*c2) +
		    intfd2d9c[i] * (b2*b2+c2*c2))/size;
	}
    case 3:
        switch (k) {
	case 3:
	    return (intfd3d3a[i] * (b0*b0+c0*c0) +
		    intfd3d3b[i] * (b0*b1+c0*c1) +
		    intfd3d3c[i] * (b1*b1+c1*c1))/size;
	case 4:
	    return (intfd3d4a[i] * (b0*b0+c0*c0) +
		    intfd3d4b[i] * (b0*b1+c0*c1) +
		    intfd3d4c[i] * (b1*b1+c1*c1))/size;
	case 5:
	    return (intfd3d5a[i] * (b0*b1+c0*c1) +
		    intfd3d5b[i] * (b1*b1+c1*c1) +
		    intfd3d5c[i] * (b0*b2+c0*c2) +
		    intfd3d5d[i] * (b1*b2+c1*c2))/size;
	case 6:
	    return (intfd3d6a[i] * (b0*b1+c0*c1) +
		    intfd3d6b[i] * (b1*b1+c1*c1) +
		    intfd3d6c[i] * (b0*b2+c0*c2) +
		    intfd3d6d[i] * (b1*b2+c1*c2))/size;
	case 7:
	    return (intfd3d7a[i] * (b0*b0+c0*c0) +
		    intfd3d7b[i] * (b0*b1+c0*c1) +
		    intfd3d7c[i] * (b0*b2+c0*c2) +
		    intfd3d7d[i] * (b1*b2+c1*c2))/size;
	case 8:
	    return (intfd3d8a[i] * (b0*b0+c0*c0) +
		    intfd3d8b[i] * (b0*b1+c0*c1) +
		    intfd3d8c[i] * (b0*b2+c0*c2) +
		    intfd3d8d[i] * (b1*b2+c1*c2))/size;
	case 9:
	    return (intfd3d9a[i] * (b0*b0+c0*c0) +
		    intfd3d9b[i] * (b0*b1+c0*c1) +
		    intfd3d9c[i] * (b1*b1+c1*c1) +
		    intfd3d9d[i] * (b0*b2+c0*c2) +
		    intfd3d9e[i] * (b1*b2+c1*c2))/size;
	}
    case 4:
        switch (k) {
	case 4:
	    return (intfd4d4a[i] * (b0*b0+c0*c0) +
		    intfd4d4b[i] * (b0*b1+c0*c1) +
		    intfd4d4c[i] * (b1*b1+c1*c1))/size;
	case 5:
	    return (intfd4d5a[i] * (b0*b1+c0*c1) +
		    intfd4d5b[i] * (b1*b1+c1*c1) +
		    intfd4d5c[i] * (b0*b2+c0*c2) +
		    intfd4d5d[i] * (b1*b2+c1*c2))/size;
	case 6:
	    return (intfd4d6a[i] * (b0*b1+c0*c1) +
		    intfd4d6b[i] * (b1*b1+c1*c1) +
		    intfd4d6c[i] * (b0*b2+c0*c2) +
		    intfd4d6d[i] * (b1*b2+c1*c2))/size;
	case 7:
	    return (intfd4d7a[i] * (b0*b0+c0*c0) +
		    intfd4d7b[i] * (b0*b1+c0*c1) +
		    intfd4d7c[i] * (b0*b2+c0*c2) +
		    intfd4d7d[i] * (b1*b2+c1*c2))/size;
	case 8:
	    return (intfd4d8a[i] * (b0*b0+c0*c0) +
		    intfd4d8b[i] * (b0*b1+c0*c1) +
		    intfd4d8c[i] * (b0*b2+c0*c2) +
		    intfd4d8d[i] * (b1*b2+c1*c2))/size;
	case 9:
	    return (intfd4d9a[i] * (b0*b0+c0*c0) +
		    intfd4d9b[i] * (b0*b1+c0*c1) +
		    intfd4d9c[i] * (b1*b1+c1*c1) +
		    intfd4d9d[i] * (b0*b2+c0*c2) +
		    intfd4d9e[i] * (b1*b2+c1*c2))/size;
	}
    case 5:
        switch (k) {
	case 5:
	    return (intfd5d5a[i] * (b1*b1+c1*c1) +
		    intfd5d5b[i] * (b1*b2+c1*c2) +
		    intfd5d5c[i] * (b2*b2+c2*c2))/size;
	case 6:
	    return (intfd5d6a[i] * (b1*b1+c1*c1) +
		    intfd5d6b[i] * (b1*b2+c1*c2) +
		    intfd5d6c[i] * (b2*b2+c2*c2))/size;
	case 7:
	    return (intfd5d7a[i] * (b0*b1+c0*c1) +
		    intfd5d7b[i] * (b0*b2+c0*c2) +
		    intfd5d7c[i] * (b1*b2+c1*c2) +
		    intfd5d7d[i] * (b2*b2+c2*c2))/size;
	case 8:
	    return (intfd5d8a[i] * (b0*b1+c0*c1) +
		    intfd5d8b[i] * (b0*b2+c0*c2) +
		    intfd5d8c[i] * (b1*b2+c1*c2) +
		    intfd5d8d[i] * (b2*b2+c2*c2))/size;
	case 9:
	    return (intfd5d9a[i] * (b0*b1+c0*c1) +
		    intfd5d9b[i] * (b1*b1+c1*c1) +
		    intfd5d9c[i] * (b0*b2+c0*c2) +
		    intfd5d9d[i] * (b1*b2+c1*c2) +
		    intfd5d9e[i] * (b2*b2+c2*c2))/size;
	}
    case 6:
        switch (k) {
	case 6:
	    return (intfd6d6a[i] * (b1*b1+c1*c1) +
		    intfd6d6b[i] * (b1*b2+c1*c2) +
		    intfd6d6c[i] * (b2*b2+c2*c2))/size;
	case 7:
	    return (intfd6d7a[i] * (b0*b1+c0*c1) +
		    intfd6d7b[i] * (b0*b2+c0*c2) +
		    intfd6d7c[i] * (b1*b2+c1*c2) +
		    intfd6d7d[i] * (b2*b2+c2*c2))/size;
	case 8:
	    return (intfd6d8a[i] * (b0*b1+c0*c1) +
		    intfd6d8b[i] * (b0*b2+c0*c2) +
		    intfd6d8c[i] * (b1*b2+c1*c2) +
		    intfd6d8d[i] * (b2*b2+c2*c2))/size;
	case 9:
	    return (intfd6d9a[i] * (b0*b1+c0*c1) +
		    intfd6d9b[i] * (b1*b1+c1*c1) +
		    intfd6d9c[i] * (b0*b2+c0*c2) +
		    intfd6d9d[i] * (b1*b2+c1*c2) +
		    intfd6d9e[i] * (b2*b2+c2*c2))/size;
	}
    case 7:
        switch (k) {
	case 7:
	    return (intfd7d7a[i] * (b0*b0+c0*c0) +
		    intfd7d7b[i] * (b0*b2+c0*c2) +
		    intfd7d7c[i] * (b2*b2+c2*c2))/size;
	case 8:
	    return (intfd7d8a[i] * (b0*b0+c0*c0) +
		    intfd7d8b[i] * (b0*b2+c0*c2) +
		    intfd7d8c[i] * (b2*b2+c2*c2))/size;
	case 9:
	    return (intfd7d9a[i] * (b0*b0+c0*c0) +
		    intfd7d9b[i] * (b0*b1+c0*c1) +
		    intfd7d9c[i] * (b0*b2+c0*c2) +
		    intfd7d9d[i] * (b1*b2+c1*c2) +
		    intfd7d9e[i] * (b2*b2+c2*c2))/size;
	}
    case 8:
        switch (k) {
	case 8:
	    return (intfd8d8a[i] * (b0*b0+c0*c0) +
		    intfd8d8b[i] * (b0*b2+c0*c2) +
		    intfd8d8c[i] * (b2*b2+c2*c2))/size;
	case 9:
	    return (intfd8d9a[i] * (b0*b0+c0*c0) +
		    intfd8d9b[i] * (b0*b1+c0*c1) +
		    intfd8d9c[i] * (b0*b2+c0*c2) +
		    intfd8d9d[i] * (b1*b2+c1*c2) +
		    intfd8d9e[i] * (b2*b2+c2*c2))/size;
	}
    case 9:
        switch (k) {
	case 9:
	    return (intfd9d9a[i] * (b0*b0+c0*c0) +
		    intfd9d9b[i] * (b0*b1+c0*c1) +
		    intfd9d9c[i] * (b1*b1+c1*c1) +
		    intfd9d9d[i] * (b0*b2+c0*c2) +
		    intfd9d9e[i] * (b1*b2+c1*c2) +
		    intfd9d9f[i] * (b2*b2+c2*c2))/size;
	}
    }
    xERROR ("Index out of range");
    return 0.0;
}

double Triangle10::IntPDD (int j, int k, const RVector &P) const
{
    int i;
    double sum1, sum2, sum3, sum4, sum5, pi;

    if (j > k) { int tmp=j; j=k; k=tmp; } // symmetry
    switch (j) {
    case 0:
        switch (k) {
	case 0:
	    for (i = 0, sum1 = 0.0; i < 10; i++)
	        sum1 += intfd0d0[i] * P[Node[i]];
	    return sum1 * (b0*b0+c0*c0)/size;
	case 1:
	    for (i = 0, sum1 = 0.0; i < 10; i++)
	        sum1 += intfd0d1[i] * P[Node[i]];
	    return sum1 * (b0*b1+c0*c1)/size;
	case 2:
	    for (i = 0, sum1 = 0.0; i < 10; i++)
	        sum1 += intfd0d2[i] * P[Node[i]];
	    return sum1 * (b0*b2+c0*c2)/size;
	case 3:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd0d3a[i] * P[Node[i]];
		sum2 += intfd0d3b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b0+c0*c0) + 
		    sum2 * (b0*b1+c0*c1))/size;
	case 4:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd0d4a[i] * P[Node[i]];
		sum2 += intfd0d4b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1))/size;
	case 5:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd0d5a[i] * P[Node[i]];
		sum2 += intfd0d5b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b0*b2+c0*c2))/size;
	case 6:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd0d6a[i] * P[Node[i]];
		sum2 += intfd0d6b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b0*b2+c0*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	       sum1 += intfd0d7a[i] * P[Node[i]];
	       sum2 += intfd0d7b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b2+c0*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	       sum1 += intfd0d8a[i] * P[Node[i]];
	       sum2 += intfd0d8b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b2+c0*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd0d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd0d9b[i] * pi;
		sum3 += intfd0d9c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2))/size;
	}
    case 1:
        switch (k) {
	case 1:
	    for (i = 0, sum1 = 0.0; i < 10; i++)
	        sum1 += intfd1d1[i] * P[Node[i]];
	    return sum1 * (b1*b1+c1*c1)/size;
	case 2:
	    for (i = 0, sum1 = 0.0; i < 10; i++)
	        sum1 += intfd1d2[i] * P[Node[i]];
	    return sum1 * (b1*b2+c1*c2)/size;
	case 3:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd1d3a[i] * P[Node[i]];
		sum2 += intfd1d3b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1))/size;
	case 4:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd1d4a[i] * P[Node[i]];
		sum2 += intfd1d4b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1))/size;
	case 5:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd1d5a[i] * P[Node[i]];
		sum2 += intfd1d5b[i] * P[Node[i]];
	    }
	    return (sum1 * (b1*b1+c1*c1) +
		    sum2 * (b1*b2+c1*c2))/size;
	case 6:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd1d6a[i] * P[Node[i]];
		sum2 += intfd1d6b[i] * P[Node[i]];
	    }
	    return (sum1 * (b1*b1+c1*c1) +
		    sum2 * (b1*b2+c1*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd1d7a[i] * P[Node[i]];
		sum2 += intfd1d7b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b2+c1*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd1d8a[i] * P[Node[i]];
		sum2 += intfd1d8b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b2+c1*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd1d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd1d9b[i] * pi;
		sum3 += intfd1d9c[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b1*b2+c1*c2))/size;
	}
    case 2:
        switch (k) {
	case 2:
	    for (i = 0, sum1 = 0.0; i < 10; i++)
	        sum1 += intfd2d2[i] * P[Node[i]];
	    return sum1 * (b2*b2+c2*c2)/size;
	case 3:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd2d3a[i] * P[Node[i]];
		sum2 += intfd2d3b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b2+c0*c2) +
		    sum2 * (b1*b2+c1*c2))/size;
	case 4:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd2d4a[i] * P[Node[i]];
		sum2 += intfd2d4b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b2+c0*c2) +
		    sum2 * (b1*b2+c1*c2))/size;
	case 5:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd2d5a[i] * P[Node[i]];
		sum2 += intfd2d5b[i] * P[Node[i]];
	    }
	    return (sum1 * (b1*b2+c1*c2) +
		    sum2 * (b2*b2+c2*c2))/size;
	case 6:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd2d6a[i] * P[Node[i]];
		sum2 += intfd2d6b[i] * P[Node[i]];
	    }
	    return (sum1 * (b1*b2+c1*c2) +
		    sum2 * (b2*b2+c2*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd2d7a[i] * P[Node[i]];
		sum2 += intfd2d7b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b2+c0*c2) +
		    sum2 * (b2*b2+c2*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2 = 0.0; i < 10; i++) {
	        sum1 += intfd2d8a[i] * P[Node[i]];
		sum2 += intfd2d8b[i] * P[Node[i]];
	    }
	    return (sum1 * (b0*b2+c0*c2) +
		    sum2 * (b2*b2+c2*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd2d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd2d9b[i] * pi;
		sum3 += intfd2d9c[i] * pi;
	    }
	    return (sum1 * (b0*b2+c0*c2) +
		    sum2 * (b1*b2+c1*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	}
    case 3:
        switch (k) {
	case 3:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd3d3a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d3b[i] * pi;
		sum3 += intfd3d3c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b1*b1+c1*c1))/size;
	case 4:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd3d4a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d4b[i] * pi;
		sum3 += intfd3d4c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b1*b1+c1*c1))/size;
	case 5:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd3d5a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d5b[i] * pi;
		sum3 += intfd3d5c[i] * pi;
		sum4 += intfd3d5d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 6:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd3d6a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d6b[i] * pi;
		sum3 += intfd3d6c[i] * pi;
		sum4 += intfd3d6d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd3d7a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d7b[i] * pi;
		sum3 += intfd3d7c[i] * pi;
		sum4 += intfd3d7d[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd3d8a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d8b[i] * pi;
		sum3 += intfd3d8c[i] * pi;
		sum4 += intfd3d8d[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd3d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd3d9b[i] * pi;
		sum3 += intfd3d9c[i] * pi;
		sum4 += intfd3d9d[i] * pi;
		sum5 += intfd3d9e[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b1*b1+c1*c1) +
		    sum4 * (b0*b2+c0*c2) +
		    sum5 * (b1*b2+c1*c2))/size;
	}
    case 4:
        switch (k) {
	case 4:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd4d4a[i] * (pi = P[Node[i]]);
		sum2 += intfd4d4b[i] * pi;
		sum3 += intfd4d4c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b1*b1+c1*c1))/size;
	case 5:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd4d5a[i] * (pi = P[Node[i]]);
		sum2 += intfd4d5b[i] * pi;
		sum3 += intfd4d5c[i] * pi;
		sum4 += intfd4d5d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 6:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd4d6a[i] * (pi = P[Node[i]]);
		sum2 += intfd4d6b[i] * pi;
		sum3 += intfd4d6c[i] * pi;
		sum4 += intfd4d6d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd4d7a[i] * (pi = P[Node[i]]);
		sum2 += intfd4d7b[i] * pi;
		sum3 += intfd4d7c[i] * pi;
		sum4 += intfd4d7d[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd4d8a[i] * (pi = P[Node[i]]);
		sum2 += intfd4d8b[i] * pi;
		sum3 += intfd4d8c[i] * pi;
		sum4 += intfd4d8d[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd4d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd4d9b[i] * pi;
		sum3 += intfd4d9c[i] * pi;
		sum4 += intfd4d9d[i] * pi;
		sum5 += intfd4d9e[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b1*b1+c1*c1) +
		    sum4 * (b0*b2+c0*c2) +
		    sum5 * (b1*b2+c1*c2))/size;
	}
    case 5:
        switch (k) {
	case 5:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd5d5a[i] * (pi = P[Node[i]]);
		sum2 += intfd5d5b[i] * pi;
		sum3 += intfd5d5c[i] * pi;
	    }
	    return (sum1 * (b1*b1+c1*c1) +
		    sum2 * (b1*b2+c1*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	case 6:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd5d6a[i] * (pi = P[Node[i]]);
		sum2 += intfd5d6b[i] * pi;
		sum3 += intfd5d6c[i] * pi;
	    }
	    return (sum1 * (b1*b1+c1*c1) +
		    sum2 * (b1*b2+c1*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd5d7a[i] * (pi = P[Node[i]]);
		sum2 += intfd5d7b[i] * pi;
		sum3 += intfd5d7c[i] * pi;
		sum4 += intfd5d7d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b1*b2+c1*c2) +
		    sum4 * (b2*b2+c2*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd5d8a[i] * (pi = P[Node[i]]);
		sum2 += intfd5d8b[i] * pi;
		sum3 += intfd5d8c[i] * pi;
		sum4 += intfd5d8d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b1*b2+c1*c2) +
		    sum4 * (b2*b2+c2*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd5d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd5d9b[i] * pi;
		sum3 += intfd5d9c[i] * pi;
		sum4 += intfd5d9d[i] * pi;
		sum5 += intfd5d9e[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2) +
		    sum5 * (b2*b2+c2*c2))/size;
	}
    case 6:
        switch (k) {
	case 6:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd6d6a[i] * (pi = P[Node[i]]);
		sum2 += intfd6d6b[i] * pi;
		sum3 += intfd6d6c[i] * pi;
	    }
	    return (sum1 * (b1*b1+c1*c1) +
		    sum2 * (b1*b2+c1*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	case 7:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd6d7a[i] * (pi = P[Node[i]]);
		sum2 += intfd6d7b[i] * pi;
		sum3 += intfd6d7c[i] * pi;
		sum4 += intfd6d7d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b1*b2+c1*c2) +
		    sum4 * (b2*b2+c2*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2=sum3=sum4 = 0.0; i < 10; i++) {
	        sum1 += intfd6d8a[i] * (pi = P[Node[i]]);
		sum2 += intfd6d8b[i] * pi;
		sum3 += intfd6d8c[i] * pi;
		sum4 += intfd6d8d[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b1*b2+c1*c2) +
		    sum4 * (b2*b2+c2*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd6d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd6d9b[i] * pi;
		sum3 += intfd6d9c[i] * pi;
		sum4 += intfd6d9d[i] * pi;
		sum5 += intfd6d9e[i] * pi;
	    }
	    return (sum1 * (b0*b1+c0*c1) +
		    sum2 * (b1*b1+c1*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2) +
		    sum5 * (b2*b2+c2*c2))/size;
	}
    case 7:
        switch (k) {
	case 7:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd7d7a[i] * (pi = P[Node[i]]);
		sum2 += intfd7d7b[i] * pi;
		sum3 += intfd7d7c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	case 8:
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd7d8a[i] * (pi = P[Node[i]]);
		sum2 += intfd7d8b[i] * pi;
		sum3 += intfd7d8c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd7d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd7d9b[i] * pi;
		sum3 += intfd7d9c[i] * pi;
		sum4 += intfd7d9d[i] * pi;
		sum5 += intfd7d9e[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2) +
		    sum5 * (b2*b2+c2*c2))/size;
	}
    case 8:
        switch (k) {
	case 8:	    
	    for (i = 0, sum1=sum2=sum3 = 0.0; i < 10; i++) {
	        sum1 += intfd8d8a[i] * (pi = P[Node[i]]);
		sum2 += intfd8d8b[i] * pi;
		sum3 += intfd8d8c[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b2+c0*c2) +
		    sum3 * (b2*b2+c2*c2))/size;
	case 9:
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd8d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd8d9b[i] * pi;
		sum3 += intfd8d9c[i] * pi;
		sum4 += intfd8d9d[i] * pi;
		sum5 += intfd8d9e[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b0*b2+c0*c2) +
		    sum4 * (b1*b2+c1*c2) +
		    sum5 * (b2*b2+c2*c2))/size;
	}
    case 9:
        switch (k) {
	case 9: {
	    double sum6 = 0.0;
	    for (i = 0, sum1=sum2=sum3=sum4=sum5 = 0.0; i < 10; i++) {
	        sum1 += intfd9d9a[i] * (pi = P[Node[i]]);
		sum2 += intfd9d9b[i] * pi;
		sum3 += intfd9d9c[i] * pi;
		sum4 += intfd9d9d[i] * pi;
		sum5 += intfd9d9e[i] * pi;
		sum6 += intfd9d9f[i] * pi;
	    }
	    return (sum1 * (b0*b0+c0*c0) +
		    sum2 * (b0*b1+c0*c1) +
		    sum3 * (b1*b1+c1*c1) +
		    sum4 * (b0*b2+c0*c2) +
		    sum5 * (b1*b2+c1*c2) +
		    sum6 * (b2*b2+c2*c2))/size;
	}
	}
    }
    xERROR ("Index out of range");
    return 0.0;
#ifdef UNDEF
    double res = 0.0;
    for (i = 0; i < 10; i++)
        res += P[Node[i]] * IntFDD (i, j, k);
    return res;
#endif
}

static const RDenseMatrix sd0_intf0ff = RDenseMatrix (10, 10,
   "357 20 0 216 -81 0 0 0 0 0 \
    20 20 0 18 18 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    216 18 0 324 -162 0 0 0 0 0 \
    -81 18 0 -162 81 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd0_intf1ff = RDenseMatrix (10, 10,
   "20 20 0 18 18 0 0 0 0 0 \
    20 357 0 -81 216 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    18 -81 0 81 -162 0 0 0 0 0 \
    18 216 0 -162 324 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd0_intf3ff = RDenseMatrix (10, 10,
   "216 18 0 324 -162 0 0 0 0 0 \
    18 -81 0 81 -162 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    324 81 0 2187 0 0 0 0 0 0 \
    -162 -162 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd0_intf4ff = RDenseMatrix (10, 10,
   "-81 18 0 -162 81 0 0 0 0 0 \
    18 216 0 -162 324 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    -162 -162 0 0 0 0 0 0 0 0 \
    81 324 0 0 2187 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix *sd0_intfff[10] = {
    &sd0_intf0ff,
    &sd0_intf1ff,
    0,
    &sd0_intf3ff,
    &sd0_intf4ff,
    0,
    0,
    0,
    0,
    0
};

static const RDenseMatrix sd1_intf1ff = RDenseMatrix (10, 10,
   "0 0 0 0 0 0 0 0 0 0 \
    0 357 20 0 0 216 -81 0 0 0 \
    0 20 20 0 0 18 18 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 216 18 0 0 324 -162 0 0 0 \
    0 -81 18 0 0 -162 81 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd1_intf2ff = RDenseMatrix (10, 10,
   "0 0 0 0 0 0 0 0 0 0 \
    0 20 20 0 0 18 18 0 0 0 \
    0 20 357 0 0 -81 216 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 18 -81 0 0 81 -162 0 0 0 \
    0 18 216 0 0 -162 324 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd1_intf5ff = RDenseMatrix (10, 10,
   "0 0 0 0 0 0 0 0 0 0 \
    0 216 18 0 0 324 -162 0 0 0 \
    0 18 -81 0 0 81 -162 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 324 81 0 0 2187 0 0 0 0 \
    0 -162 -162 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd1_intf6ff = RDenseMatrix (10, 10,
   "0 0 0 0 0 0 0 0 0 0 \
    0 -81 18 0 0 -162 81 0 0 0 \
    0 18 216 0 0 -162 324 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 -162 -162 0 0 0 0 0 0 0 \
    0 81 324 0 0 0 2187 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix *sd1_intfff[10] = {
    0,
    &sd1_intf1ff,
    &sd1_intf2ff,
    0,
    0,
    &sd1_intf5ff,
    &sd1_intf6ff,
    0,
    0,
    0
};

static const RDenseMatrix sd2_intf0ff = RDenseMatrix (10, 10,
   "357 0 20 0 0 0 0 -81 216 0 \
    0 0 0 0 0 0 0 0 0 0 \
    20 0 20 0 0 0 0 18 18 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    -81 0 18 0 0 0 0 81 -162 0 \
    216 0 18 0 0 0 0 -162 324 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd2_intf2ff = RDenseMatrix (10, 10,
   "20 0 20 0 0 0 0 18 18 0 \
    0 0 0 0 0 0 0 0 0 0 \
    20 0 357 0 0 0 0 216 -81 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    18 0 216 0 0 0 0 324 -162 0 \
    18 0 -81 0 0 0 0 -162 81 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd2_intf7ff = RDenseMatrix (10, 10,
   "-81 0 18 0 0 0 0 81 -162 0 \
    0 0 0 0 0 0 0 0 0 0 \
    18 0 216 0 0 0 0 324 -162 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    81 0 324 0 0 0 0 2187 0 0 \
    -162 0 -162 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix sd2_intf8ff = RDenseMatrix (10, 10,
   "216 0 18 0 0 0 0 -162 324 0 \
    0 0 0 0 0 0 0 0 0 0 \
    18 0 -81 0 0 0 0 -162 81 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    0 0 0 0 0 0 0 0 0 0 \
    -162 0 -162 0 0 0 0 0 0 0 \
    324 0 81 0 0 0 0 0 2187 0 \
    0 0 0 0 0 0 0 0 0 0") * (1.0/6720.0);

static const RDenseMatrix *sd2_intfff[10] = {
    &sd2_intf0ff,
    0,
    &sd2_intf2ff,
    0,
    0,
    0,
    0,
    &sd2_intf7ff,
    &sd2_intf8ff,
    0
};

double Triangle10::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    int k, kk;

    if (bndside[0]) {
	const int idx[4] = {0,1,3,4};
        double d0 = sqrt (b2*b2 + c2*c2); // length of side 0
        for (kk = 0; kk < 4; kk++) {
	    k = idx[kk];
	    res += P[Node[k]] * d0 * sd0_intfff[k]->Get(i,j);
	}
    }
    if (bndside[1]) {
        const int idx[4] = {1,2,5,6};
        double d1 = sqrt (b0*b0 + c0*c0); // length of side 1
        for (kk = 0; kk < 4; kk++) {
	    k = idx[kk];
	    res += P[Node[k]] * d1 * sd1_intfff[k]->Get(i,j);
	}
    }
    if (bndside[2]) {
        const int idx[4] = {0,2,7,8};
        double d2 = sqrt (b1*b1 + c1*c1); // length of side 2
        for (kk = 0; kk < 4; kk++) {
	    k = idx[kk];
	    res += P[Node[k]] * d2 * sd2_intfff[k]->Get(i,j);
	}
    }
    return res;
}

double Triangle10::ComputeSize (const NodeList &nlist) const
{
    return 0.5 * (a0+a1+a2);
}

int Triangle10::GetLocalSubsampleAbsc (const Point *&absc) const
{
    absc = absc_sample;
    return nsample_tot;
}

void Triangle10::ComputeIntFD (const NodeList &nlist)
{
    intfd_0.New(10,10);
    intfd_1.New(10,10);
    intfd_0(0,0) = 128.0*b0;          intfd_1(0,0) = 128.0*c0;
    intfd_0(0,1) =  38.0*b1;          intfd_1(0,1) =  38.0*c1;
    intfd_0(0,2) =  38.0*b2;          intfd_1(0,2) =  38.0*c2;
    intfd_0(0,3) = -9.0*b0+198*b1;    intfd_1(0,3) = -9.0*c0+198*c1;
    intfd_0(0,4) = 45.0*b0-72.0*b1;   intfd_1(0,4) = 45.0*c0-72.0*c1;
    intfd_0(0,5) = 45.0*(b1+b2);      intfd_1(0,5) = 45.0*(c1+c2);
    intfd_0(0,6) = 45.0*(b1+b2);      intfd_1(0,6) = 45.0*(c1+c2);
    intfd_0(0,7) = 45.0*b0-72.0*b2;   intfd_1(0,7) = 45.0*c0-72.0*c2;
    intfd_0(0,8) = -9.0*b0+198*b2;    intfd_1(0,8) = -9.0*c0+198*c2;
    intfd_0(0,9) = 54.0*(2*b0+b1+b2); intfd_1(0,9) = 54.0*(2*c0+c1+c2);
    intfd_0(1,0) = 38.0*b0;           intfd_1(1,0) = 38.0*c0;
    intfd_0(1,1) = 128.0*b1;          intfd_1(1,1) = 128.0*c1;
    intfd_0(1,2) = 38.0*b2;           intfd_1(1,2) = 38.0*c2;
    intfd_0(1,3) = -72.0*b0+45.0*b1;  intfd_1(1,3) = -72.0*c0+45.0*c1;
    intfd_0(1,4) = 198*b0-9.0*b1;     intfd_1(1,4) = 198*c0-9.0*c1;
    intfd_0(1,5) = -9.0*b1+198*b2;    intfd_1(1,5) = -9.0*c1+198*c2;
    intfd_0(1,6) = 45.0*b1-72.0*b2;   intfd_1(1,6) = 45.0*c1-72.0*c2;
    intfd_0(1,7) = 45.0*(b0+b2);      intfd_1(1,7) = 45.0*(c0+c2);
    intfd_0(1,8) = 45.0*(b0+b2);      intfd_1(1,8) = 45.0*(c0+c2);
    intfd_0(1,9) = 54.0*(b0+2*b1+b2); intfd_1(1,9) = 54.0*(c0+2*c1+c2);
    intfd_0(2,0) = 38.0*b0;           intfd_1(2,0) = 38.0*c0;
    intfd_0(2,1) = 38.0*b1;           intfd_1(2,1) = 38.0*c1;
    intfd_0(2,2) = 128.0*b2;          intfd_1(2,2) = 128.0*c2;
    intfd_0(2,3) = 45.0*(b0+b1);      intfd_1(2,3) = 45.0*(c0+c1);
    intfd_0(2,4) = 45.0*(b0+b1);      intfd_1(2,4) = 45.0*(c0+c1);
    intfd_0(2,5) = -72.0*b1+45.0*b2;  intfd_1(2,5) = -72.0*c1+45.0*c2;
    intfd_0(2,6) = 198.0*b1-9.0*b2;   intfd_1(2,6) = 198.0*c1-9.0*c2;
    intfd_0(2,7) = 198.0*b0-9.0*b2;   intfd_1(2,7) = 198.0*c0-9.0*c2;
    intfd_0(2,8) = -72.0*b0+45.0*b2;  intfd_1(2,8) = -72.0*c0+45.0*c2;
    intfd_0(2,9) = 54.0*(b0+b1+2*b2); intfd_1(2,9) = 54.0*(c0+c1+2*c2);
    intfd_0(3,0) = 207.0*b0;          intfd_1(3,0) = 207.0*c0;
    intfd_0(3,1) = -117.0*b1;         intfd_1(3,1) = -117.0*c1;
    intfd_0(3,2) = 45.0*b2;           intfd_1(3,2) = 45.0*c2;
    intfd_0(3,3) = 648.0*(b0+b1);     intfd_1(3,3) = 648.0*(c0+c1);
    intfd_0(3,4) = -243.0*b0+81.0*b1; intfd_1(3,4) = -243.0*c0+81.0*c1;
    intfd_0(3,5) = -162.0*b1-243.0*b2;intfd_1(3,5) = -162.0*c1-243.0*c2;
    intfd_0(3,6) = -81.0*b1-162.0*b2; intfd_1(3,6) = -81.0*c1-162.0*c2;
    intfd_0(3,7) = -81.0*b0-243.0*b2; intfd_1(3,7) = -81.0*c0-243.0*c2;
    intfd_0(3,8) = 324.0*b0+648.0*b2; intfd_1(3,8) = 324.0*c0+648.0*c2;
    intfd_0(3,9) = -162.0*(b0-2*b1-4*b2); intfd_1(3,9) = -162.0*(c0-2*c1-4*c2);
    intfd_0(4,0) = -117.0*b0;         intfd_1(4,0) = -117.0*c0;
    intfd_0(4,1) = 207.0*b1;          intfd_1(4,1) = 207.0*c1;
    intfd_0(4,2) = 45.0*b2;           intfd_1(4,2) = 45.0*c2;
    intfd_0(4,3) = 81.0*b0-243.0*b1;  intfd_1(4,3) = 81.0*c0-243.0*c1;
    intfd_0(4,4) = 648.0*(b0+b1);     intfd_1(4,4) = 648.0*(c0+c1);
    intfd_0(4,5) = 324.0*b1+648.0*b2; intfd_1(4,5) = 324.0*c1+648.0*c2;
    intfd_0(4,6) = -81.0*b1-243.0*b2; intfd_1(4,6) = -81.0*c1-243.0*c2;
    intfd_0(4,7) = -81.0*b0-162.0*b2; intfd_1(4,7) = -81.0*c0-162.0*c2;
    intfd_0(4,8) = -162.0*b0-243.0*b2;intfd_1(4,8) = -162.0*c0-243.0*c2;
    intfd_0(4,9) = 162*(2*b0-b1+4*b2);intfd_1(4,9) = 162*(2*c0-c1+4*c2);
    intfd_0(5,0) = 45.0*b0;           intfd_1(5,0) = 45.0*c0;
    intfd_0(5,1) = 207.0*b1;          intfd_1(5,1) = 207.0*c1;
    intfd_0(5,2) = -117.0*b2;         intfd_1(5,2) = -117.0*c2;
    intfd_0(5,3) = -243.0*b0-81.0*b1; intfd_1(5,3) = -243.0*c0-81.0*c1;
    intfd_0(5,4) = 648.0*b0+324.0*b1; intfd_1(5,4) = 648.0*c0+324.0*c1;
    intfd_0(5,5) = 648.0*(b1+b2);     intfd_1(5,5) = 648.0*(c1+c2);
    intfd_0(5,6) = -243.0*b1-81.0*b2; intfd_1(5,6) = -243.0*c1-81.0*c2;
    intfd_0(5,7) = -243.0*b0-162.0*b2;intfd_1(5,7) = -243.0*c0-162.0*c2;
    intfd_0(5,8) = -162.0*b0-81.0*b2; intfd_1(5,8) = -162.0*c0-81.0*c2;
    intfd_0(5,9) = 162.0*(4*b0-b1+2*b2); intfd_1(5,9) = 162.0*(4*c0-c1+2*c2);
    intfd_0(6,0) = 45.0*b0;           intfd_1(6,0) = 45.0*c0;
    intfd_0(6,1) = -117.0*b1;         intfd_1(6,1) = -117.0*c1;
    intfd_0(6,2) = 207.0*b2;          intfd_1(6,2) = 207.0*c2;
    intfd_0(6,3) = -162.0*b0-81.0*b1; intfd_1(6,3) = -162.0*c0-81.0*c1;
    intfd_0(6,4) = -243.0*b0-162.0*b1;intfd_1(6,4) = -243.0*c0-162.0*c1;
    intfd_0(6,5) = 81.0*b1-243.0*b2;  intfd_1(6,5) = 81.0*c1-243.0*c2;
    intfd_0(6,6) = 648.0*(b1+b2);     intfd_1(6,6) = 648.0*(c1+c2);
    intfd_0(6,7) = 648.0*b0+324.0*b2; intfd_1(6,7) = 648.0*c0+324.0*c2;
    intfd_0(6,8) = -243.0*b0-81.0*b2; intfd_1(6,8) = -243.0*c0-81.0*c2;
    intfd_0(6,9) = 162.0*(4*b0+2*b1-b2); intfd_1(6,9) = 162.0*(4*c0+2*c1-c2);
    intfd_0(7,0) = -117.0*b0;         intfd_1(7,0) = -117.0*c0;
    intfd_0(7,1) = 45.0*b1,           intfd_1(7,1) = 45.0*c1;
    intfd_0(7,2) = 207.0*b2;          intfd_1(7,2) = 207.0*c2;
    intfd_0(7,3) = -162.0*b0-243.0*b1;intfd_1(7,3) = -162.0*c0-243.0*c1;
    intfd_0(7,4) = -81.0*b0-162.0*b1; intfd_1(7,4) = -81.0*c0-162.0*c1;
    intfd_0(7,5) = -243.0*b1-81.0*b2; intfd_1(7,5) = -243.0*c1-81.0*c2;
    intfd_0(7,6) = 648.0*b1+324.0*b2; intfd_1(7,6) = 648.0*c1+324.0*c2;
    intfd_0(7,7) = 648.0*(b0+b2);     intfd_1(7,7) = 648.0*(c0+c2);
    intfd_0(7,8) = 81.0*b0-243.0*b2;  intfd_1(7,8) = 81.0*c0-243.0*c2;
    intfd_0(7,9) = 162.0*(2*b0+4*b1-b2); intfd_1(7,9) = 162.0*(2*c0+4*c1-c2);
#ifdef UNDEF
    intfd_0(8,0) = 207.0*b0;          intfd_1(8,0) = 207.0*c0;
    intfd_0(8,1) = 45.0*b1;           intfd_1(8,1) = 45.0*c1;
    intfd_0(8,2) = -117.0*b2;         intfd_1(8,2) = -117.0*c2;
    intfd_0(8,3) = 324.0*b0+648.0*b1; intfd_1(8,3) = 324.0*c0+648.0*c1;
    intfd_0(8,4) = -81.0*b0-243.0*b1; intfd_1(8,4) = -81.0*c0-243.0*c1;
    {\(-81\)\ \((2\ b2 + 
                b3)\), \(-81\)\ \((2\ c2 + c3)\)}, {\(-81\)\ \((3\ b2 + 
                2\ b3)\), \(-81\)\ \((3\ c2 + 2\ c3)\)}, {81\ \((\(-3\)\ b1 + 
                b3)\), 81\ \((\(-3\)\ c1 + c3)\)}, {648\ \((b1 + b3)\), 
          648\ \((c1 + c3)\)}, {\(-162\)\ \((b1 - 
                2\ \((2\ b2 + b3)\))\), \(-162\)\ \((c1 - 
                2\ \((2\ c2 + 
                      c3)\))\)}}, {{\(-54\)\ b1, \(-54\)\ c1}, {\(-54\)\ b2, \
\(-54\)\ c2}, {\(-54\)\ b3, \(-54\)\ c3}, {162\ \((5\ b1 + 2\ b2)\), 
          162\ \((5\ c1 + 2\ c2)\)}, {162\ \((2\ b1 + 5\ b2)\), 
          162\ \((2\ c1 + 5\ c2)\)}, {162\ \((5\ b2 + 2\ b3)\), 
          162\ \((5\ c2 + 2\ c3)\)}, {162\ \((2\ b2 + 5\ b3)\), 
          162\ \((2\ c2 + 5\ c3)\)}, {162\ \((2\ b1 + 5\ b3)\), 
          162\ \((2\ c1 + 5\ c3)\)}, {162\ \((5\ b1 + 2\ b3)\), 
          162\ \((5\ c1 + 2\ c3)\)}, {1944\ \((b1 + b2 + b3)\), 
          1944\ \((c1 + c2 + c3)\)}}}\)], "Output"]
}, Open  ]]
}, Open  ]],
#endif
}

RSymMatrix Triangle10::ComputeIntDD (const NodeList &nlist) const
{
    static RSymMatrix dd(10);
    double bc00, bc01, bc11, bc02, bc12, bc22;
    double scale = 1.0/(320.0*size), scale6, scale27;

    dd(0,0) =  68.0 * scale * (bc00 = b0*b0 + c0*c0);
    dd(1,0) =  14.0 * scale * (bc01 = b0*b1 + c0*c1);
    dd(1,1) =  68.0 * scale * (bc11 = b1*b1 + c1*c1);
    dd(2,0) =  14.0 * scale * (bc02 = b0*b2 + c0*c2);
    dd(2,1) =  14.0 * scale * (bc12 = b1*b2 + c1*c2);
    dd(2,2) =  68.0 * scale * (bc22 = b2*b2 + c2*c2);
    dd(3,0) = (scale6 = 6.0*scale) * (bc00 + 19.0*bc01);
    dd(3,1) = scale6 * (bc11 -  8.0*bc01);
    dd(3,2) =
    dd(4,2) = scale6 * (bc02 + bc12);
    dd(3,3) =
    dd(4,4) = 270.0 * scale * (bc00 + bc01 + bc11);
    dd(4,0) = scale6 * (bc00 -  8.0*bc01);
    dd(4,1) = scale6 * (bc11 + 19.0*bc01);
    dd(4,3) = -54.0 * scale * (bc00 -  2.0*bc01 + bc11);
    dd(5,0) =
    dd(6,0) = scale6 * (bc01 + bc02);
    dd(5,1) = scale6 * (bc11 + 19.0*bc12);
    dd(5,2) = scale6 * (bc22 -  8.0*bc12);
    dd(5,3) =
    dd(6,3) =
    dd(6,4) = (scale27 = -27.0 * scale) * (bc01 + bc11 + bc12 + 2.0*bc02);
    dd(5,4) = 135.0 * scale * (bc11 + bc12 + bc01 + 2.0*bc02);
    dd(5,5) =
    dd(6,6) = 270.0 * scale * (bc11 + bc12 + bc22);
    dd(6,1) = scale6 * (bc11 -  8.0*bc12);
    dd(6,2) = scale6 * (bc22 + 19.0*bc12);
    dd(6,5) = -54.0 * scale * (bc11 - 2.0*bc12 + bc22);
    dd(7,0) = scale6 * (bc00 -  8.0*bc02);
    dd(7,1) =
    dd(8,1) = scale6 * (bc01 + bc12);
    dd(7,2) = scale6 * (bc22 + 19.0*bc02);
    dd(7,3) =
    dd(7,4) =
    dd(8,4) = scale27 * (bc00 + bc01 + bc02 + 2.0*bc12);
    dd(7,5) = scale27 * (bc12 + bc22 + bc02 + 2.0*bc01);
    dd(7,6) = 135.0 * scale * (bc12 + bc22 + bc02 + 2.0*bc01);
    dd(7,7) =
    dd(8,8) = 270.0 * scale * (bc00 + bc02 + bc22);
    dd(8,0) = scale6 * (bc00 + 19.0*bc02);
    dd(8,2) = scale6 * (bc22 -  8.0*bc02);
    dd(8,3) = 135.0 * scale * (bc00 + bc01 + bc02 + 2.0*bc12);
    dd(8,5) =
    dd(8,6) = scale27 * (bc12 + bc22 + bc02 + 2.0*bc01);
    dd(8,7) = -54.0 * scale * (bc00 + bc22 - 2.0*bc02);
    // dd(9,0) =  18.0 * scale * (bc00 + bc01 + bc02); // = 0
    // dd(9,1) =  18.0 * scale * (bc01 + bc11 + bc12); // = 0
    // dd(9,2) =  18.0 * scale * (bc02 + bc12 + bc22); // = 0
    dd(9,3) =
    dd(9,6) = 162.0 * scale * (bc11 + bc12 + bc01 + 2.0*bc02);
    dd(9,4) =
    dd(9,7) = 162.0 * scale * (bc00 + bc01 + bc02 + 2.0*bc12);
    dd(9,5) =
    dd(9,8) = 162.0 * scale * (bc12 + bc22 + bc02 + 2.0*bc01);
    dd(9,9) = 648.0 * scale * (bc00 + bc01 + bc11 + bc02 + bc12 + bc22);

    return dd;
}

double Triangle10::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 10, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 2: out of range");

    double f = bndintf(sd,i);
    if (!f) return 0.0;
    switch (sd) {
    case 0: return f * hypot (c2, b2);
    case 1: return f * hypot (c0, b0);
    case 2: return f * hypot (c1, b1);
    default: xERROR("Argument 2: out of range"); return 0.0;
    }
}

double Triangle10::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 10, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 10, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 3, "Argument 3: out of range");
    
    switch (sd) {
    case 0:
	return hypot (c2, b2) * sym_bndintff_sd0(i,j);
    case 1:
	return hypot (c0, b0) * sym_bndintff_sd1(i,j);
    case 2:
	return hypot (c1, b1) * sym_bndintff_sd2(i,j);
    default:
	xERROR("Argument 3: out of range");
	return 0.0;
    }
}

RSymMatrix Triangle10::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(10);
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

#ifdef UNDEF
double Triangle10::Jacobian (const Point &loc, RDenseMatrix &J) const
{
#ifndef TRI10_STORE_COORDS
    xERROR("Requires definition of TRI10_STORE_COORDS");
    double n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y,
           n5x, n5y, n6x, n6y, n7x, n7y, n8x, n8y, n9x, n9y;
#endif
    dASSERT(J.nRows() == 2 && J.nCols() == 2, "Parameter 3 wrong dimension");
    RDenseMatrix der = LocalShapeD (loc);

    for (int i = 0; i < 2; i++) {
        J(i,0) = der(i,0)*n0x + der(i,1)*n1x + der(i,2)*n2x +
	         der(i,3)*n3x + der(i,4)*n4x + der(i,5)*n5x +
	         der(i,6)*n6x + der(i,7)*n7x + der(i,8)*n8x +
                 der(i,9)*n9x;
	J(i,1) = der(i,0)*n0y + der(i,1)*n1y + der(i,2)*n2y +
                 der(i,3)*n3y + der(i,4)*n4y + der(i,5)*n5y +
                 der(i,6)*n6y + der(i,7)*n7y + der(i,8)*n8y +
                 der(i,9)*n9y;
    }

    return J(0,0)*J(1,1) - J(1,0)*J(0,1);
}

double Triangle10::IJacobian (const Point &loc, RDenseMatrix &IJ) const
{
    double d = Jacobian (loc, IJ);
    double id = 1.0/d;
    double tmp = IJ(0,0);
    IJ(0,0) = IJ(1,1)*id;
    IJ(0,1) *= -id;
    IJ(1,0) *= -id;
    IJ(1,1) = tmp*id;
    return d;
}

double Triangle10::DetJ (const Point &loc, const NodeList *nlist) const
{
#ifndef TRI10_STORE_COORDS
    dASSERT (nlist != 0, "Node list required");
    double n0x = (*nlist)[Node[0]][0], n0y = (*nlist)[Node[0]][1];
    double n1x = (*nlist)[Node[1]][0], n1y = (*nlist)[Node[1]][1];
    double n2x = (*nlist)[Node[2]][0], n2y = (*nlist)[Node[2]][1];
    double n3x = (*nlist)[Node[3]][0], n3y = (*nlist)[Node[3]][1];
    double n4x = (*nlist)[Node[4]][0], n4y = (*nlist)[Node[4]][1];
    double n5x = (*nlist)[Node[5]][0], n5y = (*nlist)[Node[5]][1];
    double n6x = (*nlist)[Node[6]][0], n6y = (*nlist)[Node[6]][1];
    double n7x = (*nlist)[Node[7]][0], n7y = (*nlist)[Node[7]][1];
    double n8x = (*nlist)[Node[8]][0], n8y = (*nlist)[Node[8]][1];
    double n9x = (*nlist)[Node[9]][0], n9y = (*nlist)[Node[9]][1];
#endif
    dASSERT (loc.Dim() == 2, "Parameter 1 wrong dimension");
    RDenseMatrix der = LocalShapeD (loc);
    double j00, j01, j10, j11;

    j00 = der[0][0]*n0x + der[0][1]*n1x + der[0][2]*n2x +
          der[0][3]*n3x + der[0][4]*n4x + der[0][5]*n5x +
          der[0][6]*n6x + der[0][7]*n7x + der[0][8]*n8x +
          der[0][9]*n9x;
    j01 = der[0][0]*n0y + der[0][1]*n1y + der[0][2]*n2y +
          der[0][3]*n3y + der[0][4]*n4y + der[0][5]*n5y +
          der[0][6]*n6y + der[0][7]*n7y + der[0][8]*n8y +
          der[0][9]*n9y;
    j10 = der[1][0]*n0x + der[1][1]*n1x + der[1][2]*n2x +
          der[1][3]*n3x + der[1][4]*n4x + der[1][5]*n5x +
          der[1][6]*n6x + der[1][7]*n7x + der[1][8]*n8x +
          der[1][9]*n9x;
    j11 = der[1][0]*n0y + der[1][1]*n1y + der[1][2]*n2y +
          der[1][3]*n3y + der[1][4]*n4y + der[1][5]*n5y +
          der[1][6]*n6y + der[1][7]*n7y + der[1][8]*n8y +
          der[1][9]*n9y;

    return j00*j11 - j10*j01;
}

#endif
