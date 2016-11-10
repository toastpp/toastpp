// ==========================================================================
// Module libfe
// File tet4.cc
// Definition of class Tetrahedron4
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "toastdef.h"
#include "tri_qr.h"

// general utility functions

double TriangleArea (const Point &p1, const Point &p2, const Point &p3)
{
    // returns the area of a triangle given by three vertex points in space
    dASSERT(p1.Dim() == 3 && p2.Dim() == 3 && p3.Dim() == 3,
	   "Wrong vertex dimension");
    double a = p1.Dist(p2);
    double b = p2.Dist(p3);
    double c = p3.Dist(p1);
    double s = 0.5 * (a+b+c);
    return sqrt (s * (s-a) * (s-b) * (s-c));
}

// some global constants

static bool subsampling_initialised = false;
static const int nsample_lin = NSUBSAMPLE; // from toastdef.h
static const int nsample_tot = (((nsample_lin+3)*nsample_lin+2)*nsample_lin)/6;
static Point absc_sample[nsample_tot];

static const RSymMatrix sym_intff = RSymMatrix (4,
   "2 \
    1 2 \
    1 1 2 \
    1 1 1 2") * (1.0/20.0);

static const RDenseMatrix full_intff = RDenseMatrix (4, 4,
   "2 1 1 1 \
    1 2 1 1 \
    1 1 2 1 \
    1 1 1 2") * (1.0/20.0);

// below are the 4 planes of the IntFFF tensor
static const RDenseMatrix full_intf0ff = RDenseMatrix (4, 4,
   "6 2 2 2 \
    2 2 1 1 \
    2 1 2 1 \
    2 1 1 2") / 120.0;
static const RDenseMatrix full_intf1ff = RDenseMatrix (4, 4,
   "2 2 1 1 \
    2 6 2 2 \
    1 2 2 1 \
    1 2 1 2") / 120.0;
static const RDenseMatrix full_intf2ff = RDenseMatrix (4, 4,
   "2 1 2 1 \
    1 2 2 1 \
    2 2 6 2 \
    1 1 2 2") / 120.0;
static const RDenseMatrix full_intf3ff = RDenseMatrix (4, 4,
   "2 1 1 2 \
    1 2 1 2 \
    1 1 2 2 \
    2 2 2 6") / 120.0;

static const RDenseMatrix *full_intfff[4] = {
    &full_intf0ff,
    &full_intf1ff,
    &full_intf2ff,
    &full_intf3ff
};

// the same in symmetric matrix representation
static const RSymMatrix sym_intf0ff = RSymMatrix (4,
   "6 \
    2 2 \
    2 1 2 \
    2 1 1 2") * (1.0/120.0);
static const RSymMatrix sym_intf1ff = RSymMatrix (4,
   "2 \
    2 6 \
    1 2 2 \
    1 2 1 2") * (1.0/120.0);
static const RSymMatrix sym_intf2ff = RSymMatrix (4,
   "2 \
    1 2 \
    2 2 6 \
    1 1 2 2") * (1.0/120.0);
static const RSymMatrix sym_intf3ff = RSymMatrix (4,
   "2 \
    1 2 \
    1 1 2 \
    2 2 2 6") * (1.0/120.0);

Tetrahedron4::Tetrahedron4 (const Tetrahedron4 &el)
: Element_Unstructured_3D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Tetrahedron4::Copy ()
{
    return new Tetrahedron4(*this);
}

void Tetrahedron4::Initialise (const NodeList &nlist)
{
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

#ifdef TET4_STORE_INTFF
    intff.New(4);
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

int Tetrahedron4::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 4, "Side index out of range");
    dASSERT(node >= 0 && node < 3, "Node index out of range");
    static const int SN[4][3] = {{0,1,2},{0,3,1},{0,2,3},{1,3,2}};
    return SN[side][node];
}

double Tetrahedron4::SideSize (int sd, const NodeList &nlist) const
{
    return side_size[sd];
}

Point Tetrahedron4::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");

    Point loc(3);
    double scale = 1.0/(a0+a1+a2+a3);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1] + d1*glob[2]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1] + d2*glob[2]) * scale;
    loc[2] = (a3 + b3*glob[0] + c3*glob[1] + d3*glob[2]) * scale;
    return loc;
}

Point Tetrahedron4::NodeLocal (int node) const
{
    Point nloc(3);
    switch (node) {
    case 0: nloc[0] = nloc[1] = nloc[2] = 0.0; break;
    case 1: nloc[0] = 1.0; nloc[1] = nloc[2] = 0.0; break;
    case 2: nloc[0] = nloc[2] = 0.0; nloc[1] = 1.0; break;
    case 3: nloc[0] = nloc[1] = 0.0; nloc[2] = 1.0; break;
    default: xERROR("Node index out of range");
    }
    return nloc;
}

Point Tetrahedron4::SurfToLocal (int side, const Point &p) const
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

void Tetrahedron4::MapToSide (int side, Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Arg 2 wrong vector dimension");
    switch (side) {
    case 0: loc[2] = 0.0; break;
    case 1: loc[1] = 0.0; break;
    case 2: loc[0] = 0.0; break;
    case 3: { double a = (1.0-loc[0]-loc[1]-loc[2])/3.0;
              loc[0] += a, loc[1] += a, loc[2] += a;
            } break;
    default: xERROR ("Side index out of range");
    }
}

RVector Tetrahedron4::DirectionCosine (int side, RDenseMatrix &jacin)
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

const RVector &Tetrahedron4::LNormal (int side) const
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

bool Tetrahedron4::LContains (const Point &loc, bool pad) const
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

RVector Tetrahedron4::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Invalid point dimension");
    RVector fun(4);
    fun[0] = 1.0-loc[0]-loc[1]-loc[2];
    fun[1] = loc[0];
    fun[2] = loc[1];
    fun[3] = loc[2];
    return fun;
}

RDenseMatrix Tetrahedron4::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, "Invalid point dimension");
    static const RDenseMatrix der(3, 4,
       "-1 1 0 0 \
        -1 0 1 0 \
        -1 0 0 1");
    return der;
}

RVector Tetrahedron4::GlobalShapeF (const NodeList &nlist, const Point &glob)
    const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");
    RVector fun(4);
    double scale = 1.0/(6.0*size);
    fun[0] = scale * (a0 + b0*glob[0] + c0*glob[1] + d0*glob[2]);
    fun[1] = scale * (a1 + b1*glob[0] + c1*glob[1] + d1*glob[2]);
    fun[2] = scale * (a2 + b2*glob[0] + c2*glob[1] + d2*glob[2]);
    fun[3] = scale * (a3 + b3*glob[0] + c3*glob[1] + d3*glob[2]);
    return fun;
}

RDenseMatrix Tetrahedron4::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");
    RDenseMatrix der(3,4);
    double scale = 1.0/(6.0*size);
    der(0,0) = b0*scale;
    der(1,0) = c0*scale;
    der(2,0) = d0*scale;
    der(0,1) = b1*scale;
    der(1,1) = c1*scale;
    der(2,1) = d1*scale;
    der(0,2) = b2*scale;
    der(1,2) = c2*scale;
    der(2,2) = d2*scale;
    der(0,3) = b3*scale;
    der(1,3) = c3*scale;
    der(2,3) = d3*scale;
    return der;
}

double Tetrahedron4::IntF (int i) const {
    return 0.25*size;
}

RSymMatrix Tetrahedron4::IntFF () const {
#ifdef TET4_STORE_INTFF
    return intff;
#else
    return sym_intff * size;
#endif
}

double Tetrahedron4::IntFF (int i, int j) const {
    RANGE_CHECK(i >= 0 && i < 4 && j >= 0 && j < 4);
#ifdef TET4_STORE_INTFF
    return intff(i,j);
#else
    return full_intff(i,j) * size;
#endif
}

double Tetrahedron4::IntFFF (int i, int j, int k) const
{
    RANGE_CHECK(i >= 0 && i < 4 && j >= 0 && j < 4 && k >= 0 && k < 4);
    return full_intfff[i]->Get(j,k) * size;
}

RSymMatrix Tetrahedron4::IntPFF (const RVector &P) const
{
    return (sym_intf0ff * P[Node[0]] +
	    sym_intf1ff * P[Node[1]] +
	    sym_intf2ff * P[Node[2]] +
	    sym_intf3ff * P[Node[3]]) * size;
}

double Tetrahedron4::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 4 && j >= 0 && j < 4);
    return (full_intf0ff(i,j) * P[Node[0]] +
	    full_intf1ff(i,j) * P[Node[1]] +
	    full_intf2ff(i,j) * P[Node[2]] +
	    full_intf3ff(i,j) * P[Node[3]]) * size;
}

double Tetrahedron4::IntFDD (int i, int j, int k) const
{
    // Note the result is independent of i (the index for 'F')

    if (j > k) { int tmp=j; j=k; k=tmp; } // symmetry
    double scale = 1.0/(144.0*size);
    switch (j) {
    case 0:
        switch (k) {
	case 0: return (b0*b0+c0*c0+d0*d0) * scale;
	case 1: return (b0*b1+c0*c1+d0*d1) * scale;
	case 2: return (b0*b2+c0*c2+d0*d2) * scale;
	case 3: return (b0*b3+c0*c3+d0*d3) * scale;
	}
    case 1:
        switch (k) {
	case 1: return (b1*b1+c1*c1+d1*d1) * scale;
	case 2: return (b1*b2+c1*c2+d1*d2) * scale;
	case 3: return (b1*b3+c1*c3+d1*d3) * scale;
	}
    case 2:
        switch (k) {
	case 2: return (b2*b2+c2*c2+d2*d2) * scale;
	case 3: return (b2*b3+c2*c3+d2*d3) * scale;
	}
    case 3:
        switch (k) {
	case 3: return (b3*b3+c3*c3+d3*d3) * scale;
	}
    }
    xERROR("Index out of range");
    return 0;
}

double Tetrahedron4::IntPDD (int i, int j, const RVector &P) const
{
    double res = 0.0;
    for (int k = 0; k < 4; k++) res += P[Node[k]];
    return res * IntFDD (0,i,j); // IntFDD(k,i,j) is independent of k
}

RSymMatrix Tetrahedron4::IntPDD (const RVector &P) const
{
    RSymMatrix pdd(4);
    double res = 0.0;
    for (int k = 0; k < 4; k++)
	res += P[Node[k]]; // IntFDD(k,i,j) is independent of k
    for (int i = 0; i < 4; i++) {
	for (int j = 0; j <= i; j++) {
	    pdd(i,j) = res * IntFDD (0,i,j);
	}
    }
    return pdd;
}

//Mixed Elements
double Tetrahedron4::IntFd (int i, int j, int k) const
{
    // note that the result is independent of index i.
    double scale = 1.0/24.0;
    switch (j) {
    case 0:
	switch (k) {
	case 0: return b0*scale;
	case 1: return c0*scale;
	case 2: return d0*scale;
	}
	break;
    case 1:
	switch (k) {
	case 0: return b1*scale;
	case 1: return c1*scale;
	case 2: return d1*scale;
	}
	break;
    case 2:
	switch (k) {
	case 0: return b2*scale;
	case 1: return c2*scale;
	case 2: return d2*scale;
	}
	break;
    case 3:
	switch (k) {
	case 0: return b3*scale;
	case 1: return c3*scale;
	case 2: return d3*scale;
	}
	break;
    }
    xERROR("Index out of range");
    return 0.0;
}

double Tetrahedron4::IntPd (const RVector &P, int j, int k) const
{
    // since the result of Int u_i du_j/dx_k is independent of
    // index i, we can pull the integral out of the sum over nodes

    double sum = 0.0;
    for (int i = 0; i < 4; i++)
	sum += P[Node[i]];
    return sum * IntFd (0, j, k);
}

double Tetrahedron4::IntFdd (int i, int j, int k, int l, int m) const
{
//Note as above we have independence of index i;

  if (j > k) { int tmp=j; j=k; k=tmp; tmp=l; l=m; m=tmp;  } // symmetry
    double scale = 1.0/(144.0*size);

 switch (j) {
    case 0:
        switch (k) {
	case 0: { if (l > m) { int tmp=l; l=m; m=tmp; } // symmetry
		 if(l==0&&m==0) return b0*b0*scale;
		 if(l==0&&m==1) return b0*c0*scale;
 		 if(l==0&&m==2) return b0*d0*scale; 
                 if(l==1&&m==1) return c0*c0*scale;
		 if(l==1&&m==2) return c0*d0*scale;
		 if(l==2&&m==2) return d0*d0*scale;
		 break;}  

	case 1: {if(l==0&&m==0) return b0*b1*scale;
		 if(l==0&&m==1) return b1*c0*scale;
 		 if(l==0&&m==2) return b1*d0*scale; 	
		 if(l==1&&m==0) return b0*c1*scale;
		 if(l==1&&m==1) return c0*c1*scale;
		 if(l==1&&m==2) return c1*d0*scale;
		 if(l==2&&m==0) return b0*d1*scale;
		 if(l==2&&m==1) return c0*d1*scale;
		 if(l==2&&m==2) return d0*d1*scale;
		 break;}  

	case 2: {if(l==0&&m==0) return b0*b2*scale;
		 if(l==0&&m==1) return b2*c0*scale;
 		 if(l==0&&m==2) return b2*d0*scale; 	
		 if(l==1&&m==0) return b0*c2*scale;
		 if(l==1&&m==1) return c0*c2*scale;
		 if(l==1&&m==2) return c2*d0*scale;
		 if(l==2&&m==0) return b0*d2*scale;
		 if(l==2&&m==1) return c0*d2*scale;
		 if(l==2&&m==2) return d0*d2*scale;
                 break;}   

	case 3: {if(l==0&&m==0) return b0*b3*scale;
		 if(l==0&&m==1) return b3*c0*scale;
 		 if(l==0&&m==2) return b3*d0*scale; 	
		 if(l==1&&m==0) return b0*c3*scale;
		 if(l==1&&m==1) return c0*c3*scale;
		 if(l==1&&m==2) return c3*d0*scale;
		 if(l==2&&m==0) return b0*d3*scale;
		 if(l==2&&m==1) return c0*d3*scale;
		 if(l==2&&m==2) return d0*d3*scale;
		  break;} break; 
	}
    case 1:
        switch (k) {
	case 1: { if (l > m) { int tmp=l; l=m; m=tmp; } // symmetry
		 if(l==0&&m==0) return b1*b1*scale;
		 if(l==0&&m==1) return b1*c1*scale;
 		 if(l==0&&m==2) return b1*d1*scale; 	
		 if(l==1&&m==1) return c1*c1*scale;
		 if(l==1&&m==2) return c1*d1*scale;
		 if(l==2&&m==2) return d1*d1*scale;
		 break;}   

	case 2: {if(l==0&&m==0) return b1*b2*scale;
		 if(l==0&&m==1) return b2*c1*scale;
 		 if(l==0&&m==2) return b2*d1*scale; 	
		 if(l==1&&m==0) return b1*c2*scale;
		 if(l==1&&m==1) return c1*c2*scale;
		 if(l==1&&m==2) return c2*d1*scale;
		 if(l==2&&m==0) return b1*d2*scale;
		 if(l==2&&m==1) return c1*d2*scale;
		 if(l==2&&m==2) return d1*d2*scale;
		 break;}   

	case 3: {if(l==0&&m==0) return b1*b3*scale;
		 if(l==0&&m==1) return b3*c1*scale;
 		 if(l==0&&m==2) return b3*d1*scale; 	
		 if(l==1&&m==0) return b1*c3*scale;
		 if(l==1&&m==1) return c1*c3*scale;
		 if(l==1&&m==2) return c3*d1*scale;
		 if(l==2&&m==0) return b1*d3*scale;
		 if(l==2&&m==1) return c1*d3*scale;
		 if(l==2&&m==2) return d1*d3*scale;
		 break; } break;   
	}
    case 2:
        switch (k) {
	case 2:{ if (l > m) { int tmp=l; l=m; m=tmp; } // symmetry
		 if(l==0&&m==0) return b2*b2*scale;
		 if(l==0&&m==1) return b2*c2*scale;
 		 if(l==0&&m==2) return b2*d2*scale; 	
		 if(l==1&&m==1) return c2*c2*scale;
		 if(l==1&&m==2) return c2*d2*scale;
		 if(l==2&&m==2) return d2*d2*scale;
		 break;}  

	case 3: {if(l==0&&m==0) return b2*b3*scale;
		 if(l==0&&m==1) return b3*c2*scale;
 		 if(l==0&&m==2) return b3*d2*scale; 	
		 if(l==1&&m==0) return b2*c3*scale;
		 if(l==1&&m==1) return c2*c3*scale;
		 if(l==1&&m==2) return c3*d2*scale;
		 if(l==2&&m==0) return b2*d3*scale;
		 if(l==2&&m==1) return c2*d3*scale;
		 if(l==2&&m==2) return d2*d3*scale;
		 break;}   
	}
    case 3:
        switch (k) {
	case 3: { if (l > m) { int tmp=l; l=m; m=tmp; } // symmetry
		 if(l==0&&m==0) return b3*b3*scale;
		 if(l==0&&m==1) return b3*c3*scale;
 		 if(l==0&&m==2) return b3*d3*scale; 	
		 if(l==1&&m==1) return c3*c3*scale;
		 if(l==1&&m==2) return c3*d3*scale;
		 if(l==2&&m==2) return d3*d3*scale;
		 break;}  break; 
	}
    }
    xERROR("Index out of range");
    return 0;
}

double Tetrahedron4::IntPdd (const RVector &P, int j, int k, int l, int m) const
{
    double val = 0;
    for (int i = 0; i < 4; i++)
	val += IntFdd (i,j,k,l,m) * P[Node[i]];
    return val;
}

RSymMatrix Tetrahedron4::Intdd () const
{
    double dd[12];
    dd[0] = b0; dd[1] = c0; dd[2] = d0;
    dd[3] = b1; dd[4] = c1; dd[5] = d1;
    dd[6] = b2; dd[7] = c2; dd[8] = d2;
    dd[9] = b3; dd[10]= c3; dd[11]= d3;

    RSymMatrix MDD(12,12);
    int i, j;
    double scale = 1.0/(36.0*size);
    for (i = 0; i < 12; i++)
	for (j = 0; j <= i; j++)
	    MDD(i,j) = dd[i]*dd[j] * scale;


    /*for (i = 0; i < 12; i++)
	std::cout << dd[i] << "\t";
    std::cout << std::endl;
    std::cout << 36*size << std::endl;
*/
    return MDD;
}

double Tetrahedron4::Intd (int i, int k) const
{
    const double fac = 1.0/6.0;
    switch(i) {
    case 0:
	switch(k) {
	case 0: return b0*fac;
	case 1: return c0*fac;
	case 2: return d0*fac;
	default: xERROR("Invalid index"); return 0;
	}
    case 1:
	switch(k) {
	case 0: return b1*fac;
	case 1: return c1*fac;
	case 2: return d1*fac;
	default: xERROR("Invalid index"); return 0;
	}
    case 2:
	switch(k) {
	case 0: return b2*fac;
	case 1: return c2*fac;
	case 2: return d2*fac;
	default: xERROR("Invalid index"); return 0;
	}
    case 3:
	switch(k) {
	case 0: return b3*fac;
	case 1: return c3*fac;
	case 2: return d3*fac;
	default: xERROR("Invalid index"); return 0;
	}
    default: xERROR("Invalid index"); return 0;
    }
}

double Tetrahedron4::IntFfd (int i, int j, int k, int l) const
{
    //Coefficient
    double coeff=120.0; 
    if(i==j) coeff=60.0;  

    switch (k) {
    case 0:
	switch (l) {
	case 0: return b0/coeff;
	case 1: return c0/coeff;
	case 2: return d0/coeff;
	default: xERROR("Invalid index"); return 0;
	}
    case 1:
	switch (l) {
	case 0: return b1/coeff;
	case 1: return c1/coeff;
	case 2: return d1/coeff;
	default: xERROR("Invalid index"); return 0;
	}
    case 2:
	switch (l) {
	case 0: return b2/coeff;
	case 1: return c2/coeff;
	case 2: return d2/coeff;
	default: xERROR("Invalid index"); return 0;
	}
    case 3:
	switch (l) {
	case 0: return b3/coeff;
	case 1: return c3/coeff;
	case 2: return d3/coeff;
	default: xERROR("Invalid index"); return 0;
	}
    default:
	xERROR("Invalid index"); return 0;
    }
}

double Tetrahedron4::IntPfd (const RVector &P, int j, int k, int l) const
{
    double val = 0;
    for (int i = 0; i < 4; i++)
	val += IntFfd (i,j,k,l) * P[Node[i]];
    return val;
}

static const RDenseMatrix sd0_intf0ff = RDenseMatrix (4, 4,
   "6 2 2 0 \
    2 2 1 0 \
    2 1 2 0 \
    0 0 0 0") * (2.0/120.0);
static const RDenseMatrix sd0_intf1ff = RDenseMatrix (4, 4,
   "2 2 1 0 \
    2 6 2 0 \
    1 2 2 0 \
    0 0 0 0") * (2.0/120.0);
static const RDenseMatrix sd0_intf2ff = RDenseMatrix (4, 4,
   "2 1 2 0 \
    1 2 2 0 \
    2 2 6 0 \
    0 0 0 0") * (2.0/120.0);
static const RDenseMatrix sd0_intf3ff = RDenseMatrix (4, 4); // = 0
static const RDenseMatrix *sd0_intfff[4] = {
    &sd0_intf0ff,
    &sd0_intf1ff,
    &sd0_intf2ff,
    &sd0_intf3ff
};
static const RDenseMatrix sd1_intf0ff = RDenseMatrix (4, 4,
   "6 2 0 2 \
    2 2 0 1 \
    0 0 0 0 \
    2 1 0 2") * (2.0/120.0);
static const RDenseMatrix sd1_intf1ff = RDenseMatrix (4, 4,
   "2 2 0 1 \
    2 6 0 2 \
    0 0 0 0 \
    1 2 0 2") * (2.0/120.0);
static const RDenseMatrix sd1_intf2ff = RDenseMatrix (4, 4); // = 0
static const RDenseMatrix sd1_intf3ff = RDenseMatrix (4, 4,
   "2 1 0 2 \
    1 2 0 2 \
    0 0 0 0 \
    2 2 0 6") * (2.0/120.0);
static const RDenseMatrix *sd1_intfff[4] = {
    &sd1_intf0ff,
    &sd1_intf1ff,
    &sd1_intf2ff,
    &sd1_intf3ff
};
static const RDenseMatrix sd2_intf0ff = RDenseMatrix (4, 4,
   "6 0 2 2 \
    0 0 0 0 \
    2 0 2 1 \
    2 0 1 2") * (2.0/120.0);
static const RDenseMatrix sd2_intf1ff = RDenseMatrix (4, 4); // = 0
static const RDenseMatrix sd2_intf2ff = RDenseMatrix (4, 4,
   "2 0 2 1 \
    0 0 0 0 \
    2 0 6 2 \
    1 0 2 2") * (2.0/120.0);
static const RDenseMatrix sd2_intf3ff = RDenseMatrix (4, 4,
   "2 0 1 2 \
    0 0 0 0 \
    1 0 2 2 \
    2 0 2 6") * (2.0/120.0);
static const RDenseMatrix *sd2_intfff[4] = {
    &sd2_intf0ff,
    &sd2_intf1ff,
    &sd2_intf2ff,
    &sd2_intf3ff
};
static const RDenseMatrix sd3_intf0ff = RDenseMatrix (4, 4); // = 0
static const RDenseMatrix sd3_intf1ff = RDenseMatrix (4, 4,
   "0 0 0 0 \
    0 6 2 2 \
    0 2 2 1 \
    0 2 1 2") * (2.0/120.0);
static const RDenseMatrix sd3_intf2ff = RDenseMatrix (4, 4,
   "0 0 0 0 \
    0 2 2 1 \
    0 2 6 2 \
    0 1 2 2") * (2.0/120.0);
static const RDenseMatrix sd3_intf3ff = RDenseMatrix (4, 4,
   "0 0 0 0 \
    0 2 1 2 \
    0 1 2 2 \
    0 2 2 6") * (2.0/120.0);
static const RDenseMatrix *sd3_intfff[4] = {
    &sd3_intf0ff,
    &sd3_intf1ff,
    &sd3_intf2ff,
    &sd3_intf3ff
};
static const RDenseMatrix **bndintfff[4] = {
    sd0_intfff,
    sd1_intfff,
    sd2_intfff,
    sd3_intfff
};

double Tetrahedron4::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    for (int sd = 0; sd < 4; sd++) {
        if (bndside[sd]) {
	    double sres = 0.0;
	    for (int k = 0; k < 4; k++)
	        sres += P[Node[k]] * bndintfff[sd][k]->Get(i,j);
	    res += sres * side_size[sd];
	}
    }
    return res;
}

int Tetrahedron4::GetLocalSubsampleAbsc (const Point *&absc) const
{
    absc = absc_sample;
    return nsample_tot;
}

int Tetrahedron4::GetBndSubsampleAbsc (int side, const Point *&absc) const
{
    extern int Triangle_GetLocalSubsampleAbsc (const Point *&);
    return Triangle_GetLocalSubsampleAbsc (absc);
}

RDenseMatrix Tetrahedron4::StrainDisplacementMatrix (const Point& /*glob*/)
    const
{
    RDenseMatrix B(6,12);
    B(0,0) = b0;  B(0,3) = b1;  B(0,6) = b2;  B(0, 9) = b3;
    B(1,1) = c0;  B(1,4) = c1;  B(1,7) = c2;  B(1,10) = c3;
    B(2,2) = d0;  B(2,5) = d1;  B(2,8) = d2;  B(2,11) = d3;
    B(3,0) = c0;  B(3,3) = c1;  B(3,6) = c2;  B(3, 9) = c3;  
    B(3,1) = b0;  B(3,4) = b1;  B(3,7) = b2;  B(3,10) = b3;
    B(4,1) = d0;  B(4,4) = d1;  B(4,7) = d2;  B(4,10) = d3;
    B(4,2) = c0;  B(4,5) = c1;  B(4,8) = c2;  B(4,11) = c3;
    B(5,0) = d0;  B(5,3) = d1;  B(5,6) = d2;  B(5, 9) = d3;
    B(5,2) = b0;  B(5,5) = b1;  B(5,8) = b2;  B(5,11) = b3;
    return B;
}

RDenseMatrix Tetrahedron4::ElasticityStiffnessMatrix (const RDenseMatrix &D)
    const
{
    // This forms a 12x12 element stiffness matrix, built as 4x4 blocks, where
    // each block ij is a 3x3 submatrix: K_ij = B_i^T D B_j V_e \in R^{3x3}
    // Input matrix D is a 6x6 elasticity matrix for the element

    int i;
    double scale = 1.0/(6.0*size);

    // First form 12x6 matrix P = B^T D
    // (note we only sum over nozero elements of B)
    static RDenseMatrix P(12,6);
    for (i = 0; i < 6; i++) {
        P( 0,i) = (b0*D(0,i) + c0*D(3,i) + d0*D(5,i)) * scale;
	P( 1,i) = (c0*D(1,i) + b0*D(3,i) + d0*D(4,i)) * scale;
	P( 2,i) = (d0*D(2,i) + c0*D(4,i) + b0*D(5,i)) * scale;
	P( 3,i) = (b1*D(0,i) + c1*D(3,i) + d1*D(5,i)) * scale;
	P( 4,i) = (c1*D(1,i) + b1*D(3,i) + d1*D(4,i)) * scale;
	P( 5,i) = (d1*D(2,i) + c1*D(4,i) + b1*D(5,i)) * scale;
	P( 6,i) = (b2*D(0,i) + c2*D(3,i) + d2*D(5,i)) * scale;
	P( 7,i) = (c2*D(1,i) + b2*D(3,i) + d2*D(4,i)) * scale;
	P( 8,i) = (d2*D(2,i) + c2*D(4,i) + b2*D(5,i)) * scale;
	P( 9,i) = (b3*D(0,i) + c3*D(3,i) + d3*D(5,i)) * scale;
	P(10,i) = (c3*D(1,i) + b3*D(3,i) + d3*D(4,i)) * scale;
	P(11,i) = (d3*D(2,i) + c3*D(4,i) + b3*D(5,i)) * scale;
    }

    // Now form 12x12 matrix K = PB
    RDenseMatrix K(12,12);
    for (i = 0; i < 12; i++) {
        K(i, 0) = (b0*P(i,0) + c0*P(i,3) + d0*P(i,5)) * scale;
	K(i, 1) = (c0*P(i,1) + b0*P(i,3) + d0*P(i,4)) * scale;
	K(i, 2) = (d0*P(i,2) + c0*P(i,4) + b0*P(i,5)) * scale;
	K(i, 3) = (b1*P(i,0) + c1*P(i,3) + d1*P(i,5)) * scale;
	K(i, 4) = (c1*P(i,1) + b1*P(i,3) + d1*P(i,4)) * scale;
	K(i, 5) = (d1*P(i,2) + c1*P(i,4) + b1*P(i,5)) * scale;
	K(i, 6) = (b2*P(i,0) + c2*P(i,3) + d2*P(i,5)) * scale;
	K(i, 7) = (c2*P(i,1) + b2*P(i,3) + d2*P(i,4)) * scale;
	K(i, 8) = (d2*P(i,2) + c2*P(i,4) + b2*P(i,5)) * scale;
	K(i, 9) = (b3*P(i,0) + c3*P(i,3) + d3*P(i,5)) * scale;
	K(i,10) = (c3*P(i,1) + b3*P(i,3) + d3*P(i,4)) * scale;
	K(i,11) = (d3*P(i,2) + c3*P(i,4) + b3*P(i,5)) * scale;
    }

    return K * size;
}

RDenseMatrix Tetrahedron4::ElasticityStiffnessMatrix (double E, double nu)
    const
{
    // This forms a 12x12 element stiffness matrix as above, but
    // assumes an isotropic material defined by Young's modulus E and
    // Poisson's ratio nu, which leads to a sparse symmetric elasticity
    // matrix (Zienkiewicz p.94)

    double mu = E/(2.0*(1.0+nu));                 // shear modulus
    double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu)); // Lame modulus

    // For now we explicitly form the elasticity matrix D, but K should
    // be formed directly for performance

    RDenseMatrix D(6,6);
    D(0,0) = D(1,1) = D(2,2) = lambda + 2.0*mu;
    D(0,1) = D(0,2) = D(1,2) = D(1,0) = D(2,0) = D(2,1) = lambda;
    D(3,3) = D(4,4) = D(5,5) = mu;

    double r = lambda + 2.0*mu;
    double s = lambda;
    double t = mu;
    
    RDenseMatrix BTDB(12,12);
    BTDB(0, 0) = b0*b0*r + (c0*c0 + d0*d0)*t;
    BTDB(0, 1) = b0*c0*s + b0*c0*t;
    BTDB(0, 2) = b0*d0*s + b0*d0*t;
    BTDB(0, 3) = b0*b1*r + (c0*c1 + d0*d1)*t;
    BTDB(0, 4) = b0*c1*s + b1*c0*t;
    BTDB(0, 5) = b0*d1*s + b1*d0*t;
    BTDB(0, 6) = b0*b2*r + (c0*c2 + d0*d2)*t;
    BTDB(0, 7) = b0*c2*s + b2*c0*t;
    BTDB(0, 8) = b0*d2*s + b2*d0*t;
    BTDB(0, 9) = b0*b3*r + (c0*c3 + d0*d3)*t;
    BTDB(0,10) = b0*c3*s + b3*c0*t;
    BTDB(0,11) = b0*d3*s + b3*d0*t;

    BTDB(1, 1) = c0*c0*r + (b0*b0 + d0*d0)*t;
    BTDB(1, 2) = c0*d0*s + c0*d0*t;
    BTDB(1, 3) = b1*c0*s + b0*c1*t;
    BTDB(1, 4) = c0*c1*r + (b0*b1 + d0*d1)*t;
    BTDB(1, 5) = c0*d1*s + c1*d0*t;
    BTDB(1, 6) = b2*c0*s + b0*c2*t;
    BTDB(1, 7) = c0*c2*r + (b0*b2 + d0*d2)*t;
    BTDB(1, 8) = c0*d2*s + c2*d0*t;
    BTDB(1, 9) = b3*c0*s + b0*c3*t;
    BTDB(1,10) = c0*c3*r + (b0*b3 + d0*d3)*t;
    BTDB(1,11) = c0*d3*s + c3*d0*t;
    BTDB /= 36.0*size;

    RDenseMatrix BTDB2 = ElasticityStiffnessMatrix (D);
    return BTDB2;
    //return ElasticityStiffnessMatrix (D);
}

RVector Tetrahedron4::InitialStrainVector (double E, double nu,
    const RVector &e0)
{
    double mu = E/(2.0*(1.0+nu));                 // shear modulus
    double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu)); // Lame modulus

    double r = lambda + 2.0*mu;
    double s = lambda;
    double t = mu;
    
    RVector res(12);
    
    res[ 0] = b0*e0[0]*r + (b0*e0[1] + b0*e0[2])*s + (c0*e0[3] + d0*e0[5])*t;
    res[ 1] = c0*e0[1]*r + (c0*e0[0] + c0*e0[2])*s + (b0*e0[3] + d0*e0[4])*t;
    res[ 2] = d0*e0[2]*r + (d0*e0[0] + d0*e0[1])*s + (c0*e0[4] + b0*e0[5])*t;
    res[ 3] = b1*e0[0]*r + (b1*e0[1] + b1*e0[2])*s + (c1*e0[3] + d1*e0[5])*t;
    res[ 4] = c1*e0[1]*r + (c1*e0[0] + c1*e0[2])*s + (b1*e0[3] + d1*e0[4])*t;
    res[ 5] = d1*e0[2]*r + (d1*e0[0] + d1*e0[1])*s + (c1*e0[4] + b1*e0[5])*t;
    res[ 6] = b2*e0[0]*r + (b2*e0[1] + b2*e0[2])*s + (c2*e0[3] + d2*e0[5])*t;
    res[ 7] = c2*e0[1]*r + (c2*e0[0] + c2*e0[2])*s + (b2*e0[3] + d2*e0[4])*t;
    res[ 8] = d2*e0[2]*r + (d2*e0[0] + d2*e0[1])*s + (c2*e0[4] + b2*e0[5])*t;
    res[ 9] = b3*e0[0]*r + (b3*e0[1] + b3*e0[2])*s + (c3*e0[3] + d3*e0[5])*t;
    res[10] = c3*e0[1]*r + (c3*e0[0] + c3*e0[2])*s + (b3*e0[3] + d3*e0[4])*t;
    res[11] = d3*e0[2]*r + (d3*e0[0] + d3*e0[1])*s + (c3*e0[4] + b3*e0[5])*t;

    return res/6.0;
}

RVector Tetrahedron4::ThermalExpansionVector (double E, double nu,
    double alpha, double dT)
{
    RVector e0(6);
    e0[0] = e0[1] = e0[2] = alpha*dT;
    return InitialStrainVector (E, nu, e0);
}

RVector Tetrahedron4::DThermalExpansionVector (double E, double nu)
{
    // returns derivative of thermal expansion vector with respect to
    // expansion coefficient (assuming temperature change 1)

    RVector de0(6);
    de0[0] = de0[1] = de0[2] = 1.0;
    return InitialStrainVector (E, nu, de0);
}

double Tetrahedron4::ComputeSize (const NodeList &nlist) const
{
    return (1.0/6.0) * (a0+a1+a2+a3);
}

RSymMatrix Tetrahedron4::ComputeIntDD (const NodeList &nlist) const
{
    RSymMatrix dd(4);
    double scale = 1.0/(36.0*size);
    dd(0,0) = scale * (b0*b0+c0*c0+d0*d0);
    dd(1,0) = scale * (b0*b1+c0*c1+d0*d1);
    dd(1,1) = scale * (b1*b1+c1*c1+d1*d1);
    dd(2,0) = scale * (b0*b2+c0*c2+d0*d2);
    dd(2,1) = scale * (b1*b2+c1*c2+d1*d2);
    dd(2,2) = scale * (b2*b2+c2*c2+d2*d2);
    dd(3,0) = scale * (b0*b3+c0*c3+d0*d3);
    dd(3,1) = scale * (b1*b3+c1*c3+d1*d3);
    dd(3,2) = scale * (b2*b3+c2*c3+d2*d3);
    dd(3,3) = scale * (b3*b3+c3*c3+d3*d3);
    return dd;
};

static const RDenseMatrix bndintf = RDenseMatrix (4, 4,
   "1 1 1 0\
    1 1 0 1\
    1 0 1 1\
    0 1 1 1") * (2.0/6.0);

static const RSymMatrix sym_bndintff_sd0 = RSymMatrix (4,
   "2 \
    1 2 \
    1 1 2 \
    0 0 0 0") * (2.0/24.0);
// note 2 in numerator because surface area in local element is 1/2
// global_bndint = local_bndint * global_area / local_area
static const RSymMatrix sym_bndintff_sd1 = RSymMatrix (4,
   "2 \
    1 2 \
    0 0 0 \
    1 1 0 2") * (2.0/24.0);
static const RSymMatrix sym_bndintff_sd2 = RSymMatrix (4,
   "2 \
    0 0 \
    1 0 2 \
    1 0 1 2") * (2.0/24.0);
static const RSymMatrix sym_bndintff_sd3 = RSymMatrix (4,
   "0 \
    0 2 \
    0 1 2 \
    0 1 1 2") * (2.0/24.0);
// note surface area of side 3 in local element is not actually 1/2
// but integration values are not correct either -> should compensate
static const RSymMatrix *sym_bndintff[4] = {
   &sym_bndintff_sd0,
   &sym_bndintff_sd1,
   &sym_bndintff_sd2,
   &sym_bndintff_sd3 
   };


RVector Tetrahedron4::BndIntF () const
{
    RVector bf(4);
    for (int sd = 0; sd < 4; sd++)
	if (bndside[sd])
	    bf += bndintf.Row(sd) * side_size[sd];
    return bf;
}

double Tetrahedron4::BndIntF (int i) const
{
    dASSERT(i >= 0 && i < 4, "Argument 1: out of range");
    double bf = 0.0;
    for (int sd = 0; sd < 4; sd++)
	if (bndside[sd])
	    bf += bndintf(sd, i) * side_size[sd];
    return bf;
}

double Tetrahedron4::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 4, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 4, "Argument 2: out of range");
    return bndintf(sd,i) * side_size[sd];
}


double Tetrahedron4::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 4, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 4, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 4, "Argument 3: out of range");
    RSymMatrix bff;
    bff= *sym_bndintff[sd];
    return bff(i,j) * side_size[sd];
}


RSymMatrix Tetrahedron4::ComputeBndIntFF (const NodeList &nlist) const
{
    RSymMatrix bff(4);
    if (!(bndel||interfaceel)) return bff;
    //    cout << "In ComputeBndIntFF with " << (interfaceel ? "interface" : "boundary ") << " sides ";
    for (int sd = 0; sd < 4; sd++) {
      //        cout << (bndside[sd] ? "true " : "false ");
        if (bndside[sd]) {
	    bff += *sym_bndintff[sd] * side_size[sd];
	}
    }
    //    cout << endl;
    return bff;
}

RVector Tetrahedron4::BndIntFX (int side, double (*func)(const Point&),
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

double Tetrahedron4::BndIntFD (int sd, int i, int j, int k)
{
    // computes \int_s du_i/dx_j u_k ds over side sd
    dASSERT (sd >= 0 && sd < 4, "Argument 1: index out of range");
    dASSERT (j >= 0 && j < 3, "Argument 3: index out of range");
#ifdef FEM_DEBUG
    int nd;
    int nsdnd = nSideNode(sd);
    for (nd = 0; nd < nsdnd; nd++)
	if (i == SideNode (sd,nd)) break;
    dASSERT (nd < nsdnd, "Argument 2: node index not found in side");
    for (nd = 0; nd < nsdnd; nd++)
	if (k == SideNode (sd,nd)) break;
    dASSERT (nd < nsdnd, "Argument 4: node index not found in side");
#endif
    double intf = IntF(k);
    intf *= 1.0/(6.0*size);
    switch (i) {
    case 0:
	intf *= (j == 0 ? b0 : j == 1 ? c0 : d0);
	break;
    case 1:
	intf *= (j == 0 ? b1 : j == 1 ? c1 : d1);
	break;
    case 2:
	intf *= (j == 0 ? b2 : j == 1 ? c2 : d2);
	break;
    case 3:
	intf *= (j == 0 ? b3 : j == 1 ? c3 : d3);
	break;
    }
	return intf;
}

#ifdef UNDEF
double Tetrahedron4::BndIntFD (Mesh &mesh, int el2, int i, int j, int k)
{
    // todo
}
#endif

RVector Tetrahedron4::BndIntFCos (int side, const RVector &cntcos, double a,
    const NodeList &nlist) const
{
    int j, p;
    double d, f;
    RVector sum(4);
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
    
int Tetrahedron4::Intersection (const Point &p1, const Point &p2,
    Point *s, bool add_endpoints, bool boundary_only)
{
    if (boundary_only && !bndel) return 0;
    
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
    if ((!boundary_only || bndside[3]) && (dx+dy+dz)) {
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
