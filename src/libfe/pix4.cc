// ==========================================================================
// Module libfe
// File pix4.cc
// Definition of class Pixel4
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <string.h>
#include "felib.h"

using namespace std;

Pixel4::Pixel4 (const Pixel4 &el): Element_Structured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Pixel4::Copy ()
{
    return new Pixel4(*this);
}

void Pixel4::Initialise (const NodeList &nlist)
{
    static bool need_setup = true;
    x0 = nlist[Node[0]][0];
    y0 = nlist[Node[0]][1];
    if (need_setup) {
        dx   = nlist[Node[1]][0] - x0;
	dy   = nlist[Node[2]][1] - y0;
	size = dx*dy;
	ComputeIntFF();
	ComputeIntDD();
	ComputeIntFDD();
	ComputeBndIntFFF();
	need_setup = false;
    }
    Element_Structured_2D::Initialise (nlist);
}

int Pixel4::SideNode (int side, int node) const
{
    dASSERT (side >= 0 && side < 4, "Side index out of range");
    dASSERT (node >= 0 && node < 2, "Node index out of range");
    static int SN[4][2] = {{0,1},{3,2},{2,0},{1,3}};
    return SN[side][node];
}

Point Pixel4::Local (const Point &glob) const
{
    dASSERT (glob.Dim() == 2, "Invalid point dimension");

    Point loc(2);
    loc[0] = (glob[0]-x0)/dx;
    loc[1] = (glob[1]-y0)/dy;
    return loc;
}

Point Pixel4::NodeLocal (int node) const
{
    Point nloc(2);
    switch (node) {
    case 0: break;
    case 1: nloc[0] = 1.0; break;
    case 2: nloc[1] = 1.0; break;
    case 3: nloc[0] = nloc[1] = 1.0; break;
    default: xERROR ("Node index out of range");
    }
    return nloc;
}

const RVector &Pixel4::LNormal (int side) const
{
    static const RVector lnm0 = RVector (2, " 0 -1");
    static const RVector lnm1 = RVector (2, " 0  1");
    static const RVector lnm2 = RVector (2, "-1  0");
    static const RVector lnm3 = RVector (2, " 1  0");
    static const RVector *lnm[4] = { &lnm0, &lnm1, &lnm2, &lnm3 };
    dASSERT (side >= 0 && side < 4, "Argument 1 index out of range");
    return *lnm[side];
}

bool Pixel4::LContains (const Point &loc, bool pad) const
{
    dASSERT (loc.Dim() == 2, "Argument 1 invalid dimension");

    if (pad) {
        const double EPS = 1e-10;
	return (loc[0] >= -EPS && loc[0] <= dx+EPS &&
		loc[1] >= -EPS && loc[1] <= dy+EPS);
    } else {
        return (loc[0] >= 0 && loc[0] <= dx &&
		loc[1] >= 0 && loc[1] <= dy);
    }
}

bool Pixel4::GContains (const Point &glob, const NodeList&) const
{
    // this version pads the element bounding box to include points on
    // the boundary
    const double EPS = 1e-8;

    return (glob[0] >= x0-EPS && glob[0] <= x0+dx+EPS &&
	    glob[1] >= y0-EPS && glob[1] <= y0+dy+EPS);
}

RVector Pixel4::LocalShapeF (const Point &loc) const
{
    dASSERT (loc.Dim() == 2, "Argument 1 invalid dimension");
    RVector fun(4);
    double iloc0 = 1.0-loc[0], iloc1 = 1.0-loc[1];

    fun[0] = iloc0  * iloc1;
    fun[1] = loc[0] * iloc1;
    fun[2] = iloc0  * loc[1];
    fun[3] = loc[0] * loc[1];
    return fun;
}

RDenseMatrix Pixel4::LocalShapeD (const Point &loc) const
{
    dASSERT (loc.Dim() == 2, "Argument 1 invalid dimension");

    double xm = 1.0-loc[0], ym = 1.0-loc[1];
    RDenseMatrix der (2,4);
    der(0,0) = -ym;
    der(1,0) = -xm;
    der(0,1) =  ym;
    der(1,1) = -loc[0];
    der(0,2) = -loc[1];
    der(1,2) =  xm;
    der(0,3) =  loc[1];
    der(1,3) =  loc[0];
    return der;
}

void Pixel4::ComputeIntFF () const
{
    double scale = size/36.0;
    intff.New(4,4);
    intff(0,0) = 4*scale;
    intff(0,1) = 2*scale;
    intff(0,2) = 2*scale;
    intff(0,3) = scale;
    intff(1,1) = 4*scale;
    intff(1,2) = scale;
    intff(1,3) = 2*scale;
    intff(2,2) = 4*scale;
    intff(2,3) = 2*scale;
    intff(3,3) = 4*scale;
}

double Pixel4::IntFFF(int i, int j, int k) const
{
    // index set for retrieving IntFFF(i,j,k)
    static const int intfff_index[4][4][4] = {
	{{0,1,1,2},{1,1,2,2},{1,2,1,2},{2,2,2,2}},
	{{1,1,2,2},{1,0,2,1},{2,2,2,2},{2,1,2,1}},
	{{1,2,1,2},{2,2,2,2},{1,2,0,1},{2,2,1,1}},
	{{2,2,2,2},{2,1,2,1},{2,2,1,1},{2,1,1,0}}
    };
    // scaling factors for IntFFF entries
    static const double intfff_scale[3] = {9.0/144.0, 3.0/144.0, 1.0/144.0};

    RANGE_CHECK(i >= 0 && i < 4 && j >= 0 && j < 4 && k >= 0 && k < 4);
    return size * intfff_scale[intfff_index[i][j][k]];
}

RSymMatrix Pixel4::IntPFF (const RVector &P) const
{
    static RSymMatrix pff(4);
    // static const double f0 = 1.0/16.0, f1 = 1.0/48.0, f2 = 1.0/144.0;
    double p0 = P[Node[0]], p1 = P[Node[1]], p2 = P[Node[2]], p3 = P[Node[3]];
    double fac = size/144.0;

    pff(0,0) = fac * (9*p0 + 3*p1 + 3*p2 +   p3);
    pff(1,0) = fac * (3*p0 + 3*p1 +   p2 +   p3);
    pff(2,0) = fac * (3*p0 +   p1 + 3*p2 +   p3);
    pff(3,0) = fac * (  p0 +   p1 +   p2 +   p3);
    pff(1,1) = fac * (3*p0 + 9*p1 +   p2 + 3*p3);
    pff(2,1) = fac * (  p0 +   p1 +   p2 +   p3);
    pff(3,1) = fac * (  p0 + 3*p1 +   p2 + 3*p3);
    pff(2,2) = fac * (3*p0 +   p1 + 9*p2 + 3*p3);
    pff(3,2) = fac * (  p0 +   p1 + 3*p2 + 3*p3);
    pff(3,3) = fac * (  p0 + 3*p1 + 3*p2 + 9*p3);
    return pff;
}

double Pixel4::IntPFF (int i, int j, const RVector &P) const
{
    RANGE_CHECK(i >= 0 && i < 4 && j >= 0 && j < 4);
    if (i < j) { int tmp = i; i = j; j = tmp; }
    double fac = size/144.0;
    switch (i) {
    case 0:
	return fac * (9*P[Node[0]]+3*P[Node[1]]+3*P[Node[2]]+P[Node[3]]);
    case 1:
	switch (j) {
	case 0:
	    return fac * (3*P[Node[0]]+3*P[Node[1]]+P[Node[2]]+P[Node[3]]);
	case 1:
	    return fac * (3*P[Node[0]]+9*P[Node[1]]+P[Node[2]]+3*P[Node[3]]);
	}
    case 2:
	switch (j) {
	case 0:
	    return fac * (3*P[Node[0]]+P[Node[1]]+3*P[Node[2]]+P[Node[3]]);
	case 1:
	    return fac * (P[Node[0]]+P[Node[1]]+P[Node[2]]+P[Node[3]]);
	case 2:
	    return fac * (3*P[Node[0]]+P[Node[1]]+9*P[Node[2]]+3*P[Node[3]]);
	}
    case 3:
	switch (j) {
	case 0:
	    return fac * (P[Node[0]]+P[Node[1]]+P[Node[2]]+P[Node[3]]);
	case 1:
	    return fac * (P[Node[0]]+3*P[Node[1]]+P[Node[2]]+3*P[Node[3]]);
	case 2:
	    return fac * (P[Node[0]]+P[Node[1]]+3*P[Node[2]]+3*P[Node[3]]);
	case 3:
	    return fac * (P[Node[0]]+3*P[Node[1]]+3*P[Node[2]]+9*P[Node[3]]);
	}
    }
    
    return 0;
}

void Pixel4::ComputeIntDD () const
{
    double dx2 = dx*dx, dy2 = dy*dy;
    intdd.New(4,4);
    intdd(0,0) = 2*(dx2 + dy2);
    intdd(0,1) = dx2 - 2*dy2;
    intdd(0,2) = -2*dx2 + dy2;
    intdd(0,3) = -dx2 - dy2;
    intdd(1,1) = 2*(dx2 + dy2);
    intdd(1,2) = -dx2 - dy2;
    intdd(1,3) = -2*dx2 + dy2;
    intdd(2,2) = 2*(dx2 + dy2);
    intdd(2,3) = dx2 - 2*dy2;
    intdd(3,3) = 2*(dx2 + dy2);
    intdd /= 6.0*dx*dy;
}

void Pixel4::ComputeIntFDD () const
{
    int i;
    double dx2 = dx*dx, dy2 = dy*dy;

    for (i = 0; i < 4; i++)
        intfdd[i].New(4,4);

    intfdd[0](0,0) = 3*(dx2 + dy2);
    intfdd[0](0,1) = dx2 - 3*dy2;
    intfdd[0](0,2) = -3*dx2 + dy2;
    intfdd[0](0,3) = -dx2 - dy2;
    intfdd[0](1,1) = dx2 + 3*dy2;
    intfdd[0](1,2) = -dx2 - dy2;
    intfdd[0](1,3) = -dx2 + dy2;
    intfdd[0](2,2) = 3*dx2 + dy2;
    intfdd[0](2,3) = dx2 - dy2;
    intfdd[0](3,3) = dx2 + dy2;
    intfdd[1](0,0) = dx2 + 3*dy2;
    intfdd[1](0,1) = dx2 - 3*dy2;
    intfdd[1](0,2) = -dx2 + dy2;
    intfdd[1](0,3) = -dx2 - dy2;
    intfdd[1](1,1) = 3*(dx2 + dy2);
    intfdd[1](1,2) = -dx2 - dy2;
    intfdd[1](1,3) = -3*dx2 + dy2;
    intfdd[1](2,2) = dx2 + dy2;
    intfdd[1](2,3) = dx2 - dy2;
    intfdd[1](3,3) = 3*dx2 + dy2;
    intfdd[2](0,0) = 3*dx2 + dy2;
    intfdd[2](0,1) = dx2 - dy2;
    intfdd[2](0,2) = -3*dx2 + dy2;
    intfdd[2](0,3) = -dx2 - dy2;
    intfdd[2](1,1) = dx2 + dy2;
    intfdd[2](1,2) = -dx2 - dy2;
    intfdd[2](1,3) = -dx2 + dy2;
    intfdd[2](2,2) = 3*(dx2 + dy2);
    intfdd[2](2,3) = dx2 - 3*dy2;
    intfdd[2](3,3) = dx2 + 3*dy2;
    intfdd[3](0,0) = dx2 + dy2;
    intfdd[3](0,1) = dx2 - dy2;
    intfdd[3](0,2) = -dx2 + dy2;
    intfdd[3](0,3) = -dx2 - dy2;
    intfdd[3](1,1) = 3*dx2 + dy2;
    intfdd[3](1,2) = -dx2 - dy2;
    intfdd[3](1,3) = -3*dx2 + dy2;
    intfdd[3](2,2) = dx2 + 3*dy2;
    intfdd[3](2,3) = dx2 - 3*dy2;
    intfdd[3](3,3) = 3*(dx2 + dy2);

    for (i = 0; i < 4; i++)
        intfdd[i] /= 24*dx*dy;
}

RSymMatrix Pixel4::Intdd () const
{
    static RSymMatrix dd;
    if (dd.nRows() == 0) {
	dd.New (8,8);
	double dx2 = dx*dx, dy2 = dy*dy, dxy = dx*dy;
	dd(0,0) = 4*dy2;
	dd(1,0) = 3*dxy, dd(1,1) = 4*dx2;
	dd(2,0) = -4*dy2, dd(2,1) = -3*dxy, dd(2,2) = 4*dy2;
	dd(3,0) = 3*dxy,  dd(3,1) = 2*dx2,  dd(3,2) = -3*dxy, dd(3,3) = 4*dx2;
	dd(4,0) = 2*dy2,  dd(4,1) = 3*dxy,  dd(4,2) = -2*dy2, dd(4,3) = 3*dxy,
	    dd(4,4) = 4*dy2;
	dd(5,0) = -3*dxy, dd(5,1) = -4*dx2, dd(5,2) = 3*dxy,  dd(5,3) = -2*dx2,
	    dd(5,4) = -3*dxy, dd(5,5) = 4*dx2;
	dd(6,0) = -2*dy2, dd(6,1) = -3*dxy, dd(6,2) = 2*dy2,  dd(6,3) = -3*dxy,
	    dd(6,4) = -4*dy2, dd(6,5) = 3*dxy, dd(6,6) = 4*dy2;
	dd(7,0) = -3*dxy, dd(7,1) = -2*dx2, dd(7,2) = 3*dxy, dd(7,3) = -4*dx2,
	    dd(7,4) = -3*dxy, dd(7,5) = 2*dx2, dd(7,6) = 3*dxy,
	    dd(7,7) = 4*dx2;

	dd *= 1.0/(12.0*dxy);
    }
    return dd;
}

double Pixel4::IntFd (int i, int j, int k) const
{
    static RDenseMatrix fid;
    if (fid.nRows() == 0) {
	fid.New (8,4);
	fid(0,0) = -2*dy; fid(0,1) = 2*dy;  fid(0,2) = -dy;   fid(0,3) = dy;
	fid(1,0) = -2*dx; fid(1,1) = -dx;   fid(1,2) = 2*dx;  fid(1,3) = dx;
	fid(2,0) = -2*dy; fid(2,1) = 2*dy;  fid(2,2) = -dy;   fid(2,3) = dy;
	fid(3,0) = -dx;   fid(3,1) = -2*dx; fid(3,2) = dx;    fid(3,3) = 2*dx;
	fid(4,0) = -dy;   fid(4,1) = dy;    fid(4,2) = -2*dy; fid(4,3) = 2*dy;
	fid(5,0) = -2*dx; fid(5,1) = -dx;   fid(5,2) = 2*dx;  fid(5,3) = dx;
	fid(6,0) = -dy;   fid(6,1) = dy;    fid(6,2) = -2*dy; fid(6,3) = 2*dy;
	fid(7,0) = -dx;   fid(7,1) = -2*dx; fid(7,2) = dx;    fid(7,3) = 2*dx;

	fid *= 1.0/12.0;
    }
    return fid(i*2+k,j);
}

double Pixel4::IntFdd (int i, int j, int k, int l, int m) const
{
    static RSymMatrix *fidd = 0;
    if (!fidd) {
	int z;
	fidd = new RSymMatrix[4];
	for (z = 0; z < 4; z++)
	    fidd[z].New(8,8);
	RSymMatrix &m0 = fidd[0];
	RSymMatrix &m1 = fidd[1];
	RSymMatrix &m2 = fidd[2];
	RSymMatrix &m3 = fidd[3];
	
	double dx2 = dx*dx, dy2 = dy*dy, dxy = dx*dy;

	m0(0,0) = 9*dy2;
	m0(1,0) = 8*dxy,  m0(1,1) = 9*dx2;
	m0(2,0) = -9*dy2, m0(2,1) = -8*dxy, m0(2,2) = 9*dy2;
	m0(3,0) = 4*dxy,  m0(3,1) = 3*dx2, m0(3,2) = -4*dxy, m0(3,3) = 3*dx2;
	m0(4,0) = 3*dy2,  m0(4,1) = 4*dxy, m0(4,2) = -3*dy2, m0(4,3) = 2*dxy,
	    m0(4,4) = 3*dy2;
	m0(5,0) = -8*dxy, m0(5,1) = -9*dx2, m0(5,2) = 8*dxy, m0(5,3) = -3*dx2,
	    m0(5,4) = -4*dxy, m0(5,5) = 9*dx2;
	m0(6,0) = -3*dy2, m0(6,1) = -4*dxy, m0(6,2) = 3*dy2, m0(6,3) = -2*dxy,
	    m0(6,4) = -3*dy2, m0(6,5) = 4*dxy, m0(6,6) = 3*dy2;
	m0(7,0) = -4*dxy, m0(7,1) = -3*dx2, m0(7,2) = 4*dxy, m0(7,3) = -3*dx2,
	    m0(7,4) = -2*dxy, m0(7,5) = 3*dx2, m0(7,6) = 2*dxy,
	    m0(7,7) = 3*dx2;

	m1(0,0) = 9*dy2;
	m1(1,0) = 4*dxy,  m1(1,1) = 3*dx2;
	m1(2,0) = -9*dy2, m1(2,1) = -4*dxy, m1(2,2) = 9*dy2;
	m1(3,0) = 8*dxy,  m1(3,1) = 3*dx2,  m1(3,2) = -8*dxy, m1(3,3) = 9*dx2;
	m1(4,0) = 3*dy2,  m1(4,1) = 2*dxy,  m1(4,2) = -3*dy2, m1(4,3) = 4*dxy,
	    m1(4,4) = 3*dy2;
	m1(5,0) = -4*dxy, m1(5,1) = -3*dx2, m1(5,2) = 4*dxy,  m1(5,3) = -3*dx2,
	    m1(5,4) = -2*dxy, m1(5,5) = 3*dx2;
	m1(6,0) = -3*dy2, m1(6,1) = -2*dxy, m1(6,2) = 3*dy2,  m1(6,3) = -4*dxy,
	    m1(6,4) = -3*dy2, m1(6,5) = 2*dxy, m1(6,6) = 3*dy2;
	m1(7,0) = -8*dxy, m1(7,1) = -3*dx2, m1(7,2) = 8*dxy,  m1(7,3) = -9*dx2,
	    m1(7,4) = -4*dxy, m1(7,5) = 3*dx2, m1(7,6) = 4*dxy,
	    m1(7,7) = 9*dx2;

	m2(0,0) = 3*dy2;
	m2(1,0) = 4*dxy,  m2(1,1) = 9*dx2;
	m2(2,0) = -3*dy2, m2(2,1) = -4*dxy, m2(2,2) = 3*dy2;
	m2(3,0) = 2*dxy,  m2(3,1) = 3*dx2,  m2(3,2) = -2*dxy, m2(3,3) = 3*dx2;
	m2(4,0) = 3*dy2,  m2(4,1) = 8*dxy,  m2(4,2) = -3*dy2, m2(4,3) = 4*dxy,
	    m2(4,4) = 9*dy2;
	m2(5,0) = -4*dxy, m2(5,1) = -9*dx2, m2(5,2) = 4*dxy,  m2(5,3) = -3*dx2,
	    m2(5,4) = -8*dxy, m2(5,5) = 9*dx2;
	m2(6,0) = -3*dy2, m2(6,1) = -8*dxy, m2(6,2) = 3*dy2,  m2(6,3) = -4*dxy,
	    m2(6,4) = -9*dy2, m2(6,5) = 8*dxy, m2(6,6) = 9*dy2;
	m2(7,0) = -2*dxy, m2(7,1) = -3*dx2, m2(7,2) = 2*dxy,  m2(7,3) = -3*dx2,
	    m2(7,4) = -4*dxy, m2(7,5) = 3*dx2, m2(7,6) = 4*dxy,
	    m2(7,7) = 3*dx2;

	m3(0,0) = 3*dy2;
	m3(1,0) = 2*dxy,  m3(1,1) = 3*dx2;
	m3(2,0) = -3*dy2, m3(2,1) = -2*dxy, m3(2,2) = 3*dy2;
	m3(3,0) = 4*dxy,  m3(3,1) = 3*dx2,  m3(3,2) = -4*dxy, m3(3,3) = 9*dx2;
	m3(4,0) = 3*dy2,  m3(4,1) = 4*dxy,  m3(4,2) = -3*dy2, m3(4,3) = 8*dxy,
	    m3(4,4) = 9*dy2;
	m3(5,0) = -2*dxy, m3(5,1) = -3*dx2, m3(5,2) = 2*dxy,  m3(5,3) = -3*dx2,
	    m3(5,4) = -4*dxy, m3(5,5) = 3*dx2;
	m3(6,0) = -3*dy2, m3(6,1) = -4*dxy, m3(6,2) = 3*dy2,  m3(6,3) = -8*dxy,
	    m3(6,4) = -9*dy2, m3(6,5) = 4*dxy, m3(6,6) = 9*dy2;
	m3(7,0) = -4*dxy, m3(7,1) = -3*dx2, m3(7,2) = 4*dxy,  m3(7,3) = -9*dx2,
	    m3(7,4) = -8*dxy, m3(7,5) = 3*dx2, m3(7,6) = 8*dxy,
	    m3(7,7) = 9*dx2;

	for (z = 0; z < 4; z++)
	    fidd[z] /= 72.0*dx*dy;
    }
    return fidd[i](j*2+l, k*2+m);
}

double Pixel4::IntPdd (const RVector &P, int j, int k, int l, int m) const
{
    double val = 0;
    for (int i = 0; i < 4; i++)
	val += IntFdd (i,j,k,l,m) * P[Node[i]];
    return val;
}

double Pixel4::IntFfd (int i, int j, int k, int l) const
{
    static double fac[4][4][4][2] =
	{{{{-6,-6},{6,-2},{-2,6},{2,2}},
	  {{-3,-2},{3,-2},{-1,2},{1,2}},
	  {{-2,-3},{2,-1},{-2,3},{2,1}},
	  {{-1,-1},{1,-1},{-1,1},{1,1}}},
	 {{{-3,-2},{3,-2},{-1,2},{1,2}},
	  {{-6,-2},{6,-6},{-2,2},{2,6}},
	  {{-1,-1},{1,-1},{-1,1},{1,1}},
	  {{-2,-1},{2,-3},{-2,1},{2,3}}},
	 {{{-2,-3},{2,-1},{-2,3},{2,1}},
	  {{-1,-1},{1,-1},{-1,1},{1,1}},
	  {{-2,-6},{2,-2},{-6,6},{6,2}},
	  {{-1,-2},{1,-2},{-3,2},{3,2}}},
	 {{{-1,-1},{1,-1},{-1,1},{1,1}},
	  {{-2,-1},{2,-3},{-2,1},{2,3}},
	  {{-1,-2},{1,-2},{-3,2},{3,2}},
	  {{-2,-2},{2,-6},{-6,2},{6,6}}}};

    double d = (l == 0 ? dy:dx);

    return d * fac[i][j][k][l] / 72.0;
}

double Pixel4::IntPfd (const RVector &P, int j, int k, int l) const
{
    double val = 0;
    for (int i = 0; i < 4; i++)
	val += IntFfd(i,j,k,l) * P[Node[i]];
    return val;
}

double Pixel4::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 4, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 4, "Argument 2: out of range");
    switch (sd) {
    case 0: return (i==0 || i==1 ? 0.5*dx : 0.0);
    case 1: return (i==2 || i==3 ? 0.5*dx : 0.0);
    case 2: return (i==0 || i==2 ? 0.5*dy : 0.0);
    case 3: return (i==1 || i==3 ? 0.5*dy : 0.0);
    default: xERROR("Argument 2: out of range"); return 0.0;
    }
}

double Pixel4::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 4, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 4, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 4, "Argument 3: out of range");
    if (!bndintff) ComputeBndIntFF();
    return bndintff[sd](i,j);
}

RSymMatrix Pixel4::BndIntFF () const
{
    if (!bndintff) ComputeBndIntFF();
    RSymMatrix res(4,4);
    for (int sd = 0; sd < 4; sd++)
	if (bndside[sd]) res += bndintff[sd];
    return res;
}

void Pixel4::ComputeBndIntFF () const
{
    int z;
    bndintff = new RSymMatrix[4];
    for (z = 0; z < 4; z++)
	bndintff[z].New(4,4);
    // side 0: y = 0
    bndintff[0](0,0) = bndintff[0](1,1) = 2*dx;
    bndintff[0](1,0) = dx;
    // side 1: y = 1
    bndintff[1](2,2) = bndintff[1](3,3) = 2*dx;
    bndintff[1](3,2) = dx;
    // side 2: x = 0
    bndintff[2](0,0) = bndintff[2](2,2) = 2*dy;
    bndintff[2](2,0) = dy;
    // side 3: x = 1;
    bndintff[3](1,1) = bndintff[3](3,3) = 2*dy;
    bndintff[3](3,1) = dy;
    for (z = 0; z < 4; z++)
	bndintff[z] *= 1.0/6.0;
}

void Pixel4::ComputeBndIntFFF () const
{
    int i, sd;
    double scale;

    for (sd = 0; sd < 4; sd++) 
        for (i = 0; i < 4; i++) bndintfff[sd][i].New(4,4);

    scale = dx/12.0;
    bndintfff[0][0](0,0) = 3*scale;
    bndintfff[0][0](0,1) = scale;
    bndintfff[0][0](1,0) = scale;
    bndintfff[0][0](1,1) = scale;
    bndintfff[0][1](0,0) = scale;
    bndintfff[0][1](0,1) = scale;
    bndintfff[0][1](1,0) = scale;
    bndintfff[0][1](1,1) = 3*scale;

    bndintfff[1][2](2,2) = 3*scale;
    bndintfff[1][2](2,3) = scale;
    bndintfff[1][2](3,2) = scale;
    bndintfff[1][2](3,3) = scale;
    bndintfff[1][3](2,2) = scale;
    bndintfff[1][3](2,3) = scale;
    bndintfff[1][3](3,2) = scale;
    bndintfff[1][3](3,3) = 3*scale;

    scale = dy/12.0;
    bndintfff[2][0](0,0) = 3*scale;
    bndintfff[2][0](0,2) = scale;
    bndintfff[2][0](2,0) = scale;
    bndintfff[2][0](2,2) = scale;
    bndintfff[2][2](0,0) = scale;
    bndintfff[2][2](0,2) = scale;
    bndintfff[2][2](2,0) = scale;
    bndintfff[2][2](2,2) = 3*scale;

    bndintfff[3][1](1,1) = 3*scale;
    bndintfff[3][1](1,3) = scale;
    bndintfff[3][1](3,1) = scale;
    bndintfff[3][1](3,3) = scale;
    bndintfff[3][3](1,1) = scale;
    bndintfff[3][3](1,3) = scale;
    bndintfff[3][3](3,1) = scale;
    bndintfff[3][3](3,3) = 3*scale;

#ifdef UNDEF
    int j, k, nd[2];
    double area[4];
    
    area[0] = area[1] = dx;
    area[2] = area[3] = dy;

    for (sd = 0; sd < 4; sd++) {
        area[sd] /= 12.0; // scaling value;
	
	for (i = 0; i < 4; i++) bndintfff[sd][i].New(4,4);
	for (i = 0; i < 2; i++) nd[i] = SideNode (sd,i);

	for (i = 0; i < 2; i++)
  	    for (j = 0; j < 2; j++)
	        for (k = 0; k < 2; k++)
		  bndintfff[sd][nd[i]](nd[j],nd[k]) = 
		      (i==j && i==k ? 3*area[sd] : area[sd]);
    }
#endif
}

double Pixel4::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    for (int sd = 0; sd < 4; sd++) {
        if (bndside[sd]) {
	    for (int k = 0; k < 4; k++)
	        res += P[Node[k]] * bndintfff[sd][k](i,j);
	}
    }
    return res;
}

double Pixel4::Size() const
{
	return size;
}

double Pixel4::IntF (int i) const
{
	return 0.25*dx*dy;
}

double Pixel4::IntFF (int i, int j) const
{
	return intff(i,j);
}

RSymMatrix Pixel4::IntFF() const
{
	return intff;
}

RSymMatrix Pixel4::IntDD () const
{
	return intdd;
}

double Pixel4::IntDD (int i, int j) const
{
	return intdd(i,j);
}

double Pixel4::IntFDD (int i, int j, int k) const
{
	return intfdd[i](j,k);
}

double Pixel4::IntPDD (int i, int j, const RVector &P) const
{
	double res = 0.0;
    for (int k = 0; k < 4; k++) res += P[Node[k]] * intfdd[k](i,j);
    return res;
}

double Pixel4::dx = 0.0;
double Pixel4::dy = 0.0;
double Pixel4::size = 0.0;
RSymMatrix Pixel4::intff;
RSymMatrix Pixel4::intdd;
RSymMatrix Pixel4::intfdd[4];
RSymMatrix *Pixel4::bndintff = 0;
RDenseMatrix Pixel4::bndintfff[4][4];
