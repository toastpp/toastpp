// ==========================================================================
// Module libfe
// File quad4.cc
// Definition of class Quad4
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

Quad4::Quad4 (const Quad4 &el): Element_Unstructured_2D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

void Quad4::Initialise (const NodeList &nlist)
{
    n0 = nlist[Node[0]];
    len = nlist[Node[1]][0] - nlist[Node[0]][0];
}

int Quad4::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 4, "Argument 1 index out of range");
    dASSERT(node >= 0 && node < 2, "Argument 2 index out of range");
    static int SN[4][2] = {{0,1},{1,2},{2,3},{3,0}};
    return SN[side][node];
}

Point Quad4::Local (const NodeList &nlist, const Point &glob) const
{
    dASSERT(glob.Dim() == 2, "Argument 2 dimension must be 2");
    static Point loc(2);
    loc[0] = (glob[0]-n0[0])/len;
    loc[1] = (glob[1]-n0[1])/len;
    return loc;
}

Point Quad4::NodeLocal (int node) const
{
    Point n(2);
    switch (node) {
    case 1: n[0] = 1.0; break;
    case 2: n[1] = 1.0; break;
    case 3: n[0] = n[1] = 1.0; break;
    }
    return n;
}

bool Quad4::LContains (const Point &loc, bool pad) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    if (pad) {
        static const double EPS = 1e-8;
	return (loc[0]+EPS >= 0.0 && loc[0]-EPS <= 1.0 &&
		loc[1]+EPS >= 0.0 && loc[1]-EPS <= 1.0);
    } else {
        return (loc[0] >= 0.0 && loc[0] <= 1.0 &&
		loc[1] >= 0.0 && loc[1] <= 1.0);
    }
}

bool Quad4::GContains (const Point &glob, const NodeList &nlist) const
{
    dASSERT(glob.Dim() == 2, "Argument 1 invalid dimension");
    const double EPS = 1e-8;
    return (glob[0]+EPS >= n0[0] && glob[0]-EPS <= n0[0]+len &&
	    glob[1]+EPS >= n0[1] && glob[1]+EPS <= n0[1]+len);
}

RVector Quad4::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static RVector fun(4);
    fun[0] = (1.0-loc[0])*(1.0-loc[1]);
    fun[1] = loc[0]*(1.0-loc[1]);
    fun[2] = (1.0-loc[0])*loc[1];
    fun[3] = loc[0]*loc[1];
    return fun;
}

RDenseMatrix Quad4::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Argument 1 invalid dimension");
    static RDenseMatrix der(2,4);
    der(0,0) = loc[1]-1.0;
    der(1,0) = loc[0]-1.0;
    der(0,1) = 1.0-loc[1];
    der(1,1) = -loc[0];
    der(2,0) = -loc[1];
    der(2,1) = 1.0-loc[0];
    der(3,0) = loc[1];
    der(3,1) = loc[0];
    return der;
}
