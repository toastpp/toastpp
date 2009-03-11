// ==========================================================================
// Module libfe
// File cub8_reg.cc
// Definition of class Cube8_reg
// ==========================================================================

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"
#include "toastdef.h"

Cube8_reg::Cube8_reg (const Cube8_reg &el): Element3D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

void Cube8_reg::Initialise (const NodeList &nlist)
{
    dx = nlist[Node[1]][0] - nlist[Node[0]][0];
    dy = nlist[Node[2]][1] - nlist[Node[0]][1];
    dz = nlist[Node[4]][2] - nlist[Node[0]][2];
}

int Cube8_reg::SideNode (int side, int node) const
{
    dASSERT(side >= 0 && side < 6, Side index out of range);
    dASSERT(node >= 0 && node < 4, Node index out of range);
    static int SN[6][4] = {{0,1,2,3},{4,6,7,5},{0,4,5,1},
                           {2,3,7,6},{0,2,6,4},{1,5,7,3}};
    return SN[side][node];
}

Point Cube8_reg::Local (const NodeList &nlist, const Point &glob) const
{
    // this assumes that cube edges are aligned with global coordinate axes

    dASSERT(glob.Dim() == 3, Invalid point dimension);
    Point loc(3);
    double xg0 = nlist[Node[0]][0];
    double yg0 = nlist[Node[0]][1];
    double zg0 = nlist[Node[0]][2];
    double xg1 = nlist[Node[1]][0];
    double yg1 = nlist[Node[2]][1];
    double zg1 = nlist[Node[4]][2];

    loc[0] = (glob[0]-xg0)/(xg1-xg0);
    loc[1] = (glob[1]-yg0)/(yg1-yg0);
    loc[2] = (glob[2]-zg0)/(zg1-zg0);
    return loc;
}

Point Cube8_reg::NodeLocal (int node) const
{
    Point nloc(3);
    switch (node) {
    case 0: nloc[0] = nloc[1] = nloc[2] = 0.0; break;
    case 1: nloc[0] = 1.0; nloc[1] = nloc[2] = 0.0; break;
    case 2: nloc[0] = nloc[2] = 0.0; nloc[1] = 1.0; break;
    case 3: nloc[0] = nloc[1] = 1.0; nloc[2] = 0.0; break;
    case 4: nloc[0] = nloc[1] = 0.0; nloc[2] = 1.0; break;
    case 5: nloc[0] = nloc[2] = 1.0; nloc[1] = 0.0; break;
    case 6: nloc[0] = 0.0; nloc[1] = nloc[2] = 1.0; break;
    case 7: nloc[0] = nloc[1] = nloc[2] = 1.0; break;
    default: xERROR(Node index out of range);
    }
    return nloc;
}

RVector Cube8_reg::LocalShapeF (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, Invalid point dimension);
    RVector fun(8);
    double x0_term = (1.0-loc[0]);
    double x1_term = loc[0];
    double y0_term = (1.0-loc[1]);
    double y1_term = loc[1];
    double z0_term = (1.0-loc[2]);
    double z1_term = loc[2];
    fun[0] = x0_term * y0_term * z0_term;
    fun[1] = x1_term * y0_term * z0_term;
    fun[2] = x0_term * y1_term * z0_term;
    fun[3] = x1_term * y1_term * z0_term;
    fun[4] = x0_term * y0_term * z1_term;
    fun[5] = x1_term * y0_term * z1_term;
    fun[6] = x0_term * y1_term * z1_term;
    fun[7] = x1_term * y1_term * z1_term;
    return fun;
}

RDenseMatrix Cube8_reg::LocalShapeD (const Point &loc) const
{
    dASSERT(loc.Dim() == 3, Invalid point dimension);
    RDenseMatrix der(3, 8);
    double x0_term = (1.0-loc[0]);
    double x1_term = loc[0];
    double y0_term = (1.0-loc[1]);
    double y1_term = loc[1];
    double z0_term = (1.0-loc[2]);
    double z1_term = loc[2];
    der(0,0) = -x1_term * y0_term * z0_term;
    der(1,0) = -y1_term * x0_term * z0_term;
    der(2,0) = -z1_term * x0_term * y0_term;
    der(0,1) =            y0_term * z0_term;
    der(1,1) = -x1_term * y1_term * z0_term;
    der(2,1) = -x1_term * y0_term * z1_term;
    der(0,2) = -x1_term * y1_term * z0_term;
    der(1,2) =  x0_term *           z0_term;
    der(2,2) = -x0_term * y1_term * z1_term;
    der(0,3) =            y1_term * z0_term;
    der(1,3) =  x1_term *           z0_term;
    der(2,3) = -x1_term * y1_term * z1_term;
    der(0,4) = -x1_term * y0_term * z1_term;
    der(1,4) = -x0_term * y1_term * z1_term;
    der(2,4) =  x0_term * y0_term;
    der(0,5) =            y0_term * z1_term;
    der(1,5) = -x1_term * y1_term * z1_term;
    der(2,5) =  x1_term * y0_term;
    der(0,6) = -x1_term * y1_term * z1_term;
    der(1,6) =  x0_term *           z1_term;
    der(2,6) =  x0_term * y1_term;
    der(0,7) =            y1_term * z1_term;
    der(1,7) =  x1_term *           z1_term;
    der(2,7) =  x1_term * y1_term;
    return der;
}

RVector Cube8_reg::GlobalShapeF (const NodeList &nlist, const Point &glob)
    const
{
    dASSERT(glob.Dim() == 3, Invalid point dimension);
    RVector fun(8);
    double x0_term = (nlist[Node[1]][0]-glob[0])/
                     (nlist[Node[1]][0]-nlist[Node[0]][0]);
    double x1_term = 1.0-x0_term;
    double y0_term = (nlist[Node[2]][1]-glob[1])/
                     (nlist[Node[2]][1]-nlist[Node[0]][1]);
    double y1_term = 1.0-y0_term;
    double z0_term = (nlist[Node[4]][2]-glob[2])/
                     (nlist[Node[4]][2]-nlist[Node[0]][2]);
    double z1_term = 1.0-z0_term;
    fun[0] = x0_term * y0_term * z0_term;
    fun[1] = x1_term * y0_term * z0_term;
    fun[2] = x0_term * y1_term * z0_term;
    fun[3] = x1_term * y1_term * z0_term;
    fun[4] = x0_term * y0_term * z1_term;
    fun[5] = x1_term * y0_term * z1_term;
    fun[6] = x0_term * y1_term * z1_term;
    fun[7] = x1_term * y1_term * z1_term;
    return fun;
}

RDenseMatrix Cube8_reg::GlobalShapeD (const NodeList &nlist,
    const Point &glob) const
{
    dASSERT(glob.Dim() == 3, Invalid point dimension);
    RDenseMatrix der (LocalShapeD (Local (nlist, glob)));
    double scalex = 1.0/(nlist[Node[1]][0]-nlist[Node[0]][0]);
    double scaley = 1.0/(nlist[Node[2]][1]-nlist[Node[0]][1]);
    double scalez = 1.0/(nlist[Node[4]][2]-nlist[Node[0]][2]);
    for (int i = 0; i < 8; i++) {
        der(0,i) *= scalex;
	der(1,i) *= scaley;
	der(2,i) *= scalez;
    }
    return der;
}

double Cube8_reg::IntF (int i) const
{
    return 0.125*dx*dy*dz;
}

RDenseMatrix Cube8_reg::ElasticityStiffnessMatrix (double E, double nu)
    const
{
    const int nq = 0;
    int k;

    RDenseMatrix B(6,24);
    RDenseMatrix D(6,6);
    RDenseMatrix BTDB(24,24);
    
    // Numerical integration
    for (k = 0; k < nq; k++) {
    }
    

    return BTDB;
}
