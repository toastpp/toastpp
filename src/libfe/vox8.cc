// ==========================================================================
// Module libfe
// File vox8.cc
// Definition of class Voxel8
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"
#include "toastdef.h"

bool Voxel8::need_setup = true;

Voxel8::Voxel8 (const Voxel8 &el): Element_Structured_3D (el)
{
    Node = new int[nNode()];
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

Element *Voxel8::Copy ()
{
    return new Voxel8(*this);
}

void Voxel8::Initialise (const NodeList &nlist)
{
    x0 = nlist[Node[0]][0];
    y0 = nlist[Node[0]][1];
    z0 = nlist[Node[0]][2];
    if (need_setup) {
        dx    = nlist[Node[1]][0] - x0;
	dy    = nlist[Node[2]][1] - y0;
	dz    = nlist[Node[4]][2] - z0;
	size  = dx*dy*dz;
	intf  = size*0.125;
	ComputeIntFF();
	ComputeIntFFF();
	ComputeIntDD();
	ComputeIntFDD();
	ComputeBndIntF();
	ComputeBndIntFF();
	ComputeBndIntFFF();
	need_setup = false;
    }
    Element_Structured_3D::Initialise (nlist);
}

void Voxel8::PostInitialisation (const NodeList &nlist)
{
    need_setup = true;
    Element_Structured_3D::PostInitialisation (nlist);
}

int Voxel8::SideNode (int side, int node) const
{
    dASSERT (side >= 0 && side < 6, "Side index out of range");
    dASSERT (node >= 0 && node < 4, "Node index out of range");
    static int SN[6][4] = {{0,1,3,2},{7,5,4,6},{0,4,5,1},
                           {7,6,2,3},{0,2,6,4},{7,3,1,5}};
    return SN[side][node];
}

Point Voxel8::Local (const Point &glob) const
{
    dASSERT(glob.Dim() == 3, "Invalid point dimension");

    Point loc(3);
    loc[0] = (glob[0]-x0)/dx;
    loc[1] = (glob[1]-y0)/dy;
    loc[2] = (glob[2]-z0)/dz;
    return loc;
}

Point Voxel8::NodeLocal (int node) const
{
    Point nloc(3);
    switch (node) {
    case 0: break;
    case 1: nloc[0] = 1.0; break;
    case 2: nloc[1] = 1.0; break;
    case 3: nloc[0] = nloc[1] = 1.0; break;
    case 4: nloc[2] = 1.0; break;
    case 5: nloc[0] = nloc[2] = 1.0; break;
    case 6: nloc[1] = nloc[2] = 1.0; break;
    case 7: nloc[0] = nloc[1] = nloc[2] = 1.0; break;
    default: xERROR ("Node index out of range");
    }
    return nloc;
}

const RVector &Voxel8::LNormal (int side) const
{
    static const RVector lnm0 = RVector (3, "0 0 -1");
    static const RVector lnm1 = RVector (3, "0 0  1");
    static const RVector lnm2 = RVector (3, "0 -1 0");
    static const RVector lnm3 = RVector (3, "0  1 0");
    static const RVector lnm4 = RVector (3, "-1 0 0");
    static const RVector lnm5 = RVector (3, " 1 0 0");
    static const RVector *lnm[6] = {
        &lnm0, &lnm1, &lnm2, &lnm3, &lnm4, &lnm5
    };
    dASSERT (side >= 0 && side < 6, "Argument 1 index out of range");
    return *lnm[side];
}

bool Voxel8::LContains (const Point &loc, bool pad) const
{
    dASSERT (loc.Dim() == 3, "Argument 1 invalid dimension");

    if (pad) {
        const double EPS = 1e-8;
        return (loc[0] >= -EPS && loc[0] <= 1+EPS &&
		loc[1] >= -EPS && loc[1] <= 1+EPS &&
		loc[2] >= -EPS && loc[2] <= 1+EPS);
    } else {
        return (loc[0] >= 0 && loc[0] <= 1 &&
		loc[1] >= 0 && loc[1] <= 1 &&
		loc[2] >= 0 && loc[2] <= 1);
    }
}

bool Voxel8::GContains (const Point &glob, const NodeList&) const
{
    const double EPS = 1e-6;
    return (glob[0] >= x0-EPS && glob[0] <= x0+dx+EPS &&
	    glob[1] >= y0-EPS && glob[1] <= y0+dy+EPS &&
	    glob[2] >= z0-EPS && glob[2] <= z0+dz+EPS);
}

RVector Voxel8::LocalShapeF (const Point &loc) const
{
    dASSERT (loc.Dim() == 3, "Argument 1 invalid dimension");
    RVector fun(8);
    double iloc0 = 1.0-loc[0], iloc1 = 1.0-loc[1], iloc2 = 1.0-loc[2];

    fun[0] = iloc0  * iloc1  * iloc2;
    fun[1] = loc[0] * iloc1  * iloc2;
    fun[2] = iloc0  * loc[1] * iloc2;
    fun[3] = loc[0] * loc[1] * iloc2;
    fun[4] = iloc0  * iloc1  * loc[2];
    fun[5] = loc[0] * iloc1  * loc[2];
    fun[6] = iloc0  * loc[1] * loc[2];
    fun[7] = loc[0] * loc[1] * loc[2];
    return fun;
}

RDenseMatrix Voxel8::LocalShapeD (const Point &loc) const
{
    dASSERT (loc.Dim() == 3, "Argument 1 invalid dimension");

    double iloc0 = 1.0-loc[0], iloc1 = 1.0-loc[1], iloc2 = 1.0-loc[2];
    RDenseMatrix der(3,8);
    der(0,0) = -iloc1*iloc2;
    der(1,0) = -iloc0*iloc2;
    der(2,0) = -iloc0*iloc1;
    der(0,1) =  iloc1*iloc2;
    der(1,1) = -loc[0]*iloc2;
    der(2,1) = -loc[0]*iloc1;
    der(0,2) = -loc[1]*iloc2;
    der(1,2) =  iloc0*iloc2;
    der(2,2) = -iloc0*loc[1];
    der(0,3) =  loc[1]*iloc2;
    der(1,3) =  loc[0]*iloc2;
    der(2,3) = -loc[0]*loc[1];
    der(0,4) = -iloc1*loc[2];
    der(1,4) = -iloc0*loc[2];
    der(2,4) =  iloc0*iloc1;
    der(0,5) =  iloc1*loc[2];
    der(1,5) = -loc[0]*loc[2];
    der(2,5) =  loc[0]*iloc1;
    der(0,6) = -loc[1]*loc[2];
    der(1,6) =  iloc0*loc[2];
    der(2,6) =  iloc0*loc[1];
    der(0,7) =  loc[1]*loc[2];
    der(1,7) =  loc[0]*loc[2];
    der(2,7) =  loc[0]*loc[1];
    return der;
}

double Voxel8::Size() const
{
	return size;
}

double Voxel8::IntF (int i) const
{
	return intf;
}

double Voxel8::IntFF (int i, int j) const
{
	return intff(i,j);
}

RSymMatrix Voxel8::IntFF() const
{
	return intff;
}

double Voxel8::IntFFF (int i, int j, int k) const
{
	return intfff[i](j,k);
}

RSymMatrix Voxel8::IntDD () const
{
	return intdd;
}

double Voxel8::IntDD (int i, int j) const
{
	return intdd(i,j);
}

double Voxel8::IntFDD (int i, int j, int k) const
{
	return intfdd[i](j,k);
}

double Voxel8::IntPFF (int i, int j, const RVector &P) const
{
    return (intfff[0](i,j) * P[Node[0]] +
	    intfff[1](i,j) * P[Node[1]] +
	    intfff[2](i,j) * P[Node[2]] +
	    intfff[3](i,j) * P[Node[3]] +
	    intfff[4](i,j) * P[Node[4]] +
	    intfff[5](i,j) * P[Node[5]] +
	    intfff[6](i,j) * P[Node[6]] +
	    intfff[7](i,j) * P[Node[7]]);
}

RSymMatrix Voxel8::IntPFF (const RVector &P) const
{
    return (intfff[0] * P[Node[0]] +
	    intfff[1] * P[Node[1]] +
	    intfff[2] * P[Node[2]] +
	    intfff[3] * P[Node[3]] +
	    intfff[4] * P[Node[4]] +
	    intfff[5] * P[Node[5]] +
	    intfff[6] * P[Node[6]] +
	    intfff[7] * P[Node[7]]);
}

void Voxel8::ComputeIntFF () const
{
    static const RSymMatrix scale_intff = RSymMatrix (8,
      "8 \
       4 8 \
       4 2 8 \
       2 4 4 8 \
       4 2 2 1 8 \
       2 4 1 2 4 8 \
       2 1 4 2 4 2 8 \
       1 2 2 4 2 4 4 8") * (1.0/216.0);

    intff.New(8,8);
    intff = scale_intff * size;
}

void Voxel8::ComputeIntFFF () const
{
    static const RSymMatrix scale_intfff0 = RSymMatrix (8,
      "27  \
        9  9  \
        9  3  9  \
        3  3  3  3  \
        9  3  3  1  9  \
        3  3  1  1  3  3  \
        3  1  3  1  3  1  3  \
        1  1  1  1  1  1  1  1") * (1.0/1728.0);
    static const RSymMatrix scale_intfff1 = RSymMatrix (8,
      " 9  \
        9 27  \
        3  3  3  \
        3  9  3  9  \
        3  3  1  1  3  \
        3  9  1  3  3  9  \
        1  1  1  1  1  1  1  \
        1  3  1  3  1  3  1  3") * (1.0/1728.0);
    static const RSymMatrix scale_intfff2 = RSymMatrix (8,
      " 9  \
        3  3  \
        9  3 27  \
        3  3  9  9  \
        3  1  3  1  3  \
        1  1  1  1  1  1  \
        3  1  9  3  3  1  9  \
        1  1  3  3  1  1  3  3") * (1.0/1728.0);
    static const RSymMatrix scale_intfff3 = RSymMatrix (8,
      " 3  \
        3  9  \
        3  3  9  \
        3  9  9 27  \
        1  1  1  1  1  \
        1  3  1  3  1  3  \
        1  1  3  3  1  1  3  \
        1  3  3  9  1  3  3  9") * (1.0/1728.0);
    static const RSymMatrix scale_intfff4 = RSymMatrix (8,
      " 9  \
        3  3  \
        3  1  3  \
        1  1  1  1  \
        9  3  3  1 27  \
        3  3  1  1  9  9  \
        3  1  3  1  9  3  9  \
        1  1  1  1  3  3  3  3") * (1.0/1728.0);
    static const RSymMatrix scale_intfff5 = RSymMatrix (8,
      " 3  \
        3  9  \
        1  1  1  \
        1  3  1  3  \
        3  3  1  1  9  \
        3  9  1  3  9 27  \
        1  1  1  1  3  3  3  \
        1  3  1  3  3  9  3  9") * (1.0/1728.0);
    static const RSymMatrix scale_intfff6 = RSymMatrix (8,
      " 3  \
        1  1  \
        3  1  9  \
        1  1  3  3  \
        3  1  3  1  9  \
        1  1  1  1  3  3  \
        3  1  9  3  9  3 27  \
        1  1  3  3  3  3  9  9") * (1.0/1728.0);
    static const RSymMatrix scale_intfff7 = RSymMatrix (8,
      " 1  \
        1  3  \
        1  1  3  \
        1  3  3  9  \
        1  1  1  1  3  \
        1  3  1  3  3  9  \
        1  1  3  3  3  3  9  \
        1  3  3  9  3  9  9 27") * (1.0/1728.0);

    for (int i = 0; i < 8; i++)
        intfff[i].New(8,8);
    intfff[0] = scale_intfff0 * size;
    intfff[1] = scale_intfff1 * size;
    intfff[2] = scale_intfff2 * size;
    intfff[3] = scale_intfff3 * size;
    intfff[4] = scale_intfff4 * size;
    intfff[5] = scale_intfff5 * size;
    intfff[6] = scale_intfff6 * size;
    intfff[7] = scale_intfff7 * size;
}

void Voxel8::ComputeIntDD () const
{
    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;
    intdd.New(8,8);

    intdd(0,0) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(0,1) = 2*(-2*dy2*dz2 + dx2*(dy2 + dz2));
    intdd(0,2) = 2*(dy2*dz2 + dx2*(dy2 - 2*dz2));
    intdd(0,3) = -2*dy2*dz2 + dx2*(dy2 - 2*dz2);
    intdd(0,4) = 2*dy2*dz2 + dx2*(-4*dy2 + 2*dz2);
    intdd(0,5) = -2*dy2*dz2 + dx2*(-2*dy2 + dz2);
    intdd(0,6) = dy2*dz2 - 2*dx2*(dy2 + dz2);
    intdd(0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intdd(1,1) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(1,2) = -2*dy2*dz2 + dx2*(dy2 - 2*dz2);
    intdd(1,3) = 2*(dy2*dz2 + dx2*(dy2 - 2*dz2));
    intdd(1,4) = -2*dy2*dz2 + dx2*(-2*dy2 + dz2);
    intdd(1,5) = 2*dy2*dz2 + dx2*(-4*dy2 + 2*dz2);
    intdd(1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intdd(1,7) = dy2*dz2 - 2*dx2*(dy2 + dz2);
    intdd(2,2) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(2,3) = 2*(-2*dy2*dz2 + dx2*(dy2 + dz2));
    intdd(2,4) = dy2*dz2 - 2*dx2*(dy2 + dz2);
    intdd(2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intdd(2,6) = 2*dy2*dz2 + dx2*(-4*dy2 + 2*dz2);
    intdd(2,7) = -2*dy2*dz2 + dx2*(-2*dy2 + dz2);
    intdd(3,3) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intdd(3,5) = dy2*dz2 - 2*dx2*(dy2 + dz2);
    intdd(3,6) = -2*dy2*dz2 + dx2*(-2*dy2 + dz2);
    intdd(3,7) = 2*dy2*dz2 + dx2*(-4*dy2 + 2*dz2);
    intdd(4,4) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(4,5) = 2*(-2*dy2*dz2 + dx2*(dy2 + dz2));
    intdd(4,6) = 2*(dy2*dz2 + dx2*(dy2 - 2*dz2));
    intdd(4,7) = -2*dy2*dz2 + dx2*(dy2 - 2*dz2);
    intdd(5,5) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(5,6) = -2*dy2*dz2 + dx2*(dy2 - 2*dz2);
    intdd(5,7) = 2*(dy2*dz2 + dx2*(dy2 - 2*dz2));
    intdd(6,6) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd(6,7) = 2*(-2*dy2*dz2 + dx2*(dy2 + dz2));
    intdd(7,7) = 4*(dy2*dz2 + dx2*(dy2 + dz2));
    intdd /= 36*size;
}

void Voxel8::ComputeIntFDD () const
{
    int i;
    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;

    for (i = 0; i < 8; i++)
        intfdd[i].New(8,8);

    intfdd[0](0,0) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[0](0,1) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[0](0,2) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[0](0,3) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[0](0,4) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[0](0,5) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[0](0,6) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[0](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[0](1,1) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[0](1,2) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[0](1,3) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[0](1,4) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[0](1,5) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[0](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[0](1,7) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[0](2,2) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[0](2,3) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[0](2,4) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[0](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[0](2,6) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[0](2,7) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[0](3,3) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[0](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[0](3,5) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[0](3,6) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[0](3,7) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[0](4,4) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[0](4,5) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[0](4,6) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[0](4,7) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[0](5,5) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[0](5,6) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[0](5,7) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[0](6,6) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[0](6,7) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[0](7,7) = dy2*dz2 + dx2*(dy2 + dz2);

    intfdd[1](0,0) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[1](0,1) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[1](0,2) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[1](0,3) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[1](0,4) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[1](0,5) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[1](0,6) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[1](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[1](1,1) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[1](1,2) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[1](1,3) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[1](1,4) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[1](1,5) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[1](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[1](1,7) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[1](2,2) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[1](2,3) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[1](2,4) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[1](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[1](2,6) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[1](2,7) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[1](3,3) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[1](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[1](3,5) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[1](3,6) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[1](3,7) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[1](4,4) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[1](4,5) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[1](4,6) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[1](4,7) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[1](5,5) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[1](5,6) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[1](5,7) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[1](6,6) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[1](6,7) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[1](7,7) = dy2*dz2 + 3*dx2*(dy2 + dz2);

    intfdd[2](0,0) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[2](0,1) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[2](0,2) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[2](0,3) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[2](0,4) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[2](0,5) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[2](0,6) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[2](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](1,0) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[2](1,1) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[2](1,2) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[2](1,3) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[2](1,4) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[2](1,5) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[2](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](1,7) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[2](2,0) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[2](2,1) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[2](2,2) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[2](2,3) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[2](2,4) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[2](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](2,6) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[2](2,7) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[2](3,0) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[2](3,1) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[2](3,2) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[2](3,3) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[2](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](3,5) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[2](3,6) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[2](3,7) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[2](4,0) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[2](4,1) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[2](4,2) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[2](4,3) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](4,4) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[2](4,5) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[2](4,6) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[2](4,7) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[2](5,0) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[2](5,1) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[2](5,2) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](5,3) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[2](5,4) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[2](5,5) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[2](5,6) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[2](5,7) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[2](6,0) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[2](6,1) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](6,2) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[2](6,3) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[2](6,4) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[2](6,5) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[2](6,6) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[2](6,7) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[2](7,0) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[2](7,1) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[2](7,2) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[2](7,3) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[2](7,4) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[2](7,5) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[2](7,6) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[2](7,7) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);

    intfdd[3](0,0) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[3](0,1) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[3](0,2) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[3](0,3) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[3](0,4) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[3](0,5) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[3](0,6) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[3](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[3](1,1) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[3](1,2) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[3](1,3) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[3](1,4) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[3](1,5) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[3](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[3](1,7) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[3](2,2) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[3](2,3) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[3](2,4) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[3](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[3](2,6) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[3](2,7) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[3](3,3) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[3](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[3](3,5) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[3](3,6) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[3](3,7) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[3](4,4) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[3](4,5) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[3](4,6) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[3](4,7) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[3](5,5) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[3](5,6) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[3](5,7) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[3](6,6) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[3](6,7) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[3](7,7) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));

    intfdd[4](0,0) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[4](0,1) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[4](0,2) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[4](0,3) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[4](0,4) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[4](0,5) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[4](0,6) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[4](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[4](1,1) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[4](1,2) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[4](1,3) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[4](1,4) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[4](1,5) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[4](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[4](1,7) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[4](2,2) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[4](2,3) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[4](2,4) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[4](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[4](2,6) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[4](2,7) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[4](3,3) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[4](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[4](3,5) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[4](3,6) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[4](3,7) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[4](4,4) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[4](4,5) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[4](4,6) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[4](4,7) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[4](5,5) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[4](5,6) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[4](5,7) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[4](6,6) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[4](6,7) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[4](7,7) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);

    intfdd[5](0,0) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[5](0,1) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[5](0,2) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[5](0,3) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[5](0,4) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[5](0,5) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[5](0,6) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[5](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[5](1,1) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[5](1,2) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[5](1,3) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[5](1,4) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[5](1,5) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[5](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[5](1,7) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[5](2,2) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[5](2,3) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[5](2,4) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[5](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[5](2,6) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[5](2,7) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[5](3,3) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[5](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[5](3,5) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[5](3,6) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[5](3,7) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[5](4,4) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[5](4,5) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[5](4,6) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[5](4,7) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[5](5,5) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[5](5,6) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[5](5,7) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[5](6,6) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[5](6,7) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[5](7,7) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));

    intfdd[6](0,0) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[6](0,1) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[6](0,2) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[6](0,3) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[6](0,4) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[6](0,5) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[6](0,6) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[6](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[6](1,1) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[6](1,2) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[6](1,3) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[6](1,4) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[6](1,5) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[6](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[6](1,7) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[6](2,2) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[6](2,3) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[6](2,4) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[6](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[6](2,6) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[6](2,7) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[6](3,3) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[6](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[6](3,5) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[6](3,6) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[6](3,7) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[6](4,4) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[6](4,5) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[6](4,6) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[6](4,7) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[6](5,5) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[6](5,6) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[6](5,7) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[6](6,6) = 9*(dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[6](6,7) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[6](7,7) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));

    intfdd[7](0,0) = dy2*dz2 + dx2*(dy2 + dz2);
    intfdd[7](0,1) = -(dy2*dz2) + dx2*(dy2 + dz2);
    intfdd[7](0,2) = dy2*dz2 + dx2*(dy2 - dz2);
    intfdd[7](0,3) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[7](0,4) = dy2*dz2 + dx2*(-dy2 + dz2);
    intfdd[7](0,5) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[7](0,6) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[7](0,7) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[7](1,1) = dy2*dz2 + 3*dx2*(dy2 + dz2);
    intfdd[7](1,2) = -(dy2*dz2) + dx2*(dy2 - dz2);
    intfdd[7](1,3) = dy2*dz2 + 3*dx2*(dy2 - dz2);
    intfdd[7](1,4) = -(dy2*dz2) + dx2*(-dy2 + dz2);
    intfdd[7](1,5) = dy2*dz2 - 3*dx2*(dy2 - dz2);
    intfdd[7](1,6) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[7](1,7) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[7](2,2) = 3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[7](2,3) = -3*dy2*dz2 + dx2*(3*dy2 + dz2);
    intfdd[7](2,4) = dy2*dz2 - dx2*(dy2 + dz2);
    intfdd[7](2,5) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[7](2,6) = 3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[7](2,7) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[7](3,3) = 3*(dy2*dz2 + dx2*(3*dy2 + dz2));
    intfdd[7](3,4) = -(dy2*dz2) - dx2*(dy2 + dz2);
    intfdd[7](3,5) = dy2*dz2 - 3*dx2*(dy2 + dz2);
    intfdd[7](3,6) = -3*dy2*dz2 + dx2*(-3*dy2 + dz2);
    intfdd[7](3,7) = 3*(dy2*dz2 + dx2*(-3*dy2 + dz2));
    intfdd[7](4,4) = 3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[7](4,5) = -3*dy2*dz2 + dx2*(dy2 + 3*dz2);
    intfdd[7](4,6) = 3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[7](4,7) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[7](5,5) = 3*(dy2*dz2 + dx2*(dy2 + 3*dz2));
    intfdd[7](5,6) = -3*dy2*dz2 + dx2*(dy2 - 3*dz2);
    intfdd[7](5,7) = 3*(dy2*dz2 + dx2*(dy2 - 3*dz2));
    intfdd[7](6,6) = 3*(3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[7](6,7) = 3*(-3*dy2*dz2 + dx2*(dy2 + dz2));
    intfdd[7](7,7) = 9*(dy2*dz2 + dx2*(dy2 + dz2));

    for (i = 0; i < 8; i++)
        intfdd[i] /= 288*dx*dy*dz;
}

void Voxel8::ComputeBndIntF () const
{
    double scale;

    // side 0: z = 0
    scale = dx*dy/4.0;
    bndintf[0].New(8);
    bndintf[0][0] = scale;
    bndintf[0][1] = scale;
    bndintf[0][2] = scale;
    bndintf[0][3] = scale;

    // side 1: z = dz
    bndintf[1].New(8);
    bndintf[1][4] = scale;
    bndintf[1][5] = scale;
    bndintf[1][6] = scale;
    bndintf[1][7] = scale;

    // side 2: y = 0
    scale = dx*dz/4.0;
    bndintf[2].New(8);
    bndintf[2][0] = scale;
    bndintf[2][1] = scale;
    bndintf[2][4] = scale;
    bndintf[2][5] = scale;

    // side 3: y = dy
    bndintf[3].New(8);
    bndintf[3][2] = scale;
    bndintf[3][3] = scale;
    bndintf[3][6] = scale;
    bndintf[3][7] = scale;

    // side 4: x = 0
    scale = dy*dz/4.0;
    bndintf[4].New(8);
    bndintf[4][0] = scale;
    bndintf[4][2] = scale;
    bndintf[4][4] = scale;
    bndintf[4][6] = scale;

    // side 5: x = dx
    bndintf[5].New(8);
    bndintf[5][1] = scale;
    bndintf[5][3] = scale;
    bndintf[5][5] = scale;
    bndintf[5][7] = scale;
}

void Voxel8::ComputeBndIntFF () const
{
    double scale;

    // side 0: z = 0
    scale = dx*dy/36.0;
    bndintff[0].New(8,8);
    bndintff[0](0,0) = 4*scale;
    bndintff[0](0,1) = 2*scale;
    bndintff[0](0,2) = 2*scale;
    bndintff[0](0,3) = scale;
    bndintff[0](1,1) = 4*scale;
    bndintff[0](1,2) = scale;
    bndintff[0](1,3) = 2*scale;
    bndintff[0](2,2) = 4*scale;
    bndintff[0](2,3) = 2*scale;
    bndintff[0](3,3) = 4*scale;

    // side 1: z = dz
    bndintff[1].New(8,8);
    bndintff[1](4,4) = 4*scale;
    bndintff[1](4,5) = 2*scale;
    bndintff[1](4,6) = 2*scale;
    bndintff[1](4,7) = scale;
    bndintff[1](5,5) = 4*scale;
    bndintff[1](5,6) = scale;
    bndintff[1](5,7) = 2*scale;
    bndintff[1](6,6) = 4*scale;
    bndintff[1](6,7) = 2*scale;
    bndintff[1](7,7) = 4*scale;

    // side 2: y = 0
    scale = dx*dz/36.0;
    bndintff[2].New(8,8);
    bndintff[2](0,0) = 4*scale;
    bndintff[2](0,1) = 2*scale;
    bndintff[2](0,4) = 2*scale;
    bndintff[2](0,5) = scale;
    bndintff[2](1,1) = 4*scale;
    bndintff[2](1,4) = scale;
    bndintff[2](1,5) = 2*scale;
    bndintff[2](4,4) = 4*scale;
    bndintff[2](4,5) = 2*scale;
    bndintff[2](5,5) = 4*scale;

    // side 3: y = dy
    bndintff[3].New(8,8);
    bndintff[3](2,2) = 4*scale;
    bndintff[3](2,3) = 2*scale;
    bndintff[3](2,6) = 2*scale;
    bndintff[3](2,7) = scale;
    bndintff[3](3,3) = 4*scale;
    bndintff[3](3,6) = scale;
    bndintff[3](3,7) = 2*scale;
    bndintff[3](6,6) = 4*scale;
    bndintff[3](6,7) = 2*scale;
    bndintff[3](7,7) = 4*scale;

    // side 4: x = 0
    scale = dy*dz/36.0;
    bndintff[4].New(8,8);
    bndintff[4](0,0) = 4*scale;
    bndintff[4](0,2) = 2*scale;
    bndintff[4](0,4) = 2*scale;
    bndintff[4](0,6) = scale;
    bndintff[4](2,2) = 4*scale;
    bndintff[4](2,4) = scale;
    bndintff[4](2,6) = 2*scale;
    bndintff[4](4,4) = 4*scale;
    bndintff[4](4,6) = 2*scale;
    bndintff[4](6,6) = 4*scale;

    // side 5: x = dx
    bndintff[5].New(8,8);
    bndintff[5](1,1) = 4*scale;
    bndintff[5](1,3) = 2*scale;
    bndintff[5](1,5) = 2*scale;
    bndintff[5](1,7) = scale;
    bndintff[5](3,3) = 4*scale;
    bndintff[5](3,5) = scale;
    bndintff[5](3,7) = 2*scale;
    bndintff[5](5,5) = 4*scale;
    bndintff[5](5,7) = 2*scale;
    bndintff[5](7,7) = 4*scale;
}

void Voxel8::ComputeBndIntFFF () const
{
  char scales[6][8][8][8] = {
  {{{9,3,3,1,0,0,0,0},{3,3,1,1,0,0,0,0},{3,1,3,1,0,0,0,0},{1,1,1,1,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{3,3,1,1,0,0,0,0},{3,9,1,3,0,0,0,0},{1,1,1,1,0,0,0,0},{1,3,1,3,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{3,1,3,1,0,0,0,0},{1,1,1,1,0,0,0,0},{3,1,9,3,0,0,0,0},{1,1,3,3,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{1,1,1,1,0,0,0,0},{1,3,1,3,0,0,0,0},{1,1,3,3,0,0,0,0},{1,3,3,9,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}},
  
  {{{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,9,3,3,1},{0,0,0,0,3,3,1,1},{0,0,0,0,3,1,3,1},{0,0,0,0,1,1,1,1}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,3,3,1,1},{0,0,0,0,3,9,1,3},{0,0,0,0,1,1,1,1},{0,0,0,0,1,3,1,3}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,3,1,3,1},{0,0,0,0,1,1,1,1},{0,0,0,0,3,1,9,3},{0,0,0,0,1,1,3,3}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,1,1,1,1},{0,0,0,0,1,3,1,3},{0,0,0,0,1,1,3,3},{0,0,0,0,1,3,3,9}}},

  {{{9,3,0,0,3,1,0,0},{3,3,0,0,1,1,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {3,1,0,0,3,1,0,0},{1,1,0,0,1,1,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{3,3,0,0,1,1,0,0},{3,9,0,0,1,3,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {1,1,0,0,1,1,0,0},{1,3,0,0,1,3,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{3,1,0,0,3,1,0,0},{1,1,0,0,1,1,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {3,1,0,0,9,3,0,0},{1,1,0,0,3,3,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{1,1,0,0,1,1,0,0},{1,3,0,0,1,3,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {1,1,0,0,3,3,0,0},{1,3,0,0,3,9,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}},

  {{{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,9,3,0,0,3,1},{0,0,3,3,0,0,1,1},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,3,1,0,0,3,1},{0,0,1,1,0,0,1,1}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,3,3,0,0,1,1},{0,0,3,9,0,0,1,3},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,1,1,0,0,1,1},{0,0,1,3,0,0,1,3}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,3,1,0,0,3,1},{0,0,1,1,0,0,1,1},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,3,1,0,0,9,3},{0,0,1,1,0,0,3,3}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,1,1,0,0,1,1},{0,0,1,3,0,0,1,3},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,1,1,0,0,3,3},{0,0,1,3,0,0,3,9}}},

  {{{9,0,3,0,3,0,1,0},{0,0,0,0,0,0,0,0},{3,0,3,0,1,0,1,0},{0,0,0,0,0,0,0,0},
    {3,0,1,0,3,0,1,0},{0,0,0,0,0,0,0,0},{1,0,1,0,1,0,1,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{3,0,3,0,1,0,1,0},{0,0,0,0,0,0,0,0},{3,0,9,0,1,0,3,0},{0,0,0,0,0,0,0,0},
    {1,0,1,0,1,0,1,0},{0,0,0,0,0,0,0,0},{1,0,3,0,1,0,3,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{3,0,1,0,3,0,1,0},{0,0,0,0,0,0,0,0},{1,0,1,0,1,0,1,0},{0,0,0,0,0,0,0,0},
    {3,0,1,0,9,0,3,0},{0,0,0,0,0,0,0,0},{1,0,1,0,3,0,3,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{1,0,1,0,1,0,1,0},{0,0,0,0,0,0,0,0},{1,0,3,0,1,0,3,0},{0,0,0,0,0,0,0,0},
    {1,0,1,0,3,0,3,0},{0,0,0,0,0,0,0,0},{1,0,3,0,3,0,9,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}},
  
  {{{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,9,0,3,0,3,0,1},{0,0,0,0,0,0,0,0},{0,3,0,3,0,1,0,1},
    {0,0,0,0,0,0,0,0},{0,3,0,1,0,3,0,1},{0,0,0,0,0,0,0,0},{0,1,0,1,0,1,0,1}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,3,0,3,0,1,0,1},{0,0,0,0,0,0,0,0},{0,3,0,9,0,1,0,3},
    {0,0,0,0,0,0,0,0},{0,1,0,1,0,1,0,1},{0,0,0,0,0,0,0,0},{0,1,0,3,0,1,0,3}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,3,0,1,0,3,0,1},{0,0,0,0,0,0,0,0},{0,1,0,1,0,1,0,1},
    {0,0,0,0,0,0,0,0},{0,3,0,1,0,9,0,3},{0,0,0,0,0,0,0,0},{0,1,0,1,0,3,0,3}},
   {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}},
   {{0,0,0,0,0,0,0,0},{0,1,0,1,0,1,0,1},{0,0,0,0,0,0,0,0},{0,1,0,3,0,1,0,3},
    {0,0,0,0,0,0,0,0},{0,1,0,1,0,3,0,3},{0,0,0,0,0,0,0,0},{0,1,0,3,0,3,0,9}}}
  };

    int i, j, k, sd;
    double area[6], factor;

    area[0] = area[1] = dx*dy;
    area[2] = area[3] = dx*dz;
    area[4] = area[5] = dy*dz;

    for (sd = 0; sd < 6; sd++) { // loop over sides
        factor = area[sd] / 144.0;
        for (i = 0; i < 8; i++) {
	    bndintfff[sd][i].New(8,8);
	    for (j = 0; j < 8; j++)
	        for (k = 0; k < 8; k++)
		    bndintfff[sd][i](j,k) = factor*scales[sd][i][j][k];
	}
    }
}

double Voxel8::IntPDD (int i, int j, const RVector &P) const
{
    double res = 0.0;
    for (int k = 0; k < 8; k++) res += P[Node[k]] * intfdd[k](i,j);
    return res;
}

RSymMatrix Voxel8::IntPDD (const RVector &P) const
{
    RSymMatrix pdd(8);
    for (int i = 0; i < 8; i++)
        for (int j = i; j < 8; j++)
	    pdd(i,j) = IntPDD (i,j,P);
    return pdd;
}

RVector Voxel8::BndIntF () const
{
    RVector intf(8);
    for (int sd = 0; sd < 6; sd++) {
	if (bndside[sd])
	    for (int i = 0; i < 8; i++)
		intf[i] += bndintf[sd][i];
    }
    return intf;
}

double Voxel8::BndIntF (int i)
{
    dASSERT(i >= 0 && i < 8, "Argument 1: out of range");
    if (!bndel) return 0.0;
    double res = 0.0;
    for (int sd = 0; sd < 6; sd++) {
	if (bndside[sd])
	    res += bndintf[sd][i];
    }
    return res;
}

double Voxel8::SurfIntF (int i, int sd) const
{
    dASSERT(i >= 0 && i < 8, "Argument 1: out of range");
    dASSERT(sd >= 0 && sd < 6, "Argument 3: out of range");
    return bndintf[sd][i];
}

double Voxel8::BndIntFF (int i, int j)
{
    if (!bndel) return 0.0;
    double res = 0.0;
    for (int sd = 0; sd < 6; sd++) {
        if (bndside[sd])
	    res += bndintff[sd](i,j);
    }
    return res;
}

double Voxel8::SurfIntFF (int i, int j, int sd) const
{
    dASSERT(i >= 0 && i < 8, "Argument 1: out of range");
    dASSERT(j >= 0 && j < 8, "Argument 2: out of range");
    dASSERT(sd >= 0 && sd < 6, "Argument 3: out of range");
    return bndintff[sd](i,j);

}
double Voxel8::BndIntPFF (int i, int j, const RVector &P) const
{
    if (!bndel) return 0.0;
    double res = 0.0;
    for (int sd = 0; sd < 6; sd++) {
        if (bndside[sd]) {
	    for (int k = 0; k < 8; k++)
	        res += P[Node[k]] * bndintfff[sd][k](i,j);
	}
    }
    return res;
}

double Voxel8::IntFd (int i, int j, int k) const
{
    static RDenseMatrix fid;
    if (fid.nRows() == 0) {
	fid.New (24,8);
	fid(0,0) = -4*dy*dz,  fid(1,0) = -4*dx*dz,  fid(2,0) = -4*dx*dy;
	fid(0,1) =  4*dy*dz,  fid(1,1) = -2*dx*dz,  fid(2,1) = -2*dx*dy;
	fid(0,2) = -2*dy*dz,  fid(1,2) =  4*dx*dz,  fid(2,2) = -2*dx*dy;
	fid(0,3) =  2*dy*dz,  fid(1,3) =  2*dx*dz,  fid(2,3) = -(dx*dy);
	fid(0,4) = -2*dy*dz,  fid(1,4) = -2*dx*dz,  fid(2,4) =  4*dx*dy;
	fid(0,5) =  2*dy*dz,  fid(1,5) = -(dx*dz),  fid(2,5) =  2*dx*dy;
	fid(0,6) = -(dy*dz),  fid(1,6) =  2*dx*dz,  fid(2,6) =  2*dx*dy;
	fid(0,7) =  dy*dz,    fid(1,7) =  dx*dz,    fid(2,7) =  dx*dy;

	fid(3,0) = -4*dy*dz,  fid(4,0) = -2*dx*dz,  fid(5,0) = -2*dx*dy;
	fid(3,1) =  4*dy*dz,  fid(4,1) = -4*dx*dz,  fid(5,1) = -4*dx*dy;
	fid(3,2) = -2*dy*dz,  fid(4,2) =  2*dx*dz,  fid(5,2) = -(dx*dy);
	fid(3,3) =  2*dy*dz,  fid(4,3) =  4*dx*dz,  fid(5,3) = -2*dx*dy;
	fid(3,4) = -2*dy*dz,  fid(4,4) = -(dx*dz),  fid(5,4) =  2*dx*dy;
	fid(3,5) =  2*dy*dz,  fid(4,5) = -2*dx*dz,  fid(5,5) =  4*dx*dy;
	fid(3,6) = -(dy*dz),  fid(4,6) =  dx*dz,    fid(5,6) =  dx*dy;
	fid(3,7) =  dy*dz,    fid(4,7) =  2*dx*dz,  fid(5,7) =  2*dx*dy;

	fid(6,0) = -2*dy*dz,  fid(7,0) = -4*dx*dz,  fid(8,0) = -2*dx*dy;
	fid(6,1) =  2*dy*dz,  fid(7,1) = -2*dx*dz,  fid(8,1) = -(dx*dy);
	fid(6,2) = -4*dy*dz,  fid(7,2) =  4*dx*dz,  fid(8,2) = -4*dx*dy;
	fid(6,3) =  4*dy*dz,  fid(7,3) =  2*dx*dz,  fid(8,3) = -2*dx*dy;
	fid(6,4) = -(dy*dz),  fid(7,4) = -2*dx*dz,  fid(8,4) =  2*dx*dy;
	fid(6,5) =  dy*dz,    fid(7,5) = -(dx*dz),  fid(8,5) =  dx*dy;
	fid(6,6) = -2*dy*dz,  fid(7,6) =  2*dx*dz,  fid(8,6) =  4*dx*dy;
	fid(6,7) =  2*dy*dz,  fid(7,7) =  dx*dz,    fid(8,7) =  2*dx*dy;

	fid(9,0) = -2*dy*dz,  fid(10,0) = -2*dx*dz, fid(11,0) = -(dx*dy);
	fid(9,1) =  2*dy*dz,  fid(10,1) = -4*dx*dz, fid(11,1) = -2*dx*dy;
	fid(9,2) = -4*dy*dz,  fid(10,2) =  2*dx*dz, fid(11,2) = -2*dx*dy;
	fid(9,3) =  4*dy*dz,  fid(10,3) =  4*dx*dz, fid(11,3) = -4*dx*dy;
	fid(9,4) = -(dy*dz),  fid(10,4) = -(dx*dz), fid(11,4) =  dx*dy;
	fid(9,5) =  dy*dz,    fid(10,5) = -2*dx*dz, fid(11,5) =  2*dx*dy;
	fid(9,6) = -2*dy*dz,  fid(10,6) =  dx*dz,   fid(11,6) =  2*dx*dy;
	fid(9,7) =  2*dy*dz,  fid(10,7) =  2*dx*dz, fid(11,7) =  4*dx*dy;

	fid(12,0) = -2*dy*dz, fid(13,0) = -2*dx*dz, fid(14,0) = -4*dx*dy;
	fid(12,1) =  2*dy*dz, fid(13,1) = -(dx*dz), fid(14,1) = -2*dx*dy;
	fid(12,2) = -(dy*dz), fid(13,2) =  2*dx*dz, fid(14,2) = -2*dx*dy;
	fid(12,3) =  dy*dz,   fid(13,3) =  dx*dz,   fid(14,3) = -(dx*dy);
	fid(12,4) = -4*dy*dz, fid(13,4) = -4*dx*dz, fid(14,4) =  4*dx*dy;
	fid(12,5) =  4*dy*dz, fid(13,5) = -2*dx*dz, fid(14,5) =  2*dx*dy;
	fid(12,6) = -2*dy*dz, fid(13,6) =  4*dx*dz, fid(14,6) =  2*dx*dy;
	fid(12,7) =  2*dy*dz, fid(13,7) =  2*dx*dz, fid(14,7) =  dx*dy;

	fid(15,0) = -2*dy*dz, fid(16,0) = -(dx*dz), fid(17,0) = -2*dx*dy;
	fid(15,1) =  2*dy*dz, fid(16,1) = -2*dx*dz, fid(17,1) = -4*dx*dy;
	fid(15,2) = -(dy*dz), fid(16,2) =  dx*dz,   fid(17,2) = -(dx*dy);
	fid(15,3) =  dy*dz,   fid(16,3) =  2*dx*dz, fid(17,3) = -2*dx*dy;
	fid(15,4) = -4*dy*dz, fid(16,4) = -2*dx*dz, fid(17,4) =  2*dx*dy;
	fid(15,5) =  4*dy*dz, fid(16,5) = -4*dx*dz, fid(17,5) =  4*dx*dy;
	fid(15,6) = -2*dy*dz, fid(16,6) =  2*dx*dz, fid(17,6) =  dx*dy;
	fid(15,7) =  2*dy*dz, fid(16,7) =  4*dx*dz, fid(17,7) =  2*dx*dy;

	fid(18,0) = -(dy*dz), fid(19,0) = -2*dx*dz, fid(20,0) = -2*dx*dy;
	fid(18,1) =  dy*dz,   fid(19,1) = -(dx*dz), fid(20,1) = -(dx*dy);
	fid(18,2) = -2*dy*dz, fid(19,2) =  2*dx*dz, fid(20,2) = -4*dx*dy;
	fid(18,3) =  2*dy*dz, fid(19,3) =  dx*dz,   fid(20,3) = -2*dx*dy;
	fid(18,4) = -2*dy*dz, fid(19,4) = -4*dx*dz, fid(20,4) =  2*dx*dy;
	fid(18,5) =  2*dy*dz, fid(19,5) = -2*dx*dz, fid(20,5) =  dx*dy;
	fid(18,6) = -4*dy*dz, fid(19,6) =  4*dx*dz, fid(20,6) =  4*dx*dy;
	fid(18,7) =  4*dy*dz, fid(19,7) =  2*dx*dz, fid(20,7) =  2*dx*dy;

	fid(21,0) = -(dy*dz), fid(22,0) = -(dx*dz), fid(23,0) = -(dx*dy);
	fid(21,1) =  dy*dz,   fid(22,1) = -2*dx*dz, fid(23,1) = -2*dx*dy;
	fid(21,2) = -2*dy*dz, fid(22,2) =  dx*dz,   fid(23,2) = -2*dx*dy;
	fid(21,3) =  2*dy*dz, fid(22,3) =  2*dx*dz, fid(23,3) = -4*dx*dy;
	fid(21,4) = -2*dy*dz, fid(22,4) = -2*dx*dz, fid(23,4) =  dx*dy;
	fid(21,5) =  2*dy*dz, fid(22,5) = -4*dx*dz, fid(23,5) =  2*dx*dy;
	fid(21,6) = -4*dy*dz, fid(22,6) =  2*dx*dz, fid(23,6) =  2*dx*dy;
	fid(21,7) =  4*dy*dz, fid(22,7) =  4*dx*dz, fid(23,7) =  4*dx*dy;

	fid *= 1.0/72.0;
    }

    return fid(i*3+k,j);
}

RSymMatrix Voxel8::Intdd () const
{
    static RSymMatrix *MixedDD = 0;
    if (!MixedDD) {
	MixedDD = new RSymMatrix (24,24);
	RSymMatrix &MDD = *MixedDD;
	double dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz;
	double dxy = dx*dy, dxz = dx*dz, dyz = dy*dz;
	double dxyz2 = dxy*dz2, dxy2z = dxz*dy2, dx2yz = dyz*dx2;
	double dy2z2 = dy2*dz2, dx2z2 = dx2*dz2, dx2y2 = dx2*dy2;

	MDD( 0, 0) = 8*dy2z2;
	MDD( 1, 0) = 6*dxyz2; MDD( 1, 1) = 8*dx2z2;
	MDD( 2, 0) = 6*dxy2z; MDD( 2, 1) = 6*dx2yz; MDD( 2, 2) = 8*dx2y2;

	MDD( 3, 0) =-8*dy2z2; MDD( 3, 1) =-6*dxyz2; MDD( 3, 2) =-6*dxy2z;
	MDD( 4, 0) = 6*dxyz2; MDD( 4, 1) = 4*dx2z2; MDD( 4, 2) = 3*dx2yz;
	MDD( 5, 0) = 6*dxy2z; MDD( 5, 1) = 3*dx2yz; MDD( 5, 2) = 4*dx2y2;

	MDD( 3, 3) = 8*dy2z2;
	MDD( 4, 3) =-6*dxyz2; MDD( 4, 4) = 8*dx2z2;
	MDD( 5, 3) =-6*dxy2z; MDD( 5, 4) = 6*dx2yz; MDD( 5, 5) = 8*dx2y2;

	MDD( 6, 0) = 4*dy2z2; MDD( 6, 1) = 6*dxyz2; MDD( 6, 2) = 3*dxy2z;
	MDD( 7, 0) =-6*dxyz2; MDD( 7, 1) =-8*dx2z2; MDD( 7, 2) =-6*dx2yz;
	MDD( 8, 0) = 3*dxy2z; MDD( 8, 1) = 6*dx2yz; MDD( 8, 2) = 4*dx2y2;

	MDD( 6, 3) =-4*dy2z2; MDD( 6, 4) = 6*dxyz2; MDD( 6, 5) = 3*dxy2z;
	MDD( 7, 3) = 6*dxyz2; MDD( 7, 4) =-4*dx2z2; MDD( 7, 5) =-3*dx2yz;
	MDD( 8, 3) =-3*dxy2z; MDD( 8, 4) = 3*dx2yz; MDD( 8, 5) = 2*dx2y2;

	MDD( 6, 6) = 8*dy2z2;
	MDD( 7, 6) =-6*dxyz2; MDD( 7, 7) = 8*dx2z2;
	MDD( 8, 6) = 6*dxy2z; MDD( 8, 7) =-6*dx2yz; MDD( 8, 8) = 8*dx2y2;

	MDD( 9, 0) =-4*dy2z2; MDD( 9, 1) =-6*dxyz2; MDD( 9, 2) =-3*dxy2z;
	MDD(10, 0) =-6*dxyz2; MDD(10, 1) =-4*dx2z2; MDD(10, 2) =-3*dx2yz;
	MDD(11, 0) = 3*dxy2z; MDD(11, 1) = 3*dx2yz; MDD(11, 2) = 2*dx2y2;

	MDD( 9, 3) = 4*dy2z2; MDD( 9, 4) =-6*dxyz2; MDD( 9, 5) =-3*dxy2z;
	MDD(10, 3) = 6*dxyz2; MDD(10, 4) =-8*dx2z2; MDD(10, 5) =-6*dx2yz;
	MDD(11, 3) =-3*dxy2z; MDD(11, 4) = 6*dx2yz; MDD(11, 5) = 4*dx2y2;

	MDD( 9, 6) =-8*dy2z2; MDD( 9, 7) = 6*dxyz2; MDD( 9, 8) =-6*dxy2z;
	MDD(10, 6) =-6*dxyz2; MDD(10, 7) = 4*dx2z2; MDD(10, 8) =-3*dx2yz;
	MDD(11, 6) = 6*dxy2z; MDD(11, 7) =-3*dx2yz; MDD(11, 8) = 4*dx2y2;

	MDD( 9, 9) = 8*dy2z2;
	MDD(10, 9) = 6*dxyz2; MDD(10,10) = 8*dx2z2;
	MDD(11, 9) =-6*dxy2z; MDD(11,10) =-6*dx2yz; MDD(11,11) = 8*dx2y2;

	MDD(12, 0) = 4*dy2z2; MDD(12, 1) = 3*dxyz2; MDD(12, 2) = 6*dxy2z;
	MDD(13, 0) = 3*dxyz2; MDD(13, 1) = 4*dx2z2; MDD(13, 2) = 6*dx2yz;
	MDD(14, 0) =-6*dxy2z; MDD(14, 1) =-6*dx2yz; MDD(14, 2) =-8*dx2y2;

	MDD(12, 3) =-4*dy2z2; MDD(12, 4) = 3*dxyz2; MDD(12, 5) = 6*dxy2z;
	MDD(13, 3) =-3*dxyz2; MDD(13, 4) = 2*dx2z2; MDD(13, 5) = 3*dx2yz;
	MDD(14, 3) = 6*dxy2z; MDD(14, 4) =-3*dx2yz; MDD(14, 5) =-4*dx2y2;

	MDD(12, 6) = 2*dy2z2; MDD(12, 7) =-3*dxyz2; MDD(12, 8) = 3*dxy2z;
	MDD(13, 6) = 3*dxyz2; MDD(13, 7) =-4*dx2z2; MDD(13, 8) = 6*dx2yz;
	MDD(14, 6) =-3*dxy2z; MDD(14, 7) = 6*dx2yz; MDD(14, 8) =-4*dx2y2;

	MDD(12, 9) =-2*dy2z2; MDD(12,10) =-3*dxyz2; MDD(12,11) = 3*dxy2z;
	MDD(13, 9) =-3*dxyz2; MDD(13,10) =-2*dx2z2; MDD(13,11) = 3*dx2yz;
	MDD(14, 9) = 3*dxy2z; MDD(14,10) = 3*dx2yz; MDD(14,11) =-2*dx2y2;

	MDD(12,12) = 8*dy2z2;
	MDD(13,12) = 6*dxyz2; MDD(13,13) = 8*dx2z2;
	MDD(14,12) =-6*dxy2z; MDD(14,13) =-6*dx2yz; MDD(14,14) = 8*dx2y2;

	MDD(15, 0) =-4*dy2z2; MDD(15, 1) =-3*dxyz2; MDD(15, 2) =-6*dxy2z;
	MDD(16, 0) = 3*dxyz2; MDD(16, 1) = 2*dx2z2; MDD(16, 2) = 3*dx2yz;
	MDD(17, 0) =-6*dxy2z; MDD(17, 1) =-3*dx2yz; MDD(17, 2) =-4*dx2y2;

	MDD(15, 3) = 4*dy2z2; MDD(15, 4) =-3*dxyz2; MDD(15, 5) =-6*dxy2z;
	MDD(16, 3) =-3*dxyz2; MDD(16, 4) = 4*dx2z2; MDD(16, 5) = 6*dx2yz;
	MDD(17, 3) = 6*dxy2z; MDD(17, 4) =-6*dx2yz; MDD(17, 5) =-8*dx2y2;

	MDD(15, 6) =-2*dy2z2; MDD(15, 7) = 3*dxyz2; MDD(15, 8) =-3*dxy2z;
	MDD(16, 6) = 3*dxyz2; MDD(16, 7) =-2*dx2z2; MDD(16, 8) = 3*dx2yz;
	MDD(17, 6) =-3*dxy2z; MDD(17, 7) = 3*dx2yz; MDD(17, 8) =-2*dx2y2;

	MDD(15, 9) = 2*dy2z2; MDD(15,10) = 3*dxyz2; MDD(15,11) =-3*dxy2z;
	MDD(16, 9) =-3*dxyz2; MDD(16,10) =-4*dx2z2; MDD(16,11) = 6*dx2yz;
	MDD(17, 9) = 3*dxy2z; MDD(17,10) = 6*dx2yz; MDD(17,11) =-4*dx2y2;

	MDD(15,12) =-8*dy2z2; MDD(15,13) =-6*dxyz2; MDD(15,14) = 6*dxy2z;
	MDD(16,12) = 6*dxyz2; MDD(16,13) = 4*dx2z2; MDD(16,14) =-3*dx2yz;
	MDD(17,12) =-6*dxy2z; MDD(17,13) =-3*dx2yz; MDD(17,14) = 4*dx2y2;

	MDD(15,15) = 8*dy2z2;
	MDD(16,15) =-6*dxyz2; MDD(16,16) = 8*dx2z2;
	MDD(17,15) = 6*dxy2z; MDD(17,16) =-6*dx2yz; MDD(17,17) = 8*dx2y2;

	MDD(18, 0) = 2*dy2z2; MDD(18, 1) = 3*dxyz2; MDD(18, 2) = 3*dxy2z;
	MDD(19, 0) =-3*dxyz2; MDD(19, 1) =-4*dx2z2; MDD(19, 2) =-6*dx2yz;
	MDD(20, 0) =-3*dxy2z; MDD(20, 1) =-6*dx2yz; MDD(20, 2) =-4*dx2y2;

	MDD(18, 3) =-2*dy2z2; MDD(18, 4) = 3*dxyz2; MDD(18, 5) = 3*dxy2z;
	MDD(19, 3) = 3*dxyz2; MDD(19, 4) =-2*dx2z2; MDD(19, 5) =-3*dx2yz;
	MDD(20, 3) = 3*dxy2z; MDD(20, 4) =-3*dx2yz; MDD(20, 5) =-2*dx2y2;

	MDD(18, 6) = 4*dy2z2; MDD(18, 7) =-3*dxyz2; MDD(18, 8) = 6*dxy2z;
	MDD(19, 6) =-3*dxyz2; MDD(19, 7) = 4*dx2z2; MDD(19, 8) =-6*dx2yz;
	MDD(20, 6) =-6*dxy2z; MDD(20, 7) = 6*dx2yz; MDD(20, 8) =-8*dx2y2;

	MDD(18, 9) =-4*dy2z2; MDD(18,10) =-3*dxyz2; MDD(18,11) = 6*dxy2z;
	MDD(19, 9) = 3*dxyz2; MDD(19,10) = 2*dx2z2; MDD(19,11) =-3*dx2yz;
	MDD(20, 9) = 6*dxy2z; MDD(20,10) = 3*dx2yz; MDD(20,11) =-4*dx2y2;

	MDD(18,12) = 4*dy2z2; MDD(18,13) = 6*dxyz2; MDD(18,14) =-3*dxy2z;
	MDD(19,12) =-6*dxyz2; MDD(19,13) =-8*dx2z2; MDD(19,14) = 6*dx2yz;
	MDD(20,12) =-3*dxy2z; MDD(20,13) =-6*dx2yz; MDD(20,14) = 4*dx2y2;

	MDD(18,15) =-4*dy2z2; MDD(18,16) = 6*dxyz2; MDD(18,17) =-3*dxy2z;
	MDD(19,15) = 6*dxyz2; MDD(19,16) =-4*dx2z2; MDD(19,17) = 3*dx2yz;
	MDD(20,15) = 3*dxy2z; MDD(20,16) =-3*dx2yz; MDD(20,17) = 2*dx2y2;

	MDD(18,18) = 8*dy2z2;
	MDD(19,18) =-6*dxyz2; MDD(19,19) = 8*dx2z2;
	MDD(20,18) =-6*dxy2z; MDD(20,19) = 6*dx2yz; MDD(20,20) = 8*dx2y2;

	MDD(21, 0) =-2*dy2z2; MDD(21, 1) =-3*dxyz2; MDD(21, 2) =-3*dxy2z;
	MDD(22, 0) =-3*dxyz2; MDD(22, 1) =-2*dx2z2; MDD(22, 2) =-3*dx2yz;
	MDD(23, 0) =-3*dxy2z; MDD(23, 1) =-3*dx2yz; MDD(23, 2) =-2*dx2y2;

	MDD(21, 3) = 2*dy2z2; MDD(21, 4) =-3*dxyz2; MDD(21, 5) =-3*dxy2z;
	MDD(22, 3) = 3*dxyz2; MDD(22, 4) =-4*dx2z2; MDD(22, 5) =-6*dx2yz;
	MDD(23, 3) = 3*dxy2z; MDD(23, 4) =-6*dx2yz; MDD(23, 5) =-4*dx2y2;

	MDD(21, 6) =-4*dy2z2; MDD(21, 7) = 3*dxyz2; MDD(21, 8) =-6*dxy2z;
	MDD(22, 6) =-3*dxyz2; MDD(22, 7) = 2*dx2z2; MDD(22, 8) =-3*dx2yz;
	MDD(23, 6) =-6*dxy2z; MDD(23, 7) = 3*dx2yz; MDD(23, 8) =-4*dx2y2;
 
	MDD(21, 9) = 4*dy2z2; MDD(21,10) = 3*dxyz2; MDD(21,11) =-6*dxy2z;
	MDD(22, 9) = 3*dxyz2; MDD(22,10) = 4*dx2z2; MDD(22,11) =-6*dx2yz;
	MDD(23, 9) = 6*dxy2z; MDD(23,10) = 6*dx2yz; MDD(23,11) =-8*dx2y2;

	MDD(21,12) =-4*dy2z2; MDD(21,13) =-6*dxyz2; MDD(21,14) = 3*dxy2z;
	MDD(22,12) =-6*dxyz2; MDD(22,13) =-4*dx2z2; MDD(22,14) = 3*dx2yz;
	MDD(23,12) =-3*dxy2z; MDD(23,13) =-3*dx2yz; MDD(23,14) = 2*dx2y2;

	MDD(21,15) = 4*dy2z2; MDD(21,16) =-6*dxyz2; MDD(21,17) = 3*dxy2z;
	MDD(22,15) = 6*dxyz2; MDD(22,16) =-8*dx2z2; MDD(22,17) = 6*dx2yz;
	MDD(23,15) = 3*dxy2z; MDD(23,16) =-6*dx2yz; MDD(23,17) = 4*dx2y2;

	MDD(21,18) =-8*dy2z2; MDD(21,19) = 6*dxyz2; MDD(21,20) = 6*dxy2z;
	MDD(22,18) =-6*dxyz2; MDD(22,19) = 4*dx2z2; MDD(22,20) = 3*dx2yz;
	MDD(23,18) =-6*dxy2z; MDD(23,19) = 3*dx2yz; MDD(23,20) = 4*dx2y2;

	MDD(21,21) = 8*dy2z2;
	MDD(22,21) = 6*dxyz2; MDD(22,22) = 8*dx2z2;
	MDD(23,21) = 6*dxy2z; MDD(23,22) = 6*dx2yz; MDD(23,23) = 8*dx2y2;

	MDD /= 72.0*dx*dy*dz;
    }
    return *MixedDD;
}

double Voxel8::IntFdd (int i, int j, int k, int l, int m) const
{
    static RSymMatrix *fidd = 0;
    if (!fidd) {
	int z;
	fidd = new RSymMatrix[8];
	for (z = 0; z < 8; z++)
	    fidd[z].New(24,24);
	RSymMatrix &m0 = fidd[0];
	RSymMatrix &m1 = fidd[1];
	RSymMatrix &m2 = fidd[2];
	RSymMatrix &m3 = fidd[3];
	RSymMatrix &m4 = fidd[4];
	RSymMatrix &m5 = fidd[5];
	RSymMatrix &m6 = fidd[6];
	RSymMatrix &m7 = fidd[7];

	double dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz;
	double dx2y2 = dx2*dy2, dx2z2 = dx2*dz2, dy2z2 = dy2*dz2;
	double dx2yz = dx2*dy*dz, dxy2z = dx*dy2*dz, dxyz2 = dx*dy*dz2;

	m0( 0, 0) = 27*dy2z2; 
	m0( 1, 0) = 24*dxyz2;  m0( 1, 1) = 27*dx2z2; 
	m0( 2, 0) = 24*dxy2z;  m0( 2, 1) = 24*dx2yz;  m0( 2, 2) = 27*dx2y2; 
	m0( 3, 0) = -27*dy2z2; m0( 3, 1) = -24*dxyz2; m0( 3, 2) = -24*dxy2z; 
	m0( 4, 0) = 12*dxyz2;  m0( 4, 1) = 9*dx2z2;   m0( 4, 2) = 8*dx2yz; 
	m0( 5, 0) = 12*dxy2z;  m0( 5, 1) = 8*dx2yz;   m0( 5, 2) = 9*dx2y2; 
	m0( 3, 3) = 27*dy2z2; 
	m0( 4, 3) = -12*dxyz2; m0( 4, 4) = 9*dx2z2; 
	m0( 5, 3) = -12*dxy2z; m0( 5, 4) = 8*dx2yz;   m0( 5, 5) = 9*dx2y2; 
	m0( 6, 0) = 9*dy2z2;   m0( 6, 1) = 12*dxyz2;  m0( 6, 2) = 8*dxy2z; 
	m0( 7, 0) = -24*dxyz2; m0( 7, 1) = -27*dx2z2; m0( 7, 2) = -24*dx2yz; 
	m0( 8, 0) = 8*dxy2z;   m0( 8, 1) = 12*dx2yz;  m0( 8, 2) = 9*dx2y2; 
	m0( 6, 3) = -9*dy2z2;  m0( 6, 4) = 6*dxyz2;   m0( 6, 5) = 4*dxy2z; 
	m0( 7, 3) = 24*dxyz2;  m0( 7, 4) = -9*dx2z2;  m0( 7, 5) = -8*dx2yz; 
	m0( 8, 3) = -8*dxy2z;  m0( 8, 4) = 4*dx2yz;   m0( 8, 5) = 3*dx2y2; 
	m0( 6, 6) = 9*dy2z2; 
	m0( 7, 6) = -12*dxyz2; m0( 7, 7) = 27*dx2z2; 
	m0( 8, 6) = 8*dxy2z; m0( 8, 7) = -12*dx2yz; m0( 8, 8) = 9*dx2y2; 
	m0( 9, 0) = -9*dy2z2; m0( 9, 1) = -12*dxyz2; m0( 9, 2) = -8*dxy2z; 
	m0(10, 0) = -12*dxyz2; m0(10, 1) = -9*dx2z2; m0(10, 2) = -8*dx2yz; 
	m0(11, 0) = 4*dxy2z; m0(11, 1) = 4*dx2yz; m0(11, 2) = 3*dx2y2; 
	m0( 9, 3) = 9*dy2z2; m0( 9, 4) = -6*dxyz2; m0( 9, 5) = -4*dxy2z; 
	m0(10, 3) = 12*dxyz2; m0(10, 4) = -9*dx2z2; m0(10, 5) = -8*dx2yz; 
	m0(11, 3) = -4*dxy2z; m0(11, 4) = 4*dx2yz; m0(11, 5) = 3*dx2y2; 
	m0( 9, 6) = -9*dy2z2; m0( 9, 7) = 12*dxyz2; m0( 9, 8) = -8*dxy2z; 
	m0(10, 6) = -6*dxyz2; m0(10, 7) = 9*dx2z2; m0(10, 8) = -4*dx2yz; 
	m0(11, 6) = 4*dxy2z; m0(11, 7) = -4*dx2yz; m0(11, 8) = 3*dx2y2; 
	m0( 9, 9) = 9*dy2z2; 
	m0(10, 9) = 6*dxyz2; m0(10,10) = 9*dx2z2; 
	m0(11, 9) = -4*dxy2z; m0(11,10) = -4*dx2yz; m0(11,11) = 3*dx2y2; 
	m0(12, 0) = 9*dy2z2; m0(12, 1) = 8*dxyz2; m0(12, 2) = 12*dxy2z; 
	m0(13, 0) = 8*dxyz2; m0(13, 1) = 9*dx2z2; m0(13, 2) = 12*dx2yz; 
	m0(14, 0) = -24*dxy2z; m0(14, 1) = -24*dx2yz; m0(14, 2) = -27*dx2y2; 
	m0(12, 3) = -9*dy2z2; m0(12, 4) = 4*dxyz2; m0(12, 5) = 6*dxy2z; 
	m0(13, 3) = -8*dxyz2; m0(13, 4) = 3*dx2z2; m0(13, 5) = 4*dx2yz; 
	m0(14, 3) = 24*dxy2z; m0(14, 4) = -8*dx2yz; m0(14, 5) = -9*dx2y2; 
	m0(12, 6) = 3*dy2z2; m0(12, 7) = -8*dxyz2; m0(12, 8) = 4*dxy2z; 
	m0(13, 6) = 4*dxyz2; m0(13, 7) = -9*dx2z2; m0(13, 8) = 6*dx2yz; 
	m0(14, 6) = -8*dxy2z; m0(14, 7) = 24*dx2yz; m0(14, 8) = -9*dx2y2; 
	m0(12, 9) = -3*dy2z2; m0(12,10) = -4*dxyz2; m0(12,11) = 2*dxy2z; 
	m0(13, 9) = -4*dxyz2; m0(13,10) = -3*dx2z2; m0(13,11) = 2*dx2yz; 
	m0(14, 9) = 8*dxy2z; m0(14,10) = 8*dx2yz; m0(14,11) = -3*dx2y2; 
	m0(12,12) = 9*dy2z2; 
	m0(13,12) = 8*dxyz2; m0(13,13) = 9*dx2z2; 
	m0(14,12) = -12*dxy2z; m0(14,13) = -12*dx2yz; m0(14,14) = 27*dx2y2; 
	m0(15, 0) = -9*dy2z2; m0(15, 1) = -8*dxyz2; m0(15, 2) = -12*dxy2z; 
	m0(16, 0) = 4*dxyz2; m0(16, 1) = 3*dx2z2; m0(16, 2) = 4*dx2yz; 
	m0(17, 0) = -12*dxy2z; m0(17, 1) = -8*dx2yz; m0(17, 2) = -9*dx2y2; 
	m0(15, 3) = 9*dy2z2; m0(15, 4) = -4*dxyz2; m0(15, 5) = -6*dxy2z; 
	m0(16, 3) = -4*dxyz2; m0(16, 4) = 3*dx2z2; m0(16, 5) = 4*dx2yz; 
	m0(17, 3) = 12*dxy2z; m0(17, 4) = -8*dx2yz; m0(17, 5) = -9*dx2y2; 
	m0(15, 6) = -3*dy2z2; m0(15, 7) = 8*dxyz2; m0(15, 8) = -4*dxy2z; 
	m0(16, 6) = 2*dxyz2; m0(16, 7) = -3*dx2z2; m0(16, 8) = 2*dx2yz; 
	m0(17, 6) = -4*dxy2z; m0(17, 7) = 8*dx2yz; m0(17, 8) = -3*dx2y2; 
	m0(15, 9) = 3*dy2z2; m0(15,10) = 4*dxyz2; m0(15,11) = -2*dxy2z; 
	m0(16, 9) = -2*dxyz2; m0(16,10) = -3*dx2z2; m0(16,11) = 2*dx2yz; 
	m0(17, 9) = 4*dxy2z; m0(17,10) = 8*dx2yz; m0(17,11) = -3*dx2y2; 
	m0(15,12) = -9*dy2z2; m0(15,13) = -8*dxyz2; m0(15,14) = 12*dxy2z; 
	m0(16,12) = 4*dxyz2; m0(16,13) = 3*dx2z2; m0(16,14) = -4*dx2yz; 
	m0(17,12) = -6*dxy2z; m0(17,13) = -4*dx2yz; m0(17,14) = 9*dx2y2; 
	m0(15,15) = 9*dy2z2; 
	m0(16,15) = -4*dxyz2; m0(16,16) = 3*dx2z2; 
	m0(17,15) = 6*dxy2z; m0(17,16) = -4*dx2yz; m0(17,17) = 9*dx2y2; 
	m0(18, 0) = 3*dy2z2; m0(18, 1) = 4*dxyz2; m0(18, 2) = 4*dxy2z; 
	m0(19, 0) = -8*dxyz2; m0(19, 1) = -9*dx2z2; m0(19, 2) = -12*dx2yz; 
	m0(20, 0) = -8*dxy2z; m0(20, 1) = -12*dx2yz; m0(20, 2) = -9*dx2y2; 
	m0(18, 3) = -3*dy2z2; m0(18, 4) = 2*dxyz2; m0(18, 5) = 2*dxy2z; 
	m0(19, 3) = 8*dxyz2; m0(19, 4) = -3*dx2z2; m0(19, 5) = -4*dx2yz; 
	m0(20, 3) = 8*dxy2z; m0(20, 4) = -4*dx2yz; m0(20, 5) = -3*dx2y2; 
	m0(18, 6) = 3*dy2z2; m0(18, 7) = -4*dxyz2; m0(18, 8) = 4*dxy2z; 
	m0(19, 6) = -4*dxyz2; m0(19, 7) = 9*dx2z2; m0(19, 8) = -6*dx2yz; 
	m0(20, 6) = -8*dxy2z; m0(20, 7) = 12*dx2yz; m0(20, 8) = -9*dx2y2; 
	m0(18, 9) = -3*dy2z2; m0(18,10) = -2*dxyz2; m0(18,11) = 2*dxy2z; 
	m0(19, 9) = 4*dxyz2; m0(19,10) = 3*dx2z2; m0(19,11) = -2*dx2yz; 
	m0(20, 9) = 8*dxy2z; m0(20,10) = 4*dx2yz; m0(20,11) = -3*dx2y2; 
	m0(18,12) = 3*dy2z2; m0(18,13) = 4*dxyz2; m0(18,14) = -4*dxy2z; 
	m0(19,12) = -8*dxyz2; m0(19,13) = -9*dx2z2; m0(19,14) = 12*dx2yz; 
	m0(20,12) = -4*dxy2z; m0(20,13) = -6*dx2yz; m0(20,14) = 9*dx2y2; 
	m0(18,15) = -3*dy2z2; m0(18,16) = 2*dxyz2; m0(18,17) = -2*dxy2z; 
	m0(19,15) = 8*dxyz2; m0(19,16) = -3*dx2z2; m0(19,17) = 4*dx2yz; 
	m0(20,15) = 4*dxy2z; m0(20,16) = -2*dx2yz; m0(20,17) = 3*dx2y2; 
	m0(18,18) = 3*dy2z2; 
	m0(19,18) = -4*dxyz2; m0(19,19) = 9*dx2z2; 
	m0(20,18) = -4*dxy2z; m0(20,19) = 6*dx2yz; m0(20,20) = 9*dx2y2; 
	m0(21, 0) = -3*dy2z2; m0(21, 1) = -4*dxyz2; m0(21, 2) = -4*dxy2z; 
	m0(22, 0) = -4*dxyz2; m0(22, 1) = -3*dx2z2; m0(22, 2) = -4*dx2yz; 
	m0(23, 0) = -4*dxy2z; m0(23, 1) = -4*dx2yz; m0(23, 2) = -3*dx2y2; 
	m0(21, 3) = 3*dy2z2; m0(21, 4) = -2*dxyz2; m0(21, 5) = -2*dxy2z; 
	m0(22, 3) = 4*dxyz2; m0(22, 4) = -3*dx2z2; m0(22, 5) = -4*dx2yz; 
	m0(23, 3) = 4*dxy2z; m0(23, 4) = -4*dx2yz; m0(23, 5) = -3*dx2y2; 
	m0(21, 6) = -3*dy2z2; m0(21, 7) = 4*dxyz2; m0(21, 8) = -4*dxy2z; 
	m0(22, 6) = -2*dxyz2; m0(22, 7) = 3*dx2z2; m0(22, 8) = -2*dx2yz; 
	m0(23, 6) = -4*dxy2z; m0(23, 7) = 4*dx2yz; m0(23, 8) = -3*dx2y2; 
	m0(21, 9) = 3*dy2z2; m0(21,10) = 2*dxyz2; m0(21,11) = -2*dxy2z; 
	m0(22, 9) = 2*dxyz2; m0(22,10) = 3*dx2z2; m0(22,11) = -2*dx2yz; 
	m0(23, 9) = 4*dxy2z; m0(23,10) = 4*dx2yz; m0(23,11) = -3*dx2y2; 
	m0(21,12) = -3*dy2z2; m0(21,13) = -4*dxyz2; m0(21,14) = 4*dxy2z; 
	m0(22,12) = -4*dxyz2; m0(22,13) = -3*dx2z2; m0(22,14) = 4*dx2yz; 
	m0(23,12) = -2*dxy2z; m0(23,13) = -2*dx2yz; m0(23,14) = 3*dx2y2; 
	m0(21,15) = 3*dy2z2; m0(21,16) = -2*dxyz2; m0(21,17) = 2*dxy2z; 
	m0(22,15) = 4*dxyz2; m0(22,16) = -3*dx2z2; m0(22,17) = 4*dx2yz; 
	m0(23,15) = 2*dxy2z; m0(23,16) = -2*dx2yz; m0(23,17) = 3*dx2y2; 
	m0(21,18) = -3*dy2z2; m0(21,19) = 4*dxyz2; m0(21,20) = 4*dxy2z; 
	m0(22,18) = -2*dxyz2; m0(22,19) = 3*dx2z2; m0(22,20) = 2*dx2yz; 
	m0(23,18) = -2*dxy2z; m0(23,19) = 2*dx2yz; m0(23,20) = 3*dx2y2; 
	m0(21,21) = 3*dy2z2; 
	m0(22,21) = 2*dxyz2; m0(22,22) = 3*dx2z2; 
	m0(23,21) = 2*dxy2z; m0(23,22) = 2*dx2yz; m0(23,23) = 3*dx2y2; 
	m1( 0, 0) = 27*dy2z2; 
	m1( 1, 0) = 12*dxyz2; m1( 1, 1) = 9*dx2z2; 
	m1( 2, 0) = 12*dxy2z; m1( 2, 1) = 8*dx2yz; m1( 2, 2) = 9*dx2y2; 
	m1( 3, 0) = -27*dy2z2; m1( 3, 1) = -12*dxyz2; m1( 3, 2) = -12*dxy2z; 
	m1( 4, 0) = 24*dxyz2; m1( 4, 1) = 9*dx2z2; m1( 4, 2) = 8*dx2yz; 
	m1( 5, 0) = 24*dxy2z; m1( 5, 1) = 8*dx2yz; m1( 5, 2) = 9*dx2y2; 
	m1( 3, 3) = 27*dy2z2; 
	m1( 4, 3) = -24*dxyz2; m1( 4, 4) = 27*dx2z2; 
	m1( 5, 3) = -24*dxy2z; m1( 5, 4) = 24*dx2yz; m1( 5, 5) = 27*dx2y2; 
	m1( 6, 0) = 9*dy2z2; m1( 6, 1) = 6*dxyz2; m1( 6, 2) = 4*dxy2z; 
	m1( 7, 0) = -12*dxyz2; m1( 7, 1) = -9*dx2z2; m1( 7, 2) = -8*dx2yz; 
	m1( 8, 0) = 4*dxy2z; m1( 8, 1) = 4*dx2yz; m1( 8, 2) = 3*dx2y2; 
	m1( 6, 3) = -9*dy2z2; m1( 6, 4) = 12*dxyz2; m1( 6, 5) = 8*dxy2z; 
	m1( 7, 3) = 12*dxyz2; m1( 7, 4) = -9*dx2z2; m1( 7, 5) = -8*dx2yz; 
	m1( 8, 3) = -4*dxy2z; m1( 8, 4) = 4*dx2yz; m1( 8, 5) = 3*dx2y2; 
	m1( 6, 6) = 9*dy2z2; 
	m1( 7, 6) = -6*dxyz2; m1( 7, 7) = 9*dx2z2; 
	m1( 8, 6) = 4*dxy2z; m1( 8, 7) = -4*dx2yz; m1( 8, 8) = 3*dx2y2; 
	m1( 9, 0) = -9*dy2z2; m1( 9, 1) = -6*dxyz2; m1( 9, 2) = -4*dxy2z; 
	m1(10, 0) = -24*dxyz2; m1(10, 1) = -9*dx2z2; m1(10, 2) = -8*dx2yz; 
	m1(11, 0) = 8*dxy2z; m1(11, 1) = 4*dx2yz; m1(11, 2) = 3*dx2y2; 
	m1( 9, 3) = 9*dy2z2; m1( 9, 4) = -12*dxyz2; m1( 9, 5) = -8*dxy2z; 
	m1(10, 3) = 24*dxyz2; m1(10, 4) = -27*dx2z2; m1(10, 5) = -24*dx2yz; 
	m1(11, 3) = -8*dxy2z; m1(11, 4) = 12*dx2yz; m1(11, 5) = 9*dx2y2; 
	m1( 9, 6) = -9*dy2z2; m1( 9, 7) = 6*dxyz2; m1( 9, 8) = -4*dxy2z; 
	m1(10, 6) = -12*dxyz2; m1(10, 7) = 9*dx2z2; m1(10, 8) = -4*dx2yz; 
	m1(11, 6) = 8*dxy2z; m1(11, 7) = -4*dx2yz; m1(11, 8) = 3*dx2y2; 
	m1( 9, 9) = 9*dy2z2; 
	m1(10, 9) = 12*dxyz2; m1(10,10) = 27*dx2z2; 
	m1(11, 9) = -8*dxy2z; m1(11,10) = -12*dx2yz; m1(11,11) = 9*dx2y2; 
	m1(12, 0) = 9*dy2z2; m1(12, 1) = 4*dxyz2; m1(12, 2) = 6*dxy2z; 
	m1(13, 0) = 4*dxyz2; m1(13, 1) = 3*dx2z2; m1(13, 2) = 4*dx2yz; 
	m1(14, 0) = -12*dxy2z; m1(14, 1) = -8*dx2yz; m1(14, 2) = -9*dx2y2; 
	m1(12, 3) = -9*dy2z2; m1(12, 4) = 8*dxyz2; m1(12, 5) = 12*dxy2z; 
	m1(13, 3) = -4*dxyz2; m1(13, 4) = 3*dx2z2; m1(13, 5) = 4*dx2yz; 
	m1(14, 3) = 12*dxy2z; m1(14, 4) = -8*dx2yz; m1(14, 5) = -9*dx2y2; 
	m1(12, 6) = 3*dy2z2; m1(12, 7) = -4*dxyz2; m1(12, 8) = 2*dxy2z; 
	m1(13, 6) = 2*dxyz2; m1(13, 7) = -3*dx2z2; m1(13, 8) = 2*dx2yz; 
	m1(14, 6) = -4*dxy2z; m1(14, 7) = 8*dx2yz; m1(14, 8) = -3*dx2y2; 
	m1(12, 9) = -3*dy2z2; m1(12,10) = -8*dxyz2; m1(12,11) = 4*dxy2z; 
	m1(13, 9) = -2*dxyz2; m1(13,10) = -3*dx2z2; m1(13,11) = 2*dx2yz; 
	m1(14, 9) = 4*dxy2z; m1(14,10) = 8*dx2yz; m1(14,11) = -3*dx2y2; 
	m1(12,12) = 9*dy2z2; 
	m1(13,12) = 4*dxyz2; m1(13,13) = 3*dx2z2; 
	m1(14,12) = -6*dxy2z; m1(14,13) = -4*dx2yz; m1(14,14) = 9*dx2y2; 
	m1(15, 0) = -9*dy2z2; m1(15, 1) = -4*dxyz2; m1(15, 2) = -6*dxy2z; 
	m1(16, 0) = 8*dxyz2; m1(16, 1) = 3*dx2z2; m1(16, 2) = 4*dx2yz; 
	m1(17, 0) = -24*dxy2z; m1(17, 1) = -8*dx2yz; m1(17, 2) = -9*dx2y2; 
	m1(15, 3) = 9*dy2z2; m1(15, 4) = -8*dxyz2; m1(15, 5) = -12*dxy2z; 
	m1(16, 3) = -8*dxyz2; m1(16, 4) = 9*dx2z2; m1(16, 5) = 12*dx2yz; 
	m1(17, 3) = 24*dxy2z; m1(17, 4) = -24*dx2yz; m1(17, 5) = -27*dx2y2; 
	m1(15, 6) = -3*dy2z2; m1(15, 7) = 4*dxyz2; m1(15, 8) = -2*dxy2z; 
	m1(16, 6) = 4*dxyz2; m1(16, 7) = -3*dx2z2; m1(16, 8) = 2*dx2yz; 
	m1(17, 6) = -8*dxy2z; m1(17, 7) = 8*dx2yz; m1(17, 8) = -3*dx2y2; 
	m1(15, 9) = 3*dy2z2; m1(15,10) = 8*dxyz2; m1(15,11) = -4*dxy2z; 
	m1(16, 9) = -4*dxyz2; m1(16,10) = -9*dx2z2; m1(16,11) = 6*dx2yz; 
	m1(17, 9) = 8*dxy2z; m1(17,10) = 24*dx2yz; m1(17,11) = -9*dx2y2; 
	m1(15,12) = -9*dy2z2; m1(15,13) = -4*dxyz2; m1(15,14) = 6*dxy2z; 
	m1(16,12) = 8*dxyz2; m1(16,13) = 3*dx2z2; m1(16,14) = -4*dx2yz; 
	m1(17,12) = -12*dxy2z; m1(17,13) = -4*dx2yz; m1(17,14) = 9*dx2y2; 
	m1(15,15) = 9*dy2z2; 
	m1(16,15) = -8*dxyz2; m1(16,16) = 9*dx2z2; 
	m1(17,15) = 12*dxy2z; m1(17,16) = -12*dx2yz; m1(17,17) = 27*dx2y2; 
	m1(18, 0) = 3*dy2z2; m1(18, 1) = 2*dxyz2; m1(18, 2) = 2*dxy2z; 
	m1(19, 0) = -4*dxyz2; m1(19, 1) = -3*dx2z2; m1(19, 2) = -4*dx2yz; 
	m1(20, 0) = -4*dxy2z; m1(20, 1) = -4*dx2yz; m1(20, 2) = -3*dx2y2; 
	m1(18, 3) = -3*dy2z2; m1(18, 4) = 4*dxyz2; m1(18, 5) = 4*dxy2z; 
	m1(19, 3) = 4*dxyz2; m1(19, 4) = -3*dx2z2; m1(19, 5) = -4*dx2yz; 
	m1(20, 3) = 4*dxy2z; m1(20, 4) = -4*dx2yz; m1(20, 5) = -3*dx2y2; 
	m1(18, 6) = 3*dy2z2; m1(18, 7) = -2*dxyz2; m1(18, 8) = 2*dxy2z; 
	m1(19, 6) = -2*dxyz2; m1(19, 7) = 3*dx2z2; m1(19, 8) = -2*dx2yz; 
	m1(20, 6) = -4*dxy2z; m1(20, 7) = 4*dx2yz; m1(20, 8) = -3*dx2y2; 
	m1(18, 9) = -3*dy2z2; m1(18,10) = -4*dxyz2; m1(18,11) = 4*dxy2z; 
	m1(19, 9) = 2*dxyz2; m1(19,10) = 3*dx2z2; m1(19,11) = -2*dx2yz; 
	m1(20, 9) = 4*dxy2z; m1(20,10) = 4*dx2yz; m1(20,11) = -3*dx2y2; 
	m1(18,12) = 3*dy2z2; m1(18,13) = 2*dxyz2; m1(18,14) = -2*dxy2z; 
	m1(19,12) = -4*dxyz2; m1(19,13) = -3*dx2z2; m1(19,14) = 4*dx2yz; 
	m1(20,12) = -2*dxy2z; m1(20,13) = -2*dx2yz; m1(20,14) = 3*dx2y2; 
	m1(18,15) = -3*dy2z2; m1(18,16) = 4*dxyz2; m1(18,17) = -4*dxy2z; 
	m1(19,15) = 4*dxyz2; m1(19,16) = -3*dx2z2; m1(19,17) = 4*dx2yz; 
	m1(20,15) = 2*dxy2z; m1(20,16) = -2*dx2yz; m1(20,17) = 3*dx2y2; 
	m1(18,18) = 3*dy2z2; 
	m1(19,18) = -2*dxyz2; m1(19,19) = 3*dx2z2; 
	m1(20,18) = -2*dxy2z; m1(20,19) = 2*dx2yz; m1(20,20) = 3*dx2y2; 
	m1(21, 0) = -3*dy2z2; m1(21, 1) = -2*dxyz2; m1(21, 2) = -2*dxy2z; 
	m1(22, 0) = -8*dxyz2; m1(22, 1) = -3*dx2z2; m1(22, 2) = -4*dx2yz; 
	m1(23, 0) = -8*dxy2z; m1(23, 1) = -4*dx2yz; m1(23, 2) = -3*dx2y2; 
	m1(21, 3) = 3*dy2z2; m1(21, 4) = -4*dxyz2; m1(21, 5) = -4*dxy2z; 
	m1(22, 3) = 8*dxyz2; m1(22, 4) = -9*dx2z2; m1(22, 5) = -12*dx2yz; 
	m1(23, 3) = 8*dxy2z; m1(23, 4) = -12*dx2yz; m1(23, 5) = -9*dx2y2; 
	m1(21, 6) = -3*dy2z2; m1(21, 7) = 2*dxyz2; m1(21, 8) = -2*dxy2z; 
	m1(22, 6) = -4*dxyz2; m1(22, 7) = 3*dx2z2; m1(22, 8) = -2*dx2yz; 
	m1(23, 6) = -8*dxy2z; m1(23, 7) = 4*dx2yz; m1(23, 8) = -3*dx2y2; 
	m1(21, 9) = 3*dy2z2; m1(21,10) = 4*dxyz2; m1(21,11) = -4*dxy2z; 
	m1(22, 9) = 4*dxyz2; m1(22,10) = 9*dx2z2; m1(22,11) = -6*dx2yz; 
	m1(23, 9) = 8*dxy2z; m1(23,10) = 12*dx2yz; m1(23,11) = -9*dx2y2; 
	m1(21,12) = -3*dy2z2; m1(21,13) = -2*dxyz2; m1(21,14) = 2*dxy2z; 
	m1(22,12) = -8*dxyz2; m1(22,13) = -3*dx2z2; m1(22,14) = 4*dx2yz; 
	m1(23,12) = -4*dxy2z; m1(23,13) = -2*dx2yz; m1(23,14) = 3*dx2y2; 
	m1(21,15) = 3*dy2z2; m1(21,16) = -4*dxyz2; m1(21,17) = 4*dxy2z; 
	m1(22,15) = 8*dxyz2; m1(22,16) = -9*dx2z2; m1(22,17) = 12*dx2yz; 
	m1(23,15) = 4*dxy2z; m1(23,16) = -6*dx2yz; m1(23,17) = 9*dx2y2; 
	m1(21,18) = -3*dy2z2; m1(21,19) = 2*dxyz2; m1(21,20) = 2*dxy2z; 
	m1(22,18) = -4*dxyz2; m1(22,19) = 3*dx2z2; m1(22,20) = 2*dx2yz; 
	m1(23,18) = -4*dxy2z; m1(23,19) = 2*dx2yz; m1(23,20) = 3*dx2y2; 
	m1(21,21) = 3*dy2z2; 
	m1(22,21) = 4*dxyz2; m1(22,22) = 9*dx2z2; 
	m1(23,21) = 4*dxy2z; m1(23,22) = 6*dx2yz; m1(23,23) = 9*dx2y2; 
	m2( 0, 0) = 9*dy2z2; 
	m2( 1, 0) = 12*dxyz2; m2( 1, 1) = 27*dx2z2; 
	m2( 2, 0) = 8*dxy2z; m2( 2, 1) = 12*dx2yz; m2( 2, 2) = 9*dx2y2; 
	m2( 3, 0) = -9*dy2z2; m2( 3, 1) = -12*dxyz2; m2( 3, 2) = -8*dxy2z; 
	m2( 4, 0) = 6*dxyz2; m2( 4, 1) = 9*dx2z2; m2( 4, 2) = 4*dx2yz; 
	m2( 5, 0) = 4*dxy2z; m2( 5, 1) = 4*dx2yz; m2( 5, 2) = 3*dx2y2; 
	m2( 3, 3) = 9*dy2z2; 
	m2( 4, 3) = -6*dxyz2; m2( 4, 4) = 9*dx2z2; 
	m2( 5, 3) = -4*dxy2z; m2( 5, 4) = 4*dx2yz; m2( 5, 5) = 3*dx2y2; 
	m2( 6, 0) = 9*dy2z2; m2( 6, 1) = 24*dxyz2; m2( 6, 2) = 8*dxy2z; 
	m2( 7, 0) = -12*dxyz2; m2( 7, 1) = -27*dx2z2; m2( 7, 2) = -12*dx2yz; 
	m2( 8, 0) = 8*dxy2z; m2( 8, 1) = 24*dx2yz; m2( 8, 2) = 9*dx2y2; 
	m2( 6, 3) = -9*dy2z2; m2( 6, 4) = 12*dxyz2; m2( 6, 5) = 4*dxy2z; 
	m2( 7, 3) = 12*dxyz2; m2( 7, 4) = -9*dx2z2; m2( 7, 5) = -4*dx2yz; 
	m2( 8, 3) = -8*dxy2z; m2( 8, 4) = 8*dx2yz; m2( 8, 5) = 3*dx2y2; 
	m2( 6, 6) = 27*dy2z2; 
	m2( 7, 6) = -24*dxyz2; m2( 7, 7) = 27*dx2z2; 
	m2( 8, 6) = 24*dxy2z; m2( 8, 7) = -24*dx2yz; m2( 8, 8) = 27*dx2y2; 
	m2( 9, 0) = -9*dy2z2; m2( 9, 1) = -24*dxyz2; m2( 9, 2) = -8*dxy2z; 
	m2(10, 0) = -6*dxyz2; m2(10, 1) = -9*dx2z2; m2(10, 2) = -4*dx2yz; 
	m2(11, 0) = 4*dxy2z; m2(11, 1) = 8*dx2yz; m2(11, 2) = 3*dx2y2; 
	m2( 9, 3) = 9*dy2z2; m2( 9, 4) = -12*dxyz2; m2( 9, 5) = -4*dxy2z; 
	m2(10, 3) = 6*dxyz2; m2(10, 4) = -9*dx2z2; m2(10, 5) = -4*dx2yz; 
	m2(11, 3) = -4*dxy2z; m2(11, 4) = 8*dx2yz; m2(11, 5) = 3*dx2y2; 
	m2( 9, 6) = -27*dy2z2; m2( 9, 7) = 24*dxyz2; m2( 9, 8) = -24*dxy2z; 
	m2(10, 6) = -12*dxyz2; m2(10, 7) = 9*dx2z2; m2(10, 8) = -8*dx2yz; 
	m2(11, 6) = 12*dxy2z; m2(11, 7) = -8*dx2yz; m2(11, 8) = 9*dx2y2; 
	m2( 9, 9) = 27*dy2z2; 
	m2(10, 9) = 12*dxyz2; m2(10,10) = 9*dx2z2; 
	m2(11, 9) = -12*dxy2z; m2(11,10) = -8*dx2yz; m2(11,11) = 9*dx2y2; 
	m2(12, 0) = 3*dy2z2; m2(12, 1) = 4*dxyz2; m2(12, 2) = 4*dxy2z; 
	m2(13, 0) = 4*dxyz2; m2(13, 1) = 9*dx2z2; m2(13, 2) = 6*dx2yz; 
	m2(14, 0) = -8*dxy2z; m2(14, 1) = -12*dx2yz; m2(14, 2) = -9*dx2y2; 
	m2(12, 3) = -3*dy2z2; m2(12, 4) = 2*dxyz2; m2(12, 5) = 2*dxy2z; 
	m2(13, 3) = -4*dxyz2; m2(13, 4) = 3*dx2z2; m2(13, 5) = 2*dx2yz; 
	m2(14, 3) = 8*dxy2z; m2(14, 4) = -4*dx2yz; m2(14, 5) = -3*dx2y2; 
	m2(12, 6) = 3*dy2z2; m2(12, 7) = -4*dxyz2; m2(12, 8) = 4*dxy2z; 
	m2(13, 6) = 8*dxyz2; m2(13, 7) = -9*dx2z2; m2(13, 8) = 12*dx2yz; 
	m2(14, 6) = -8*dxy2z; m2(14, 7) = 12*dx2yz; m2(14, 8) = -9*dx2y2; 
	m2(12, 9) = -3*dy2z2; m2(12,10) = -2*dxyz2; m2(12,11) = 2*dxy2z; 
	m2(13, 9) = -8*dxyz2; m2(13,10) = -3*dx2z2; m2(13,11) = 4*dx2yz; 
	m2(14, 9) = 8*dxy2z; m2(14,10) = 4*dx2yz; m2(14,11) = -3*dx2y2; 
	m2(12,12) = 3*dy2z2; 
	m2(13,12) = 4*dxyz2; m2(13,13) = 9*dx2z2; 
	m2(14,12) = -4*dxy2z; m2(14,13) = -6*dx2yz; m2(14,14) = 9*dx2y2; 
	m2(15, 0) = -3*dy2z2; m2(15, 1) = -4*dxyz2; m2(15, 2) = -4*dxy2z; 
	m2(16, 0) = 2*dxyz2; m2(16, 1) = 3*dx2z2; m2(16, 2) = 2*dx2yz; 
	m2(17, 0) = -4*dxy2z; m2(17, 1) = -4*dx2yz; m2(17, 2) = -3*dx2y2; 
	m2(15, 3) = 3*dy2z2; m2(15, 4) = -2*dxyz2; m2(15, 5) = -2*dxy2z; 
	m2(16, 3) = -2*dxyz2; m2(16, 4) = 3*dx2z2; m2(16, 5) = 2*dx2yz; 
	m2(17, 3) = 4*dxy2z; m2(17, 4) = -4*dx2yz; m2(17, 5) = -3*dx2y2; 
	m2(15, 6) = -3*dy2z2; m2(15, 7) = 4*dxyz2; m2(15, 8) = -4*dxy2z; 
	m2(16, 6) = 4*dxyz2; m2(16, 7) = -3*dx2z2; m2(16, 8) = 4*dx2yz; 
	m2(17, 6) = -4*dxy2z; m2(17, 7) = 4*dx2yz; m2(17, 8) = -3*dx2y2; 
	m2(15, 9) = 3*dy2z2; m2(15,10) = 2*dxyz2; m2(15,11) = -2*dxy2z; 
	m2(16, 9) = -4*dxyz2; m2(16,10) = -3*dx2z2; m2(16,11) = 4*dx2yz; 
	m2(17, 9) = 4*dxy2z; m2(17,10) = 4*dx2yz; m2(17,11) = -3*dx2y2; 
	m2(15,12) = -3*dy2z2; m2(15,13) = -4*dxyz2; m2(15,14) = 4*dxy2z; 
	m2(16,12) = 2*dxyz2; m2(16,13) = 3*dx2z2; m2(16,14) = -2*dx2yz; 
	m2(17,12) = -2*dxy2z; m2(17,13) = -2*dx2yz; m2(17,14) = 3*dx2y2; 
	m2(15,15) = 3*dy2z2; 
	m2(16,15) = -2*dxyz2; m2(16,16) = 3*dx2z2; 
	m2(17,15) = 2*dxy2z; m2(17,16) = -2*dx2yz; m2(17,17) = 3*dx2y2; 
	m2(18, 0) = 3*dy2z2; m2(18, 1) = 8*dxyz2; m2(18, 2) = 4*dxy2z; 
	m2(19, 0) = -4*dxyz2; m2(19, 1) = -9*dx2z2; m2(19, 2) = -6*dx2yz; 
	m2(20, 0) = -8*dxy2z; m2(20, 1) = -24*dx2yz; m2(20, 2) = -9*dx2y2; 
	m2(18, 3) = -3*dy2z2; m2(18, 4) = 4*dxyz2; m2(18, 5) = 2*dxy2z; 
	m2(19, 3) = 4*dxyz2; m2(19, 4) = -3*dx2z2; m2(19, 5) = -2*dx2yz; 
	m2(20, 3) = 8*dxy2z; m2(20, 4) = -8*dx2yz; m2(20, 5) = -3*dx2y2; 
	m2(18, 6) = 9*dy2z2; m2(18, 7) = -8*dxyz2; m2(18, 8) = 12*dxy2z; 
	m2(19, 6) = -8*dxyz2; m2(19, 7) = 9*dx2z2; m2(19, 8) = -12*dx2yz; 
	m2(20, 6) = -24*dxy2z; m2(20, 7) = 24*dx2yz; m2(20, 8) = -27*dx2y2; 
	m2(18, 9) = -9*dy2z2; m2(18,10) = -4*dxyz2; m2(18,11) = 6*dxy2z; 
	m2(19, 9) = 8*dxyz2; m2(19,10) = 3*dx2z2; m2(19,11) = -4*dx2yz; 
	m2(20, 9) = 24*dxy2z; m2(20,10) = 8*dx2yz; m2(20,11) = -9*dx2y2; 
	m2(18,12) = 3*dy2z2; m2(18,13) = 8*dxyz2; m2(18,14) = -4*dxy2z; 
	m2(19,12) = -4*dxyz2; m2(19,13) = -9*dx2z2; m2(19,14) = 6*dx2yz; 
	m2(20,12) = -4*dxy2z; m2(20,13) = -12*dx2yz; m2(20,14) = 9*dx2y2; 
	m2(18,15) = -3*dy2z2; m2(18,16) = 4*dxyz2; m2(18,17) = -2*dxy2z; 
	m2(19,15) = 4*dxyz2; m2(19,16) = -3*dx2z2; m2(19,17) = 2*dx2yz; 
	m2(20,15) = 4*dxy2z; m2(20,16) = -4*dx2yz; m2(20,17) = 3*dx2y2; 
	m2(18,18) = 9*dy2z2; 
	m2(19,18) = -8*dxyz2; m2(19,19) = 9*dx2z2; 
	m2(20,18) = -12*dxy2z; m2(20,19) = 12*dx2yz; m2(20,20) = 27*dx2y2; 
	m2(21, 0) = -3*dy2z2; m2(21, 1) = -8*dxyz2; m2(21, 2) = -4*dxy2z; 
	m2(22, 0) = -2*dxyz2; m2(22, 1) = -3*dx2z2; m2(22, 2) = -2*dx2yz; 
	m2(23, 0) = -4*dxy2z; m2(23, 1) = -8*dx2yz; m2(23, 2) = -3*dx2y2; 
	m2(21, 3) = 3*dy2z2; m2(21, 4) = -4*dxyz2; m2(21, 5) = -2*dxy2z; 
	m2(22, 3) = 2*dxyz2; m2(22, 4) = -3*dx2z2; m2(22, 5) = -2*dx2yz; 
	m2(23, 3) = 4*dxy2z; m2(23, 4) = -8*dx2yz; m2(23, 5) = -3*dx2y2; 
	m2(21, 6) = -9*dy2z2; m2(21, 7) = 8*dxyz2; m2(21, 8) = -12*dxy2z; 
	m2(22, 6) = -4*dxyz2; m2(22, 7) = 3*dx2z2; m2(22, 8) = -4*dx2yz; 
	m2(23, 6) = -12*dxy2z; m2(23, 7) = 8*dx2yz; m2(23, 8) = -9*dx2y2; 
	m2(21, 9) = 9*dy2z2; m2(21,10) = 4*dxyz2; m2(21,11) = -6*dxy2z; 
	m2(22, 9) = 4*dxyz2; m2(22,10) = 3*dx2z2; m2(22,11) = -4*dx2yz; 
	m2(23, 9) = 12*dxy2z; m2(23,10) = 8*dx2yz; m2(23,11) = -9*dx2y2; 
	m2(21,12) = -3*dy2z2; m2(21,13) = -8*dxyz2; m2(21,14) = 4*dxy2z; 
	m2(22,12) = -2*dxyz2; m2(22,13) = -3*dx2z2; m2(22,14) = 2*dx2yz; 
	m2(23,12) = -2*dxy2z; m2(23,13) = -4*dx2yz; m2(23,14) = 3*dx2y2; 
	m2(21,15) = 3*dy2z2; m2(21,16) = -4*dxyz2; m2(21,17) = 2*dxy2z; 
	m2(22,15) = 2*dxyz2; m2(22,16) = -3*dx2z2; m2(22,17) = 2*dx2yz; 
	m2(23,15) = 2*dxy2z; m2(23,16) = -4*dx2yz; m2(23,17) = 3*dx2y2; 
	m2(21,18) = -9*dy2z2; m2(21,19) = 8*dxyz2; m2(21,20) = 12*dxy2z; 
	m2(22,18) = -4*dxyz2; m2(22,19) = 3*dx2z2; m2(22,20) = 4*dx2yz; 
	m2(23,18) = -6*dxy2z; m2(23,19) = 4*dx2yz; m2(23,20) = 9*dx2y2; 
	m2(21,21) = 9*dy2z2; 
	m2(22,21) = 4*dxyz2; m2(22,22) = 3*dx2z2; 
	m2(23,21) = 6*dxy2z; m2(23,22) = 4*dx2yz; m2(23,23) = 9*dx2y2; 
	m3( 0, 0) = 9*dy2z2; 
	m3( 1, 0) = 6*dxyz2; m3( 1, 1) = 9*dx2z2; 
	m3( 2, 0) = 4*dxy2z; m3( 2, 1) = 4*dx2yz; m3( 2, 2) = 3*dx2y2; 
	m3( 3, 0) = -9*dy2z2; m3( 3, 1) = -6*dxyz2; m3( 3, 2) = -4*dxy2z; 
	m3( 4, 0) = 12*dxyz2; m3( 4, 1) = 9*dx2z2; m3( 4, 2) = 4*dx2yz; 
	m3( 5, 0) = 8*dxy2z; m3( 5, 1) = 4*dx2yz; m3( 5, 2) = 3*dx2y2; 
	m3( 3, 3) = 9*dy2z2; 
	m3( 4, 3) = -12*dxyz2; m3( 4, 4) = 27*dx2z2; 
	m3( 5, 3) = -8*dxy2z; m3( 5, 4) = 12*dx2yz; m3( 5, 5) = 9*dx2y2; 
	m3( 6, 0) = 9*dy2z2; m3( 6, 1) = 12*dxyz2; m3( 6, 2) = 4*dxy2z; 
	m3( 7, 0) = -6*dxyz2; m3( 7, 1) = -9*dx2z2; m3( 7, 2) = -4*dx2yz; 
	m3( 8, 0) = 4*dxy2z; m3( 8, 1) = 8*dx2yz; m3( 8, 2) = 3*dx2y2; 
	m3( 6, 3) = -9*dy2z2; m3( 6, 4) = 24*dxyz2; m3( 6, 5) = 8*dxy2z; 
	m3( 7, 3) = 6*dxyz2; m3( 7, 4) = -9*dx2z2; m3( 7, 5) = -4*dx2yz; 
	m3( 8, 3) = -4*dxy2z; m3( 8, 4) = 8*dx2yz; m3( 8, 5) = 3*dx2y2; 
	m3( 6, 6) = 27*dy2z2; 
	m3( 7, 6) = -12*dxyz2; m3( 7, 7) = 9*dx2z2; 
	m3( 8, 6) = 12*dxy2z; m3( 8, 7) = -8*dx2yz; m3( 8, 8) = 9*dx2y2; 
	m3( 9, 0) = -9*dy2z2; m3( 9, 1) = -12*dxyz2; m3( 9, 2) = -4*dxy2z; 
	m3(10, 0) = -12*dxyz2; m3(10, 1) = -9*dx2z2; m3(10, 2) = -4*dx2yz; 
	m3(11, 0) = 8*dxy2z; m3(11, 1) = 8*dx2yz; m3(11, 2) = 3*dx2y2; 
	m3( 9, 3) = 9*dy2z2; m3( 9, 4) = -24*dxyz2; m3( 9, 5) = -8*dxy2z; 
	m3(10, 3) = 12*dxyz2; m3(10, 4) = -27*dx2z2; m3(10, 5) = -12*dx2yz; 
	m3(11, 3) = -8*dxy2z; m3(11, 4) = 24*dx2yz; m3(11, 5) = 9*dx2y2; 
	m3( 9, 6) = -27*dy2z2; m3( 9, 7) = 12*dxyz2; m3( 9, 8) = -12*dxy2z; 
	m3(10, 6) = -24*dxyz2; m3(10, 7) = 9*dx2z2; m3(10, 8) = -8*dx2yz; 
	m3(11, 6) = 24*dxy2z; m3(11, 7) = -8*dx2yz; m3(11, 8) = 9*dx2y2; 
	m3( 9, 9) = 27*dy2z2; 
	m3(10, 9) = 24*dxyz2; m3(10,10) = 27*dx2z2; 
	m3(11, 9) = -24*dxy2z; m3(11,10) = -24*dx2yz; m3(11,11) = 27*dx2y2; 
	m3(12, 0) = 3*dy2z2; m3(12, 1) = 2*dxyz2; m3(12, 2) = 2*dxy2z; 
	m3(13, 0) = 2*dxyz2; m3(13, 1) = 3*dx2z2; m3(13, 2) = 2*dx2yz; 
	m3(14, 0) = -4*dxy2z; m3(14, 1) = -4*dx2yz; m3(14, 2) = -3*dx2y2; 
	m3(12, 3) = -3*dy2z2; m3(12, 4) = 4*dxyz2; m3(12, 5) = 4*dxy2z; 
	m3(13, 3) = -2*dxyz2; m3(13, 4) = 3*dx2z2; m3(13, 5) = 2*dx2yz; 
	m3(14, 3) = 4*dxy2z; m3(14, 4) = -4*dx2yz; m3(14, 5) = -3*dx2y2; 
	m3(12, 6) = 3*dy2z2; m3(12, 7) = -2*dxyz2; m3(12, 8) = 2*dxy2z; 
	m3(13, 6) = 4*dxyz2; m3(13, 7) = -3*dx2z2; m3(13, 8) = 4*dx2yz; 
	m3(14, 6) = -4*dxy2z; m3(14, 7) = 4*dx2yz; m3(14, 8) = -3*dx2y2; 
	m3(12, 9) = -3*dy2z2; m3(12,10) = -4*dxyz2; m3(12,11) = 4*dxy2z; 
	m3(13, 9) = -4*dxyz2; m3(13,10) = -3*dx2z2; m3(13,11) = 4*dx2yz; 
	m3(14, 9) = 4*dxy2z; m3(14,10) = 4*dx2yz; m3(14,11) = -3*dx2y2; 
	m3(12,12) = 3*dy2z2; 
	m3(13,12) = 2*dxyz2; m3(13,13) = 3*dx2z2; 
	m3(14,12) = -2*dxy2z; m3(14,13) = -2*dx2yz; m3(14,14) = 3*dx2y2; 
	m3(15, 0) = -3*dy2z2; m3(15, 1) = -2*dxyz2; m3(15, 2) = -2*dxy2z; 
	m3(16, 0) = 4*dxyz2; m3(16, 1) = 3*dx2z2; m3(16, 2) = 2*dx2yz; 
	m3(17, 0) = -8*dxy2z; m3(17, 1) = -4*dx2yz; m3(17, 2) = -3*dx2y2; 
	m3(15, 3) = 3*dy2z2; m3(15, 4) = -4*dxyz2; m3(15, 5) = -4*dxy2z; 
	m3(16, 3) = -4*dxyz2; m3(16, 4) = 9*dx2z2; m3(16, 5) = 6*dx2yz; 
	m3(17, 3) = 8*dxy2z; m3(17, 4) = -12*dx2yz; m3(17, 5) = -9*dx2y2; 
	m3(15, 6) = -3*dy2z2; m3(15, 7) = 2*dxyz2; m3(15, 8) = -2*dxy2z; 
	m3(16, 6) = 8*dxyz2; m3(16, 7) = -3*dx2z2; m3(16, 8) = 4*dx2yz; 
	m3(17, 6) = -8*dxy2z; m3(17, 7) = 4*dx2yz; m3(17, 8) = -3*dx2y2; 
	m3(15, 9) = 3*dy2z2; m3(15,10) = 4*dxyz2; m3(15,11) = -4*dxy2z; 
	m3(16, 9) = -8*dxyz2; m3(16,10) = -9*dx2z2; m3(16,11) = 12*dx2yz; 
	m3(17, 9) = 8*dxy2z; m3(17,10) = 12*dx2yz; m3(17,11) = -9*dx2y2; 
	m3(15,12) = -3*dy2z2; m3(15,13) = -2*dxyz2; m3(15,14) = 2*dxy2z; 
	m3(16,12) = 4*dxyz2; m3(16,13) = 3*dx2z2; m3(16,14) = -2*dx2yz; 
	m3(17,12) = -4*dxy2z; m3(17,13) = -2*dx2yz; m3(17,14) = 3*dx2y2; 
	m3(15,15) = 3*dy2z2; 
	m3(16,15) = -4*dxyz2; m3(16,16) = 9*dx2z2; 
	m3(17,15) = 4*dxy2z; m3(17,16) = -6*dx2yz; m3(17,17) = 9*dx2y2; 
	m3(18, 0) = 3*dy2z2; m3(18, 1) = 4*dxyz2; m3(18, 2) = 2*dxy2z; 
	m3(19, 0) = -2*dxyz2; m3(19, 1) = -3*dx2z2; m3(19, 2) = -2*dx2yz; 
	m3(20, 0) = -4*dxy2z; m3(20, 1) = -8*dx2yz; m3(20, 2) = -3*dx2y2; 
	m3(18, 3) = -3*dy2z2; m3(18, 4) = 8*dxyz2; m3(18, 5) = 4*dxy2z; 
	m3(19, 3) = 2*dxyz2; m3(19, 4) = -3*dx2z2; m3(19, 5) = -2*dx2yz; 
	m3(20, 3) = 4*dxy2z; m3(20, 4) = -8*dx2yz; m3(20, 5) = -3*dx2y2; 
	m3(18, 6) = 9*dy2z2; m3(18, 7) = -4*dxyz2; m3(18, 8) = 6*dxy2z; 
	m3(19, 6) = -4*dxyz2; m3(19, 7) = 3*dx2z2; m3(19, 8) = -4*dx2yz; 
	m3(20, 6) = -12*dxy2z; m3(20, 7) = 8*dx2yz; m3(20, 8) = -9*dx2y2; 
	m3(18, 9) = -9*dy2z2; m3(18,10) = -8*dxyz2; m3(18,11) = 12*dxy2z; 
	m3(19, 9) = 4*dxyz2; m3(19,10) = 3*dx2z2; m3(19,11) = -4*dx2yz; 
	m3(20, 9) = 12*dxy2z; m3(20,10) = 8*dx2yz; m3(20,11) = -9*dx2y2; 
	m3(18,12) = 3*dy2z2; m3(18,13) = 4*dxyz2; m3(18,14) = -2*dxy2z; 
	m3(19,12) = -2*dxyz2; m3(19,13) = -3*dx2z2; m3(19,14) = 2*dx2yz; 
	m3(20,12) = -2*dxy2z; m3(20,13) = -4*dx2yz; m3(20,14) = 3*dx2y2; 
	m3(18,15) = -3*dy2z2; m3(18,16) = 8*dxyz2; m3(18,17) = -4*dxy2z; 
	m3(19,15) = 2*dxyz2; m3(19,16) = -3*dx2z2; m3(19,17) = 2*dx2yz; 
	m3(20,15) = 2*dxy2z; m3(20,16) = -4*dx2yz; m3(20,17) = 3*dx2y2; 
	m3(18,18) = 9*dy2z2; 
	m3(19,18) = -4*dxyz2; m3(19,19) = 3*dx2z2; 
	m3(20,18) = -6*dxy2z; m3(20,19) = 4*dx2yz; m3(20,20) = 9*dx2y2; 
	m3(21, 0) = -3*dy2z2; m3(21, 1) = -4*dxyz2; m3(21, 2) = -2*dxy2z; 
	m3(22, 0) = -4*dxyz2; m3(22, 1) = -3*dx2z2; m3(22, 2) = -2*dx2yz; 
	m3(23, 0) = -8*dxy2z; m3(23, 1) = -8*dx2yz; m3(23, 2) = -3*dx2y2; 
	m3(21, 3) = 3*dy2z2; m3(21, 4) = -8*dxyz2; m3(21, 5) = -4*dxy2z; 
	m3(22, 3) = 4*dxyz2; m3(22, 4) = -9*dx2z2; m3(22, 5) = -6*dx2yz; 
	m3(23, 3) = 8*dxy2z; m3(23, 4) = -24*dx2yz; m3(23, 5) = -9*dx2y2; 
	m3(21, 6) = -9*dy2z2; m3(21, 7) = 4*dxyz2; m3(21, 8) = -6*dxy2z; 
	m3(22, 6) = -8*dxyz2; m3(22, 7) = 3*dx2z2; m3(22, 8) = -4*dx2yz; 
	m3(23, 6) = -24*dxy2z; m3(23, 7) = 8*dx2yz; m3(23, 8) = -9*dx2y2; 
	m3(21, 9) = 9*dy2z2; m3(21,10) = 8*dxyz2; m3(21,11) = -12*dxy2z; 
	m3(22, 9) = 8*dxyz2; m3(22,10) = 9*dx2z2; m3(22,11) = -12*dx2yz; 
	m3(23, 9) = 24*dxy2z; m3(23,10) = 24*dx2yz; m3(23,11) = -27*dx2y2; 
	m3(21,12) = -3*dy2z2; m3(21,13) = -4*dxyz2; m3(21,14) = 2*dxy2z; 
	m3(22,12) = -4*dxyz2; m3(22,13) = -3*dx2z2; m3(22,14) = 2*dx2yz; 
	m3(23,12) = -4*dxy2z; m3(23,13) = -4*dx2yz; m3(23,14) = 3*dx2y2; 
	m3(21,15) = 3*dy2z2; m3(21,16) = -8*dxyz2; m3(21,17) = 4*dxy2z; 
	m3(22,15) = 4*dxyz2; m3(22,16) = -9*dx2z2; m3(22,17) = 6*dx2yz; 
	m3(23,15) = 4*dxy2z; m3(23,16) = -12*dx2yz; m3(23,17) = 9*dx2y2; 
	m3(21,18) = -9*dy2z2; m3(21,19) = 4*dxyz2; m3(21,20) = 6*dxy2z; 
	m3(22,18) = -8*dxyz2; m3(22,19) = 3*dx2z2; m3(22,20) = 4*dx2yz; 
	m3(23,18) = -12*dxy2z; m3(23,19) = 4*dx2yz; m3(23,20) = 9*dx2y2; 
	m3(21,21) = 9*dy2z2; 
	m3(22,21) = 8*dxyz2; m3(22,22) = 9*dx2z2; 
	m3(23,21) = 12*dxy2z; m3(23,22) = 12*dx2yz; m3(23,23) = 27*dx2y2; 
	m4( 0, 0) = 9*dy2z2; 
	m4( 1, 0) = 8*dxyz2; m4( 1, 1) = 9*dx2z2; 
	m4( 2, 0) = 12*dxy2z; m4( 2, 1) = 12*dx2yz; m4( 2, 2) = 27*dx2y2; 
	m4( 3, 0) = -9*dy2z2; m4( 3, 1) = -8*dxyz2; m4( 3, 2) = -12*dxy2z; 
	m4( 4, 0) = 4*dxyz2; m4( 4, 1) = 3*dx2z2; m4( 4, 2) = 4*dx2yz; 
	m4( 5, 0) = 6*dxy2z; m4( 5, 1) = 4*dx2yz; m4( 5, 2) = 9*dx2y2; 
	m4( 3, 3) = 9*dy2z2; 
	m4( 4, 3) = -4*dxyz2; m4( 4, 4) = 3*dx2z2; 
	m4( 5, 3) = -6*dxy2z; m4( 5, 4) = 4*dx2yz; m4( 5, 5) = 9*dx2y2; 
	m4( 6, 0) = 3*dy2z2; m4( 6, 1) = 4*dxyz2; m4( 6, 2) = 4*dxy2z; 
	m4( 7, 0) = -8*dxyz2; m4( 7, 1) = -9*dx2z2; m4( 7, 2) = -12*dx2yz; 
	m4( 8, 0) = 4*dxy2z; m4( 8, 1) = 6*dx2yz; m4( 8, 2) = 9*dx2y2; 
	m4( 6, 3) = -3*dy2z2; m4( 6, 4) = 2*dxyz2; m4( 6, 5) = 2*dxy2z; 
	m4( 7, 3) = 8*dxyz2; m4( 7, 4) = -3*dx2z2; m4( 7, 5) = -4*dx2yz; 
	m4( 8, 3) = -4*dxy2z; m4( 8, 4) = 2*dx2yz; m4( 8, 5) = 3*dx2y2; 
	m4( 6, 6) = 3*dy2z2; 
	m4( 7, 6) = -4*dxyz2; m4( 7, 7) = 9*dx2z2; 
	m4( 8, 6) = 4*dxy2z; m4( 8, 7) = -6*dx2yz; m4( 8, 8) = 9*dx2y2; 
	m4( 9, 0) = -3*dy2z2; m4( 9, 1) = -4*dxyz2; m4( 9, 2) = -4*dxy2z; 
	m4(10, 0) = -4*dxyz2; m4(10, 1) = -3*dx2z2; m4(10, 2) = -4*dx2yz; 
	m4(11, 0) = 2*dxy2z; m4(11, 1) = 2*dx2yz; m4(11, 2) = 3*dx2y2; 
	m4( 9, 3) = 3*dy2z2; m4( 9, 4) = -2*dxyz2; m4( 9, 5) = -2*dxy2z; 
	m4(10, 3) = 4*dxyz2; m4(10, 4) = -3*dx2z2; m4(10, 5) = -4*dx2yz; 
	m4(11, 3) = -2*dxy2z; m4(11, 4) = 2*dx2yz; m4(11, 5) = 3*dx2y2; 
	m4( 9, 6) = -3*dy2z2; m4( 9, 7) = 4*dxyz2; m4( 9, 8) = -4*dxy2z; 
	m4(10, 6) = -2*dxyz2; m4(10, 7) = 3*dx2z2; m4(10, 8) = -2*dx2yz; 
	m4(11, 6) = 2*dxy2z; m4(11, 7) = -2*dx2yz; m4(11, 8) = 3*dx2y2; 
	m4( 9, 9) = 3*dy2z2; 
	m4(10, 9) = 2*dxyz2; m4(10,10) = 3*dx2z2; 
	m4(11, 9) = -2*dxy2z; m4(11,10) = -2*dx2yz; m4(11,11) = 3*dx2y2; 
	m4(12, 0) = 9*dy2z2; m4(12, 1) = 8*dxyz2; m4(12, 2) = 24*dxy2z; 
	m4(13, 0) = 8*dxyz2; m4(13, 1) = 9*dx2z2; m4(13, 2) = 24*dx2yz; 
	m4(14, 0) = -12*dxy2z; m4(14, 1) = -12*dx2yz; m4(14, 2) = -27*dx2y2; 
	m4(12, 3) = -9*dy2z2; m4(12, 4) = 4*dxyz2; m4(12, 5) = 12*dxy2z; 
	m4(13, 3) = -8*dxyz2; m4(13, 4) = 3*dx2z2; m4(13, 5) = 8*dx2yz; 
	m4(14, 3) = 12*dxy2z; m4(14, 4) = -4*dx2yz; m4(14, 5) = -9*dx2y2; 
	m4(12, 6) = 3*dy2z2; m4(12, 7) = -8*dxyz2; m4(12, 8) = 8*dxy2z; 
	m4(13, 6) = 4*dxyz2; m4(13, 7) = -9*dx2z2; m4(13, 8) = 12*dx2yz; 
	m4(14, 6) = -4*dxy2z; m4(14, 7) = 12*dx2yz; m4(14, 8) = -9*dx2y2; 
	m4(12, 9) = -3*dy2z2; m4(12,10) = -4*dxyz2; m4(12,11) = 4*dxy2z; 
	m4(13, 9) = -4*dxyz2; m4(13,10) = -3*dx2z2; m4(13,11) = 4*dx2yz; 
	m4(14, 9) = 4*dxy2z; m4(14,10) = 4*dx2yz; m4(14,11) = -3*dx2y2; 
	m4(12,12) = 27*dy2z2; 
	m4(13,12) = 24*dxyz2; m4(13,13) = 27*dx2z2; 
	m4(14,12) = -24*dxy2z; m4(14,13) = -24*dx2yz; m4(14,14) = 27*dx2y2; 
	m4(15, 0) = -9*dy2z2; m4(15, 1) = -8*dxyz2; m4(15, 2) = -24*dxy2z; 
	m4(16, 0) = 4*dxyz2; m4(16, 1) = 3*dx2z2; m4(16, 2) = 8*dx2yz; 
	m4(17, 0) = -6*dxy2z; m4(17, 1) = -4*dx2yz; m4(17, 2) = -9*dx2y2; 
	m4(15, 3) = 9*dy2z2; m4(15, 4) = -4*dxyz2; m4(15, 5) = -12*dxy2z; 
	m4(16, 3) = -4*dxyz2; m4(16, 4) = 3*dx2z2; m4(16, 5) = 8*dx2yz; 
	m4(17, 3) = 6*dxy2z; m4(17, 4) = -4*dx2yz; m4(17, 5) = -9*dx2y2; 
	m4(15, 6) = -3*dy2z2; m4(15, 7) = 8*dxyz2; m4(15, 8) = -8*dxy2z; 
	m4(16, 6) = 2*dxyz2; m4(16, 7) = -3*dx2z2; m4(16, 8) = 4*dx2yz; 
	m4(17, 6) = -2*dxy2z; m4(17, 7) = 4*dx2yz; m4(17, 8) = -3*dx2y2; 
	m4(15, 9) = 3*dy2z2; m4(15,10) = 4*dxyz2; m4(15,11) = -4*dxy2z; 
	m4(16, 9) = -2*dxyz2; m4(16,10) = -3*dx2z2; m4(16,11) = 4*dx2yz; 
	m4(17, 9) = 2*dxy2z; m4(17,10) = 4*dx2yz; m4(17,11) = -3*dx2y2; 
	m4(15,12) = -27*dy2z2; m4(15,13) = -24*dxyz2; m4(15,14) = 24*dxy2z; 
	m4(16,12) = 12*dxyz2; m4(16,13) = 9*dx2z2; m4(16,14) = -8*dx2yz; 
	m4(17,12) = -12*dxy2z; m4(17,13) = -8*dx2yz; m4(17,14) = 9*dx2y2; 
	m4(15,15) = 27*dy2z2; 
	m4(16,15) = -12*dxyz2; m4(16,16) = 9*dx2z2; 
	m4(17,15) = 12*dxy2z; m4(17,16) = -8*dx2yz; m4(17,17) = 9*dx2y2; 
	m4(18, 0) = 3*dy2z2; m4(18, 1) = 4*dxyz2; m4(18, 2) = 8*dxy2z; 
	m4(19, 0) = -8*dxyz2; m4(19, 1) = -9*dx2z2; m4(19, 2) = -24*dx2yz; 
	m4(20, 0) = -4*dxy2z; m4(20, 1) = -6*dx2yz; m4(20, 2) = -9*dx2y2; 
	m4(18, 3) = -3*dy2z2; m4(18, 4) = 2*dxyz2; m4(18, 5) = 4*dxy2z; 
	m4(19, 3) = 8*dxyz2; m4(19, 4) = -3*dx2z2; m4(19, 5) = -8*dx2yz; 
	m4(20, 3) = 4*dxy2z; m4(20, 4) = -2*dx2yz; m4(20, 5) = -3*dx2y2; 
	m4(18, 6) = 3*dy2z2; m4(18, 7) = -4*dxyz2; m4(18, 8) = 8*dxy2z; 
	m4(19, 6) = -4*dxyz2; m4(19, 7) = 9*dx2z2; m4(19, 8) = -12*dx2yz; 
	m4(20, 6) = -4*dxy2z; m4(20, 7) = 6*dx2yz; m4(20, 8) = -9*dx2y2; 
	m4(18, 9) = -3*dy2z2; m4(18,10) = -2*dxyz2; m4(18,11) = 4*dxy2z; 
	m4(19, 9) = 4*dxyz2; m4(19,10) = 3*dx2z2; m4(19,11) = -4*dx2yz; 
	m4(20, 9) = 4*dxy2z; m4(20,10) = 2*dx2yz; m4(20,11) = -3*dx2y2; 
	m4(18,12) = 9*dy2z2; m4(18,13) = 12*dxyz2; m4(18,14) = -8*dxy2z; 
	m4(19,12) = -24*dxyz2; m4(19,13) = -27*dx2z2; m4(19,14) = 24*dx2yz; 
	m4(20,12) = -8*dxy2z; m4(20,13) = -12*dx2yz; m4(20,14) = 9*dx2y2; 
	m4(18,15) = -9*dy2z2; m4(18,16) = 6*dxyz2; m4(18,17) = -4*dxy2z; 
	m4(19,15) = 24*dxyz2; m4(19,16) = -9*dx2z2; m4(19,17) = 8*dx2yz; 
	m4(20,15) = 8*dxy2z; m4(20,16) = -4*dx2yz; m4(20,17) = 3*dx2y2; 
	m4(18,18) = 9*dy2z2; 
	m4(19,18) = -12*dxyz2; m4(19,19) = 27*dx2z2; 
	m4(20,18) = -8*dxy2z; m4(20,19) = 12*dx2yz; m4(20,20) = 9*dx2y2; 
	m4(21, 0) = -3*dy2z2; m4(21, 1) = -4*dxyz2; m4(21, 2) = -8*dxy2z; 
	m4(22, 0) = -4*dxyz2; m4(22, 1) = -3*dx2z2; m4(22, 2) = -8*dx2yz; 
	m4(23, 0) = -2*dxy2z; m4(23, 1) = -2*dx2yz; m4(23, 2) = -3*dx2y2; 
	m4(21, 3) = 3*dy2z2; m4(21, 4) = -2*dxyz2; m4(21, 5) = -4*dxy2z; 
	m4(22, 3) = 4*dxyz2; m4(22, 4) = -3*dx2z2; m4(22, 5) = -8*dx2yz; 
	m4(23, 3) = 2*dxy2z; m4(23, 4) = -2*dx2yz; m4(23, 5) = -3*dx2y2; 
	m4(21, 6) = -3*dy2z2; m4(21, 7) = 4*dxyz2; m4(21, 8) = -8*dxy2z; 
	m4(22, 6) = -2*dxyz2; m4(22, 7) = 3*dx2z2; m4(22, 8) = -4*dx2yz; 
	m4(23, 6) = -2*dxy2z; m4(23, 7) = 2*dx2yz; m4(23, 8) = -3*dx2y2; 
	m4(21, 9) = 3*dy2z2; m4(21,10) = 2*dxyz2; m4(21,11) = -4*dxy2z; 
	m4(22, 9) = 2*dxyz2; m4(22,10) = 3*dx2z2; m4(22,11) = -4*dx2yz; 
	m4(23, 9) = 2*dxy2z; m4(23,10) = 2*dx2yz; m4(23,11) = -3*dx2y2; 
	m4(21,12) = -9*dy2z2; m4(21,13) = -12*dxyz2; m4(21,14) = 8*dxy2z; 
	m4(22,12) = -12*dxyz2; m4(22,13) = -9*dx2z2; m4(22,14) = 8*dx2yz; 
	m4(23,12) = -4*dxy2z; m4(23,13) = -4*dx2yz; m4(23,14) = 3*dx2y2; 
	m4(21,15) = 9*dy2z2; m4(21,16) = -6*dxyz2; m4(21,17) = 4*dxy2z; 
	m4(22,15) = 12*dxyz2; m4(22,16) = -9*dx2z2; m4(22,17) = 8*dx2yz; 
	m4(23,15) = 4*dxy2z; m4(23,16) = -4*dx2yz; m4(23,17) = 3*dx2y2; 
	m4(21,18) = -9*dy2z2; m4(21,19) = 12*dxyz2; m4(21,20) = 8*dxy2z; 
	m4(22,18) = -6*dxyz2; m4(22,19) = 9*dx2z2; m4(22,20) = 4*dx2yz; 
	m4(23,18) = -4*dxy2z; m4(23,19) = 4*dx2yz; m4(23,20) = 3*dx2y2; 
	m4(21,21) = 9*dy2z2; 
	m4(22,21) = 6*dxyz2; m4(22,22) = 9*dx2z2; 
	m4(23,21) = 4*dxy2z; m4(23,22) = 4*dx2yz; m4(23,23) = 3*dx2y2; 
	m5( 0, 0) = 9*dy2z2; 
	m5( 1, 0) = 4*dxyz2; m5( 1, 1) = 3*dx2z2; 
	m5( 2, 0) = 6*dxy2z; m5( 2, 1) = 4*dx2yz; m5( 2, 2) = 9*dx2y2; 
	m5( 3, 0) = -9*dy2z2; m5( 3, 1) = -4*dxyz2; m5( 3, 2) = -6*dxy2z; 
	m5( 4, 0) = 8*dxyz2; m5( 4, 1) = 3*dx2z2; m5( 4, 2) = 4*dx2yz; 
	m5( 5, 0) = 12*dxy2z; m5( 5, 1) = 4*dx2yz; m5( 5, 2) = 9*dx2y2; 
	m5( 3, 3) = 9*dy2z2; 
	m5( 4, 3) = -8*dxyz2; m5( 4, 4) = 9*dx2z2; 
	m5( 5, 3) = -12*dxy2z; m5( 5, 4) = 12*dx2yz; m5( 5, 5) = 27*dx2y2; 
	m5( 6, 0) = 3*dy2z2; m5( 6, 1) = 2*dxyz2; m5( 6, 2) = 2*dxy2z; 
	m5( 7, 0) = -4*dxyz2; m5( 7, 1) = -3*dx2z2; m5( 7, 2) = -4*dx2yz; 
	m5( 8, 0) = 2*dxy2z; m5( 8, 1) = 2*dx2yz; m5( 8, 2) = 3*dx2y2; 
	m5( 6, 3) = -3*dy2z2; m5( 6, 4) = 4*dxyz2; m5( 6, 5) = 4*dxy2z; 
	m5( 7, 3) = 4*dxyz2; m5( 7, 4) = -3*dx2z2; m5( 7, 5) = -4*dx2yz; 
	m5( 8, 3) = -2*dxy2z; m5( 8, 4) = 2*dx2yz; m5( 8, 5) = 3*dx2y2; 
	m5( 6, 6) = 3*dy2z2; 
	m5( 7, 6) = -2*dxyz2; m5( 7, 7) = 3*dx2z2; 
	m5( 8, 6) = 2*dxy2z; m5( 8, 7) = -2*dx2yz; m5( 8, 8) = 3*dx2y2; 
	m5( 9, 0) = -3*dy2z2; m5( 9, 1) = -2*dxyz2; m5( 9, 2) = -2*dxy2z; 
	m5(10, 0) = -8*dxyz2; m5(10, 1) = -3*dx2z2; m5(10, 2) = -4*dx2yz; 
	m5(11, 0) = 4*dxy2z; m5(11, 1) = 2*dx2yz; m5(11, 2) = 3*dx2y2; 
	m5( 9, 3) = 3*dy2z2; m5( 9, 4) = -4*dxyz2; m5( 9, 5) = -4*dxy2z; 
	m5(10, 3) = 8*dxyz2; m5(10, 4) = -9*dx2z2; m5(10, 5) = -12*dx2yz; 
	m5(11, 3) = -4*dxy2z; m5(11, 4) = 6*dx2yz; m5(11, 5) = 9*dx2y2; 
	m5( 9, 6) = -3*dy2z2; m5( 9, 7) = 2*dxyz2; m5( 9, 8) = -2*dxy2z; 
	m5(10, 6) = -4*dxyz2; m5(10, 7) = 3*dx2z2; m5(10, 8) = -2*dx2yz; 
	m5(11, 6) = 4*dxy2z; m5(11, 7) = -2*dx2yz; m5(11, 8) = 3*dx2y2; 
	m5( 9, 9) = 3*dy2z2; 
	m5(10, 9) = 4*dxyz2; m5(10,10) = 9*dx2z2; 
	m5(11, 9) = -4*dxy2z; m5(11,10) = -6*dx2yz; m5(11,11) = 9*dx2y2; 
	m5(12, 0) = 9*dy2z2; m5(12, 1) = 4*dxyz2; m5(12, 2) = 12*dxy2z; 
	m5(13, 0) = 4*dxyz2; m5(13, 1) = 3*dx2z2; m5(13, 2) = 8*dx2yz; 
	m5(14, 0) = -6*dxy2z; m5(14, 1) = -4*dx2yz; m5(14, 2) = -9*dx2y2; 
	m5(12, 3) = -9*dy2z2; m5(12, 4) = 8*dxyz2; m5(12, 5) = 24*dxy2z; 
	m5(13, 3) = -4*dxyz2; m5(13, 4) = 3*dx2z2; m5(13, 5) = 8*dx2yz; 
	m5(14, 3) = 6*dxy2z; m5(14, 4) = -4*dx2yz; m5(14, 5) = -9*dx2y2; 
	m5(12, 6) = 3*dy2z2; m5(12, 7) = -4*dxyz2; m5(12, 8) = 4*dxy2z; 
	m5(13, 6) = 2*dxyz2; m5(13, 7) = -3*dx2z2; m5(13, 8) = 4*dx2yz; 
	m5(14, 6) = -2*dxy2z; m5(14, 7) = 4*dx2yz; m5(14, 8) = -3*dx2y2; 
	m5(12, 9) = -3*dy2z2; m5(12,10) = -8*dxyz2; m5(12,11) = 8*dxy2z; 
	m5(13, 9) = -2*dxyz2; m5(13,10) = -3*dx2z2; m5(13,11) = 4*dx2yz; 
	m5(14, 9) = 2*dxy2z; m5(14,10) = 4*dx2yz; m5(14,11) = -3*dx2y2; 
	m5(12,12) = 27*dy2z2; 
	m5(13,12) = 12*dxyz2; m5(13,13) = 9*dx2z2; 
	m5(14,12) = -12*dxy2z; m5(14,13) = -8*dx2yz; m5(14,14) = 9*dx2y2; 
	m5(15, 0) = -9*dy2z2; m5(15, 1) = -4*dxyz2; m5(15, 2) = -12*dxy2z; 
	m5(16, 0) = 8*dxyz2; m5(16, 1) = 3*dx2z2; m5(16, 2) = 8*dx2yz; 
	m5(17, 0) = -12*dxy2z; m5(17, 1) = -4*dx2yz; m5(17, 2) = -9*dx2y2; 
	m5(15, 3) = 9*dy2z2; m5(15, 4) = -8*dxyz2; m5(15, 5) = -24*dxy2z; 
	m5(16, 3) = -8*dxyz2; m5(16, 4) = 9*dx2z2; m5(16, 5) = 24*dx2yz; 
	m5(17, 3) = 12*dxy2z; m5(17, 4) = -12*dx2yz; m5(17, 5) = -27*dx2y2; 
	m5(15, 6) = -3*dy2z2; m5(15, 7) = 4*dxyz2; m5(15, 8) = -4*dxy2z; 
	m5(16, 6) = 4*dxyz2; m5(16, 7) = -3*dx2z2; m5(16, 8) = 4*dx2yz; 
	m5(17, 6) = -4*dxy2z; m5(17, 7) = 4*dx2yz; m5(17, 8) = -3*dx2y2; 
	m5(15, 9) = 3*dy2z2; m5(15,10) = 8*dxyz2; m5(15,11) = -8*dxy2z; 
	m5(16, 9) = -4*dxyz2; m5(16,10) = -9*dx2z2; m5(16,11) = 12*dx2yz; 
	m5(17, 9) = 4*dxy2z; m5(17,10) = 12*dx2yz; m5(17,11) = -9*dx2y2; 
	m5(15,12) = -27*dy2z2; m5(15,13) = -12*dxyz2; m5(15,14) = 12*dxy2z; 
	m5(16,12) = 24*dxyz2; m5(16,13) = 9*dx2z2; m5(16,14) = -8*dx2yz; 
	m5(17,12) = -24*dxy2z; m5(17,13) = -8*dx2yz; m5(17,14) = 9*dx2y2; 
	m5(15,15) = 27*dy2z2; 
	m5(16,15) = -24*dxyz2; m5(16,16) = 27*dx2z2; 
	m5(17,15) = 24*dxy2z; m5(17,16) = -24*dx2yz; m5(17,17) = 27*dx2y2; 
	m5(18, 0) = 3*dy2z2; m5(18, 1) = 2*dxyz2; m5(18, 2) = 4*dxy2z; 
	m5(19, 0) = -4*dxyz2; m5(19, 1) = -3*dx2z2; m5(19, 2) = -8*dx2yz; 
	m5(20, 0) = -2*dxy2z; m5(20, 1) = -2*dx2yz; m5(20, 2) = -3*dx2y2; 
	m5(18, 3) = -3*dy2z2; m5(18, 4) = 4*dxyz2; m5(18, 5) = 8*dxy2z; 
	m5(19, 3) = 4*dxyz2; m5(19, 4) = -3*dx2z2; m5(19, 5) = -8*dx2yz; 
	m5(20, 3) = 2*dxy2z; m5(20, 4) = -2*dx2yz; m5(20, 5) = -3*dx2y2; 
	m5(18, 6) = 3*dy2z2; m5(18, 7) = -2*dxyz2; m5(18, 8) = 4*dxy2z; 
	m5(19, 6) = -2*dxyz2; m5(19, 7) = 3*dx2z2; m5(19, 8) = -4*dx2yz; 
	m5(20, 6) = -2*dxy2z; m5(20, 7) = 2*dx2yz; m5(20, 8) = -3*dx2y2; 
	m5(18, 9) = -3*dy2z2; m5(18,10) = -4*dxyz2; m5(18,11) = 8*dxy2z; 
	m5(19, 9) = 2*dxyz2; m5(19,10) = 3*dx2z2; m5(19,11) = -4*dx2yz; 
	m5(20, 9) = 2*dxy2z; m5(20,10) = 2*dx2yz; m5(20,11) = -3*dx2y2; 
	m5(18,12) = 9*dy2z2; m5(18,13) = 6*dxyz2; m5(18,14) = -4*dxy2z; 
	m5(19,12) = -12*dxyz2; m5(19,13) = -9*dx2z2; m5(19,14) = 8*dx2yz; 
	m5(20,12) = -4*dxy2z; m5(20,13) = -4*dx2yz; m5(20,14) = 3*dx2y2; 
	m5(18,15) = -9*dy2z2; m5(18,16) = 12*dxyz2; m5(18,17) = -8*dxy2z; 
	m5(19,15) = 12*dxyz2; m5(19,16) = -9*dx2z2; m5(19,17) = 8*dx2yz; 
	m5(20,15) = 4*dxy2z; m5(20,16) = -4*dx2yz; m5(20,17) = 3*dx2y2; 
	m5(18,18) = 9*dy2z2; 
	m5(19,18) = -6*dxyz2; m5(19,19) = 9*dx2z2; 
	m5(20,18) = -4*dxy2z; m5(20,19) = 4*dx2yz; m5(20,20) = 3*dx2y2; 
	m5(21, 0) = -3*dy2z2; m5(21, 1) = -2*dxyz2; m5(21, 2) = -4*dxy2z; 
	m5(22, 0) = -8*dxyz2; m5(22, 1) = -3*dx2z2; m5(22, 2) = -8*dx2yz; 
	m5(23, 0) = -4*dxy2z; m5(23, 1) = -2*dx2yz; m5(23, 2) = -3*dx2y2; 
	m5(21, 3) = 3*dy2z2; m5(21, 4) = -4*dxyz2; m5(21, 5) = -8*dxy2z; 
	m5(22, 3) = 8*dxyz2; m5(22, 4) = -9*dx2z2; m5(22, 5) = -24*dx2yz; 
	m5(23, 3) = 4*dxy2z; m5(23, 4) = -6*dx2yz; m5(23, 5) = -9*dx2y2; 
	m5(21, 6) = -3*dy2z2; m5(21, 7) = 2*dxyz2; m5(21, 8) = -4*dxy2z; 
	m5(22, 6) = -4*dxyz2; m5(22, 7) = 3*dx2z2; m5(22, 8) = -4*dx2yz; 
	m5(23, 6) = -4*dxy2z; m5(23, 7) = 2*dx2yz; m5(23, 8) = -3*dx2y2; 
	m5(21, 9) = 3*dy2z2; m5(21,10) = 4*dxyz2; m5(21,11) = -8*dxy2z; 
	m5(22, 9) = 4*dxyz2; m5(22,10) = 9*dx2z2; m5(22,11) = -12*dx2yz; 
	m5(23, 9) = 4*dxy2z; m5(23,10) = 6*dx2yz; m5(23,11) = -9*dx2y2; 
	m5(21,12) = -9*dy2z2; m5(21,13) = -6*dxyz2; m5(21,14) = 4*dxy2z; 
	m5(22,12) = -24*dxyz2; m5(22,13) = -9*dx2z2; m5(22,14) = 8*dx2yz; 
	m5(23,12) = -8*dxy2z; m5(23,13) = -4*dx2yz; m5(23,14) = 3*dx2y2; 
	m5(21,15) = 9*dy2z2; m5(21,16) = -12*dxyz2; m5(21,17) = 8*dxy2z; 
	m5(22,15) = 24*dxyz2; m5(22,16) = -27*dx2z2; m5(22,17) = 24*dx2yz; 
	m5(23,15) = 8*dxy2z; m5(23,16) = -12*dx2yz; m5(23,17) = 9*dx2y2; 
	m5(21,18) = -9*dy2z2; m5(21,19) = 6*dxyz2; m5(21,20) = 4*dxy2z; 
	m5(22,18) = -12*dxyz2; m5(22,19) = 9*dx2z2; m5(22,20) = 4*dx2yz; 
	m5(23,18) = -8*dxy2z; m5(23,19) = 4*dx2yz; m5(23,20) = 3*dx2y2; 
	m5(21,21) = 9*dy2z2; 
	m5(22,21) = 12*dxyz2; m5(22,22) = 27*dx2z2; 
	m5(23,21) = 8*dxy2z; m5(23,22) = 12*dx2yz; m5(23,23) = 9*dx2y2; 
	m6( 0, 0) = 3*dy2z2; 
	m6( 1, 0) = 4*dxyz2; m6( 1, 1) = 9*dx2z2; 
	m6( 2, 0) = 4*dxy2z; m6( 2, 1) = 6*dx2yz; m6( 2, 2) = 9*dx2y2; 
	m6( 3, 0) = -3*dy2z2; m6( 3, 1) = -4*dxyz2; m6( 3, 2) = -4*dxy2z; 
	m6( 4, 0) = 2*dxyz2; m6( 4, 1) = 3*dx2z2; m6( 4, 2) = 2*dx2yz; 
	m6( 5, 0) = 2*dxy2z; m6( 5, 1) = 2*dx2yz; m6( 5, 2) = 3*dx2y2; 
	m6( 3, 3) = 3*dy2z2; 
	m6( 4, 3) = -2*dxyz2; m6( 4, 4) = 3*dx2z2; 
	m6( 5, 3) = -2*dxy2z; m6( 5, 4) = 2*dx2yz; m6( 5, 5) = 3*dx2y2; 
	m6( 6, 0) = 3*dy2z2; m6( 6, 1) = 8*dxyz2; m6( 6, 2) = 4*dxy2z; 
	m6( 7, 0) = -4*dxyz2; m6( 7, 1) = -9*dx2z2; m6( 7, 2) = -6*dx2yz; 
	m6( 8, 0) = 4*dxy2z; m6( 8, 1) = 12*dx2yz; m6( 8, 2) = 9*dx2y2; 
	m6( 6, 3) = -3*dy2z2; m6( 6, 4) = 4*dxyz2; m6( 6, 5) = 2*dxy2z; 
	m6( 7, 3) = 4*dxyz2; m6( 7, 4) = -3*dx2z2; m6( 7, 5) = -2*dx2yz; 
	m6( 8, 3) = -4*dxy2z; m6( 8, 4) = 4*dx2yz; m6( 8, 5) = 3*dx2y2; 
	m6( 6, 6) = 9*dy2z2; 
	m6( 7, 6) = -8*dxyz2; m6( 7, 7) = 9*dx2z2; 
	m6( 8, 6) = 12*dxy2z; m6( 8, 7) = -12*dx2yz; m6( 8, 8) = 27*dx2y2; 
	m6( 9, 0) = -3*dy2z2; m6( 9, 1) = -8*dxyz2; m6( 9, 2) = -4*dxy2z; 
	m6(10, 0) = -2*dxyz2; m6(10, 1) = -3*dx2z2; m6(10, 2) = -2*dx2yz; 
	m6(11, 0) = 2*dxy2z; m6(11, 1) = 4*dx2yz; m6(11, 2) = 3*dx2y2; 
	m6( 9, 3) = 3*dy2z2; m6( 9, 4) = -4*dxyz2; m6( 9, 5) = -2*dxy2z; 
	m6(10, 3) = 2*dxyz2; m6(10, 4) = -3*dx2z2; m6(10, 5) = -2*dx2yz; 
	m6(11, 3) = -2*dxy2z; m6(11, 4) = 4*dx2yz; m6(11, 5) = 3*dx2y2; 
	m6( 9, 6) = -9*dy2z2; m6( 9, 7) = 8*dxyz2; m6( 9, 8) = -12*dxy2z; 
	m6(10, 6) = -4*dxyz2; m6(10, 7) = 3*dx2z2; m6(10, 8) = -4*dx2yz; 
	m6(11, 6) = 6*dxy2z; m6(11, 7) = -4*dx2yz; m6(11, 8) = 9*dx2y2; 
	m6( 9, 9) = 9*dy2z2; 
	m6(10, 9) = 4*dxyz2; m6(10,10) = 3*dx2z2; 
	m6(11, 9) = -6*dxy2z; m6(11,10) = -4*dx2yz; m6(11,11) = 9*dx2y2; 
	m6(12, 0) = 3*dy2z2; m6(12, 1) = 4*dxyz2; m6(12, 2) = 8*dxy2z; 
	m6(13, 0) = 4*dxyz2; m6(13, 1) = 9*dx2z2; m6(13, 2) = 12*dx2yz; 
	m6(14, 0) = -4*dxy2z; m6(14, 1) = -6*dx2yz; m6(14, 2) = -9*dx2y2; 
	m6(12, 3) = -3*dy2z2; m6(12, 4) = 2*dxyz2; m6(12, 5) = 4*dxy2z; 
	m6(13, 3) = -4*dxyz2; m6(13, 4) = 3*dx2z2; m6(13, 5) = 4*dx2yz; 
	m6(14, 3) = 4*dxy2z; m6(14, 4) = -2*dx2yz; m6(14, 5) = -3*dx2y2; 
	m6(12, 6) = 3*dy2z2; m6(12, 7) = -4*dxyz2; m6(12, 8) = 8*dxy2z; 
	m6(13, 6) = 8*dxyz2; m6(13, 7) = -9*dx2z2; m6(13, 8) = 24*dx2yz; 
	m6(14, 6) = -4*dxy2z; m6(14, 7) = 6*dx2yz; m6(14, 8) = -9*dx2y2; 
	m6(12, 9) = -3*dy2z2; m6(12,10) = -2*dxyz2; m6(12,11) = 4*dxy2z; 
	m6(13, 9) = -8*dxyz2; m6(13,10) = -3*dx2z2; m6(13,11) = 8*dx2yz; 
	m6(14, 9) = 4*dxy2z; m6(14,10) = 2*dx2yz; m6(14,11) = -3*dx2y2; 
	m6(12,12) = 9*dy2z2; 
	m6(13,12) = 12*dxyz2; m6(13,13) = 27*dx2z2; 
	m6(14,12) = -8*dxy2z; m6(14,13) = -12*dx2yz; m6(14,14) = 9*dx2y2; 
	m6(15, 0) = -3*dy2z2; m6(15, 1) = -4*dxyz2; m6(15, 2) = -8*dxy2z; 
	m6(16, 0) = 2*dxyz2; m6(16, 1) = 3*dx2z2; m6(16, 2) = 4*dx2yz; 
	m6(17, 0) = -2*dxy2z; m6(17, 1) = -2*dx2yz; m6(17, 2) = -3*dx2y2; 
	m6(15, 3) = 3*dy2z2; m6(15, 4) = -2*dxyz2; m6(15, 5) = -4*dxy2z; 
	m6(16, 3) = -2*dxyz2; m6(16, 4) = 3*dx2z2; m6(16, 5) = 4*dx2yz; 
	m6(17, 3) = 2*dxy2z; m6(17, 4) = -2*dx2yz; m6(17, 5) = -3*dx2y2; 
	m6(15, 6) = -3*dy2z2; m6(15, 7) = 4*dxyz2; m6(15, 8) = -8*dxy2z; 
	m6(16, 6) = 4*dxyz2; m6(16, 7) = -3*dx2z2; m6(16, 8) = 8*dx2yz; 
	m6(17, 6) = -2*dxy2z; m6(17, 7) = 2*dx2yz; m6(17, 8) = -3*dx2y2; 
	m6(15, 9) = 3*dy2z2; m6(15,10) = 2*dxyz2; m6(15,11) = -4*dxy2z; 
	m6(16, 9) = -4*dxyz2; m6(16,10) = -3*dx2z2; m6(16,11) = 8*dx2yz; 
	m6(17, 9) = 2*dxy2z; m6(17,10) = 2*dx2yz; m6(17,11) = -3*dx2y2; 
	m6(15,12) = -9*dy2z2; m6(15,13) = -12*dxyz2; m6(15,14) = 8*dxy2z; 
	m6(16,12) = 6*dxyz2; m6(16,13) = 9*dx2z2; m6(16,14) = -4*dx2yz; 
	m6(17,12) = -4*dxy2z; m6(17,13) = -4*dx2yz; m6(17,14) = 3*dx2y2; 
	m6(15,15) = 9*dy2z2; 
	m6(16,15) = -6*dxyz2; m6(16,16) = 9*dx2z2; 
	m6(17,15) = 4*dxy2z; m6(17,16) = -4*dx2yz; m6(17,17) = 3*dx2y2; 
	m6(18, 0) = 3*dy2z2; m6(18, 1) = 8*dxyz2; m6(18, 2) = 8*dxy2z; 
	m6(19, 0) = -4*dxyz2; m6(19, 1) = -9*dx2z2; m6(19, 2) = -12*dx2yz; 
	m6(20, 0) = -4*dxy2z; m6(20, 1) = -12*dx2yz; m6(20, 2) = -9*dx2y2; 
	m6(18, 3) = -3*dy2z2; m6(18, 4) = 4*dxyz2; m6(18, 5) = 4*dxy2z; 
	m6(19, 3) = 4*dxyz2; m6(19, 4) = -3*dx2z2; m6(19, 5) = -4*dx2yz; 
	m6(20, 3) = 4*dxy2z; m6(20, 4) = -4*dx2yz; m6(20, 5) = -3*dx2y2; 
	m6(18, 6) = 9*dy2z2; m6(18, 7) = -8*dxyz2; m6(18, 8) = 24*dxy2z; 
	m6(19, 6) = -8*dxyz2; m6(19, 7) = 9*dx2z2; m6(19, 8) = -24*dx2yz; 
	m6(20, 6) = -12*dxy2z; m6(20, 7) = 12*dx2yz; m6(20, 8) = -27*dx2y2; 
	m6(18, 9) = -9*dy2z2; m6(18,10) = -4*dxyz2; m6(18,11) = 12*dxy2z; 
	m6(19, 9) = 8*dxyz2; m6(19,10) = 3*dx2z2; m6(19,11) = -8*dx2yz; 
	m6(20, 9) = 12*dxy2z; m6(20,10) = 4*dx2yz; m6(20,11) = -9*dx2y2; 
	m6(18,12) = 9*dy2z2; m6(18,13) = 24*dxyz2; m6(18,14) = -8*dxy2z; 
	m6(19,12) = -12*dxyz2; m6(19,13) = -27*dx2z2; m6(19,14) = 12*dx2yz; 
	m6(20,12) = -8*dxy2z; m6(20,13) = -24*dx2yz; m6(20,14) = 9*dx2y2; 
	m6(18,15) = -9*dy2z2; m6(18,16) = 12*dxyz2; m6(18,17) = -4*dxy2z; 
	m6(19,15) = 12*dxyz2; m6(19,16) = -9*dx2z2; m6(19,17) = 4*dx2yz; 
	m6(20,15) = 8*dxy2z; m6(20,16) = -8*dx2yz; m6(20,17) = 3*dx2y2; 
	m6(18,18) = 27*dy2z2; 
	m6(19,18) = -24*dxyz2; m6(19,19) = 27*dx2z2; 
	m6(20,18) = -24*dxy2z; m6(20,19) = 24*dx2yz; m6(20,20) = 27*dx2y2; 
	m6(21, 0) = -3*dy2z2; m6(21, 1) = -8*dxyz2; m6(21, 2) = -8*dxy2z; 
	m6(22, 0) = -2*dxyz2; m6(22, 1) = -3*dx2z2; m6(22, 2) = -4*dx2yz; 
	m6(23, 0) = -2*dxy2z; m6(23, 1) = -4*dx2yz; m6(23, 2) = -3*dx2y2; 
	m6(21, 3) = 3*dy2z2; m6(21, 4) = -4*dxyz2; m6(21, 5) = -4*dxy2z; 
	m6(22, 3) = 2*dxyz2; m6(22, 4) = -3*dx2z2; m6(22, 5) = -4*dx2yz; 
	m6(23, 3) = 2*dxy2z; m6(23, 4) = -4*dx2yz; m6(23, 5) = -3*dx2y2; 
	m6(21, 6) = -9*dy2z2; m6(21, 7) = 8*dxyz2; m6(21, 8) = -24*dxy2z; 
	m6(22, 6) = -4*dxyz2; m6(22, 7) = 3*dx2z2; m6(22, 8) = -8*dx2yz; 
	m6(23, 6) = -6*dxy2z; m6(23, 7) = 4*dx2yz; m6(23, 8) = -9*dx2y2; 
	m6(21, 9) = 9*dy2z2; m6(21,10) = 4*dxyz2; m6(21,11) = -12*dxy2z; 
	m6(22, 9) = 4*dxyz2; m6(22,10) = 3*dx2z2; m6(22,11) = -8*dx2yz; 
	m6(23, 9) = 6*dxy2z; m6(23,10) = 4*dx2yz; m6(23,11) = -9*dx2y2; 
	m6(21,12) = -9*dy2z2; m6(21,13) = -24*dxyz2; m6(21,14) = 8*dxy2z; 
	m6(22,12) = -6*dxyz2; m6(22,13) = -9*dx2z2; m6(22,14) = 4*dx2yz; 
	m6(23,12) = -4*dxy2z; m6(23,13) = -8*dx2yz; m6(23,14) = 3*dx2y2; 
	m6(21,15) = 9*dy2z2; m6(21,16) = -12*dxyz2; m6(21,17) = 4*dxy2z; 
	m6(22,15) = 6*dxyz2; m6(22,16) = -9*dx2z2; m6(22,17) = 4*dx2yz; 
	m6(23,15) = 4*dxy2z; m6(23,16) = -8*dx2yz; m6(23,17) = 3*dx2y2; 
	m6(21,18) = -27*dy2z2; m6(21,19) = 24*dxyz2; m6(21,20) = 24*dxy2z; 
	m6(22,18) = -12*dxyz2; m6(22,19) = 9*dx2z2; m6(22,20) = 8*dx2yz; 
	m6(23,18) = -12*dxy2z; m6(23,19) = 8*dx2yz; m6(23,20) = 9*dx2y2; 
	m6(21,21) = 27*dy2z2; 
	m6(22,21) = 12*dxyz2; m6(22,22) = 9*dx2z2; 
	m6(23,21) = 12*dxy2z; m6(23,22) = 8*dx2yz; m6(23,23) = 9*dx2y2; 
	m7( 0, 0) = 3*dy2z2; 
	m7( 1, 0) = 2*dxyz2; m7( 1, 1) = 3*dx2z2; 
	m7( 2, 0) = 2*dxy2z; m7( 2, 1) = 2*dx2yz; m7( 2, 2) = 3*dx2y2; 
	m7( 3, 0) = -3*dy2z2; m7( 3, 1) = -2*dxyz2; m7( 3, 2) = -2*dxy2z; 
	m7( 4, 0) = 4*dxyz2; m7( 4, 1) = 3*dx2z2; m7( 4, 2) = 2*dx2yz; 
	m7( 5, 0) = 4*dxy2z; m7( 5, 1) = 2*dx2yz; m7( 5, 2) = 3*dx2y2; 
	m7( 3, 3) = 3*dy2z2; 
	m7( 4, 3) = -4*dxyz2; m7( 4, 4) = 9*dx2z2; 
	m7( 5, 3) = -4*dxy2z; m7( 5, 4) = 6*dx2yz; m7( 5, 5) = 9*dx2y2; 
	m7( 6, 0) = 3*dy2z2; m7( 6, 1) = 4*dxyz2; m7( 6, 2) = 2*dxy2z; 
	m7( 7, 0) = -2*dxyz2; m7( 7, 1) = -3*dx2z2; m7( 7, 2) = -2*dx2yz; 
	m7( 8, 0) = 2*dxy2z; m7( 8, 1) = 4*dx2yz; m7( 8, 2) = 3*dx2y2; 
	m7( 6, 3) = -3*dy2z2; m7( 6, 4) = 8*dxyz2; m7( 6, 5) = 4*dxy2z; 
	m7( 7, 3) = 2*dxyz2; m7( 7, 4) = -3*dx2z2; m7( 7, 5) = -2*dx2yz; 
	m7( 8, 3) = -2*dxy2z; m7( 8, 4) = 4*dx2yz; m7( 8, 5) = 3*dx2y2; 
	m7( 6, 6) = 9*dy2z2; 
	m7( 7, 6) = -4*dxyz2; m7( 7, 7) = 3*dx2z2; 
	m7( 8, 6) = 6*dxy2z; m7( 8, 7) = -4*dx2yz; m7( 8, 8) = 9*dx2y2; 
	m7( 9, 0) = -3*dy2z2; m7( 9, 1) = -4*dxyz2; m7( 9, 2) = -2*dxy2z; 
	m7(10, 0) = -4*dxyz2; m7(10, 1) = -3*dx2z2; m7(10, 2) = -2*dx2yz; 
	m7(11, 0) = 4*dxy2z; m7(11, 1) = 4*dx2yz; m7(11, 2) = 3*dx2y2; 
	m7( 9, 3) = 3*dy2z2; m7( 9, 4) = -8*dxyz2; m7( 9, 5) = -4*dxy2z; 
	m7(10, 3) = 4*dxyz2; m7(10, 4) = -9*dx2z2; m7(10, 5) = -6*dx2yz; 
	m7(11, 3) = -4*dxy2z; m7(11, 4) = 12*dx2yz; m7(11, 5) = 9*dx2y2; 
	m7( 9, 6) = -9*dy2z2; m7( 9, 7) = 4*dxyz2; m7( 9, 8) = -6*dxy2z; 
	m7(10, 6) = -8*dxyz2; m7(10, 7) = 3*dx2z2; m7(10, 8) = -4*dx2yz; 
	m7(11, 6) = 12*dxy2z; m7(11, 7) = -4*dx2yz; m7(11, 8) = 9*dx2y2; 
	m7( 9, 9) = 9*dy2z2; 
	m7(10, 9) = 8*dxyz2; m7(10,10) = 9*dx2z2; 
	m7(11, 9) = -12*dxy2z; m7(11,10) = -12*dx2yz; m7(11,11) = 27*dx2y2; 
	m7(12, 0) = 3*dy2z2; m7(12, 1) = 2*dxyz2; m7(12, 2) = 4*dxy2z; 
	m7(13, 0) = 2*dxyz2; m7(13, 1) = 3*dx2z2; m7(13, 2) = 4*dx2yz; 
	m7(14, 0) = -2*dxy2z; m7(14, 1) = -2*dx2yz; m7(14, 2) = -3*dx2y2; 
	m7(12, 3) = -3*dy2z2; m7(12, 4) = 4*dxyz2; m7(12, 5) = 8*dxy2z; 
	m7(13, 3) = -2*dxyz2; m7(13, 4) = 3*dx2z2; m7(13, 5) = 4*dx2yz; 
	m7(14, 3) = 2*dxy2z; m7(14, 4) = -2*dx2yz; m7(14, 5) = -3*dx2y2; 
	m7(12, 6) = 3*dy2z2; m7(12, 7) = -2*dxyz2; m7(12, 8) = 4*dxy2z; 
	m7(13, 6) = 4*dxyz2; m7(13, 7) = -3*dx2z2; m7(13, 8) = 8*dx2yz; 
	m7(14, 6) = -2*dxy2z; m7(14, 7) = 2*dx2yz; m7(14, 8) = -3*dx2y2; 
	m7(12, 9) = -3*dy2z2; m7(12,10) = -4*dxyz2; m7(12,11) = 8*dxy2z; 
	m7(13, 9) = -4*dxyz2; m7(13,10) = -3*dx2z2; m7(13,11) = 8*dx2yz; 
	m7(14, 9) = 2*dxy2z; m7(14,10) = 2*dx2yz; m7(14,11) = -3*dx2y2; 
	m7(12,12) = 9*dy2z2; 
	m7(13,12) = 6*dxyz2; m7(13,13) = 9*dx2z2; 
	m7(14,12) = -4*dxy2z; m7(14,13) = -4*dx2yz; m7(14,14) = 3*dx2y2; 
	m7(15, 0) = -3*dy2z2; m7(15, 1) = -2*dxyz2; m7(15, 2) = -4*dxy2z; 
	m7(16, 0) = 4*dxyz2; m7(16, 1) = 3*dx2z2; m7(16, 2) = 4*dx2yz; 
	m7(17, 0) = -4*dxy2z; m7(17, 1) = -2*dx2yz; m7(17, 2) = -3*dx2y2; 
	m7(15, 3) = 3*dy2z2; m7(15, 4) = -4*dxyz2; m7(15, 5) = -8*dxy2z; 
	m7(16, 3) = -4*dxyz2; m7(16, 4) = 9*dx2z2; m7(16, 5) = 12*dx2yz; 
	m7(17, 3) = 4*dxy2z; m7(17, 4) = -6*dx2yz; m7(17, 5) = -9*dx2y2; 
	m7(15, 6) = -3*dy2z2; m7(15, 7) = 2*dxyz2; m7(15, 8) = -4*dxy2z; 
	m7(16, 6) = 8*dxyz2; m7(16, 7) = -3*dx2z2; m7(16, 8) = 8*dx2yz; 
	m7(17, 6) = -4*dxy2z; m7(17, 7) = 2*dx2yz; m7(17, 8) = -3*dx2y2; 
	m7(15, 9) = 3*dy2z2; m7(15,10) = 4*dxyz2; m7(15,11) = -8*dxy2z; 
	m7(16, 9) = -8*dxyz2; m7(16,10) = -9*dx2z2; m7(16,11) = 24*dx2yz; 
	m7(17, 9) = 4*dxy2z; m7(17,10) = 6*dx2yz; m7(17,11) = -9*dx2y2; 
	m7(15,12) = -9*dy2z2; m7(15,13) = -6*dxyz2; m7(15,14) = 4*dxy2z; 
	m7(16,12) = 12*dxyz2; m7(16,13) = 9*dx2z2; m7(16,14) = -4*dx2yz; 
	m7(17,12) = -8*dxy2z; m7(17,13) = -4*dx2yz; m7(17,14) = 3*dx2y2; 
	m7(15,15) = 9*dy2z2; 
	m7(16,15) = -12*dxyz2; m7(16,16) = 27*dx2z2; 
	m7(17,15) = 8*dxy2z; m7(17,16) = -12*dx2yz; m7(17,17) = 9*dx2y2; 
	m7(18, 0) = 3*dy2z2; m7(18, 1) = 4*dxyz2; m7(18, 2) = 4*dxy2z; 
	m7(19, 0) = -2*dxyz2; m7(19, 1) = -3*dx2z2; m7(19, 2) = -4*dx2yz; 
	m7(20, 0) = -2*dxy2z; m7(20, 1) = -4*dx2yz; m7(20, 2) = -3*dx2y2; 
	m7(18, 3) = -3*dy2z2; m7(18, 4) = 8*dxyz2; m7(18, 5) = 8*dxy2z; 
	m7(19, 3) = 2*dxyz2; m7(19, 4) = -3*dx2z2; m7(19, 5) = -4*dx2yz; 
	m7(20, 3) = 2*dxy2z; m7(20, 4) = -4*dx2yz; m7(20, 5) = -3*dx2y2; 
	m7(18, 6) = 9*dy2z2; m7(18, 7) = -4*dxyz2; m7(18, 8) = 12*dxy2z; 
	m7(19, 6) = -4*dxyz2; m7(19, 7) = 3*dx2z2; m7(19, 8) = -8*dx2yz; 
	m7(20, 6) = -6*dxy2z; m7(20, 7) = 4*dx2yz; m7(20, 8) = -9*dx2y2; 
	m7(18, 9) = -9*dy2z2; m7(18,10) = -8*dxyz2; m7(18,11) = 24*dxy2z; 
	m7(19, 9) = 4*dxyz2; m7(19,10) = 3*dx2z2; m7(19,11) = -8*dx2yz; 
	m7(20, 9) = 6*dxy2z; m7(20,10) = 4*dx2yz; m7(20,11) = -9*dx2y2; 
	m7(18,12) = 9*dy2z2; m7(18,13) = 12*dxyz2; m7(18,14) = -4*dxy2z; 
	m7(19,12) = -6*dxyz2; m7(19,13) = -9*dx2z2; m7(19,14) = 4*dx2yz; 
	m7(20,12) = -4*dxy2z; m7(20,13) = -8*dx2yz; m7(20,14) = 3*dx2y2; 
	m7(18,15) = -9*dy2z2; m7(18,16) = 24*dxyz2; m7(18,17) = -8*dxy2z; 
	m7(19,15) = 6*dxyz2; m7(19,16) = -9*dx2z2; m7(19,17) = 4*dx2yz; 
	m7(20,15) = 4*dxy2z; m7(20,16) = -8*dx2yz; m7(20,17) = 3*dx2y2; 
	m7(18,18) = 27*dy2z2; 
	m7(19,18) = -12*dxyz2; m7(19,19) = 9*dx2z2; 
	m7(20,18) = -12*dxy2z; m7(20,19) = 8*dx2yz; m7(20,20) = 9*dx2y2; 
	m7(21, 0) = -3*dy2z2; m7(21, 1) = -4*dxyz2; m7(21, 2) = -4*dxy2z; 
	m7(22, 0) = -4*dxyz2; m7(22, 1) = -3*dx2z2; m7(22, 2) = -4*dx2yz; 
	m7(23, 0) = -4*dxy2z; m7(23, 1) = -4*dx2yz; m7(23, 2) = -3*dx2y2; 
	m7(21, 3) = 3*dy2z2; m7(21, 4) = -8*dxyz2; m7(21, 5) = -8*dxy2z; 
	m7(22, 3) = 4*dxyz2; m7(22, 4) = -9*dx2z2; m7(22, 5) = -12*dx2yz; 
	m7(23, 3) = 4*dxy2z; m7(23, 4) = -12*dx2yz; m7(23, 5) = -9*dx2y2; 
	m7(21, 6) = -9*dy2z2; m7(21, 7) = 4*dxyz2; m7(21, 8) = -12*dxy2z; 
	m7(22, 6) = -8*dxyz2; m7(22, 7) = 3*dx2z2; m7(22, 8) = -8*dx2yz; 
	m7(23, 6) = -12*dxy2z; m7(23, 7) = 4*dx2yz; m7(23, 8) = -9*dx2y2; 
	m7(21, 9) = 9*dy2z2; m7(21,10) = 8*dxyz2; m7(21,11) = -24*dxy2z; 
	m7(22, 9) = 8*dxyz2; m7(22,10) = 9*dx2z2; m7(22,11) = -24*dx2yz; 
	m7(23, 9) = 12*dxy2z; m7(23,10) = 12*dx2yz; m7(23,11) = -27*dx2y2; 
	m7(21,12) = -9*dy2z2; m7(21,13) = -12*dxyz2; m7(21,14) = 4*dxy2z; 
	m7(22,12) = -12*dxyz2; m7(22,13) = -9*dx2z2; m7(22,14) = 4*dx2yz; 
	m7(23,12) = -8*dxy2z; m7(23,13) = -8*dx2yz; m7(23,14) = 3*dx2y2; 
	m7(21,15) = 9*dy2z2; m7(21,16) = -24*dxyz2; m7(21,17) = 8*dxy2z; 
	m7(22,15) = 12*dxyz2; m7(22,16) = -27*dx2z2; m7(22,17) = 12*dx2yz; 
	m7(23,15) = 8*dxy2z; m7(23,16) = -24*dx2yz; m7(23,17) = 9*dx2y2; 
	m7(21,18) = -27*dy2z2; m7(21,19) = 12*dxyz2; m7(21,20) = 12*dxy2z; 
	m7(22,18) = -24*dxyz2; m7(22,19) = 9*dx2z2; m7(22,20) = 8*dx2yz; 
	m7(23,18) = -24*dxy2z; m7(23,19) = 8*dx2yz; m7(23,20) = 9*dx2y2; 
	m7(21,21) = 27*dy2z2; 
	m7(22,21) = 24*dxyz2; m7(22,22) = 27*dx2z2; 
	m7(23,21) = 24*dxy2z; m7(23,22) = 24*dx2yz; m7(23,23) = 27*dx2y2; 
	for (z = 0; z < 8; z++)
	    fidd[z] /= 864.0*dx*dy*dz;
    }
    return fidd[i](j*3+l, k*3+m);
}

double Voxel8::IntPdd (const RVector &P, int j, int k, int l, int m) const
{
    double val = 0;
    for (int i = 0; i < 8; i++)
	val += IntFdd (i,j,k,l,m) * P[Node[i]];
    return val;
}

double Voxel8::IntFfd (int i, int j, int k, int l) const
{
static RVector coeff(864,"-18 -18 -18 18 -6 -6 -6 18 -6 6 6 -2 -6 -6 18 6 -2 6 -2 6 6 2 2 2 -9 -6 -6 9 -6 -6 -3 6 -2 3 6 -2 -3 -2 6 3 -2 6 -1 2 2 1 2 2 -6 -9 -6 6 -3 -2 -6 9 -6 6 3 -2 -2 -3 6 2 -1 2 -2 3 6 2 1 2 -3 -3 -2 3 -3 -2 -3 3 -2 3 3 -2 -1 -1 2 1 -1 2 -1 1 2 1 1 2 -6 -6 -9 6 -2 -3 -2 6 -3 2 2 -1 -6 -6 9 6 -2 3 -2 6 3 2 2 1 -3 -2 -3 3 -2 -3 -1 2 -1 1 2 -1 -3 -2 3 3 -2 3 -1 2 1 1 2 1 -2 -3 -3 2 -1 -1 -2 3 -3 2 1 -1 -2 -3 3 2 -1 1 -2 3 3 2 1 1 -1 -1 -1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 -18 -6 -6 18 -18 -18 -6 6 -2 6 18 -6 -6 -2 6 6 -6 18 -2 2 2 2 6 6 -3 -3 -2 3 -3 -2 -3 3 -2 3 3 -2 -1 -1 2 1 -1 2 -1 1 2 1 1 2 -6 -3 -2 6 -9 -6 -6 3 -2 6 9 -6 -2 -1 2 2 -3 6 -2 1 2 2 3 6 -3 -2 -3 3 -2 -3 -1 2 -1 1 2 -1 -3 -2 3 3 -2 3 -1 2 1 1 2 1 -6 -2 -3 6 -6 -9 -2 2 -1 2 6 -3 -6 -2 3 6 -6 9 -2 2 1 2 6 3 -1 -1 -1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 -2 -1 -1 2 -3 -3 -2 1 -1 2 3 -3 -2 -1 1 2 -3 3 -2 1 1 2 3 3 -6 -18 -6 6 -6 -2 -18 18 -18 18 6 -6 -2 -6 6 2 -2 2 -6 6 18 6 2 6 -3 -6 -2 3 -6 -2 -9 6 -6 9 6 -6 -1 -2 2 1 -2 2 -3 2 6 3 2 6 -2 -3 -3 2 -1 -1 -2 3 -3 2 1 -1 -2 -3 3 2 -1 1 -2 3 3 2 1 1 -1 -1 -1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 -2 -6 -3 2 -2 -1 -6 6 -9 6 2 -3 -2 -6 3 2 -2 1 -6 6 9 6 2 3 -1 -2 -1 1 -2 -1 -3 2 -3 3 2 -3 -1 -2 1 1 -2 1 -3 2 3 3 2 3 -6 -6 -2 6 -18 -6 -18 6 -6 18 18 -18 -2 -2 2 2 -6 6 -6 2 6 6 6 18 -1 -1 -1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 -2 -1 -1 2 -3 -3 -2 1 -1 2 3 -3 -2 -1 1 2 -3 3 -2 1 1 2 3 3 -1 -2 -1 1 -2 -1 -3 2 -3 3 2 -3 -1 -2 1 1 -2 1 -3 2 3 3 2 3 -2 -2 -1 2 -6 -3 -6 2 -3 6 6 -9 -2 -2 1 2 -6 3 -6 2 3 6 6 9 -6 -6 -18 6 -2 -6 -2 6 -6 2 2 -2 -18 -18 18 18 -6 6 -6 18 6 6 6 2 -3 -2 -6 3 -2 -6 -1 2 -2 1 2 -2 -9 -6 6 9 -6 6 -3 6 2 3 6 2 -2 -3 -6 2 -1 -2 -2 3 -6 2 1 -2 -6 -9 6 6 -3 2 -6 9 6 6 3 2 -1 -1 -2 1 -1 -2 -1 1 -2 1 1 -2 -3 -3 2 3 -3 2 -3 3 2 3 3 2 -6 -2 -6 6 -6 -18 -2 2 -2 2 6 -6 -18 -6 6 18 -18 18 -6 6 2 6 18 6 -1 -1 -2 1 -1 -2 -1 1 -2 1 1 -2 -3 -3 2 3 -3 2 -3 3 2 3 3 2 -2 -1 -2 2 -3 -6 -2 1 -2 2 3 -6 -6 -3 2 6 -9 6 -6 3 2 6 9 6 -2 -6 -6 2 -2 -2 -6 6 -18 6 2 -6 -6 -18 6 6 -6 2 -18 18 18 18 6 6 -1 -2 -2 1 -2 -2 -3 2 -6 3 2 -6 -3 -6 2 3 -6 2 -9 6 6 9 6 6 -2 -2 -2 2 -6 -6 -6 2 -6 6 6 -18 -6 -6 2 6 -18 6 -18 6 6 18 18 18 "); 

    int ii=i; int jj=j; 
    if(i<j) {ii=j; jj=i;}
    switch (l) {
    case 0: return (coeff[(17*jj-jj*jj)*12+(24*(ii-jj))+(k*3)+0]/864.0*dy*dz); //look up correct coeff;
    case 1: return (coeff[(17*jj-jj*jj)*12+(24*(ii-jj))+(k*3)+1]/864.0*dx*dz); //look up correct coeff;
    case 2: return (coeff[(17*jj-jj*jj)*12+(24*(ii-jj))+(k*3)+2]/864.0*dx*dy); //look up correct coeff;
    default: xERROR("Invalid coefficient"); return 0;
    }
}

double Voxel8::IntPfd (const RVector &P, int j, int k, int l) const
{
    double val = 0;
    for (int i = 0; i < 8; i++)
	val += IntFfd(i,j,k,l) * P[Node[i]];
    return val;
}

RDenseMatrix Voxel8::StrainDisplacementMatrix (const Point &glob) const
{
    RDenseMatrix B(6,24);
    RDenseMatrix der = GlobalShapeD (glob);
    for (int i = 0; i < 8; i++) {
        B(0,i*3  ) = der(0,i);
	B(1,i*3+1) = der(1,i);
	B(2,i*3+2) = der(2,i);
	B(3,i*3  ) = der(1,i);
	B(3,i*3+1) = der(0,i);
	B(4,i*3+1) = der(2,i);
	B(4,i*3+2) = der(1,i);
	B(5,i*3  ) = der(2,i);
	B(5,i*3+2) = der(0,i);
    }
    return B;
}

RDenseMatrix Voxel8::ElasticityStiffnessMatrix (double modulus, double pratio)
    const
{
    static RDenseMatrix K(24,24);
    static bool need_setup = true;

    if (need_setup) {
        double dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz;
	double e1 = modulus*(1.0-pratio)/((1.0+pratio)*(1.0-2.0*pratio));
	double e2 = modulus*pratio/((1.0+pratio)*(1.0-2.0*pratio));
	double e3 = modulus/((2.0*(1.0+pratio)));

	K( 0, 0) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 0, 1) = 6*dx*dy*dz2*(e2 + e3);
	K( 0, 2) = 6*dx*dy2*dz*(e2 + e3);
	K( 0, 3) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K( 0, 4) = 6*dx*dy*dz2*(e2 - e3);
	K( 0, 5) = 6*dx*dy2*dz*(e2 - e3);
	K( 0, 6) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 0, 7) = -6*dx*dy*dz2*(e2 - e3);
	K( 0, 8) = 3*dx*dy2*dz*(e2 + e3);
	K( 0, 9) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K( 0,10) = -6*dx*dy*dz2*(e2 + e3);
	K( 0,11) = 3*dx*dy2*dz*(e2 - e3);
	K( 0,12) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K( 0,13) = 3*dx*dy*dz2*(e2 + e3);
	K( 0,14) = -6*dx*dy2*dz*(e2 - e3);
	K( 0,15) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K( 0,16) = 3*dx*dy*dz2*(e2 - e3);
	K( 0,17) = -6*dx*dy2*dz*(e2 + e3);
	K( 0,18) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K( 0,19) = -3*dx*dy*dz2*(e2 - e3);
	K( 0,20) = -3*dx*dy2*dz*(e2 - e3);
	K( 0,21) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 0,22) = -3*dx*dy*dz2*(e2 + e3);
	K( 0,23) = -3*dx*dy2*dz*(e2 + e3);

	K( 1, 0) = 6*dx*dy*dz2*(e2 + e3);
	K( 1, 1) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 1, 2) = 6*dx2*dy*dz*(e2 + e3);
	K( 1, 3) = -6*dx*dy*dz2*(e2 - e3);
	K( 1, 4) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 1, 5) = 3*dx2*dy*dz*(e2 + e3);
	K( 1, 6) = 6*dx*dy*dz2*(e2 - e3);
	K( 1, 7) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K( 1, 8) = 6*dx2*dy*dz*(e2 - e3);
	K( 1, 9) = -6*dx*dy*dz2*(e2 + e3);
	K( 1,10) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K( 1,11) = 3*dx2*dy*dz*(e2 - e3);
	K( 1,12) = 3*dx*dy*dz2*(e2 + e3);
	K( 1,13) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K( 1,14) = -6*dx2*dy*dz*(e2 - e3);
	K( 1,15) = -3*dx*dy*dz2*(e2 - e3);
	K( 1,16) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K( 1,17) = -3*dx2*dy*dz*(e2 - e3);
	K( 1,18) = 3*dx*dy*dz2*(e2 - e3);
	K( 1,19) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K( 1,20) = -6*dx2*dy*dz*(e2 + e3);
	K( 1,21) = -3*dx*dy*dz2*(e2 + e3);
	K( 1,22) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 1,23) = -3*dx2*dy*dz*(e2 + e3);

	K( 2, 0) = 6*dx*dy2*dz*(e2 + e3);
	K( 2, 1) = 6*dx2*dy*dz*(e2 + e3);
	K( 2, 2) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 2, 3) = -6*dx*dy2*dz*(e2 - e3);
	K( 2, 4) = 3*dx2*dy*dz*(e2 + e3);
	K( 2, 5) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 2, 6) = 3*dx*dy2*dz*(e2 + e3);
	K( 2, 7) = -6*dx2*dy*dz*(e2 - e3);
	K( 2, 8) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K( 2, 9) = -3*dx*dy2*dz*(e2 - e3);
	K( 2,10) = -3*dx2*dy*dz*(e2 - e3);
	K( 2,11) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K( 2,12) = 6*dx*dy2*dz*(e2 - e3);
	K( 2,13) = 6*dx2*dy*dz*(e2 - e3);
	K( 2,14) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K( 2,15) = -6*dx*dy2*dz*(e2 + e3);
	K( 2,16) = 3*dx2*dy*dz*(e2 - e3);
	K( 2,17) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K( 2,18) = 3*dx*dy2*dz*(e2 - e3);
	K( 2,19) = -6*dx2*dy*dz*(e2 + e3);
	K( 2,20) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K( 2,21) = -3*dx*dy2*dz*(e2 + e3);
	K( 2,22) = -3*dx2*dy*dz*(e2 + e3);
	K( 2,23) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));

	K( 3, 0) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K( 3, 1) = -6*dx*dy*dz2*(e2 - e3);
	K( 3, 2) = -6*dx*dy2*dz*(e2 - e3);
	K( 3, 3) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 3, 4) = -6*dx*dy*dz2*(e2 + e3);
	K( 3, 5) = -6*dx*dy2*dz*(e2 + e3);
	K( 3, 6) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K( 3, 7) = 6*dx*dy*dz2*(e2 + e3);
	K( 3, 8) = -3*dx*dy2*dz*(e2 - e3);
	K( 3, 9) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 3,10) = 6*dx*dy*dz2*(e2 - e3);
	K( 3,11) = -3*dx*dy2*dz*(e2 + e3);
	K( 3,12) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K( 3,13) = -3*dx*dy*dz2*(e2 - e3);
	K( 3,14) = 6*dx*dy2*dz*(e2 + e3);
	K( 3,15) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K( 3,16) = -3*dx*dy*dz2*(e2 + e3);
	K( 3,17) = 6*dx*dy2*dz*(e2 - e3);
	K( 3,18) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 3,19) = 3*dx*dy*dz2*(e2 + e3);
	K( 3,20) = 3*dx*dy2*dz*(e2 + e3);
	K( 3,21) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K( 3,22) = 3*dx*dy*dz2*(e2 - e3);
	K( 3,23) = 3*dx*dy2*dz*(e2 - e3);

	K( 4, 0) = 6*dx*dy*dz2*(e2 - e3);
	K( 4, 1) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 4, 2) = 3*dx2*dy*dz*(e2 + e3);
	K( 4, 3) = -6*dx*dy*dz2*(e2 + e3);
	K( 4, 4) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 4, 5) = 6*dx2*dy*dz*(e2 + e3);
	K( 4, 6) = 6*dx*dy*dz2*(e2 + e3);
	K( 4, 7) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K( 4, 8) = 3*dx2*dy*dz*(e2 - e3);
	K( 4, 9) = -6*dx*dy*dz2*(e2 - e3);
	K( 4,10) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K( 4,11) = 6*dx2*dy*dz*(e2 - e3);
	K( 4,12) = 3*dx*dy*dz2*(e2 - e3);
	K( 4,13) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K( 4,14) = -3*dx2*dy*dz*(e2 - e3);
	K( 4,15) = -3*dx*dy*dz2*(e2 + e3);
	K( 4,16) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K( 4,17) = -6*dx2*dy*dz*(e2 - e3);
	K( 4,18) = 3*dx*dy*dz2*(e2 + e3);
	K( 4,19) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 4,20) = -3*dx2*dy*dz*(e2 + e3);
	K( 4,21) = -3*dx*dy*dz2*(e2 - e3);
	K( 4,22) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K( 4,23) = -6*dx2*dy*dz*(e2 + e3);

	K( 5, 0) = 6*dx*dy2*dz*(e2 - e3);
	K( 5, 1) = 3*dx2*dy*dz*(e2 + e3);
	K( 5, 2) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 5, 3) = -6*dx*dy2*dz*(e2 + e3);
	K( 5, 4) = 6*dx2*dy*dz*(e2 + e3);
	K( 5, 5) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 5, 6) = 3*dx*dy2*dz*(e2 - e3);
	K( 5, 7) = -3*dx2*dy*dz*(e2 - e3);
	K( 5, 8) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K( 5, 9) = -3*dx*dy2*dz*(e2 + e3);
	K( 5,10) = -6*dx2*dy*dz*(e2 - e3);
	K( 5,11) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K( 5,12) = 6*dx*dy2*dz*(e2 + e3);
	K( 5,13) = 3*dx2*dy*dz*(e2 - e3);
	K( 5,14) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K( 5,15) = 6*dx*dy2*dz*(-e2 + e3);
	K( 5,16) = 6*dx2*dy*dz*(e2 - e3);
	K( 5,17) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K( 5,18) = 3*dx*dy2*dz*(e2 + e3);
	K( 5,19) = -3*dx2*dy*dz*(e2 + e3);
	K( 5,20) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 5,21) = 3*dx*dy2*dz*(-e2 + e3);
	K( 5,22) = -6*dx2*dy*dz*(e2 + e3);
	K( 5,23) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);

	K( 6, 0) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 6, 1) = 6*dx*dy*dz2*(e2 - e3);
	K( 6, 2) = 3*dx*dy2*dz*(e2 + e3);
	K( 6, 3) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K( 6, 4) = 6*dx*dy*dz2*(e2 + e3);
	K( 6, 5) = 3*dx*dy2*dz*(e2 - e3);
	K( 6, 6) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 6, 7) = -6*dx*dy*dz2*(e2 + e3);
	K( 6, 8) = 6*dx*dy2*dz*(e2 + e3);
	K( 6, 9) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K( 6,10) = -6*dx*dy*dz2*(e2 - e3);
	K( 6,11) = 6*dx*dy2*dz*(e2 - e3);
	K( 6,12) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K( 6,13) = 3*dx*dy*dz2*(e2 - e3);
	K( 6,14) = -3*dx*dy2*dz*(e2 - e3);
	K( 6,15) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 6,16) = 3*dx*dy*dz2*(e2 + e3);
	K( 6,17) = -3*dx*dy2*dz*(e2 + e3);
	K( 6,18) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K( 6,19) = -3*dx*dy*dz2*(e2 + e3);
	K( 6,20) = -6*dx*dy2*dz*(e2 - e3);
	K( 6,21) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K( 6,22) = -3*dx*dy*dz2*(e2 - e3);
	K( 6,23) = -6*dx*dy2*dz*(e2 + e3);

	K( 7, 0) = -6*dx*dy*dz2*(e2 - e3);
	K( 7, 1) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K( 7, 2) = -6*dx2*dy*dz*(e2 - e3);
	K( 7, 3) = 6*dx*dy*dz2*(e2 + e3);
	K( 7, 4) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K( 7, 5) = -3*dx2*dy*dz*(e2 - e3);
	K( 7, 6) = -6*dx*dy*dz2*(e2 + e3);
	K( 7, 7) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 7, 8) = -6*dx2*dy*dz*(e2 + e3);
	K( 7, 9) = 6*dx*dy*dz2*(e2 - e3);
	K( 7,10) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 7,11) = -3*dx2*dy*dz*(e2 + e3);
	K( 7,12) = -3*dx*dy*dz2*(e2 - e3);
	K( 7,13) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K( 7,14) = 6*dx2*dy*dz*(e2 + e3);
	K( 7,15) = 3*dx*dy*dz2*(e2 + e3);
	K( 7,16) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K( 7,17) = 3*dx2*dy*dz*(e2 + e3);
	K( 7,18) = -3*dx*dy*dz2*(e2 + e3);
	K( 7,19) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K( 7,20) = 6*dx2*dy*dz*(e2 - e3);
	K( 7,21) = 3*dx*dy*dz2*(e2 - e3);
	K( 7,22) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K( 7,23) = 3*dx2*dy*dz*(e2 - e3);

	K( 8, 0) = 3*dx*dy2*dz*(e2 + e3);
	K( 8, 1) = 6*dx2*dy*dz*(e2 - e3);
	K( 8, 2) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K( 8, 3) = -3*dx*dy2*dz*(e2 - e3);
	K( 8, 4) = 3*dx2*dy*dz*(e2 - e3);
	K( 8, 5) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K( 8, 6) = 6*dx*dy2*dz*(e2 + e3);
	K( 8, 7) = -6*dx2*dy*dz*(e2 + e3);
	K( 8, 8) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 8, 9) = -6*dx*dy2*dz*(e2 - e3);
	K( 8,10) = -3*dx2*dy*dz*(e2 + e3);
	K( 8,11) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 8,12) = 3*dx*dy2*dz*(e2 - e3);
	K( 8,13) = 6*dx2*dy*dz*(e2 + e3);
	K( 8,14) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K( 8,15) = -3*dx*dy2*dz*(e2 + e3);
	K( 8,16) = 3*dx2*dy*dz*(e2 + e3);
	K( 8,17) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K( 8,18) = 6*dx*dy2*dz*(e2 - e3);
	K( 8,19) = 6*dx2*dy*dz*(-e2 + e3);
	K( 8,20) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K( 8,21) = -6*dx*dy2*dz*(e2 + e3);
	K( 8,22) = 3*dx2*dy*dz*(-e2 + e3);
	K( 8,23) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);

	K( 9, 0) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K( 9, 1) = -6*dx*dy*dz2*(e2 + e3);
	K( 9, 2) = -3*dx*dy2*dz*(e2 - e3);
	K( 9, 3) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 9, 4) = -6*dx*dy*dz2*(e2 - e3);
	K( 9, 5) = -3*dx*dy2*dz*(e2 + e3);
	K( 9, 6) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K( 9, 7) = 6*dx*dy*dz2*(e2 - e3);
	K( 9, 8) = -6*dx*dy2*dz*(e2 - e3);
	K( 9, 9) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 9,10) = 6*dx*dy*dz2*(e2 + e3);
	K( 9,11) = -6*dx*dy2*dz*(e2 + e3);
	K( 9,12) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K( 9,13) = -3*dx*dy*dz2*(e2 + e3);
	K( 9,14) = 3*dx*dy2*dz*(e2 + e3);
	K( 9,15) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K( 9,16) = -3*dx*dy*dz2*(e2 - e3);
	K( 9,17) = 3*dx*dy2*dz*(e2 - e3);
	K( 9,18) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K( 9,19) = 3*dx*dy*dz2*(e2 - e3);
	K( 9,20) = 6*dx*dy2*dz*(e2 + e3);
	K( 9,21) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K( 9,22) = 3*dx*dy*dz2*(e2 + e3);
	K( 9,23) = 6*dx*dy2*dz*(e2 - e3);

	K(10, 0) = -6*dx*dy*dz2*(e2 + e3);
	K(10, 1) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K(10, 2) = -3*dx2*dy*dz*(e2 - e3);
	K(10, 3) = 6*dx*dy*dz2*(e2 - e3);
	K(10, 4) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K(10, 5) = -6*dx2*dy*dz*(e2 - e3);
	K(10, 6) = -6*dx*dy*dz2*(e2 - e3);
	K(10, 7) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(10, 8) = -3*dx2*dy*dz*(e2 + e3);
	K(10, 9) = 6*dx*dy*dz2*(e2 + e3);
	K(10,10) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(10,11) = -6*dx2*dy*dz*(e2 + e3);
	K(10,12) = -3*dx*dy*dz2*(e2 + e3);
	K(10,13) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(10,14) = 3*dx2*dy*dz*(e2 + e3);
	K(10,15) = 3*dx*dy*dz2*(e2 - e3);
	K(10,16) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K(10,17) = 6*dx2*dy*dz*(e2 + e3);
	K(10,18) = -3*dx*dy*dz2*(e2 - e3);
	K(10,19) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K(10,20) = 3*dx2*dy*dz*(e2 - e3);
	K(10,21) = 3*dx*dy*dz2*(e2 + e3);
	K(10,22) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K(10,23) = 6*dx2*dy*dz*(e2 - e3);

	K(11, 0) = 3*dx*dy2*dz*(e2 - e3);
	K(11, 1) = 3*dx2*dy*dz*(e2 - e3);
	K(11, 2) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K(11, 3) = -3*dx*dy2*dz*(e2 + e3);
	K(11, 4) = 6*dx2*dy*dz*(e2 - e3);
	K(11, 5) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K(11, 6) = 6*dx*dy2*dz*(e2 - e3);
	K(11, 7) = -3*dx2*dy*dz*(e2 + e3);
	K(11, 8) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(11, 9) = -6*dx*dy2*dz*(e2 + e3);
	K(11,10) = -6*dx2*dy*dz*(e2 + e3);
	K(11,11) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(11,12) = 3*dx*dy2*dz*(e2 + e3);
	K(11,13) = 3*dx2*dy*dz*(e2 + e3);
	K(11,14) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(11,15) = 3*dx*dy2*dz*(-e2 + e3);
	K(11,16) = 6*dx2*dy*dz*(e2 + e3);
	K(11,17) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K(11,18) = 6*dx*dy2*dz*(e2 + e3);
	K(11,19) = 3*dx2*dy*dz*(-e2 + e3);
	K(11,20) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K(11,21) = 6*dx*dy2*dz*(-e2 + e3);
	K(11,22) = 6*dx2*dy*dz*(-e2 + e3);
	K(11,23) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));

	K(12, 0) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K(12, 1) = 3*dx*dy*dz2*(e2 + e3);
	K(12, 2) = 6*dx*dy2*dz*(e2 - e3);
	K(12, 3) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K(12, 4) = 3*dx*dy*dz2*(e2 - e3);
	K(12, 5) = 6*dx*dy2*dz*(e2 + e3);
	K(12, 6) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K(12, 7) = -3*dx*dy*dz2*(e2 - e3);
	K(12, 8) = 3*dx*dy2*dz*(e2 - e3);
	K(12, 9) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(12,10) = -3*dx*dy*dz2*(e2 + e3);
	K(12,11) = 3*dx*dy2*dz*(e2 + e3);
	K(12,12) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(12,13) = 6*dx*dy*dz2*(e2 + e3);
	K(12,14) = -6*dx*dy2*dz*(e2 + e3);
	K(12,15) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K(12,16) = 6*dx*dy*dz2*(e2 - e3);
	K(12,17) = -6*dx*dy2*dz*(e2 - e3);
	K(12,18) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(12,19) = -6*dx*dy*dz2*(e2 - e3);
	K(12,20) = -3*dx*dy2*dz*(e2 + e3);
	K(12,21) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K(12,22) = -6*dx*dy*dz2*(e2 + e3);
	K(12,23) = -3*dx*dy2*dz*(e2 - e3);

	K(13, 0) = 3*dx*dy*dz2*(e2 + e3);
	K(13, 1) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K(13, 2) = 6*dx2*dy*dz*(e2 - e3);
	K(13, 3) = -3*dx*dy*dz2*(e2 - e3);
	K(13, 4) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K(13, 5) = 3*dx2*dy*dz*(e2 - e3);
	K(13, 6) = 3*dx*dy*dz2*(e2 - e3);
	K(13, 7) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K(13, 8) = 6*dx2*dy*dz*(e2 + e3);
	K(13, 9) = -3*dx*dy*dz2*(e2 + e3);
	K(13,10) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(13,11) = 3*dx2*dy*dz*(e2 + e3);
	K(13,12) = 6*dx*dy*dz2*(e2 + e3);
	K(13,13) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(13,14) = -6*dx2*dy*dz*(e2 + e3);
	K(13,15) = -6*dx*dy*dz2*(e2 - e3);
	K(13,16) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(13,17) = -3*dx2*dy*dz*(e2 + e3);
	K(13,18) = 6*dx*dy*dz2*(e2 - e3);
	K(13,19) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K(13,20) = -6*dx2*dy*dz*(e2 - e3);
	K(13,21) = -6*dx*dy*dz2*(e2 + e3);
	K(13,22) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K(13,23) = -3*dx2*dy*dz*(e2 - e3);

	K(14, 0) = -6*dx*dy2*dz*(e2 - e3);
	K(14, 1) = -6*dx2*dy*dz*(e2 - e3);
	K(14, 2) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K(14, 3) = 6*dx*dy2*dz*(e2 + e3);
	K(14, 4) = -3*dx2*dy*dz*(e2 - e3);
	K(14, 5) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K(14, 6) = -3*dx*dy2*dz*(e2 - e3);
	K(14, 7) = 6*dx2*dy*dz*(e2 + e3);
	K(14, 8) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K(14, 9) = 3*dx*dy2*dz*(e2 + e3);
	K(14,10) = 3*dx2*dy*dz*(e2 + e3);
	K(14,11) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(14,12) = -6*dx*dy2*dz*(e2 + e3);
	K(14,13) = -6*dx2*dy*dz*(e2 + e3);
	K(14,14) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(14,15) = 6*dx*dy2*dz*(e2 - e3);
	K(14,16) = -3*dx2*dy*dz*(e2 + e3);
	K(14,17) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(14,18) = -3*dx*dy2*dz*(e2 + e3);
	K(14,19) = 6*dx2*dy*dz*(e2 - e3);
	K(14,20) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K(14,21) = 3*dx*dy2*dz*(e2 - e3);
	K(14,22) = 3*dx2*dy*dz*(e2 - e3);
	K(14,23) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);

	K(15, 0) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K(15, 1) = -3*dx*dy*dz2*(e2 - e3);
	K(15, 2) = -6*dx*dy2*dz*(e2 + e3);
	K(15, 3) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K(15, 4) = -3*dx*dy*dz2*(e2 + e3);
	K(15, 5) = 6*dx*dy2*dz*(-e2 + e3);
	K(15, 6) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(15, 7) = 3*dx*dy*dz2*(e2 + e3);
	K(15, 8) = -3*dx*dy2*dz*(e2 + e3);
	K(15, 9) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K(15,10) = 3*dx*dy*dz2*(e2 - e3);
	K(15,11) = 3*dx*dy2*dz*(-e2 + e3);
	K(15,12) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K(15,13) = -6*dx*dy*dz2*(e2 - e3);
	K(15,14) = 6*dx*dy2*dz*(e2 - e3);
	K(15,15) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(15,16) = -6*dx*dy*dz2*(e2 + e3);
	K(15,17) = 6*dx*dy2*dz*(e2 + e3);
	K(15,18) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K(15,19) = 6*dx*dy*dz2*(e2 + e3);
	K(15,20) = 3*dx*dy2*dz*(e2 - e3);
	K(15,21) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(15,22) = 6*dx*dy*dz2*(e2 - e3);
	K(15,23) = 3*dx*dy2*dz*(e2 + e3);

	K(16, 0) = 3*dx*dy*dz2*(e2 - e3);
	K(16, 1) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K(16, 2) = 3*dx2*dy*dz*(e2 - e3);
	K(16, 3) = -3*dx*dy*dz2*(e2 + e3);
	K(16, 4) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K(16, 5) = 6*dx2*dy*dz*(e2 - e3);
	K(16, 6) = 3*dx*dy*dz2*(e2 + e3);
	K(16, 7) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(16, 8) = 3*dx2*dy*dz*(e2 + e3);
	K(16, 9) = -3*dx*dy*dz2*(e2 - e3);
	K(16,10) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K(16,11) = 6*dx2*dy*dz*(e2 + e3);
	K(16,12) = 6*dx*dy*dz2*(e2 - e3);
	K(16,13) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(16,14) = -3*dx2*dy*dz*(e2 + e3);
	K(16,15) = -6*dx*dy*dz2*(e2 + e3);
	K(16,16) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(16,17) = -6*dx2*dy*dz*(e2 + e3);
	K(16,18) = 6*dx*dy*dz2*(e2 + e3);
	K(16,19) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K(16,20) = -3*dx2*dy*dz*(e2 - e3);
	K(16,21) = -6*dx*dy*dz2*(e2 - e3);
	K(16,22) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K(16,23) = -6*dx2*dy*dz*(e2 - e3);

	K(17, 0) = -6*dx*dy2*dz*(e2 + e3);
	K(17, 1) = -3*dx2*dy*dz*(e2 - e3);
	K(17, 2) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K(17, 3) = 6*dx*dy2*dz*(e2 - e3);
	K(17, 4) = -6*dx2*dy*dz*(e2 - e3);
	K(17, 5) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K(17, 6) = -3*dx*dy2*dz*(e2 + e3);
	K(17, 7) = 3*dx2*dy*dz*(e2 + e3);
	K(17, 8) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(17, 9) = 3*dx*dy2*dz*(e2 - e3);
	K(17,10) = 6*dx2*dy*dz*(e2 + e3);
	K(17,11) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K(17,12) = -6*dx*dy2*dz*(e2 - e3);
	K(17,13) = -3*dx2*dy*dz*(e2 + e3);
	K(17,14) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(17,15) = 6*dx*dy2*dz*(e2 + e3);
	K(17,16) = -6*dx2*dy*dz*(e2 + e3);
	K(17,17) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(17,18) = -3*dx*dy2*dz*(e2 - e3);
	K(17,19) = 3*dx2*dy*dz*(e2 - e3);
	K(17,20) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K(17,21) = 3*dx*dy2*dz*(e2 + e3);
	K(17,22) = 6*dx2*dy*dz*(e2 - e3);
	K(17,23) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));

	K(18, 0) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K(18, 1) = 3*dx*dy*dz2*(e2 - e3);
	K(18, 2) = 3*dx*dy2*dz*(e2 - e3);
	K(18, 3) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(18, 4) = 3*dx*dy*dz2*(e2 + e3);
	K(18, 5) = 3*dx*dy2*dz*(e2 + e3);
	K(18, 6) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K(18, 7) = -3*dx*dy*dz2*(e2 + e3);
	K(18, 8) = 6*dx*dy2*dz*(e2 - e3);
	K(18, 9) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K(18,10) = -3*dx*dy*dz2*(e2 - e3);
	K(18,11) = 6*dx*dy2*dz*(e2 + e3);
	K(18,12) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(18,13) = 6*dx*dy*dz2*(e2 - e3);
	K(18,14) = -3*dx*dy2*dz*(e2 + e3);
	K(18,15) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K(18,16) = 6*dx*dy*dz2*(e2 + e3);
	K(18,17) = -3*dx*dy2*dz*(e2 - e3);
	K(18,18) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(18,19) = -6*dx*dy*dz2*(e2 + e3);
	K(18,20) = -6*dx*dy2*dz*(e2 + e3);
	K(18,21) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K(18,22) = -6*dx*dy*dz2*(e2 - e3);
	K(18,23) = -6*dx*dy2*dz*(e2 - e3);

	K(19, 0) = -3*dx*dy*dz2*(e2 - e3);
	K(19, 1) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K(19, 2) = -6*dx2*dy*dz*(e2 + e3);
	K(19, 3) = 3*dx*dy*dz2*(e2 + e3);
	K(19, 4) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(19, 5) = -3*dx2*dy*dz*(e2 + e3);
	K(19, 6) = -3*dx*dy*dz2*(e2 + e3);
	K(19, 7) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K(19, 8) = 6*dx2*dy*dz*(-e2 + e3);
	K(19, 9) = 3*dx*dy*dz2*(e2 - e3);
	K(19,10) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K(19,11) = 3*dx2*dy*dz*(-e2 + e3);
	K(19,12) = -6*dx*dy*dz2*(e2 - e3);
	K(19,13) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K(19,14) = 6*dx2*dy*dz*(e2 - e3);
	K(19,15) = 6*dx*dy*dz2*(e2 + e3);
	K(19,16) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K(19,17) = 3*dx2*dy*dz*(e2 - e3);
	K(19,18) = -6*dx*dy*dz2*(e2 + e3);
	K(19,19) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(19,20) = 6*dx2*dy*dz*(e2 + e3);
	K(19,21) = 6*dx*dy*dz2*(e2 - e3);
	K(19,22) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(19,23) = 3*dx2*dy*dz*(e2 + e3);

	K(20, 0) = -3*dx*dy2*dz*(e2 - e3);
	K(20, 1) = -6*dx2*dy*dz*(e2 + e3);
	K(20, 2) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K(20, 3) = 3*dx*dy2*dz*(e2 + e3);
	K(20, 4) = -3*dx2*dy*dz*(e2 + e3);
	K(20, 5) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(20, 6) = -6*dx*dy2*dz*(e2 - e3);
	K(20, 7) = 6*dx2*dy*dz*(e2 - e3);
	K(20, 8) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K(20, 9) = 6*dx*dy2*dz*(e2 + e3);
	K(20,10) = 3*dx2*dy*dz*(e2 - e3);
	K(20,11) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K(20,12) = -3*dx*dy2*dz*(e2 + e3);
	K(20,13) = -6*dx2*dy*dz*(e2 - e3);
	K(20,14) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K(20,15) = 3*dx*dy2*dz*(e2 - e3);
	K(20,16) = -3*dx2*dy*dz*(e2 - e3);
	K(20,17) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K(20,18) = -6*dx*dy2*dz*(e2 + e3);
	K(20,19) = 6*dx2*dy*dz*(e2 + e3);
	K(20,20) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(20,21) = 6*dx*dy2*dz*(e2 - e3);
	K(20,22) = 3*dx2*dy*dz*(e2 + e3);
	K(20,23) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));

	K(21, 0) = -2*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(21, 1) = -3*dx*dy*dz2*(e2 + e3);
	K(21, 2) = -3*dx*dy2*dz*(e2 + e3);
	K(21, 3) = -4*dx2*dz2*e3 + 2*dy2*(dz2*e1 - 2*dx2*e3);
	K(21, 4) = -3*dx*dy*dz2*(e2 - e3);
	K(21, 5) = 3*dx*dy2*dz*(-e2 + e3);
	K(21, 6) = 2*dx2*dz2*e3 - 4*dy2*(dz2*e1 + dx2*e3);
	K(21, 7) = 3*dx*dy*dz2*(e2 - e3);
	K(21, 8) = -6*dx*dy2*dz*(e2 + e3);
	K(21, 9) = 4*(dx2*dz2*e3 + dy2*(dz2*e1 - 2*dx2*e3));
	K(21,10) = 3*dx*dy*dz2*(e2 + e3);
	K(21,11) = 6*dx*dy2*dz*(-e2 + e3);
	K(21,12) = -4*dx2*dz2*e3 + dy2*(-4*dz2*e1 + 2*dx2*e3);
	K(21,13) = -6*dx*dy*dz2*(e2 + e3);
	K(21,14) = 3*dx*dy2*dz*(e2 - e3);
	K(21,15) = 4*(-2*dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(21,16) = -6*dx*dy*dz2*(e2 - e3);
	K(21,17) = 3*dx*dy2*dz*(e2 + e3);
	K(21,18) = 4*(dx2*dz2*e3 + dy2*(-2*dz2*e1 + dx2*e3));
	K(21,19) = 6*dx*dy*dz2*(e2 - e3);
	K(21,20) = 6*dx*dy2*dz*(e2 - e3);
	K(21,21) = 8*(dx2*dz2*e3 + dy2*(dz2*e1 + dx2*e3));
	K(21,22) = 6*dx*dy*dz2*(e2 + e3);
	K(21,23) = 6*dx*dy2*dz*(e2 + e3);

	K(22, 0) = -3*dx*dy*dz2*(e2 + e3);
	K(22, 1) = -2*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(22, 2) = -3*dx2*dy*dz*(e2 + e3);
	K(22, 3) = 3*dx*dy*dz2*(e2 - e3);
	K(22, 4) = 2*dy2*dz2*e3 - 4*dx2*(dz2*e1 + dy2*e3);
	K(22, 5) = -6*dx2*dy*dz*(e2 + e3);
	K(22, 6) = -3*dx*dy*dz2*(e2 - e3);
	K(22, 7) = -4*dy2*dz2*e3 + 2*dx2*(dz2*e1 - 2*dy2*e3);
	K(22, 8) = 3*dx2*dy*dz*(-e2 + e3);
	K(22, 9) = 3*dx*dy*dz2*(e2 + e3);
	K(22,10) = 4*(dy2*dz2*e3 + dx2*(dz2*e1 - 2*dy2*e3));
	K(22,11) = 6*dx2*dy*dz*(-e2 + e3);
	K(22,12) = -6*dx*dy*dz2*(e2 + e3);
	K(22,13) = -4*dy2*dz2*e3 + dx2*(-4*dz2*e1 + 2*dy2*e3);
	K(22,14) = 3*dx2*dy*dz*(e2 - e3);
	K(22,15) = 6*dx*dy*dz2*(e2 - e3);
	K(22,16) = 4*(dy2*dz2*e3 + dx2*(-2*dz2*e1 + dy2*e3));
	K(22,17) = 6*dx2*dy*dz*(e2 - e3);
	K(22,18) = -6*dx*dy*dz2*(e2 - e3);
	K(22,19) = 4*(-2*dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(22,20) = 3*dx2*dy*dz*(e2 + e3);
	K(22,21) = 6*dx*dy*dz2*(e2 + e3);
	K(22,22) = 8*(dy2*dz2*e3 + dx2*(dz2*e1 + dy2*e3));
	K(22,23) = 6*dx2*dy*dz*(e2 + e3);

	K(23, 0) = -3*dx*dy2*dz*(e2 + e3);
	K(23, 1) = -3*dx2*dy*dz*(e2 + e3);
	K(23, 2) = -2*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(23, 3) = 3*dx*dy2*dz*(e2 - e3);
	K(23, 4) = -6*dx2*dy*dz*(e2 + e3);
	K(23, 5) = 2*dy2*dz2*e3 - 4*dx2*(dy2*e1 + dz2*e3);
	K(23, 6) = -6*dx*dy2*dz*(e2 + e3);
	K(23, 7) = 3*dx2*dy*dz*(e2 - e3);
	K(23, 8) = -4*dy2*dz2*e3 + dx2*(-4*dy2*e1 + 2*dz2*e3);
	K(23, 9) = 6*dx*dy2*dz*(e2 - e3);
	K(23,10) = 6*dx2*dy*dz*(e2 - e3);
	K(23,11) = 4*(dy2*dz2*e3 + dx2*(-2*dy2*e1 + dz2*e3));
	K(23,12) = -3*dx*dy2*dz*(e2 - e3);
	K(23,13) = -3*dx2*dy*dz*(e2 - e3);
	K(23,14) = -4*dy2*dz2*e3 + 2*dx2*(dy2*e1 - 2*dz2*e3);
	K(23,15) = 3*dx*dy2*dz*(e2 + e3);
	K(23,16) = -6*dx2*dy*dz*(e2 - e3);
	K(23,17) = 4*(dy2*dz2*e3 + dx2*(dy2*e1 - 2*dz2*e3));
	K(23,18) = -6*dx*dy2*dz*(e2 - e3);
	K(23,19) = 3*dx2*dy*dz*(e2 + e3);
	K(23,20) = 4*(-2*dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));
	K(23,21) = 6*dx*dy2*dz*(e2 + e3);
	K(23,22) = 6*dx2*dy*dz*(e2 + e3);
	K(23,23) = 8*(dy2*dz2*e3 + dx2*(dy2*e1 + dz2*e3));

        K /= 72.0*dx*dy*dz;
	need_setup = false;
    }
    return K;
}

int Voxel8::Intersection (const Point &p1, const Point &p2, Point *s,
    bool add_endpoints, bool boundary_only)
{
    dASSERT(p1.Dim() == 3 && p2.Dim() == 3, "Points must be 3D.");

    const double EPS = 1e-12;
    Point d = p2-p1;
    double a, rx, ry, rz;
    Point p(3);
    int n = 0;
    
    // check surface z=0
    if (d[2] && (!boundary_only || bndside[0])) {
	a = -p1[2]/d[2];
	rx = p1[0] + a*d[0];
	ry = p1[1] + a*d[1];
	if (rx > -EPS && rx < 1.0+EPS && ry > -EPS && ry < 1.0+EPS) {
	    p[0] = rx, p[1] = ry, p[2] = 0.0;
	    s[n++] = p;
	}
    }

    // check surface z=1
    if (d[2] && (!boundary_only || bndside[1])) {
	a = (1.0-p1[2])/d[2];
	rx = p1[0] + a*d[0];
	ry = p1[1] + a*d[1];
	if (rx > -EPS && rx < 1.0+EPS && ry > -EPS && ry < 1.0+EPS) {
	    p[0] = rx, p[1] = ry, p[2] = 0.0;
	    s[n++] = p;
	}
    }
    
    // check surface y=0
    if (d[1] && (!boundary_only || bndside[2])) {
	a = -p1[1]/d[1];
	rx = p1[0] + a*d[0];
	rz = p1[2] + a*d[2];
	if (rx > -EPS && rx < 1.0+EPS && rz > -EPS && rz < 1.0+EPS) {
	    p[0] = rx, p[1] = 0.0, p[2] = rz;
	    s[n++] = p;
	}
    }

    // check surface y=1
    if (d[1] && (!boundary_only || bndside[3])) {
	a = (1.0-p1[1])/d[1];
	rx = p1[0] + a*d[0];
	rz = p1[2] + a*d[2];
	if (rx > -EPS && rx < 1.0+EPS && rz > -EPS && rz < 1.0+EPS) {
	    p[0] = rx, p[1] = 0.0, p[2] = rz;
	    s[n++] = p;
	}	
    }

    // check surface x=0
    if (d[0] && (!boundary_only || bndside[4])) {
	a = -p1[0]/d[0];
	ry = p1[1] + a*d[1];
	rz = p1[2] + a*d[2];
	if (ry > -EPS && ry < 1.0+EPS && rz > -EPS && rz < 1.0+EPS) {
	    p[0] = 0.0, p[1] = ry, p[2] = rz;
	    s[n++] = p;
	}
    }

    // check surface x=1
    if (d[0] && (!boundary_only || bndside[5])) {
	a = (1.0-p1[0])/d[0];
	ry = p1[1] + a*d[1];
	rz = p1[2] + a*d[2];
	if (ry > -EPS && ry < 1.0+EPS && rz > -EPS && rz < 1.0+EPS) {
	    p[0] = 0.0, p[1] = ry, p[2] = rz;
	    s[n++] = p;
	}
    }
    return n;
}

double Voxel8::dx = 0.0;
double Voxel8::dy = 0.0;
double Voxel8::dz = 0.0;
double Voxel8::size = 0.0;
double Voxel8::intf = 0.0;
RSymMatrix Voxel8::intff;
RSymMatrix Voxel8::intfff[8];
RSymMatrix Voxel8::intdd;
RSymMatrix Voxel8::intfdd[8];
RVector Voxel8::bndintf[6];
RSymMatrix Voxel8::bndintff[6];
RDenseMatrix Voxel8::bndintfff[6][8];
