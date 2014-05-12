// ==========================================================================
// Module libfe
// File wdg18inf.cc
// Definition of class Wedge18inf
// ==========================================================================

#include <math.h>
#include <stdlib.h>
#include <mathlib.h>
#include <string.h>
#include "felib.h"

using namespace std;

void Wedge18inf::Initialise (const NodeList &nlist)
{
    z0 = -1.0; z1 = 0.0; // defaults (just use local coords)
}

RDenseMatrix Wedge18inf::LocalShapeD_Z (double loc) const
{
    RDenseMatrix der(1,3);
    der(0,0) = loc-0.5;
    der(0,1) = -2.0*loc;
    der(0,2) = loc+0.5;
    return der;
}

// Jacobian for z-axis only

double Wedge18inf::JacobianZ (double loc, RDenseMatrix &J) const
{
    RDenseMatrix der = LocalShapeD_Z(loc);
    J(0,0) = der(0,0)*z0 + der(0,1)*z1;
    return J(0,0); // note that J is a 1x1 matrix
}

RSymMatrix Wedge18inf::IntFF () const
{
    return RSymMatrix();
}

istream &Wedge18inf::operator>> (istream &i)
{
    char cbuf[200];
    int v, n;

    for (n = 0; n < 6; n++) {
        i >> v;
	Node[n] = v-1;
    }
    i >> z0 >> z1; // read the z-coordinates of nodes in infinite direction
    i.getline (cbuf, 200); // read to eoln to skip comments
    return i;
}

ostream &Wedge18inf::operator<< (ostream &o) const
{
    o << (char)(Type()-1+'a');
    for (int n = 0; n < nNode(); n++) o << " " << Node[n]+1;
    o << " " << z0 << " " << z1 << endl;
    return o;
}
