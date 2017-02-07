// ==========================================================================
// Module libfe
// File element.cc
// Definition of class Element
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "felib.h"

using namespace std;

Element::Element ()
{
    Node = 0;
    bndside = 0;
    sdnbhr = 0;
    sdnbhridx = 0;
    region = -1;
    subdivdata = 0;
    bndel = false;
    interfaceel = false;
}

Element::Element (const Element& el)   // copy constructor
{
    Node = 0;	// this must be done by derived class
    bndside = 0;
    sdnbhr = 0;
    sdnbhridx = 0;
    region = -1;
    subdivdata = 0;
    bndel = false;
    interfaceel = false;
}

Element::~Element ()
{
    if (bndside) delete []bndside;
    if (sdnbhr)  delete []sdnbhr;
    if (sdnbhridx) delete []sdnbhridx;
    if (subdivdata) delete subdivdata;
}

void Element::Initialise (const NodeList& nlist)
{
    // set boundary side flags
    //*** logic changed 17-09-06 to include Internal interfaces
    //*** this is open to improvement, certainly 
    if (!bndside) bndside = new bool[nSide()];
    bndel   = false;
    interfaceel = false;
    for (int side = 0; side < nSide(); side++) {
	bndside[side] = false;
#ifdef UNDEF
        bool internal = false;
	bool boundary = false;
        bndside[side] = true;
        for (int n = 0; n < nSideNode (side); n++) {
	    int nd = Node[SideNode (side, n)]; 
	    
	    bool thisboundary = (nlist[nd].isBnd() ? true : false);
	    bool thisinternal = (nlist[nd].isInternalInterface() ? true:false);
	    // the following logic is terrible! should be exclusive or
	    if (!(thisboundary || thisinternal)){bndside[side] = false; break;}
	    if (boundary && thisinternal || internal && thisboundary)
		{bndside[side] = false; break; }
	    boundary = thisboundary;
	    internal = thisinternal;
	}
	if (bndside[side]) 
	    if(internal) {interfaceel = true;}
	    else  {bndel = true;}
#endif
    }
}

RVector Element::IntF () const
{
    RVector intf(nNode());
    for (int i = 0; i < nNode(); i++)
	intf[i] = IntF(i);
    return intf;
}

RVector Element::SurfIntF (int sd) const
{
    dASSERT(sd >= 0 && sd < nSide(), "Argument 1: out of range");
    RVector sf(nNode());
    for (int i = 0; i < nNode(); i++)
	sf[i] = SurfIntF (i, sd);
    return sf;
}

double Element::BndIntF (int i) const
{
    dASSERT(i >= 0 && i < nNode(), "Argument 1: out of range");
    double bf = 0.0;
    for (int sd = 0; sd < nSide(); sd++)
	if (bndside[sd])
	    bf += SurfIntF (i, sd);
    return bf;
}

RVector Element::BndIntF () const
{
    RVector bf(nNode());
    for (int sd = 0; sd < nSide(); sd++)
	if (bndside[sd])
	    bf += SurfIntF (sd);
    return bf;
}

RSymMatrix Element::SurfIntFF (int sd) const
{
    dASSERT(sd >= 0 && sd < nSide(), "Argument 1: out of range");
    RSymMatrix sff(nNode());
    for (int i = 0; i < nNode(); i++)
	for (int j = 0; j <= i; j++)
	    sff(i,j) = SurfIntFF (i, j, sd);
    return sff;
}

double Element::BndIntFF (int i, int j) const
{
    dASSERT(i >= 0 && i < nNode(), "Argument 1: out of range");
    dASSERT(j >= 0 && j < nNode(), "Argument 2: out of range");
    double bff = 0.0;
    if (bndel) {
	for (int sd = 0; sd < nSide(); sd++)
	    if (bndside[sd])
		bff += SurfIntFF (i, j, sd);
    }
    return bff;
}

RSymMatrix Element::BndIntFF () const
{
    RSymMatrix bff(nNode());
    if (bndel) {
	for (int sd = 0; sd < nSide(); sd++)
	    if (bndside[sd])
		bff += SurfIntFF (sd);
    }
    return bff;
}

void Element::operator= (const Element& el)
{
    dASSERT(Type() == el.Type(), "Assignment of incompatible element types.");
    for (int i = 0; i < nNode(); i++) Node[i] = el.Node[i];
}

void Element::MapToSide (int side, Point &loc) const
{
    ERROR_UNDEF;
}

Point Element::Global (const NodeList &nlist, const Point& loc) const
{
    dASSERT(loc.Dim() == Dimension(), "Wrong point dimension");

    int i, j;
    Point glob(loc.Dim());
    RVector u = LocalShapeF (loc);

    for (i = 0; i < loc.Dim(); i++)
	for (j = 0; j < nNode(); j++)
	    glob[i] += u[j] * nlist[Node[j]][i];
    return glob;
};

bool Element::GContains (const Point& glob, const NodeList& nlist)
    const
{
    return LContains (Local (nlist, glob));
}

bool Element::IsNode (int node)
{
    int i, nnode = nNode();
    for (i = 0; i < nnode; i++) if (Node[i] == node) return true;
    return false;
}

/*
bool Element::IsSide (int node1, int node2)
{
    bool found1 = false, found2 = false;
    for (int i = 0; i < nNode(); i++) {
	if (Node[i] == node1) {
	    if (found2) return true;
	    else found1 = true;
	} else if (Node[i] == node2) {
	    if (found1) return true;
	    else found2 = true;
	}
    }
    return false;
}
*/

int Element::IsSide (int nn, int *nd)
{
    int sd, nsn, i, j, n;
    bool found_node, found_side;
    int nnd = nNode();

    // check whether all nodes in 'nd' belong to the element
    for (i = 0; i < nn; i++) {
        for (found_node = false, j = 0; j < nnd; j++) {
	    if (nd[i] == Node[j]) { found_node = true; break; }
	}
	if (!found_node) return -1;
    }

    // check whether all nodes belong to a single side of the element
    for (sd = 0; sd < nSide(); sd++) {
        nsn = nSideNode (sd);
	for (found_side = true, i = 0; i < nsn; i++) {
	    n = Node[SideNode (sd, i)];
	    for (found_node = false, j = 0; j < nn; j++)
	        if (n == nd[j]) { found_node = true; break; }
	    if (!found_node) { found_side = false; break; }
	}
	if (found_side) return sd;
    }
    return -1;
}

bool Element::IsSideNode (int side, int node)
{
    dASSERT(side >= 0 && side < nSide(), "Invalid value for argument side.");
    for (int i = 0; i < nSideNode (side); i++)
	if (Node[SideNode (side, i)] == node) return true;
    return false;
}

Point Element::SideCentre (int side) const
{
    dASSERT(side >= 0 && side < nSide(), "Invalid value for argument side.");
    Point cnt = NodeLocal (SideNode (side, 0));
    int nnode = nSideNode (side);
    for (int i = 1; i < nnode; i++) cnt += NodeLocal (SideNode (side, i));
    return cnt / (double)nnode;
}

Element *Element::SideNeighbour (int side) const
{
    dASSERT(side >= 0 && side < nSide(), "Invalid value for argument side.");
    dASSERT(sdnbhr, "Neighbour support not initialised");
    return sdnbhr[side];
}

int Element::SideNeighbourIndex (int side) const
{
    dASSERT(side >= 0 && side < nSide(), "Invalid value for argument side.");
    dASSERT(sdnbhridx, "Neighbour support not initialised");
    return sdnbhridx[side];
}

RDenseMatrix Element::Elgeom (const NodeList& nlist) const
{
    RDenseMatrix tmp(nNode(), Dimension());
    for (int i = 0; i < nNode(); i++)
	for (int j = 0; j < Dimension(); j++)
	    tmp(i,j) = nlist[Node[i]][j];
    return tmp;
}

int Element::GetSubsampleFD (int &n, double *&wght, Point *&absc,
    RVector *&F, RDenseMatrix *&D, const NodeList &nlist) const
{
    const Point *absc_loc;
    int ntot = GetLocalSubsampleAbsc (absc_loc);
    if (!ntot) return 0; // subsampling not supported
    double intot = 1.0/(double)ntot;
    if (n < ntot) {
        if (n) {
	    delete []absc;
	    delete []wght;
	    delete []F;
	    delete []D;
	}
	n = ntot;
	absc = new Point[n];
	wght = new double[n];
	F    = new RVector[n];
	D    = new RDenseMatrix[n];
    }
    for (int i = 0; i < ntot; i++) {
        absc[i] = Global (nlist, absc_loc[i]);
	wght[i] = Size()*intot;
        F[i]    = GlobalShapeF (nlist, absc[i]);
	D[i]    = GlobalShapeD (nlist, absc[i]);
    }
    return ntot;
}

int Element::GetBndSubsampleFD (int side, int &n, double *&wght, Point *&absc,
    RVector *&F, RDenseMatrix *&D, const NodeList &nlist) const
{
    const Point *babsc;
    int ntot = GetBndSubsampleAbsc (side, babsc);
    if (!ntot) return 0; // subsampling not supported
    double intot = 1.0/(double)ntot;
    if (n < ntot) {
        if (n) {
	    delete []absc;
	    delete []wght;
	    delete []F;
	    delete []D;
	}
	n = ntot;
	absc = new Point[n];
	wght = new double[n];
	F    = new RVector[n];
	D    = new RDenseMatrix[n];
    }
    for (int i = 0; i < ntot; i++) {
	Point loc = SurfToLocal (side, babsc[i]);
        absc[i]   = Global (nlist, loc);
	wght[i]   = SideSize(side,nlist)*intot;
        F[i]      = LocalShapeF (loc);
	D[i]      = LocalShapeD (loc);
    }
    return ntot;
}

// WARNING: dodgy function: does not differentiate between boundary types;
// assumes wrongly that each side with two boundary nodes is a boundary side

int Element::BndSideList (const NodeList &nlist, int *list)
{
    int i, j, bnd, recno=0;

    for (i = 0; i < nSide(); i++) {
	for (bnd = true, j = 0; j < nSideNode(i); j++)
	    if (!nlist[Node[SideNode(i,j)]].isBnd()) { bnd = false; break; }
	if (bnd) list[recno++] = i;
    }
    return recno;
}

void Element::InitNeighbourSupport ()
{
    if (!sdnbhr) sdnbhr = new Element*[nSide()];
    if (!sdnbhridx) sdnbhridx = new int[nSide()];
    for (int i = 0; i < nSide(); i++) {
	sdnbhr[i] = NULL;
	sdnbhridx[i] = -1;
    }
}

void Element::InitSubdivisionSupport ()
{
    if (!subdivdata) subdivdata = new ElementSubdivisionData;
    subdivdata->level = 0;
    subdivdata->sibling = NULL;
}

int Element::SubdivisionLevel () const
{
    return (subdivdata ? subdivdata->level : 0);
}

RDenseMatrix Element::GlobalShapeD (const NodeList &nlist, const Point &glob)
    const
{
    int d = Dimension();
    Point loc (Local(nlist, glob));
    RDenseMatrix IJ(d,d);
    IJacobian(loc, &nlist, IJ);
    return IJ * LocalShapeD (loc);
}

double Element::Jacobian (const Point &loc, const NodeList *nlist,
    RDenseMatrix &J) const
{
    RDenseMatrix der = LocalShapeD (loc);
    int i, j, k, n, dim = Dimension();
    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++) {
	    double v = 0.0;
	    for (k = 0; k < nNode(); k++) {
		n = Node[k];
		v += der(i,k)*(*nlist)[n][j];
	    }
	    J(i,j) = v;
	}
    }
    if (dim == 2) return J(0,0)*J(1,1) - J(1,0)*J(0,1);
    else          return J(0,0) * (J(1,1)*J(2,2) - J(2,1)*J(1,2)) -
		         J(0,1) * (J(1,0)*J(2,2) - J(2,0)*J(1,2)) +
		         J(0,2) * (J(1,0)*J(2,1) - J(2,0)*J(1,1));
}

double Element::IJacobian (const Point &loc, const NodeList *nlist,
    RDenseMatrix &IJ) const
{
    int dim = Dimension();
    RDenseMatrix J(dim,dim);
    double d = Jacobian (loc, nlist, J);
    double id = 1.0/d;
    if (dim == 2) {
	IJ(0,0) = J(1,1)*id;
	IJ(0,1) = J(0,1)*(-id);
	IJ(1,0) = J(1,0)*(-id);
	IJ(1,1) = J(0,0)*id;
    } else {
	IJ(0,0) = ( J(1,1)*J(2,2) - J(2,1)*J(1,2)) * id;
	IJ(1,0) = (-J(1,0)*J(2,2) + J(2,0)*J(1,2)) * id;
	IJ(2,0) = ( J(1,0)*J(2,1) - J(2,0)*J(1,1)) * id;
	IJ(0,1) = (-J(0,1)*J(2,2) + J(2,1)*J(0,2)) * id;
	IJ(1,1) = ( J(0,0)*J(2,2) - J(2,0)*J(0,2)) * id;
	IJ(2,1) = (-J(0,0)*J(2,1) + J(2,0)*J(0,1)) * id;
	IJ(0,2) = ( J(0,1)*J(1,2) - J(1,1)*J(0,2)) * id;
	IJ(1,2) = (-J(0,0)*J(1,2) + J(1,0)*J(0,2)) * id;
	IJ(2,2) = ( J(0,0)*J(1,1) - J(1,0)*J(0,1)) * id;
    }
    return d;
}

istream &operator>> (istream& i, Element &el)
{
    char cbuf[200];
    int v, n;

    for (n = 0; n < el.nNode(); n++) {
	i >> v;
	el.Node[n] = v-1;
    }
    i.getline (cbuf, 200);       // read to eoln to skip comments
    //    cout << cbuf;
    // check for element region
    el.region = -1;
    for(int k = 0;cbuf[k] !='\n';k++)
    {
      if (cbuf[k] == 'R') {
	  sscanf (cbuf+k+1, "%d", &el.region);
	  //el.region = cbuf[k+1]-'0';
          break;
      }
    };
    return i;
}

ostream &operator<< (ostream& o, const Element &el)
{
    o << (char)(el.Type()-1+'a');
    for (int n = 0; n < el.nNode(); n++) o << " " << el.Node[n]+1;
    if (el.region >= 0) o << " R" << el.region;
    o << endl;
    return o;
}

RVector Element::BndIntFX (int, double (*)(const Point&),
    const NodeList &) const
{
    ERROR_UNDEF;
    return RVector(); // dummy
}

RVector Element::BndIntFCos (int, const Surface*, const RVector&, double,
    double, const NodeList&) const
{
    ERROR_UNDEF;
    return RVector(); // dummy
}

RVector Element::BndIntFDelta (int, const Surface*, const RVector&,
    const NodeList&) const
{
    ERROR_UNDEF;
    return RVector(); // dummy
}

RDenseMatrix Element::ElasticStrainDisplacement (const RVector &loc,
    const RDenseMatrix &gder) const
{
    // Returns strain-displacement matrix for 2D plane elasticity
    // at local position loc, given global shape function derivatives gder
    // See NAG finel routine B2C2
  
    // only implemented for 2D yet

    if (Dimension() == 2) {
        RDenseMatrix B(3,2*nNode());
	for (int i = 0; i < nNode(); i++) {
	    B(0,2*i) = B(2,2*i+1) = gder(0,i);
	    B(1,2*i+1) = B(2,2*i) = gder(1,i);
	}
	return B;
    } else {
        ERROR_UNDEF;
	return RDenseMatrix();
    }
}

int Element::GlobalIntersection (const NodeList &nlist,
    const Point &p1, const Point &p2, Point *s,
    bool add_endpoints, bool boundary_only)
{
    return Intersection (Local(nlist, p1), Local(nlist, p2), s,
        add_endpoints, boundary_only);
}

// ==========================================================================
// class Element_Unstructured

void Element_Unstructured::Initialise (const NodeList &nlist)
{
    Element::Initialise (nlist);
    size = ComputeSize (nlist);
#ifdef FEM_DEBUG
    if (size <= 0)
	LOGOUT("Warning: element size not positive");
#endif
    intdd.Zero (nNode());
    intbff.Zero (nNode());

    // set element bounding box
    bbmin.New (Dimension());
    bbmax.New (Dimension());
    for (int i = 0; i < Dimension(); i++) {
        bbmin[i] = 1e10, bbmax[i] = -1e10;
	for (int j = 0; j < nNode(); j++) {
	    double val = nlist[Node[j]][i];
	    if (val < bbmin[i]) bbmin[i] = val;
	    if (val > bbmax[i]) bbmax[i] = val;
	}
    }
}    

void Element_Unstructured::PostInitialisation (const NodeList &nlist)
{
    intdd = ComputeIntDD (nlist);
    intbff = ComputeBndIntFF (nlist);
    Element::PostInitialisation (nlist);
}

void Element_Unstructured::operator= (const Element_Unstructured &el)
{
    Element::operator= (el);
    size = el.size;
}

bool Element_Unstructured::GContains (const Point& glob, const NodeList& nlist)
    const
{
    const double EPS = 1e-8;
    
    // check bounding box
    for (int i = 0; i < Dimension(); i++)
        if (glob[i] < bbmin[i]-EPS || glob[i] > bbmax[i]+EPS) return false;

    return Element::GContains (glob, nlist);
}

RDenseMatrix Element_Unstructured_3D::IsotropicElasticityMatrix (double E,
    double nu) const
{
    double mu = E/(2.0*(1.0+nu));                 // shear modulus
    double lambda = nu*E/((1.0+nu)*(1.0-2.0*nu)); // Lame modulus

    RDenseMatrix D(6,6);
    D(0,0) = D(1,1) = D(2,2) = lambda + 2.0*mu;
    D(0,1) = D(0,2) = D(1,2) = D(1,0) = D(2,0) = D(2,1) = lambda;
    D(3,3) = D(4,4) = D(5,5) = mu;

    return D;
}
