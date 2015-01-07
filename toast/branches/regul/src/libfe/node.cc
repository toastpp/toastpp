// ==========================================================================
// Module libfe
// File node.cc
// Definition of class Node
// ==========================================================================

#define FELIB_IMPLEMENTATION
#define __NODE_CC

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

// boundary type id: none, Robin, internal, Dirichlet, -, Neumann, -, any
char BndId[9] = {'N','R','I','B','X','M','?','A','V'};

Node::Node (): Point ()
{
    bndtp = BND_NONE;
    region = -1;
}

Node::Node (int dim, char _bndtp): Point (dim)
{
    bndtp = _bndtp;
    region = -1;
}

Node::Node (const Node &nd): Point (nd)
{
    bndtp  = nd.bndtp;
    region = nd.region;
}

void Node::Copy (const Node &nd)
{
    if (size != nd.Dim()) New (nd.Dim());
    *this = nd;
}

Node &Node::operator= (const Node& nd)
{
    if (size != nd.Dim()) New (nd.Dim());
    Point::operator=(nd);
    bndtp  = nd.bndtp;
    region = nd.region;
    return *this;
}

Node &Node::operator= (const Point &pt)
{
    if (size != pt.Dim()) New (pt.Dim());
    Point::operator=(pt);
    bndtp = 0;
    region = 0;
    return *this;
}

bool Node::operator== (const Node& nd) const
{
    if (size != nd.Dim()) return false;
    for (int i = 0; i < size; i++) if (data[i] != nd.data[i]) return false;
    return true;
}

bool Node::operator!= (const Node& nd) const
{
    return !(*this == nd);
}

FELIB void Swap (Node &n1, Node &n2)
{
    char tmp_bndtp = n1.bndtp;
    n1.bndtp = n2.bndtp;
    n2.bndtp = tmp_bndtp;

    int tmp_region = n1.region;
    n1.region = n2.region;
    n2.region = tmp_region;

    ::Swap ((Point&)n1, (Point&)n2);
}

FELIB double Dist (const Node &n1, const Node &n2)
{
    dASSERT(n1.Dim() == n2.Dim(), "Node dimensions do not match");
    double diff, sum = 0.0;
    for (int i = 0; i < n1.Dim(); i++) {
        diff = n1[i] - n2[i];
	sum += diff*diff;
    }
    return sqrt (sum);
}

#ifdef UNDEF
istream& operator>> (istream& is, Node& nd)
{
    double crd[3];
    char bndid, tp, c;
    int dim, i;

    is >> bndid >> c;
    dASSERT(c == '[', "Parse error reading node");

    dim = 0;
    while (dim < 3 && is >> crd[dim]) dim++;
    is.clear();
    is >> c;
    dASSERT (dim >= 2 && c == ']', "Parse error reading node.");

    nd.New (dim);
    for (i = 0; i < dim; i++) nd[i] = crd[i];

    for (nd.bndtp = tp = 0; tp < 9; tp++)
	if (bndid == BndId[tp]) nd.bndtp = tp;

    is.get(c);
    if (c == 'R') is >> nd.region;
    else nd.region = -1;

    return is;
}
#endif

FELIB istream& operator>> (istream& is, Node& nd)
{
    double crd[3];
    char cbuf[256], bndid, tp, c;
    int dim, i;

    is.getline (cbuf, 256);
    std::istringstream iss (cbuf);

    iss >> bndid >> c;
    dASSERT(c == '[', "Parse error reading node");

    dim = 0;
    while (dim < 3 && iss >> crd[dim]) dim++;
    iss.clear();
    iss >> c;
    dASSERT (dim >= 2 && c == ']', "Parse error reading node.");

    nd.New (dim);
    for (i = 0; i < dim; i++) nd[i] = crd[i];

    for (nd.bndtp = tp = 0; tp < 9; tp++)
	if (bndid == BndId[tp]) nd.bndtp = tp;

    iss.get(c);
    if (c == 'R') iss >> nd.region;
    else nd.region = -1;

    nd.phi = -1.0;
    if (nd.bndtp == BND_DIRICHLET) iss >> nd.phi;

    return is;
}

FELIB ostream& operator<< (ostream& os, Node& nd)
{
    os << BndId[nd.bndtp] << (RVector)nd;
    if (nd.region >= 0) os << 'R' << nd.region;
    return os;
}

