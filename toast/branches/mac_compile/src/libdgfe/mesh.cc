// ==========================================================================
// Module libfe
// File mesh.cc
// Definition of class Mesh
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <string.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;
using namespace toast;

#define MESH_DIRICHLET BND_DIRICHLET
#define MESH_ROBIN     BND_ROBIN
#define MESH_EXTRAPOL  BND_INTERNAL
#define MESH_ODD       BND_NONE

#define MAXNODSD 10	// maximum number of nodes per element side
#define MAXNEIGHBOUR 200 // maximum number of node neighbours

typedef struct tag_bndrec bndrec;
struct tag_bndrec {
    int el, sd;
    int nd[MAXNODSD];
    bndrec *next;
};

Mesh::Mesh ()
{
    // elk = 0; elc = 0; ela = 0;
    nbhrs = 0;  // i.e. undefined
    fullsize = 0.0; // i.e. undefined
    lastel_found = 0;
    priv_nbnd = 0;
    trap_load_error = true;
    IndexBnd2Node = 0;
    IndexNode2Bnd = 0;
    is_set_up = false;
    boundary = 0;
    bnd_param = 0;
}

Mesh::~Mesh ()
{
    //if (elk) delete []elk;
    //if (elc) delete []elc;
    //if (ela) delete []ela;
    if (nbhrs) {
	for (int el = 0; el < elen(); el++) delete []nbhrs[el];
	delete []nbhrs;
    }
    if (IndexBnd2Node) delete []IndexBnd2Node;
    if (IndexNode2Bnd) delete []IndexNode2Bnd;
    if (boundary) delete boundary;
    if (bnd_param) delete []bnd_param;
}

#ifdef TOAST_PARALLEL
void Mesh::Setup_engine (void *arg, int el0, int el1) {
    Mesh *mesh = (Mesh*)arg;
    ElementList &elist = mesh->elist;
    NodeList &nlist = mesh->nlist;
    for (int el = el0; el < el1; el++)
        elist[el]->Initialise (nlist);
}

void Mesh::PostSetup_engine (void *arg, int el0, int el1) {
    Mesh *mesh = (Mesh*)arg;
    ElementList &elist = mesh->elist;
    NodeList &nlist = mesh->nlist;
    for (int el = el0; el < el1; el++)
        elist[el]->PostInitialisation ();
}
#endif

void Mesh::Setup()
{
    int i, j, el;
    priv_ilen = 0;
    for (i = 0; i < nlist.Len(); i++)
	if (nlist[i].BndTp() != BND_DIRICHLET) priv_ilen++;
    priv_nbnd = nlist.NumberOf (BND_ANY);

#ifndef TOAST_PARALLEL
    for (el = 0; el < elist.Len(); el++)
	elist[el]->Initialise (nlist);
    for (el = 0; el < elist.Len(); el++)
	elist[el]->PostInitialisation ();
#else
    int grain = elist.Len()/(8*Task::GetThreadCount());
    g_tpool->ProcessSequence (Mesh::Setup_engine, this, 0, elist.Len(), grain);
    g_tpool->ProcessSequence (Mesh::PostSetup_engine, this, 0, elist.Len(),
			      grain);
#endif // TOAST_PARALLEL

    fullsize = CalcFullSize ();

    // Node <-> Boundary mapping indices
    if (IndexBnd2Node) delete []IndexBnd2Node;
    IndexBnd2Node = new int[priv_nbnd];
    if (IndexNode2Bnd) delete []IndexNode2Bnd;
    IndexNode2Bnd = new int[nlen()];

    for (i = j = 0; i < nlen(); i++) {
        if (nlist[i].isBnd()) {
	    IndexNode2Bnd[i] = j;
	    IndexBnd2Node[j] = i;
	    j++;
	} else {
	    IndexNode2Bnd[i] = -1;
	}
    }

    if (boundary) {
        if (bnd_param) delete []bnd_param;
	bnd_param = new RVector[priv_nbnd];
	for (i = 0; i < priv_nbnd; i++) {
	    bnd_param[i].New (boundary->ParamDim());
	    boundary->Point2Param (nlist[IndexBnd2Node[i]], bnd_param[i]);
	}
    }

    //SetupElementMatrices ();
    is_set_up = true;
}

void Mesh::Copy (const Mesh &mesh)
{
    nlist = mesh.nlist;
    elist = mesh.elist;
    plist = mesh.plist;
    Setup ();
}

RDenseMatrix Mesh::ElGeom (int el) const
{
    dASSERT(el >= 0 && el < elist.Len(), Element index out of range.);
    return elist[el]->Elgeom (nlist);
}

double Mesh::ElSize (int el) const
{
    dASSERT (el >= 0 && el < elist.Len(), Element index out of range.);
    return elist[el]->Size ();
}

Point Mesh::ElCentre (int el) const
{
    int i, nd, dim = nlist[elist[el]->Node[0]].Dim();
    Point cnt(dim);
    for (i = 0; i < elist[el]->nNode(); i++) {
	nd = elist[el]->Node[i];
	dASSERT (nlist[nd].Dim() == dim, Inconsistent node dimensions.);
	cnt += nlist[nd];
    }
    cnt /= elist[el]->nNode();
    return cnt;
}

Point Mesh::ElSideCentre (int el, int sd) const
{
    int i, nd, dim = nlist[elist[el]->Node[0]].Dim();
    Point cnt (dim);
    for (i = 0; i < elist[el]->nSideNode (sd); i++) {
	nd = elist[el]->Node[elist[el]->SideNode (sd, i)];
	dASSERT (nlist[nd].Dim() == dim, Inconsistent node dimensions.);
	cnt += nlist[nd];
    }
    cnt /= elist[el]->nSideNode (sd);
    return cnt;
}

double Mesh::ElSideSize (int el, int sd) const
{
    return elist[el]->SideSize (sd, nlist);
}

RVector Mesh::ElDirectionCosine (int el, int sd, Point *p) const
{
    Point pdummy;
    RDenseMatrix egeom = ElGeom (el);
    if (!p) {
	pdummy.New(Dimension());
	p = &pdummy;
    }
    RDenseMatrix lder  = elist[el]->LocalShapeD (*p);
    RDenseMatrix jacin = inverse (lder * egeom);
    RVector dcos  = elist[el]->DirectionCosine (sd, jacin);
    return dcos;
}

double Mesh::CalcFullSize () const
{
    double size = 0.0;
    for (int el = 0; el < elist.Len(); el++)
	size += elist[el]->Size ();
    return size;
}

double Mesh::FullSize () const
{
    if (fullsize > 0.0) return fullsize;
    else return (fullsize = CalcFullSize());
}

int Mesh::ElFind (const Point &pt) const
{
    int el, d1, d2, i;

    // we start the search from the previously found element, hoping
    // that the new point is close to the old one (and assuming that
    // elements are sorted according to geometric proximity)
    if (elist[lastel_found]->GContains (pt, nlist)) {
	return lastel_found;
    }
    d1 = lastel_found; d2 = elen()-1-lastel_found;
    if (d1 <= d2) {
        for (i = 1; i <= d1; i++) {
	    el = lastel_found-i;
	    if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
	    el = lastel_found+i;
	    if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
	}
	for (el = 2*d1+1; el < elen(); el++) {
	    if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
	}
	return -1;
    } else {
        for (i = 1; i <= d2; i++) {
	    el = lastel_found-i;
	    if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
	    el = lastel_found+i;
	    if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
	}
	for (el = 2*d1-elen(); el >= 0; el--) {
	    if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
	}
	return -1;
    }
    /*
    for (el = lastel_found; el < elen(); el++)
	if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
    for (el = lastel_found-1; el >= 0; el--)
	if (elist[el]->GContains (pt, nlist)) return (lastel_found = el);
    return -1;
    */
}

int Mesh::NdFind (const Point &pt, double &dist) const
{
    int nd, i;
    double d;
    dist = 1e20;
    for (i = 0; i < nlen(); i++)
        if ((d = pt.Dist (nlist[i])) < dist)
	    nd = i, dist = d;
    return nd;
}

double Mesh::ElDist (int el1, int el2) const
{
    return length (ElCentre(el1) - ElCentre(el2));
}

void Mesh::Reorder (int *nindex, int *iindex)
{
  NodeList ntmp;
  ParameterList ptmp;
  int nd, el;

  ntmp.New (nlist.Len()); // reorder nodes
  ptmp.New (plist.Len());

  for (nd=0; nd<nlist.Len(); nd++) {
    ntmp[nd].New (nlist[nindex[nd]].Dim());
    ntmp[nd] = nlist[nindex[nd]];
    ptmp[nd] = plist[nindex[nd]];
  }
  for (nd=0; nd<nlist.Len(); nd++) {
      nlist[nd] = ntmp[nd];
      plist[nd] = ptmp[nd];
  }
  ntmp.New (0);

  for (el=0; el<elist.Len(); el++)
    for (nd=0; nd<elist[el]->nNode(); nd++)
      elist[el]->Node[nd] = iindex[elist[el]->Node[nd]];
}

void Mesh::ScaleMesh (double scale)
{
    for (int i = 0; i < nlen(); i++) nlist[i] *= scale;
    if (boundary) boundary->Scale (scale);
    if (is_set_up) Setup (); // update settings
}

void Mesh::ScaleMesh (const RVector &scale)
{
    int i, j, dim = Dimension();
    dASSERT(scale.Dim() == Dimension(), Wrong vector dimension);
    for (i = 0; i < nlen(); i++) {
        for (j = 0; j < dim; j++) nlist[i][j] *= scale[j];
    }
    if (boundary) boundary->Scale (scale);
    if (is_set_up) Setup (); // update settings
}

double Mesh::ParamAverage (const ParameterType prmtp) const
{
    // instead of just summing up nodal values it would be better to
    // actually integrate

    double sum = 0.0;
    for (int i = 0; i < plist.Len(); i++)
	sum += plist[i].Param (prmtp);
    return sum / (double)plist.Len();
}

int Mesh::MaxNodeDiff (int mode) const
{
    dASSERT(mode == BW_TOTAL || mode == BW_INTERNAL || mode == BW_AUTO,
	Unknown mode id.);
    int nlen = nlist.Len(), elen = elist.Len();
    int el, n, node, nodmin, nodmax, diff = 0;
    if (mode == BW_AUTO)
	mode = (BndType==MESH_ROBIN ? BW_TOTAL : BW_INTERNAL);
    for (el = 0; el < elen; el++) {
	nodmin = nlen, nodmax = 0;
	for (n = 0; n < elist[el]->nNode(); n++) {
	    node = elist[el]->Node[n];
	    if (mode == BW_TOTAL || nlist[node].BndTp() != BND_DIRICHLET) {
		if (node > nodmax) nodmax = node;
		if (node < nodmin) nodmin = node;
	    }
	}
	if (nodmax-nodmin > diff) diff = nodmax-nodmin;
    }
    return diff;
}

Point Mesh::NeighbourBarycentre (int node)
{
    int i, j, k, el, elen = elist.Len();
    int nnb, nb[MAXNEIGHBOUR];

    // find node neighbours
    for (el = nnb = 0; el < elen; el++) {
	Element *pel = elist[el];
	for (i = 0; i < pel->nNode(); i++)
	    if (pel->Node[i] == node) break;
	if (i < pel->nNode()) { // found one
	    for (j = 0; j < pel->nNode(); j++) {
		if (pel->Node[j] == node) continue;
		for (k = 0; k < nnb; k++)
		    if (nb[k] == pel->Node[j]) break; // already in list
		if (k == nnb) nb[nnb++] = pel->Node[j]; // add to list
	    }
	}
    }

    // put node into barycentre of neighbours
    Point bc(Dimension());
    for (i = 0; i < nnb; i++) bc += nlist[nb[i]];
    bc /= nnb;
    return bc;
}

void Mesh::SparseRowStructure (int *&rowptr, int *&colidx, int &nzero) const
{
    // M.S. 1.10.99: Check that this works for higher-order element types
    // (TRI6 and TET10)

    typedef struct {
	int n1, n2;
    } IPair;
    IPair *pair, rra;
    int el, nn, nd1, nd2, n1, n2, i, j, l, ir, indx1, indx2, ri, ci, pi;
    int p = 0, npair = 0;

    // calc number of node pairs
    // this assumes that within an element, each node is neighbour to
    // every other node
    for (el = 0; el < elist.Len(); el++) {
	nn = elist[el]->nNode();
	npair += nn*nn;
    }
    pair = new IPair[npair];

    // collect node pairs
    for (el = 0; el < elist.Len(); el++) {
	nn = elist[el]->nNode();
	for (nd1 = 0; nd1 < nn; nd1++) {
	    n1 = elist[el]->Node[nd1];
	    for (nd2 = 0; nd2 < nn; nd2++) {
		n2 = elist[el]->Node[nd2];
		pair[p].n1 = n1;
		pair[p].n2 = n2;
		dASSERT(p < npair, Something went wrong ...);
		p++;
	    }
	}
    }

    // sort node pairs for n1 (heapsort)
    l = (npair >> 1) + 1;
    ir = npair;
    for (;;) {
	if (l > 1) {
	    l--;
	    rra.n1 = pair[l-1].n1, rra.n2 = pair[l-1].n2;
	} else {
	    rra.n1 = pair[ir-1].n1, rra.n2 = pair[ir-1].n2;
	    pair[ir-1].n1 = pair[0].n1, pair[ir-1].n2 = pair[0].n2;
	    if (--ir == 1) {
		pair[0].n1 = rra.n1, pair[0].n2 = rra.n2;
		break;
	    }
	}
	i = l, j = l << 1;
	while (j <= ir) {
	    if (j<ir && pair[j-1].n1 < pair[j].n1) j++;
	    if (rra.n1 < pair[j-1].n1) {
		pair[i-1].n1 = pair[j-1].n1, pair[i-1].n2 = pair[j-1].n2;
		i = j;
		j <<= 1;
	    } else j = ir+1;
	}
	pair[i-1].n1 = rra.n1, pair[i-1].n2 = rra.n2;
    }

    // sort for n2
    indx1 = indx2 = 0;
    while (indx1 < npair) {
	while (indx2 < npair && pair[indx2].n1 == pair[indx1].n1) indx2++;
	nn = indx2-indx1;
	if (nn == 1) goto done;	// nothing to do
	l = (nn >> 1) + 1;
	ir = nn;
	for (;;) {
	    if (l > 1) {
		l--;
		rra.n1 = pair[l-1+indx1].n1, rra.n2 = pair[l-1+indx1].n2;
	    } else {
		rra.n1 = pair[ir-1+indx1].n1, rra.n2 = pair[ir-1+indx1].n2;
		pair[ir-1+indx1].n1 = pair[indx1].n1;
		pair[ir-1+indx1].n2 = pair[indx1].n2;
		if (--ir == 1) {
		    pair[indx1].n1 = rra.n1, pair[indx1].n2 = rra.n2;
		    break;
		}
	    }
	    i = l, j = l << 1;
	    while (j <= ir) {
		if (j < ir && pair[j-1+indx1].n2 < pair[j+indx1].n2) j++;
		if (rra.n2 < pair[j-1+indx1].n2) {
		    pair[i-1+indx1].n1 = pair[j-1+indx1].n1;
		    pair[i-1+indx1].n2 = pair[j-1+indx1].n2;
		    i = j;
		    j <<= 1;
		} else j = ir+1;
	    }
	    pair[i-1+indx1].n1 = rra.n1, pair[i-1+indx1].n2 = rra.n2;
	}
	done:
	indx1 = indx2;
    }

    // mark duplicates
    indx1 = 0;
    nzero = npair;
    for (i = 1; i < npair; i++) {
	if (pair[i].n1 == pair[indx1].n1 && pair[i].n2 == pair[indx1].n2)
	    pair[i].n1 = -1, nzero--;
	else indx1 = i;
    }
    //if (pair[0].n1 == 0 && pair[0].n2 == 0) pair[0].n1 = -1, nzero--;

    colidx = new int[nzero];
    rowptr = new int[nlen()+1];
    
    for (i = pi = ri = ci = 0; i < npair; i++) {
	if (pair[i].n1 < 0) continue;
	if ((i == 0) || (pair[i].n1 > pair[pi].n1)) {
	    dASSERT(ri < nlen()+1, ri index out of range);
	    rowptr[ri++] = ci;
	}
	pi = i;
	dASSERT(ci < nzero, ci index out of range);
	colidx[ci++] = pair[i].n2;
    }
    dASSERT(ci == nzero, ci index out of sync);
    dASSERT(ri == nlen(), ri index out of sync);
    rowptr[ri] = nzero;
    delete []pair;
}

void Mesh::NeighbourCount (int *plist, int nnode, bool include_self) const
{
    dASSERT(nnode <= nlen(), Paramter out of range);
    int i, nzero, *rowptr, *colidx;
    SparseRowStructure (rowptr, colidx, nzero);

    for (i = 0; i < nnode; i++)
	plist[i] = rowptr[i+1]-rowptr[i];

    if (!include_self)  // subtract diagonal elements
        for (i = 0; i < nnode; i++)
	    plist[i]--;

    delete []rowptr;
    delete []colidx;
}

void Mesh::SysMatrixStructure (int *_nz, int **_row_ind, int **_col_ind)
{
    typedef struct {
	int n1, n2;
    } IPair;
    IPair *pair, rra;
    int el, nd1, nd2, i, j, n, n1, n2, nn, l, ir, indx1, indx2, nterm;
    int p = 0, npair = 0;

    int nnode = nlist.Len();

    // calc number of node pairs
    for (el = 0; el < elist.Len(); el++) {
	for (nterm = 0, n = 0; n < elist[el]->nNode(); n++)
	    if (elist[el]->Node[n] < nnode) nterm++;
	for (n = 1; n < nterm; n++) npair += n;
    }
    pair = new IPair[npair];

    // collect node pairs
    for (el = 0; el < elist.Len(); el++)
	for (nd1 = 0; nd1 < elist[el]->nNode()-1; nd1++) {
	    n1 = elist[el]->Node[nd1];
	    if (n1 >= nnode) continue;
	    for (nd2 = nd1+1; nd2 < elist[el]->nNode(); nd2++) {
		n2 = elist[el]->Node[nd2];
		if (n2 >= nnode) continue;
		if (n1 < n2) pair[p].n1 = n1, pair[p].n2 = n2;
		else         pair[p].n1 = n2, pair[p].n2 = n1;
		dASSERT(p < npair, Something went wrong ...);
		p++;
	    }
	}

    // sort node pairs for n1 (heapsort)
    l = (npair >> 1) + 1;
    ir = npair;
    for (;;) {
	if (l > 1) {
	    l--;
	    rra.n1 = pair[l-1].n1, rra.n2 = pair[l-1].n2;
	} else {
	    rra.n1 = pair[ir-1].n1, rra.n2 = pair[ir-1].n2;
	    pair[ir-1].n1 = pair[0].n1, pair[ir-1].n2 = pair[0].n2;
	    if (--ir == 1) {
		pair[0].n1 = rra.n1, pair[0].n2 = rra.n2;
		break;
	    }
	}
	i = l, j = l << 1;
	while (j <= ir) {
	    if (j<ir && pair[j-1].n1 < pair[j].n1) j++;
	    if (rra.n1 < pair[j-1].n1) {
		pair[i-1].n1 = pair[j-1].n1, pair[i-1].n2 = pair[j-1].n2;
		i = j;
		j <<= 1;
	    } else j = ir+1;
	}
	pair[i-1].n1 = rra.n1, pair[i-1].n2 = rra.n2;
    }

    // sort for n2
    indx1 = indx2 = 0;
    while (indx1 < npair) {
	while (indx2 < npair && pair[indx2].n1 == pair[indx1].n1) indx2++;
	nn = indx2-indx1;
	if (nn == 1) goto done;	// nothing to do
	l = (nn >> 1) + 1;
	ir = nn;
	for (;;) {
	    if (l > 1) {
		l--;
		rra.n1 = pair[l-1+indx1].n1, rra.n2 = pair[l-1+indx1].n2;
	    } else {
		rra.n1 = pair[ir-1+indx1].n1, rra.n2 = pair[ir-1+indx1].n2;
		pair[ir-1+indx1].n1 = pair[indx1].n1;
		pair[ir-1+indx1].n2 = pair[indx1].n2;
		if (--ir == 1) {
		    pair[indx1].n1 = rra.n1, pair[indx1].n2 = rra.n2;
		    break;
		}
	    }
	    i = l, j = l << 1;
	    while (j <= ir) {
		if (j < ir && pair[j-1+indx1].n2 < pair[j+indx1].n2) j++;
		if (rra.n2 < pair[j-1+indx1].n2) {
		    pair[i-1+indx1].n1 = pair[j-1+indx1].n1;
		    pair[i-1+indx1].n2 = pair[j-1+indx1].n2;
		    i = j;
		    j <<= 1;
		} else j = ir+1;
	    }
	    pair[i-1+indx1].n1 = rra.n1, pair[i-1+indx1].n2 = rra.n2;
	}
	done:
	indx1 = indx2;
    }

    // mark duplicates
    indx1 = 0;
    for (i = 1; i < npair; i++) {
	if (pair[i].n1 == pair[indx1].n1 && pair[i].n2 == pair[indx1].n2)
	    pair[i].n1 = -1;
	else indx1 = i;
    }
    if (pair[0].n1 == 0 && pair[0].n2 == 0) pair[0].n1 = -1;

    int nz;


    for (nz = nnode, i = 0; i < npair; i++) {
	if (pair[i].n1 < 0) continue;
	nz += 2;
    }
    int *row_ind = new int[nz];
    int *col_ind = new int[nz];

    for (i = 0; i < nnode; i++) {  // diagonal elements
	row_ind[i] = i;
	col_ind[i] = i;
    }

    for (i = 0, j = nnode; i < npair; i++) {
	if (pair[i].n1 < 0) continue;
	row_ind[j] = pair[i].n1;  col_ind[j] = pair[i].n2;  j++;
	row_ind[j] = pair[i].n2;  col_ind[j] = pair[i].n1;  j++;
    }

    *_nz = nz;
    *_row_ind = row_ind;
    *_col_ind = col_ind;

    delete []pair;
}

int Mesh::FindBoundarySegment (const Point &p, int *n1, int *n2,
    double *dist1, double *dist2) const
{
    int el, nd, vertex, bnd;
    int closest_el = -1;
    int ndlist[2];
    bool vertical = (fabs (p[0]) < 1e-10);
    double m, d, dmin=1e10;
    double x[2], y[2];

    if (!vertical) m = p[1] / p[0];
    for (el = 0; el < elist.Len(); el++) {
	for (bnd = vertex = 0; vertex < elist[el]->nNode() && bnd<2; vertex++) {
	    nd = elist[el]->Node[vertex];
	    if (nlist[nd].BndTp() == BndType) {
		x[bnd] = nlist[nd][0];
		y[bnd] = nlist[nd][1];
		ndlist[bnd] = nd;
		bnd++;
	    }
	}
	if (BndType == MESH_EXTRAPOL) {
	    for (vertex = 0; vertex < elist[el]->nNode(); vertex++)
		if (nlist[elist[el]->Node[vertex]].BndTp() == BND_DIRICHLET ||
		    nlist[elist[el]->Node[vertex]].BndTp() == XLAYER_INTERNAL)
		    bnd=0;   // wrong side of the physical boundary
	}
	if (bnd == 2) {  // boundary element
	    if (vertical ? x[0]*x[1]<=0.0 : (m*x[0]-y[0])*(m*x[1]-y[1])<=0.0) {
		d = hypot (0.5*(x[0]+x[1])-p[0], 0.5*(y[0]+y[1])-p[1]);
		if (d<dmin) {
		    dmin = d;
		    *n1 = ndlist[0];
		    *n2 = ndlist[1];
		    *dist1 = p.Dist (nlist[*n1]);
		    *dist2 = p.Dist (nlist[*n2]);
		    closest_el = el;
		}
	    }
	}
    }
    return closest_el;
}

double Mesh::BoundaryDistance (int node) const
{
    if (nlist[node].isBnd()) return 0.0;
    // if the node is itself a boundary node, then its distance from the
    // boundary is zero

    Point p = nlist[node];
    double d, dmin = 1e10;

    for (int i = 0; i < nlist.Len(); i++)
	if (nlist[i].isBnd()) {
	    d = p.Dist (nlist[i]);
	    if (d < dmin) dmin = d;
	}
    return dmin;
}

void Mesh::SetBoundary (const Surface &_boundary)
{
    if (boundary) delete boundary;
    boundary = _boundary.Clone();
    if (is_set_up) {
        if (bnd_param) delete []bnd_param;
	bnd_param = new RVector[priv_nbnd];
	for (int i = 0; i < priv_nbnd; i++) {
	    bnd_param[i].New (boundary->ParamDim());
	    boundary->Point2Param (nlist[IndexBnd2Node[i]], bnd_param[i]);
	}
    }
}

int Mesh::BoundaryList (int **bndellist, int **bndsdlist) const
{
    const int MAXSDND = 10;
    int i, j, el, sd, n, nnd, nsd, ndi, found, listsize = 0;
    int *nsdnd, **sdnd, nnsdnd = 0, nsdndi, *sdndi;
    bndrec *bndlist = NULL, *nptr, *pptr;
    bool removed, valid, isbndsd;

    for (el = elist.Len()-1; el >= 0; el--) {
	Element *pel = elist[el];
	nnd = pel->nNode();
	nsd = pel->nSide();
	if (nsd > nnsdnd) {
	    if (nnsdnd) {
		delete []nsdnd;
		for (sd = 0; sd < nnsdnd; sd++) delete []sdnd[sd];
		delete []sdnd;
	    }
	    nsdnd = new int[nnsdnd = nsd];
	    sdnd = new int*[nsd];
	    for (sd = 0; sd < nsd; sd++) sdnd[sd] = new int[MAXSDND];
	}
	for (sd = 0; sd < nsd; sd++) {
	    nsdnd[sd] = pel->nSideNode(sd);
	    for (n = 0; n < nsdnd[sd]; n++)
		sdnd[sd][n] = pel->Node[pel->SideNode(sd,n)];
	}

	for (sd = 0; sd < nsd; sd++) {
	    nsdndi = nsdnd[sd];
	    sdndi = sdnd[sd];
	    for (isbndsd = true, n = 0; n < nsdndi; n++)
		if (!nlist[sdndi[n]].BndTp()) {
		    isbndsd = false; break;
		}
	    if (!isbndsd) continue;

	    if (BndType == MESH_EXTRAPOL) {
		for (isbndsd = true, n = 0; n < nsdndi; n++)
		    if (nlist[sdndi[n]].BndTp() != BND_INTERNAL) {
			    isbndsd=false; break;
			}
		if (!isbndsd) continue;
		for (valid = true, i = 0; i < nnd; i++) {
		    n = pel->Node[i];
		    if (nlist[n].BndTp() != BND_INTERNAL && nlist[n].BndTp()
			!= BND_NONE) {
			    valid = false; break;
			}
		}
		if (!valid) continue;
	    }

	    // if combination is already in the list, delete (no real boundary)
	    removed = false;
	    for (pptr=NULL, nptr=bndlist; nptr; pptr=nptr, nptr=nptr->next) {
		for (found = i = 0; i < nsdndi; i++) {
		    ndi = nptr->nd[i];
		    for (j = 0; j < nsdndi; j++)
			if (ndi == sdndi[j]) found++;
		}
		if (found == nsdndi) {
		    if (pptr) pptr->next = nptr->next;
		    else bndlist = nptr->next;
		    delete nptr;
		    removed = true;
		    listsize--;
		    break;
		}
	    }

	    // otherwise, add entry to list
	    if (!removed) {
		nptr = new bndrec;
		nptr->el = el;
		nptr->sd = sd;
		for (i = 0; i < nsdndi; i++)
		    nptr->nd[i] = sdndi[i];
		nptr->next = bndlist;
		bndlist = nptr;
		listsize++;
	    }
	}
    }
    if (nnsdnd) {
	delete []nsdnd;
	for (sd = 0; sd < nnsdnd; sd++) delete []sdnd[sd];
	delete []sdnd;
    }

    // create final list
    *bndellist = new int[listsize];
    *bndsdlist = new int[listsize];
    for (nptr = bndlist, i = 0; i < listsize; i++) {
	(*bndellist)[i] = nptr->el;
	(*bndsdlist)[i] = nptr->sd;
	pptr = nptr; nptr = nptr->next;
	delete pptr;
    }
    return listsize;
}

bool Mesh::PullToBoundary (const Point &p, Point &pshift, int &element,
    int &side, double sub) const
{
    // MS 25.11.99 This is a new version which works by mapping p into
    // local coordinates of eligible elements. The old version had problems
    // when applied to coarse meshes

    const int MAXSIDE = 20;
    int i, n, nmin, el, sd, nbndsd, bndsd[MAXSIDE];
    double d, dmin;

    // 1. Find boundary node closest to p
    for (n = 0, dmin = 1e10; n < nlen(); n++)
	if (nlist[n].isBnd() && (d = p.Dist (nlist[n])) < dmin)
	    dmin = d, nmin = n;

    // 2. Loop over elements containing nmin
    for (el = 0; el < elen(); el++) {
	if (!elist[el]->IsNode (nmin)) continue;
	nbndsd = elist[el]->BndSideList (nlist, bndsd);
	if (!nbndsd) continue; // element contains no boundary sides

	// 2.1 map to local coordinates
	Point ploc = elist[el]->Local (nlist, p);

	// 2.2 project to boundary sides
	for (i = 0; i < nbndsd; i++) {
	    sd = bndsd[i];
	    const RVector lnm = elist[el]->LNormal (sd);
	    Point nloc = elist[el]->Local (nlist,
		nlist[elist[el]->Node[elist[el]->SideNode(sd, 0)]]);
	    double pdist = (ploc - nloc) & lnm;
	    Point pnew = ploc - (lnm*pdist);
	    if (elist[el]->LContains (pnew)) {
	        RDenseMatrix egeom = ElGeom (el);
	        pshift = elist[el]->Global (nlist, pnew);
	        if (sub) {
		    RDenseMatrix lder  = elist[el]->LocalShapeD (pnew);
		    RDenseMatrix jacin = inverse (lder * egeom);
		    RVector dcos  = elist[el]->DirectionCosine (sd, jacin);
		    double ofs = (sub > 0 ? sub : 1.0/plist[nmin].Mus());
		    pshift -= dcos*ofs;
		    element = ElFind (pshift);
	        } else {
		    element = el;
		}
		side = sd;
		if (element >= 0) return true;
	    }
	}
    }

    return PullToBoundary_old (p, pshift, element, side, sub);
    // If this didn't work, fall back to the old method as a last resort
}

bool Mesh::PullToBoundary_old (const Point &p, Point &pshift, int &element,
    int &side, double sub) const
{
    const double EPS = 1e-10;
    int i, n, nmin, el, sd;
    double d, dmin, ofs, wght, wghtmax = 0.0;
    bool containsmin, isbnd;
    Point pnew (p);
    RVector normal (nlist[0].Dim());

    // 1. find boundary node closest to p
    for (n = 0, dmin = 1e10; n < nlen(); n++)
	if (nlist[n].isBnd() && (d = p.Dist (nlist[n])) < dmin)
	    dmin = d, nmin = n;

    // 2. generate a local normal direction
    for (el = 0; el < elen(); el++) {
	if (!elist[el]->IsNode (nmin)) continue;
	for (sd = 0; sd < elist[el]->nSide(); sd++) {
	    containsmin = false, isbnd = true;
	    for (i = 0; i < elist[el]->nSideNode (sd); i++) {
		n = elist[el]->Node[elist[el]->SideNode (sd, i)];
		if (n == nmin) containsmin = true;
		else if (!nlist[n].isBnd()) { isbnd = false; break; }
	    }
	    if (!containsmin || !isbnd) continue;

	    // normal for this side
	    RDenseMatrix egeom = ElGeom (el);
	    Point lcnt = elist[el]->SideCentre (sd);
	    Point gcnt = elist[el]->Global (nlist, lcnt);
	    RDenseMatrix lder  = elist[el]->LocalShapeD (lcnt);
	    RDenseMatrix jacin = inverse (lder * egeom);
	    RVector dcos  = elist[el]->DirectionCosine (sd, jacin);

	    d = p.Dist (ElSideCentre (el, sd));
	    if (d < EPS) d = EPS;
	    wght = 1.0/d;
	    if (wght > wghtmax) {
		element = el, side = sd;
		wghtmax = wght;
		if (!sub) {
		    double pdist = (p - gcnt) & dcos; // distance p to plane
		    pnew = p - (dcos * pdist);        // proj. of p on plane
		}
	    }
	    if (sub) {
		normal += dcos * wght;
	    }
	}
    }
    if (sub) {
	normal /= length (normal);
	ofs = (sub > 0 ? sub : 1.0/plist[nmin].Mus());
	pnew = p - (normal * ofs);
	element = ElFind (pnew);
    }
    pshift = pnew;
    return (element >= 0);
}

// ***********************************************
/*
void Mesh::SetupElementMatrices (void)
{
    int i, blen, *bndel, *bndsd, el, sd, elen = elist.Len();
    blen = BoundaryList (&bndel, &bndsd);

    if (elk) delete []elk;
    if (elc) delete []elc;
    if (ela) delete []ela;
    elk = new SymMatrix[elen];   // Int[Fi Fj] over el
    elc = new SymMatrix[elen];   // Int[Di Dj] over el

    for (el = 0; el < elen; el++) {
	int nodel = elist[el]->nNode();
	elk[el].New (nodel);
	elc[el].New (nodel);
	elist[el]->FTF_DTD (ElGeom (el), elc[el], elk[el]);
    }
    if (BndType == MESH_ROBIN) {     // need surface integrals
	ela = new SymMatrix[elen];   // Int[Fi Fj] over boundary sides of el
	for (el = 0; el < elen; el++) {
	    int nodel = elist[el]->nNode();
	    ela[el].New (nodel);
	}
	for (i = 0; i < blen; i++) {
	    el = bndel[i];
	    sd = bndsd[i];
	    ela[el] += elist[el]->FTF_bnd (ElGeom (el), sd);
	}
    }
    delete []bndel;
    delete []bndsd;
}
*/
// ***********************************************

void Mesh::SetupNeighbourList ()
{
    int el, side;
    if (nbhrs) return;   // list exists already
    nbhrs = new int*[elen()];
    for (el = 0; el < elen(); el++) {
	nbhrs[el] = new int[elist[el]->nSide()];
	for (side = 0; side < elist[el]->nSide(); side++)
	    nbhrs[el][side] = elist.EdgeAdjacentElement (el, side);
    }
}

void Mesh::NodeNeighbourList (int **_nnbhrs, int ***_nbhrs)
{
    const int chunksize = 64;
    int i, j, k, el, nodel, ni, nj, nn = nlen();

    // allocate memory for list
    int *nnbhrs  = new int[nn];
    int **nbhrs  = new int*[nn];
    int *bufsize = new int[nn];
    for (i = 0; i < nn; i++) {
        nbhrs[i] = new int[bufsize[i] = chunksize];
	nnbhrs[i] = 0;
    }

    // go through mesh and generate list
    for (el = 0; el < elen(); el++) {
        nodel = elist[el]->nNode();
	for (i = 0; i < nodel; i++) {
	    ni = elist[el]->Node[i];
	    for (j = 0; j < nodel; j++) {
	        if (j == i) continue;   // node isn't its own neighbour
		nj = elist[el]->Node[j];
		for (k = 0; k < nnbhrs[ni] && nbhrs[ni][k] != nj; k++);
		if (k == nnbhrs[ni]) {  // not already present
		    if (nnbhrs[ni] == bufsize[ni]) { // reallocate
		        int *tmp = new int[bufsize[ni]+chunksize];
			memcpy (tmp, nbhrs[ni], nnbhrs[ni]*sizeof(int));
			delete []nbhrs[ni];
			nbhrs[ni] = tmp;
			bufsize[ni] += chunksize;
			cerr << "Increasing neighbour buffer to "
			     << bufsize[ni] << endl;
		    }
		    nbhrs[ni][nnbhrs[ni]++] = nj;
		}
	    }
	}
    }

    // shrink list to required size
    for (i = 0; i < nn; i++) {
        int *tmp = new int[nnbhrs[i]];
	for (j = 0; j < nnbhrs[i]; j++) tmp[j] = nbhrs[i][j];
	delete []nbhrs[i];
	nbhrs[i] = tmp;
    }

    *_nnbhrs = nnbhrs;
    *_nbhrs = nbhrs;
    delete []bufsize;
}

// ***********************************************
// nodal centre of gravity

double Mesh::Size (Point *centre) const
{
    int i, d, dim = nlist[0].Dim();
    double size = 0.0;
    Point mx(dim), mn(dim);
    for (d = 0; d < dim; d++) {
	mx[d] = -1e10, mn[d] = 1e10;
	for (i = 0; i < nlist.Len(); i++) {
	    mx[d] = ::max (mx[d], nlist[i][d]);
	    mn[d] = ::min (mn[d], nlist[i][d]);
	}
	size = ::max (size, mx[d]-mn[d]);
    }
    if (centre) *centre = (mx+mn) * 0.5;
    return 0.5 * size;
}

// ***********************************************

void Mesh::BoundingBox (Point &mmin, Point &mmax) const
{
    int d, n, dim = nlist[0].Dim();
    mmin.New (dim);
    mmax.New (dim);
    for (d = 0; d < dim; d++) {
	mmin[d] = 1e6, mmax[d] = -1e6;
    }
    for (n = 0; n < nlist.Len(); n++) {
	for (d = 0; d < dim; d++) {
	    if (nlist[n][d] < mmin[d]) mmin[d] = nlist[n][d];
	    if (nlist[n][d] > mmax[d]) mmax[d] = nlist[n][d];
	}
    }
}

Point Mesh::BndIntersect (const Point &pt1, const Point &pt2)
{
    // Calculates the intersection of a straight line, given by points pt1
    // and pt2, with the mesh surface. The boundary must be bracketed by
    // pt1 and pt2 such that pt1 is inside, pt2 outside the mesh
    // This just performs a binary search, and assumes that there is only a
    // single intersection between pt1 and pt2.
    // Works in 2D and 3D, but is fairly inefficient

    xASSERT (ElFind (pt1) >= 0, First point must be inside mesh);
    xASSERT (ElFind (pt2) <  0, Second point must be outside mesh);

    const double acc = 1e-6;
    int i;
    RVector p1 = pt1;
    RVector p2 = pt2;
    RVector pm = (p1+p2) * 0.5;
    Point m (pm.Dim());
    double dist = length (p2-p1);
    while (dist > acc) {
        for (i = 0; i < pm.Dim(); i++) m[i] = pm[i];
        bool inside = (ElFind (m) >= 0);
	if (inside) p1 = pm;
	else        p2 = pm;
	dist *= 0.5;
	pm = (p1+p2) * 0.5;
    }
    for (i = 0; i < p1.Dim(); i++) m[i] = p1[i];
    // we use p1 to ensure that the returned point is inside the mesh
    return m;
}

// ***********************************************

#ifdef UNDEF
Point Mesh::BndIntersect (const Point &pt1, const Point &pt2)
{
    // this only works in 2D

    const double EPS = 1e-8;
    int n, nl, nr, nhi, nlo;
    double dx, dy, dir, dir2, a1, b1, a2, b2;
    double dirl = -100.0, dirr = 100.0, dirhi = -100.0, dirlo = 100.0;
    Point sc(2);
    dx = pt2[0]-pt1[0];
    dy = pt2[1]-pt1[1];
    a1 = dy/dx;
    b1 = pt1[1] - a1*pt1[0];
    dir = atan2 (dy, dx);
    for (n = 0; n < nlist.Len(); n++)
	if ((BndType == MESH_EXTRAPOL && nlist[n].BndTp() == BND_INTERNAL) ||
		(BndType != MESH_EXTRAPOL && nlist[n].isBnd())) {
	    dx = nlist[n][0]-pt1[0];
	    dy = nlist[n][1]-pt1[1];
	    dir2 = atan2 (dy, dx);
	    if (dir2 > dir) {
		if (dir2 < dirr) dirr = dir2, nr = n;
	    } else {
		if (dir2 > dirl) dirl = dir2, nl = n;
	    }
	    if (dir2 > dirhi) dirhi = dir2, nhi = n;
	    if (dir2 < dirlo) dirlo = dir2, nlo = n;
	}
    if (dirr > 99.0) dirr = dirlo, nr = nlo;	// take care of phase wrapping
    if (dirl < -99.0) dirl = dirhi, nl = nhi;
    dx = nlist[nr][0]-nlist[nl][0];
    dy = nlist[nr][1]-nlist[nl][1];

    // now intersection of two lines
    if (fabs (dx) < EPS) { // boundary segment vertical
	sc[0] = nlist[nr][0];
    } else {
	a2 = dy/dx;
	b2 = nlist[nl][1] - a2*nlist[nl][0];	
	sc[0] = -(b2-b1)/(a2-a1);
    }
    sc[1] = a1*sc[0] + b1;
    return sc;
}
#endif

// ***********************************************

int Mesh::CheckConsistency () const
{
    int el, nd, nnd, vx, vxx, err, dim, d, i, errorid=0;
    char *checklist;
    const double eps = 1e-8;

    dim = nlist[0].Dim();
    for (i = 1; i < nlist.Len(); i++)
	if (nlist[i].Dim() != dim) {
	    errorid = 1;
	    dERROR(Inconsistent node dimensions.);
	}

    for (el = 0; el < elist.Len(); el++)
	for (vx = 0; vx < elist[el]->nNode(); vx++)
	    if (elist[el]->Node[vx] < 0 || elist[el]->Node[vx] >= nlist.Len()) {
		errorid = 2;
		dERROR(Node reference out of range.);
	    }

    for (el = 0; el < elist.Len(); el++)
	for (vx = 0; vx < elist[el]->nNode()-1; vx++)
	    for (vxx = vx+1; vxx < elist[el]->nNode(); vxx++)
		if (elist[el]->Node[vx] == elist[el]->Node[vxx]) {
		    errorid = 3;
		    dERROR(Dublicate node reference in single element.);
		}

    checklist = new char[nlist.Len()];
    for (nd = 0; nd < nlist.Len(); nd++) checklist[nd] = 0;
    for (el = 0; el < elist.Len(); el++)
	for (vx = 0; vx < elist[el]->nNode(); vx++)
	    checklist[elist[el]->Node[vx]] = 1;
    for (err = nd = 0; nd < nlist.Len(); nd++)
	if (!checklist[nd]) err = 1;
    delete checklist;
    if (err) {
	errorid = 4;
	dERROR(Node list contains nodes not used by any element.);
    }

    for (nd = 0; nd < nlist.Len()-1; nd++)
	for (nnd = nd+1; nnd < nlist.Len(); nnd++) {
	    int ident = 1;
	    for (d = 0; d < dim; d++)
		if (fabs (nlist[nd][d] - nlist[nnd][d]) > eps)
		    ident = 0;
	    if (ident) {
		errorid = 5;
		dERROR(Nodes with identical coordinates found.);
	    }
	}
    return errorid;
}

// ***********************************************

typedef struct {
    int n;
    int *nd;
} SideRec;

// local MarkBoundary helper function: comparison function for qsort
int MarkBoundary_qsort_comp (const void *arg1, const void *arg2)
{
    int k;
    const SideRec *sd1 = (SideRec*)arg1;
    const SideRec *sd2 = (SideRec*)arg2;
    int nnd = ::min (sd1->n, sd2->n);
    for (k = 0; k < nnd; k++) {
        if (sd1->nd[k] < sd2->nd[k]) return -1;
        else if (sd1->nd[k] > sd2->nd[k]) return 1;
    }
    return 0; // should not happen
}

void Mesh::MarkBoundary ()
{
    int i, j, k, slen;
    int nel = elen();
    Element *pel;

    // first we generate a list of sides
    for (i = slen = 0; i < nel; i++)
        slen += elist[i]->nSide();

    SideRec *sd = new SideRec[slen];

    for (i = slen = 0; i < nel; i++) {
        pel = elist[i];
	for (j = 0; j < pel->nSide(); j++) {
	    sd[slen].n  = pel->nSideNode(j);
	    sd[slen].nd = new int[sd[slen].n];
	    for (k = 0; k < sd[slen].n; k++)
	        sd[slen].nd[k] = pel->Node[pel->SideNode(j,k)];
	    slen++;
	}
    }

    // order each side for nodes
    for (i = 0; i < slen; i++) {
        for (j = 0; j < sd[i].n-1; j++)
	    for (k = j+1; k < sd[i].n; k++)
	        if (sd[i].nd[j] > sd[i].nd[k]) {
		    int itmp = sd[i].nd[j];
		    sd[i].nd[j] = sd[i].nd[k];
		    sd[i].nd[k] = itmp;
		}
    }

    // order side list
    qsort (sd, slen, sizeof(SideRec), MarkBoundary_qsort_comp);
		
    // now all non-boundary sides should be listed as consecutive pairs,
    // so go through the list and look for duplicates

    int bndsd = 0, bndnd = 0;
    for (i = 0; i < nlen(); i++)
        nlist[i].SetBndTp (BND_NONE); // reset boundary flags

    for (i = 0; i < slen; i++) {
        bool isbnd = true;
        int n = sd[i].n;
	if (i) { // check left neighbour
	    if (n == sd[i-1].n) {
	        for (j = 0; j < n; j++)
		    if (sd[i].nd[j] != sd[i-1].nd[j]) break;
		if (j == n) isbnd = false;
	    }
	}
	if (i < slen-1 && isbnd) { // check right neighbour
	    if (n == sd[i+1].n) {
	        for (j = 0; j < n; j++)
		    if (sd[i].nd[j] != sd[i+1].nd[j]) break;
		if (j == n) isbnd = false;
	    }
	}
	if (isbnd) {
	    bndsd++;
	    for (j = 0; j < n; j++)
	        nlist[sd[i].nd[j]].SetBndTp (BND_DIRICHLET);
	}
    }

    // count boundary nodes
    for (i = 0; i < nlen(); i++)
        if (nlist[i].isBnd()) bndnd++;

    // cleanup
    for (i = 0; i < slen; i++)
	delete []sd[i].nd;
    delete []sd;
}

// ***********************************************

void Mesh::ResetCoeff_region (ParameterType prmtp, double val, int region)
{
    for (int i = 0; i < nlist.Len(); i++)
	if (nlist[i].Region() == region)
	    plist[i].SetParam (prmtp, val);
}

void Mesh::ResetCoeff_sqrt (ParameterType prmtp, double cnt, double bnd)
{
    double rad, rad0, val;
    Point pcnt(nlist[0].Dim());
    Size (&pcnt);

    for (int i = 0; i < nlist.Len(); i++) {
	Point sec = BndIntersect (pcnt, nlist[i]);
	rad0 = pcnt.Dist (sec);
	rad  = pcnt.Dist (nlist[i]);
	val  = sqrt (cnt*cnt - (cnt*cnt-bnd*bnd)*rad/rad0);
	plist[i].SetParam (prmtp, val);
    }
}

FELIB istream& operator>> (istream& i, Mesh& mesh)
{
    char cbuf[200];
    bool bDirichlet = false, bInternal = false, bRobin = false;
    do {
	i.getline (cbuf, 200);
    } while (i.good() && strncmp (cbuf, "MeshData 5.0", 12));
    if (!i.good()) {
        if (mesh.trap_load_error)
	    xERROR(Mesh file not found or invalid);
	return i;
    }
    i >> mesh.nlist >> mesh.elist >> mesh.plist;
    //cout << "Parameter list type " << (mesh.plist.pltype() == PLIST_BY_ELEMENT ? "by element" : "by node (default) " ) <<endl;
    if( mesh.plist.pltype() == PLIST_BY_ELEMENT) {
      xASSERT(mesh.elist.Len() == mesh.plist.Len(),
	    Element list and parameter list sizes differ.);
    }
    else {
      xASSERT(mesh.nlist.Len() == mesh.plist.Len(),
	    Node list and parameter list sizes differ.);
    }
    for (int nd = 0; nd < mesh.nlist.Len(); nd++)
	switch (mesh.nlist[nd].BndTp()) {
	    case BND_DIRICHLET: bDirichlet=true; break;
	    case BND_INTERNAL:  bInternal=true;  break;
	    case BND_ROBIN:     bRobin=true;     break;
	}
    if (bDirichlet) {
	if (bRobin) mesh.BndType = MESH_ODD;
	else mesh.BndType = (bInternal ? MESH_EXTRAPOL : MESH_DIRICHLET);
    } else if (bRobin) {
	mesh.BndType = (bInternal ? MESH_ODD : MESH_ROBIN);
    } else mesh.BndType = MESH_ODD;

    // read mesh surface, if exists
    if (mesh.boundary) delete mesh.boundary;
    mesh.boundary = ReadSurface (i);

    mesh.lastel_found = 0;
    return i;
}

FELIB ostream& operator<< (ostream& o, Mesh& mesh)
{
    o << "MeshData 5.0\n\n";
    o << mesh.nlist << endl << mesh.elist << endl << mesh.plist;
    if (mesh.Boundary()) {
        Surface *bnd = mesh.Boundary();
        o << endl << (*bnd) << endl;
    }
    return o;
}

void Mesh::put (ostream &os, ParameterType p1, ParameterType p2,
		ParameterType p3)
{
    plist.SetOutputParameterTypes (p1, p2, p3);
    os << "MeshData 5.0\n\n";
    os << nlist << endl << elist << endl << plist;
    if (Boundary()) {
        Surface *bnd = Boundary();
        os << endl << (*bnd) << endl;
    }
}

// =========================================================================
// Nonmember functions
// =========================================================================

// Add a component to element matrix 'M', given 'mesh' and 'el'
// Element matrix type is defined by 'mode' (see mesh.h)
// nodal or element coefficients are given by 'coeff'

void AddToElMatrix (const Mesh &mesh, int el, CGenericSparseMatrix &M,
    const RVector *coeff, int mode)
{
    int i, j, is, js, nnode;
    complex entry;

    nnode = mesh.elist[el]->nNode();
    for (i = 0; i < nnode; i++) {
	is = mesh.elist[el]->Node[i];
	for (j = 0; j < nnode; j++) {
	    js = mesh.elist[el]->Node[j];
	    switch (mode) {
	    case ASSEMBLE_FF:
		entry.re = mesh.elist[el]->IntFF (i, j);
		break;
	    case ASSEMBLE_DD:
		entry.re = mesh.elist[el]->IntDD (i, j);
		break;
	    case ASSEMBLE_PFF:
		entry.re = mesh.elist[el]->IntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PDD:
		entry.re = mesh.elist[el]->IntPDD (i, j, *coeff);
		break;
	    case ASSEMBLE_BNDPFF:
		entry.re = mesh.elist[el]->BndIntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PFF_EL:
		entry.re = mesh.elist[el]->IntFF (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_PDD_EL:
		entry.re = mesh.elist[el]->IntDD (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_BNDPFF_EL:
		entry.re = mesh.elist[el]->BndIntFF (i, j) * (*coeff)[el];
		break;
	    }
	    M.Add (is, js, entry);
	}
    }
}

// Add a component to system matrix 'M', given 'mesh'
// Element matrix type is defined by 'mode' (see mesh.h)
// nodal coefficients are given by 'coeff'

void AddToSysMatrix (const Mesh &mesh, RGenericSparseMatrix &M,
    const RVector *coeff, int mode)
{
    int i, j, is, js, el, nnode;
    double entry;

    for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	for (i = 0; i < nnode; i++) {
	    is = mesh.elist[el]->Node[i];
	    for (j = 0; j < nnode; j++) {
	        js = mesh.elist[el]->Node[j];
	        switch (mode) {
		case ASSEMBLE_FF:
		    entry = mesh.elist[el]->IntFF (i, j);
		    break;
		case ASSEMBLE_DD:
		    entry = mesh.elist[el]->IntDD (i, j);
		    break;
		case ASSEMBLE_PFF:
		    entry = mesh.elist[el]->IntPFF (i, j, *coeff);
		    break;
		case ASSEMBLE_PDD:
		    entry = mesh.elist[el]->IntPDD (i, j, *coeff);
		    break;
		case ASSEMBLE_BNDPFF:
		    entry = mesh.elist[el]->BndIntPFF (i, j, *coeff);
		    break;
		case ASSEMBLE_PFF_EL:
		    entry = mesh.elist[el]->IntFF (i, j) * (*coeff)[el];
		    break;
		case ASSEMBLE_PDD_EL:
		    entry = mesh.elist[el]->IntDD (i, j) * (*coeff)[el];
		    break;
		case ASSEMBLE_BNDPFF_EL:
		    entry = mesh.elist[el]->BndIntFF (i, j) * (*coeff)[el];
		    break;
		}
		M.Add (is, js, entry);
		//M(is, js) += entry;
	    }
	}
    }
}

void AddToSysMatrix (const Mesh &mesh, CGenericSparseMatrix &M,
    const RVector *coeff, int mode)
{
    int el;
    complex entry;

    for (el = 0; el < mesh.elen(); el++) {
	AddToElMatrix (mesh, el, M, coeff, mode);
    }
}

void AddToSysMatrix (const Mesh &mesh, CGenericSparseMatrix &M,
    const double coeff, int mode)
{
    int i, j, is, js, el, nnode;
    complex entry;

    for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	for (i = 0; i < nnode; i++) {
	    is = mesh.elist[el]->Node[i];
	    for (j = 0; j < nnode; j++) {
	        js = mesh.elist[el]->Node[j];
	        switch (mode) {
		case ASSEMBLE_CFF:
		    entry.re = mesh.elist[el]->IntFF (i, j) * coeff;
		    break;
		case ASSEMBLE_CDD:
		    entry.re = mesh.elist[el]->IntDD (i, j) * coeff;
		    break;
		case ASSEMBLE_iCFF:
		    entry.im = mesh.elist[el]->IntFF (i, j) * coeff;
		    break;
		case ASSEMBLE_iCDD:
		    entry.im = mesh.elist[el]->IntDD (i, j) * coeff;
		    break;
		}
		M.Add (is, js, entry);
		//M(is, js) += entry;
	    }
	}
    }  
}

void AddToSysMatrix_elasticity (const Mesh &mesh, RGenericSparseMatrix &M,
    const RVector &modulus, const RVector &pratio)
{
    int i, j, k, m, is, js, id, jd, el, nnode, *node;
    int dim = mesh.Dimension();
    Element *pel;

    for (el = 0; el < mesh.elen(); el++) {
        pel   = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
        RDenseMatrix K = pel->ElasticityStiffnessMatrix (
	    modulus[el], pratio[el]);
	for (i = 0; i < nnode; i++) {
  	    is = node[i]*dim; // row block offset into system matrix
	    id = i*dim;       // row block offset into element matrix
	    for (j = 0; j < nnode; j++) {
	        js = node[j]*dim;  // column block offset into system matrix
		jd = j*dim;        // column block offset into element matrix

		// copy each dim x dim block from element to system matrix
		for (k = 0; k < dim; k++)
		    for (m = 0; m < dim; m++)
			M.Add (is+k, js+m, K(id+k,jd+m));
		        //M(is+k,js+m) += K(id+k,jd+m);
	    }
	}
    }
}

void AddToRHS_thermal_expansion (const Mesh &mesh, RVector &rhs,
    const RVector &modulus, const RVector &pratio,
    const RVector &thermal_expansion, double deltaT)
{
    int i, d, el, nnode, *node;
    int dim = mesh.Dimension();
    Element *pel;
    RVector eps0(dim*2);
    Point dummy(dim);

    for (el = 0; el < mesh.elen(); el++) {
	pel = mesh.elist[el];
	nnode = pel->nNode();
	node = pel->Node;
	double size = pel->Size();
	RVector fe = pel->ThermalExpansionVector (modulus[el], pratio[el], thermal_expansion[el], deltaT);

	//for (d = 0; d < dim; d++) eps0[d] = thermal_expansion[el]*deltaT;

	//RVector fe = ATx (pel->StrainDisplacementMatrix (dummy),
	//  Ax (pel->IsotropicElasticityMatrix (modulus[el], pratio[el]), eps0));

	for (i = 0; i < nnode; i++) {
	    for (d = 0; d < 3; d++)
		rhs[node[i]*3+d] += fe[i*3+d];
	}
    }
}

void AddToRHS_elasticity (const Mesh &mesh, RVector &rhs, const RVector *coeff,
    int mode)
{
    int i, j, is, js, el, d, nnode, *node;
    int dim = mesh.Dimension();
    double f, ff;
    Element *pel;

    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    is = pel->Node[i];
	    switch (mode) {
	    case RHS_P:
	        f = pel->IntF(i);
		for (d = 0; d < dim; d++)
		    rhs[is*3+d] += (*coeff)[is*3+d] * f;
		break;
	    default:
	        for (j = 0; j < nnode; j++) {
		    js = pel->Node[j];
		    switch (mode) {
		    case RHS_PF:
		        ff = pel->IntFF(i,j);
			for (d = 0; d < dim; d++)
			    rhs[js*3+d] += (*coeff)[is*3+d] * ff;
			break;
		    case RHS_BNDPF:
		        ff = pel->BndIntFF(i,j);
			for (d = 0; d < dim; d++)
			    rhs[js*3+d] += (*coeff)[is*3+d] * ff;
			break;
		    default:
		        xERROR(Assembly mode not supported);
			return;
		    }
		}
	    }
	}
    }
}

int MakeNodalFreedomArray (int *&nf, int nlen, int dofnod, bool *rest)
{
    int i, j, idx;
    nf = new int[nlen*dofnod];
    for (i = idx = 0; i < nlen; i++)
        for (j = 0; j < dofnod; j++)
	    nf[i*dofnod+j] = (rest[i*dofnod+j] ? -1 : idx++);
    return idx;
}

void AddElasticStrainDisplacementToSysMatrix (const Mesh &mesh,
    RGenericSparseMatrix &M, double E, double nu, int *nf)
{
    // EXPERIMENTAL
    // This does not implement analytic integration over elements yet

    const int dofnod = 2; // for now

    int i, j, is, js, ii, jj, ip, jp, el, nq, dim = mesh.Dimension();
    const double *wght;
    const Point *absc;
    RDenseMatrix D;
    if (dim == 2) {
        D.Zero (3,3);
	D(0,0) = D(1,1) = 1.0;
	D(0,1) = D(1,0) = nu/(1.0-nu);
	D(2,2) = (1.0-2.0*nu)/(2.0*(1.0-nu));
	D = D * (E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu)));
    } else {
        xERROR(Not implemented);
    }
    
    for (el = 0; el < mesh.elen(); el++) {
        Element *pel = mesh.elist[el];
	int nnode = pel->nNode();
	int dofel = nnode*dofnod;
	RDenseMatrix egeom = mesh.ElGeom (el);
	RDenseMatrix elk(dofel, dofel);
	nq = pel->QuadRule (2, &wght, &absc); // check order!
	for (i = 0; i < nq; i++) {
	    RDenseMatrix lder = pel->LocalShapeD (absc[i]);
	    RDenseMatrix jac = lder*egeom;
	    RDenseMatrix jacin (dim,dim);
	    double d = det (jac, &jacin);
	    RDenseMatrix gder = jacin*lder;
	    RDenseMatrix B = pel->ElasticStrainDisplacement (absc[i], gder);
	    elk = elk + (transpose (B) * (D*B)) * (wght[i]*fabs(d));
	}
	for (i = 0; i < nnode; i++) {
	    is = pel->Node[i];
	    for (ii = 0; ii < dofnod; ii++) {
	        if ((ip = nf[is*dofnod+ii]) < 0) continue; // restricted node
	        for (j = 0; j < nnode; j++) {
		    js = pel->Node[j];
		    for (jj = 0; jj < dofnod; jj++) {
		        if ((jp = nf[js*dofnod+jj]) < 0) continue;
			M.Add (ip, jp, elk(i*dofnod+ii,j*dofnod+jj));
			//M(ip,jp) += elk(i*dofnod+ii,j*dofnod+jj);
		    }
		}
	    }
	}
    }
}

// Generate a sparse matrix which allows mapping of a solution from
// mesh nodal base to a regular pixel grid

RGenericSparseMatrix *GridMapMatrix (const Mesh &mesh, const IVector &gdim,
    const Point *bbmin, const Point *bbmax, const int *elref)
{
    const double EPS = 1e-8;
    int i, j, ix, iy, iz, idx, d, dim = mesh.Dimension();
    int nx = gdim[0], ny = gdim[1], nz = (dim > 2 ? gdim[2]:1);
    int rlen = nx*ny*nz;
    Point p(dim), loc(dim);
    bool local_elref;

    dASSERT(gdim.Dim() == dim, Parameter 3 wrong dimension);

    // bounding box of grid region  
    Point mesh_bbmin, mesh_bbmax;
    if (!bbmin || !bbmax) {
        mesh.BoundingBox (mesh_bbmin, mesh_bbmax);
	if (!bbmin) bbmin = &mesh_bbmin;
	if (!bbmax) bbmax = &mesh_bbmax;
    }
    dASSERT(bbmin->Dim() == dim, Parameter 4 wrong dimension);
    dASSERT(bbmax->Dim() == dim, Parameter 5 wrong dimension);

    // grid size and spacing
    RVector W(dim), dW(dim);
    for (d = 0; d < dim; d++)
        dW[d] = (W[d] = ((*bbmax)[d]-(*bbmin)[d]))/(double)(gdim[d]-1);

    // allocate element index list
    local_elref = (elref == 0);
    if (local_elref)
        elref = GenerateElementPixelRef (mesh, gdim, bbmin, bbmax);
#ifdef UNDEF
    int *node, nnode, el;
    int rlen = nx*ny*nz;
    int *elref = new int[rlen];
    for (i = 0; i < rlen; i++) elref[i] = -1;

    // associate elements with grid points
    Point elbbmin(dim), elbbmax(dim);
    for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	node  = mesh.elist[el]->Node;

        // compute element bounding box
        for (d = 0; d < dim; d++) elbbmin[d] = 1e10, elbbmax[d] = -1e10;
	for (i = 0; i < nnode; i++) {
	    Point &nd = mesh.nlist[node[i]];
	    for (d = 0; d < dim; d++) {
	        if (nd[d] < elbbmin[d]) elbbmin[d] = nd[d];
		if (nd[d] > elbbmax[d]) elbbmax[d] = nd[d];
	    }
	}
	for (d = 0; d < dim; d++) { // inflate a bit to avoid edge problems
	    elbbmin[d] -= EPS;
	    elbbmax[d] += EPS;
	}

	// search all grid points inside the element bounding box
	int rmin[3], rmax[3];
	for (d = 0; d < dim; d++) {
	    if (elbbmin[d] > (*bbmax)[d]) goto no_overlap;
	    if (elbbmax[d] < (*bbmin)[d]) goto no_overlap;
	    if (elbbmin[d] <= (*bbmin)[d] || dW[d] == 0) rmin[d] = 0;
	    else rmin[d] = (int)((elbbmin[d]-(*bbmin)[d])/dW[d]);
	    if (elbbmax[d] >= (*bbmax)[d] || dW[d] == 0) rmax[d] = gdim[d];
	    else rmax[d] = (int)((elbbmax[d]-(*bbmin)[d])/dW[d] + 1.0);
	}
	if (dim < 3) rmin[2] = 0, rmax[2] = 1;
	for (iz = rmin[2]; iz < rmax[2]; iz++) {
	    if (dim > 2) p[2] = iz*dW[2] + (*bbmin)[2];
	    for (iy = rmin[1]; iy < rmax[1]; iy++) {
	        p[1] = iy*dW[1] + (*bbmin)[1];
		for (ix = rmin[0]; ix < rmax[0]; ix++) {
		    idx = ix + nx*(iy + ny*iz);
		    if (elref[idx] < 0) {
		        p[0] = ix*dW[0] + (*bbmin)[0];
			if (mesh.elist[el]->GContains (p, mesh.nlist))
			    elref[idx] = el;
		    }
		}
	    }
	}
    no_overlap:;
    }

#endif
    // get map matrix structure
    int nzero = 0;
    for (i = 0; i < rlen; i++)
        if (elref[i] >= 0) nzero += mesh.elist[elref[i]]->nNode();
    int *rowptr = new int[rlen+1];
    int *colidx = new int[nzero];
    double *data = new double[nzero];
    rowptr[0] = 0;

    // populate the mapping matrix
    for (iz = i = idx = 0; iz < nz; iz++) {
        if (dim > 2) p[2] = iz*dW[2] + (*bbmin)[2];
	for (iy = 0; iy < ny; iy++) {
	    p[1] = iy*dW[1] + (*bbmin)[1];
	    for (ix = 0; ix < nx; ix++) {
	        if (elref[i] >= 0) {
		    p[0] = ix*dW[0] + (*bbmin)[0];
		    Element *pel = mesh.elist[elref[i]];
		    RVector fun = pel->GlobalShapeF (mesh.nlist, p);
		    for (j = 0; j < fun.Dim(); j++) {
		        colidx[idx] = pel->Node[j];
			data[idx] = fun[j];
			idx++;
		    }
		}
		rowptr[++i] = idx;
	    }
	}
    }
    if (local_elref) delete []elref;

    RCompRowMatrix *B = new RCompRowMatrix (rlen, mesh.nlen(),
        rowptr, colidx, data);
    delete []rowptr;
    delete []colidx;
    delete []data;
    return B;
}

#ifdef UNDEF // not finished yet!
RGenericSparseMatrix *GridMapMatrix (const Mesh &mesh, const IVector &gdim,
    const Point *bbmin, const Point *bbmax)
{
    const double EPS = 1e-8;
    const int nsub = 10; // sub-sampling density

    int i, j, ix, iy, iz, idx, d, dim = mesh.Dimension();
    int nx = gdim[0], ny = gdim[1], nz = (dim > 2 ? gdim[2]:1);
    int rlen = nx*ny*nz;
    Point p(dim), loc(dim);

    dASSERT(gdim.Dim() == dim, Parameter 3 wrong dimension);

    // bounding box of grid region  
    Point mesh_bbmin, mesh_bbmax;
    if (!bbmin || !bbmax) {
        mesh.BoundingBox (mesh_bbmin, mesh_bbmax);
	if (!bbmin) bbmin = &mesh_bbmin;
	if (!bbmax) bbmax = &mesh_bbmax;
    }
    dASSERT(bbmin->Dim() == dim, Parameter 4 wrong dimension);
    dASSERT(bbmax->Dim() == dim, Parameter 5 wrong dimension);

    IVector gsdim(dim); // subsampled grid
    for (i = 0; i < dim; i++) gsdim[i] = gdim[i]*nsub-1;
    int nsx = gsdim[0], nsy = gsdim[1], nsz = (dim > 2 ? gsdim[2]:1);
    int rslen = nsx*nsy*nsz;

    // grid size and spacing
    RVector W(dim), dW(dim);
    for (d = 0; d < dim; d++)
        dW[d] = (W[d] = ((*bbmax)[d]-(*bbmin)[d]))/(double)(gsdim[d]-1);

    // allocate element index list
    int *elref = GenerateElementPixelRef (mesh, gsdim, bbmin, bbmax);

    // get map matrix structure
    int nzero = 0;
    for (i = 0; i < rslen; i++)
        if (elref[i] >= 0) nzero += mesh.elist[elref[i]]->nNode();
    int *rowptr = new int[rlen+1];
    int *colidx = new int[nzero];
    double *data = new double[nzero];
    rowptr[0] = 0;

    // populate the mapping matrix
    for (iz = i = idx = 0; iz < nsz; iz++) {
        if (dim > 2) p[2] = iz*dW[2] + (*bbmin)[2];
	for (iy = 0; iy < nsy; iy++) {
	    p[1] = iy*dW[1] + (*bbmin)[1];
	    for (ix = 0; ix < nsx; ix++) {
	        if (elref[i] >= 0) {
		    p[0] = ix*dW[0] + (*bbmin)[0];
		    Element *pel = mesh.elist[elref[i]];
		    RVector fun = pel->GlobalShapeF (mesh.nlist, p);
		    RVector gfun = GridShapeF (ix,nsx,nsub, iy,nsy,nsub, iz,nsz,nsub);

		    for (j = 0; j < fun.Dim(); j++) {
		        colidx[idx] = pel->Node[j];
			data[idx] = fun[j];
			idx++;
		    }
		}
		rowptr[++i] = idx;
	    }
	}
    }
    delete []elref;

    RCompRowMatrix *B = new RCompRowMatrix (rlen, mesh.nlen(),
        rowptr, colidx, data);
    delete []rowptr;
    delete []colidx;
    delete []data;
    return B;
}
#endif
   
RGenericSparseMatrix *NodeMapMatrix (const Mesh &mesh, const IVector &gdim,
    const Point *bbmin, const Point *bbmax, const int *elref)
{
    const double EPS = 1e-8;
    int i, n, d, i0[3], idx, idx0, dim = mesh.Dimension();
    int nzero, *rowptr, *colidx;
    int nx = gdim[0], ny = gdim[1], nz = (dim > 2 ? gdim[2] : 1);
    int pix;
    int rlen = nx*ny*nz;
    double *data, d0[3];
    Point p(dim);
    bool local_elref;

    dASSERT(gdim.Dim() == dim, Parameter 2 wrong dimension);

    // bounding box of grid region  
    Point mesh_bbmin, mesh_bbmax;
    if (!bbmin || !bbmax) {
        mesh.BoundingBox (mesh_bbmin, mesh_bbmax);
	if (!bbmin) bbmin = &mesh_bbmin;
	if (!bbmax) bbmax = &mesh_bbmax;
    }
    dASSERT(bbmin->Dim() == dim, Parameter 4 wrong dimension);
    dASSERT(bbmax->Dim() == dim, Parameter 5 wrong dimension);

    // grid size and spacing
    RVector W(dim), dW(dim);
    for (d = 0; d < dim; d++)
        dW[d] = (W[d] = ((*bbmax)[d]-(*bbmin)[d]))/(double)(gdim[d]-1);

    // allocate element index list
    local_elref = (elref == 0);
    if (local_elref)
        elref = GenerateElementPixelRef (mesh, gdim, bbmin, bbmax);
#ifdef UNDEF
    int *node, nnode, ix, iy, iz, el;
    // create mask for pixel outside mesh support
    Point elbbmin(dim), elbbmax(dim);
    int *elref = new int[rlen];
    for (i = 0; i < rlen; i++) elref[i] = -1;
    for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	node  = mesh.elist[el]->Node;

        // compute element bounding box
        for (d = 0; d < dim; d++) elbbmin[d] = 1e10, elbbmax[d] = -1e10;
	for (i = 0; i < nnode; i++) {
	    Point &nd = mesh.nlist[node[i]];
	    for (d = 0; d < dim; d++) {
	        if (nd[d] < elbbmin[d]) elbbmin[d] = nd[d];
		if (nd[d] > elbbmax[d]) elbbmax[d] = nd[d];
	    }
	}
	for (d = 0; d < dim; d++) { // inflate a bit to avoid edge problems
	    elbbmin[d] -= EPS;
	    elbbmax[d] += EPS;
	}

	// search all grid points inside the element bounding box
	int rmin[3], rmax[3];
	for (d = 0; d < dim; d++) {
	    if (elbbmin[d] > (*bbmax)[d]) goto no_overlap;
	    if (elbbmax[d] < (*bbmin)[d]) goto no_overlap;
	    if (elbbmin[d] <= (*bbmin)[d] || dW[d] == 0) rmin[d] = 0;
	    else rmin[d] = (int)((elbbmin[d]-(*bbmin)[d])/dW[d]);
	    if (elbbmax[d] >= (*bbmax)[d] || dW[d] == 0) rmax[d] = gdim[d];
	    else rmax[d] = (int)((elbbmax[d]-(*bbmin)[d])/dW[d] + 1.0);
	}
	if (dim < 3) rmin[2] = 0, rmax[2] = 1;
	for (iz = rmin[2]; iz < rmax[2]; iz++) {
	    if (dim > 2) p[2] = iz*dW[2] + (*bbmin)[2];
	    for (iy = rmin[1]; iy < rmax[1]; iy++) {
	        p[1] = iy*dW[1] + (*bbmin)[1];
		for (ix = rmin[0]; ix < rmax[0]; ix++) {
		    idx = ix + nx*(iy + ny*iz);
		    if (elref[idx] < 0) {
		        p[0] = ix*dW[0] + (*bbmin)[0];
			if (mesh.elist[el]->GContains (p, mesh.nlist))
			    elref[idx] = el;
		    }
		}
	    }
	}
    no_overlap:;
    }
#endif
    for (nzero = mesh.nlen(), d = 0; d < dim; d++) nzero *= 2;
    // assuming linear interpolation between nearest pixels

    rowptr = new int[mesh.nlen()+1];
    colidx = new int[nzero];
    data   = new double[nzero];
    rowptr[0] = 0;

    // build map matrix structure
    for (n = idx = 0; n < mesh.nlen(); n++) {
        Point &p = mesh.nlist[n];
	for (d = 0; d < dim; d++) {
	    if (p[d] >= (*bbmin)[d] && p[d] <= (*bbmax)[d]) {
	        i0[d] = (int)((p[d]-(*bbmin)[d])/dW[d]);
		if      (i0[d] < 0)          i0[d] = 0;
		else if (i0[d] >= gdim[d]-1) i0[d] = gdim[d]-2;
		d0[d] = (p[d] - (i0[d]*dW[d] + (*bbmin)[d]))/dW[d];
		if      (d0[d] < 0.0) d0[d] = 0.0;
		else if (d0[d] > 1.0) d0[d] = 1.0;
	    } else goto nosupport;
	}
	if (dim < 3) {  // 2D case
	    idx0 = idx;
	    if (elref[pix = i0[0] + nx*i0[1]] >= 0) {
	        colidx[idx] = pix;
		data[idx]   = (1.0-d0[0]) * (1.0-d0[1]);
		idx++;
	    }
	    if (elref[pix = (i0[0]+1) + nx*i0[1]] >= 0) {
	        colidx[idx] = pix;
		data[idx]   = d0[0] * (1.0-d0[1]);
		idx++;
	    }
	    if (elref[pix = i0[0] + nx*(i0[1]+1)] >= 0) {
	        colidx[idx] = pix;
		data[idx]   = (1.0-d0[0]) * d0[1];
		idx++;
	    }
	    if (elref[pix = (i0[0]+1) + nx*(i0[1]+1)] >= 0) {
	        colidx[idx] = pix;
		data[idx]   = d0[0] * d0[1];
		idx++;
	    }
	    if (idx == idx0) {
	        cerr << "Node without pixel support found" << endl;
	    }
	    if (idx > idx0 && idx-idx0 < 4) {
	        // incomplete support: re-normalise pixel weights
	        double weight = 0.0;
		for (i = idx0; i < idx; i++) weight += data[i];
		if (weight > EPS)
		    for (i = idx0; i < idx; i++)
		        data[i] /= weight;
		else // last resort - equal weights
		    for (i = idx0; i < idx; i++)
		        data[i] = 1.0/(double)(idx-idx0);
	    }
	} else {  // 3D case
	    idx0 = idx;
	    if (elref[pix = i0[0] + nx*(i0[1] + ny*i0[2])] >= 0) {
	        colidx[idx] = pix;
		data[idx] = (1.0-d0[0]) * (1.0-d0[1]) * (1.0-d0[2]);
		idx++;
	    }
	    if (elref[pix = (i0[0]+1) + nx*(i0[1] + ny*i0[2])] >= 0) {
	        colidx[idx] = pix;
		data[idx] = d0[0] * (1.0-d0[1]) * (1.0-d0[2]);
		idx++;
	    }
	    if (elref[pix = i0[0] + nx*((i0[1]+1) + ny*i0[2])] >= 0) {
	        colidx[idx] = pix;
		data[idx] = (1.0-d0[0]) * d0[1] * (1.0-d0[2]);
		idx++;
	    }
	    if (elref[pix =(i0[0]+1) + nx*((i0[1]+1) + ny*i0[2])] >= 0) {
	        colidx[idx] = pix;
		data[idx] = d0[0] * d0[1] * (1.0-d0[2]);
		idx++;
	    }
	    if (elref[pix = i0[0] + nx*(i0[1] + ny*(i0[2]+1))] >= 0) {
	        colidx[idx] = pix;
		data[idx] = (1.0-d0[0]) * (1.0-d0[1]) * d0[2];
		idx++;
	    }
	    if (elref[pix = (i0[0]+1) + nx*(i0[1] + ny*(i0[2]+1))] >= 0) {
	        colidx[idx] = pix;
		data[idx] = d0[0] * (1.0-d0[1]) * d0[2];
		idx++;
	    }
	    if (elref[pix = i0[0] + nx*((i0[1]+1) + ny*(i0[2]+1))] >= 0) {
	        colidx[idx] = pix;
		data[idx] = (1.0-d0[0]) * d0[1] * d0[2];
		idx++;
	    }
	    if (elref[pix = (i0[0]+1) + nx*((i0[1]+1) + ny*(i0[2]+1))] >= 0){
 	        colidx[idx] = pix;
		data[idx] = d0[0] * d0[1] * d0[2];
		idx++;
	    }
	    if (idx > idx0 && idx-idx0 < 8) {
	        // incomplete support: re-normalise pixel weights
	        double weight = 0.0;
		for (i = idx0; i < idx; i++) weight += data[i];
		if (weight > EPS)
		    for (i = idx0; i < idx; i++)
		        data[i] /= weight;
		else // last resort - equal weights
		    for (i = idx0; i < idx; i++)
		        data[i] = 1.0/(double)(idx-idx0);
	    }
	}
    nosupport:
	rowptr[n+1] = idx;
    }
    if (local_elref) delete []elref;

    RCompRowMatrix *BI = new RCompRowMatrix (mesh.nlen(), rlen,
        rowptr, colidx, data);
    delete []rowptr;
    delete []colidx;
    delete []data;
    return BI;
}

RGenericSparseMatrix *NodeMapMatrix2 (const Mesh &mesh, const IVector &gdim,
    const Point *bbmin, const Point *bbmax, const int *elref)
{
    int dim = mesh.Dimension();
    int nx = gdim[0], ny = gdim[1], nz = (dim > 2 ? gdim[2] : 1);
    int rlen = nx*ny*nz;
    Point p(dim);

    int i, j, ix, iy, iz, n, el, nzero;
    int elen = mesh.elen();
    int nlen = mesh.nlen();
    int *egno = new int[elen];
    int *ngno = new int[nlen];
    double *scal = new double[nlen];

    for (el = 0; el < elen; el++)
	egno[el] = 0;
    for (n = 0; n < nlen; n++) {
	ngno[n] = 0;
	scal[n] = 0.0;
    }
    for (i = 0; i < rlen; i++)
	if (elref[i] >= 0) egno[elref[i]]++;
    for (el = 0; el < elen; el++) {
	Element *pel = mesh.elist[el];
	int nn = pel->nNode();
	int *nd = pel->Node;
	for (i = 0; i < nn; i++)
	    ngno[nd[i]] += egno[el];
    }
    int *rowptr = new int[nlen+1];
    int *rowp = new int[nlen];
    rowptr[0] = 0;
    for (i = 0; i < nlen; i++)
	rowptr[i+1] = rowptr[i]+ngno[i];
    memcpy (rowp, rowptr, nlen*sizeof(int));
    nzero = rowptr[nlen];
    int *colidx = new int[nzero];
    double *data = new double[nzero];

    double dx = ((*bbmax)[0]-(*bbmin)[0])/(double)(nx-1);
    double dy = ((*bbmax)[1]-(*bbmin)[1])/(double)(ny-1);
    double dz = (dim > 2 ? ((*bbmax)[2]-(*bbmin)[2])/(double)(nz-1) : 0);

    for (iz = i = 0; iz < nz; iz++) {
	if (dim > 2) p[2] = (*bbmin)[2] + (double)iz*dz;
	for (iy = 0; iy < ny; iy++) {
	    p[1] = (*bbmin)[1] + (double)iy*dy;
	    for (ix = 0; ix < nx; ix++) {
		p[0] = (*bbmin)[0] + (double)ix*dx;
		el = elref[i];
		if (el >= 0) {
		    Element *pel = mesh.elist[el];
		    if (!pel->GContains(p, mesh.nlist))
			cerr << "Inconsistent element reference!" << endl;
		    RVector fun = pel->GlobalShapeF (mesh.nlist, p);
		    for (j = 0; j < pel->nNode(); j++) {
			n = pel->Node[j];
			scal[n] += fun[j];
			data[rowp[n]] = fun[j];
			colidx[rowp[n]] = i;
			rowp[n]++;
		    }
		}
		i++;
	    }
	}
    }
    for (n = 0; n < nlen; n++) {
	double s = scal[n];
	if (!s) cerr << "no support for node " << n << endl;
	for (i = rowptr[n]; i < rowptr[n+1]; i++) {
	    if (s) data[i] /= s;
	    else data[i] = 0;
	}
    }

    RCompRowMatrix *BI = new RCompRowMatrix (nlen, rlen,
        rowptr, colidx, data);
    delete []rowptr;
    delete []colidx;
    delete []data;
    delete []egno;
    delete []ngno;
    delete []scal;
    delete []rowp;
    return BI;
}

int *GenerateElementPixelRef (const Mesh &mesh, const IVector &gdim,
    const Point *bbmin, const Point *bbmax)
{
    const double EPS = 1e-6;
    int i, n, ix, iy, iz, idx, el, nnode, *node;
    int d, dim = mesh.Dimension();
    Point p(dim);

    // Sanity checks
    for (i = 0; i < dim; i++) {
	xASSERT(gdim[i] >= 2, Invalid grid dimension);
	xASSERT((*bbmax)[i] > (*bbmin)[i], Invalid bounding box);
    }

    int nx = gdim[0], ny = gdim[1], nz = (dim > 2 ? gdim[2] : 1);
    int rlen = nx*ny*nz;

    int *elref = new int[rlen];
    for (i = 0; i < rlen; i++) elref[i] = -1;

    // grid size and spacing
    RVector dW(dim);
    for (d = 0; d < dim; d++)
        dW[d] = ((*bbmax)[d]-(*bbmin)[d])/(double)(gdim[d]-1);

    // associate elements with grid points
    Point elbbmin(dim), elbbmax(dim);
    for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	node  = mesh.elist[el]->Node;

        // compute element bounding box
        for (d = 0; d < dim; d++) elbbmin[d] = 1e10, elbbmax[d] = -1e10;
	for (i = 0; i < nnode; i++) {
	    Point &nd = mesh.nlist[node[i]];
	    for (d = 0; d < dim; d++) {
	        if (nd[d] < elbbmin[d]) elbbmin[d] = nd[d];
		if (nd[d] > elbbmax[d]) elbbmax[d] = nd[d];
	    }
	}
	for (d = 0; d < dim; d++) { // inflate a bit to avoid edge problems
	    elbbmin[d] -= EPS;
	    elbbmax[d] += EPS;
	}

	// search all grid points inside the element bounding box
	int rmin[3], rmax[3];
	for (d = 0; d < dim; d++) {
	    if (elbbmin[d] > (*bbmax)[d]) goto no_overlap;
	    if (elbbmax[d] < (*bbmin)[d]) goto no_overlap;
	    rmin[d] = std::max (0, (int)ceil ((elbbmin[d]-(*bbmin)[d])/dW[d]));
	    rmax[d] = std::min (gdim[d], (int)floor((elbbmax[d]-(*bbmin)[d])/dW[d])+1);
	}
	if (dim < 3) rmin[2] = 0, rmax[2] = 1;
	for (iz = rmin[2]; iz < rmax[2]; iz++) {
	    if (dim > 2) p[2] = iz*dW[2] + (*bbmin)[2];
	    for (iy = rmin[1]; iy < rmax[1]; iy++) {
	        p[1] = iy*dW[1] + (*bbmin)[1];
		for (ix = rmin[0]; ix < rmax[0]; ix++) {
		    idx = ix + nx*(iy + ny*iz);
		    if (elref[idx] < 0) {
		        p[0] = ix*dW[0] + (*bbmin)[0];
			if (mesh.elist[el]->GContains (p, mesh.nlist))
			    elref[idx] = el;
		    }
		}
	    }
	}
    no_overlap:;
    }
    return elref;
}

void SubsampleLinPixel (const RVector &gf, RVector &bf, const IVector &gdim,
    const IVector &bdim, const int *elref)
{
    int ix, iy, iz, gidx;
    int gx0, gx1, gy0, gy1, gz0, gz1, xx, yy, zz;
    int gx = gdim[0], gy = gdim[1], gz = (gdim.Dim() > 2 ? gdim[2]:1);
    int bx = bdim[0], by = bdim[1], bz = (bdim.Dim() > 2 ? bdim[2]:1);
    double nbx = bx-1.0, nby = by-1.0, nbz = bz-1.0;
    double x, y, z, x0, x1, y0, y1, z0, z1;
    double sum, wgt, dx, dy, dz, idyz, idy, idz, t;
    double dbx = 1.0/nbx, dby = 1.0/nby, dbz;
    double dgx = 1.0/(gx-1.0), dgy = 1.0/(gy-1.0), dgz;
    if (bz > 1) {
        dbz = 1.0/nbz;
	dgz = 1.0/(gz-1.0);
    }

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz ? (iz-1)*dbz:0.0);
	    z1 = (iz < bz-1 ? (iz+1)*dbz:1.0);
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy ? (iy-1)*dby:0.0);
	    y1 = (iy < by-1 ? (iy+1)*dby:1.0);
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x  = ix*dbx;
		x0 = (ix ? (ix-1)*dbx:0.0);
		x1 = (ix < bx-1 ? (ix+1)*dbx:1.0);
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));
		sum = wgt = 0.0;

		if (bz == 1) { // 2D case

		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 1.0) continue;
			idy = 1.0-dy;
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    // mesh support check
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 1.0) continue;
			    t   = (1.0-dx)*idy;
			    wgt += t;
			    sum += t*gf[gidx];
			}
		    }
		    if (wgt) bf[ix+iy*bx] = sum/wgt;
		    else bf[ix+iy*bx] = 0.0; // rather, pixel should be skipped

		} else {  // 3D case

		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 1.0) continue;
			idz = 1.0-dz;
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 1.0) continue;
			    idyz = (1.0-dy)*idz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
				if (elref && elref[gidx] < 0) continue;
				// mesh support check
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 1.0) continue;
				t   = (1.0-dx)*idyz;
				wgt += t;
				sum += t*gf[gidx];
			    }
			}
		    }
		    if (wgt) bf[ix+(iy+iz*by)*bx] = sum/wgt;
		    else bf[ix+(iy+iz*by)*bx] = 0.0;

		}
	    }
	}
    }
}

RGenericSparseMatrix *Grid2LinPixMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref)
{
    int i, ix, iy, iz, gidx, idx, idx0, rowentry, nzero;
    int gx0, gx1, gy0, gy1, gz0, gz1, xx, yy, zz;
    int gx = gdim[0], gy = gdim[1], gz = (gdim.Dim() > 2 ? gdim[2]:1);
    int bx = bdim[0], by = bdim[1], bz = (bdim.Dim() > 2 ? bdim[2]:1);
    int glen = gx*gy*gz, blen = bx*by*bz;
    //bool nop = true;
    //for (i = 0; i < gdim.Dim(); i++)
    //	if (gdim[i] != bdim[i]) { nop = false; break; }

    //if (nop) { // same grid dimensions -> return unit matrix
    //	RCompRowMatrix *B2 = new RCompRowMatrix;
    //	B2->Identity (blen);
    //	return B2;
    //}

    double nbx = bx-1.0, nby = by-1.0, nbz = bz-1.0;
    double x, y, z, x0, x1, y0, y1, z0, z1;
    double wgt, dx, dy, dz, idyz, idy, idz, t;
    double dbx = 1.0/nbx, dby = 1.0/nby, dbz;
    double dgx = 1.0/(gx-1.0), dgy = 1.0/(gy-1.0), dgz;
    if (bz > 1) {
        dbz = 1.0/nbz;
	dgz = 1.0/(gz-1.0);
    }

    // pass 1: find sparsity structure of map matrix
    int *rowptr = new int[blen+1];
    rowptr[idx=0] = 0;

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz ? (iz-1)*dbz:0.0);
	    z1 = (iz < bz-1 ? (iz+1)*dbz:1.0);
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy ? (iy-1)*dby:0.0);
	    y1 = (iy < by-1 ? (iy+1)*dby:1.0);
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x  = ix*dbx;
		x0 = (ix ? (ix-1)*dbx:0.0);
		x1 = (ix < bx-1 ? (ix+1)*dbx:1.0);
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));
		rowentry = 0;

		if (bz == 1) { // 2D case

		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 1.0) continue;
			idy = 1.0-dy;
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    // mesh support check
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 1.0) continue;
			    rowentry++;
			}
		    }

		} else {  // 3D case

		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 1.0) continue;
			idz = 1.0-dz;
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 1.0) continue;
			    idyz = (1.0-dy)*idz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
				if (elref && elref[gidx] < 0) continue;
				// mesh support check
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 1.0) continue;
				rowentry++;
			    }
			}
		    }
		}
		rowptr[idx+1] = rowptr[idx] + rowentry;
		idx++;
	    }
	}
    }
    nzero = rowptr[blen];
    cerr << "Grid2LinPixMatrix: found " << nzero << " nonzeros" << endl;

    // pass 2: construct map matrix
    int *colidx = new int[nzero];
    double *val = new double[nzero];

    // loop over coarse pixel grid
    for (iz = idx = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz ? (iz-1)*dbz:0.0);
	    z1 = (iz < bz-1 ? (iz+1)*dbz:1.0);
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy ? (iy-1)*dby:0.0);
	    y1 = (iy < by-1 ? (iy+1)*dby:1.0);
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x  = ix*dbx;
		x0 = (ix ? (ix-1)*dbx:0.0);
		x1 = (ix < bx-1 ? (ix+1)*dbx:1.0);
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));
		idx0 = idx;
		wgt = 0.0;

		if (bz == 1) { // 2D case

		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 1.0) continue;
			idy = 1.0-dy;
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    // mesh support check
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 1.0) continue;
			    t   = (1.0-dx)*idy;
			    colidx[idx] = gidx;
			    val[idx]    = t;
			    wgt += t;
			    idx++;
			}
		    }

		} else {  // 3D case

		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 1.0) continue;
			idz = 1.0-dz;
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 1.0) continue;
			    idyz = (1.0-dy)*idz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
				if (elref && elref[gidx] < 0) continue;
				// mesh support check
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 1.0) continue;
				t   = (1.0-dx)*idyz;
				colidx[idx] = gidx;
				val[idx]    = t;
				wgt += t;
				idx++;
			    }
			}
		    }
		}

		for (i = idx0; i < idx; i++)
		    if (wgt) val[i] /= wgt;
		    else val[i] = 0.0;
	    }
	}
    }
    RCompRowMatrix *B2 = new RCompRowMatrix (blen, glen, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;
    return B2;
}

double spline (double x)
{
    if (fabs (x) < 2) {
	return 1.0/12.0 * 
	    (pow (fabs(x-2), 3.0) - 
	     4.0 * pow(fabs(x-1), 3.0) +
	     6.0 * pow(fabs(x), 3.0) -
	     4.0 * pow(fabs(x+1), 3.0) +
	     pow (fabs(x+2), 3.0));
    } else {
	return 0.0;
    }
}

RGenericSparseMatrix *CubicPix2GridMatrix (const IVector &bdim,
    const IVector &gdim, const int *elref)
{
    // return transformation matrix to convert an image represented by
    // coefficients of cubic spline basis expansion into a bitmap grid

    bool is3d = (bdim.Dim() > 2);

    int bx = bdim[0], by = bdim[1], bz = (is3d ? bdim[2]:1);
    int gx = gdim[0], gy = gdim[1], gz = (is3d ? gdim[2]:1);
    int nb = bx*by*bz, ng = gx*gy*gz;

    double dx = (gx-1.0)/(bx-3.0);
    double dy = (gy-1.0)/(by-3.0);
    double dz = (is3d ? (gz-1.0)/(bz-3.0): 1.0);

    double x0, y0, z0, sx, sy, sz, v;
    int i, ix, iy, iz, xm, xp, ym, yp, zm, zp, x, y, z;
    IVector rowentry(ng);

    // pass 1: generate matrix fill-in structure
    for (iz = 0; iz < bz; iz++) {
	if (is3d) {
	    z0 = (iz-1)*dz;
	    zm = ::max (0, (int)((iz-3)*dz));
	    zp = ::min (gz-1, (int)((iz+1)*dz));
	} else {
	    zm = 0, zp = 0;
	}
	for (iy = 0; iy < by; iy++) {
	    y0 = (iy-1.0)*dy;
	    ym = ::max (0, (int)((iy-3)*dy));
	    yp = ::min (gy-1, (int)((iy+1)*dy));
	    for (ix = 0; ix < bx; ix++) {
		x0 = (ix-1.0)*dx; // pixel cnt in bitmap coords
		xm = ::max (0, (int)((ix-3)*dx));
		xp = ::min (gx-1, (int)((ix+1)*dx));

		for (z = zm; z <= zp; z++) {
		    sz = (is3d ? spline ((z-z0)/dz) : 1.0);
		    for (y = ym; y <= yp; y++) {
			sy = spline((y-y0)/dy);
			for (x = xm; x <= xp; x++) {
			    sx = spline((x-x0)/dx);
			    v = sx*sy*sz;
			    if (v) rowentry[x+y*gx+z*gx*gy]++;
			}
		    }
		}
	    }
	}
    }

    int nzero;
    int *rowptr = new int[ng+1];
    rowptr[0] = 0;
    for (i = 0; i < ng; i++) rowptr[i+1] = rowptr[i]+rowentry[i];
    nzero = rowptr[ng];
    int *colidx = new int[nzero];
    double *val = new double[nzero];

    // pass 2: generate matrix
    rowentry.Clear();
    for (iz = 0; iz < bz; iz++) {
	if (is3d) {
	    z0 = (iz-1)*dz;
	    zm = ::max (0, (int)((iz-3)*dz));
	    zp = ::min (gz-1, (int)((iz+1)*dz));
	} else {
	    zm = 0, zp = 0;
	}
	for (iy = 0; iy < by; iy++) {
	    y0 = (iy-1.0)*dy;
	    ym = ::max (0, (int)((iy-3)*dy));
	    yp = ::min (gy-1, (int)((iy+1)*dy));
	    for (ix = 0; ix < bx; ix++) {
		x0 = (ix-1.0)*dx; // pixel cnt in bitmap coords
		xm = ::max (0, (int)((ix-3)*dx));
		xp = ::min (gx-1, (int)((ix+1)*dx));

		for (z = zm; z <= zp; z++) {
		    sz = (is3d ? spline ((z-z0)/dz) : 1.0);
		    for (y = ym; y <= yp; y++) {
			sy = spline((y-y0)/dy);
			for (x = xm; x <= xp; x++) {
			    sx = spline((x-x0)/dx);
			    v = sx*sy*sz;
			    if (v) {
				int row = x+y*gx+z*gx*gy;
				int idx = rowptr[row] + rowentry[row];
				colidx[idx] = ix + iy*bx + iz*bx*by;
				val[idx] = v;
				rowentry[row]++;
			    }
			}
		    }
		}
	    }
	}
    }

    RCompRowMatrix *R = new RCompRowMatrix (ng, nb, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;
    return R;
}

#ifdef UNDEF
RGenericSparseMatrix *Grid2CubicPixMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref)
{
}
#endif

void RasterLinPixel (RVector &gf, const RVector &bf, const IVector &gdim,
    const IVector &bdim, const int *elref)
{
    int ix, iy, iz, gidx;
    int gx0, gx1, gy0, gy1, gz0, gz1, xx, yy, zz;
    int gx = gdim[0], gy = gdim[1], gz = (gdim.Dim() > 2 ? gdim[2]:1);
    int bx = bdim[0], by = bdim[1], bz = (bdim.Dim() > 2 ? bdim[2]:1);
    double nbx = bx-1.0, nby = by-1.0, nbz = bz-1.0;
    double x, y, z, x0, x1, y0, y1, z0, z1;
    double val, dx, dy, dz, idyz, idy, idz, t;
    double dbx = 1.0/nbx, dby = 1.0/nby, dbz;
    double dgx = 1.0/(gx-1.0), dgy = 1.0/(gy-1.0), dgz;
    if (bz > 1) {
        dbz = 1.0/nbz;
	dgz = 1.0/(gz-1.0);
    }
    gf.Clear(); // zero fine pixel grid

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz ? (iz-1)*dbz:0.0);
	    z1 = (iz < bz-1 ? (iz+1)*dbz:1.0);
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy ? (iy-1)*dby:0.0);
	    y1 = (iy < by-1 ? (iy+1)*dby:1.0);
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x  = ix*dbx;
		x0 = (ix ? (ix-1)*dbx:0.0);
		x1 = (ix < bx-1 ? (ix+1)*dbx:1.0);
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));

		if (bz == 1) { // 2D case

		    val = bf[ix+iy*bx];
		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 1.0) continue;
			idy = 1.0-dy;
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 1.0) continue;
			    t   = (1.0-dx)*idy;
			    gf[gidx] += t*val;
			}
		    }

		} else { // 3D case

		    val = bf[ix+(iy+iz*by)*bx];
		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 1.0) continue;
			idz = 1.0-dz;
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 1.0) continue;
			    idyz = (1.0-dy)*idz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
				if (elref && elref[gidx] < 0) continue;
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 1.0) continue;
				t   = (1.0-dx)*idyz;
				gf[gidx] += t*val;
			    }
			}
		    }
		}
	    }
	}
    }
}

RGenericSparseMatrix *LinPix2GridMatrix (const IVector &gdim,
    const IVector &bdim, const int *elref)
{
    int i, ix, iy, iz, gidx, nzero, idx, c;
    int gx0, gx1, gy0, gy1, gz0, gz1, xx, yy, zz;
    int gx = gdim[0], gy = gdim[1], gz = (gdim.Dim() > 2 ? gdim[2]:1);
    int bx = bdim[0], by = bdim[1], bz = (bdim.Dim() > 2 ? bdim[2]:1);
    int glen = gx*gy*gz, blen = bx*by*bz;
    //bool nop = true;
    //for (i = 0; i < gdim.Dim(); i++)
    //	if (gdim[i] != bdim[i]) { nop = false; break; }
    //
    //if (nop) { // same grid dimensions -> return unit matrix
    //	RCompRowMatrix *IB2 = new RCompRowMatrix;
    //	IB2->Identity (blen);
    //	return IB2;
    //}

    double nbx = bx-1.0, nby = by-1.0, nbz = bz-1.0;
    double x, y, z, x0, x1, y0, y1, z0, z1;
    double dx, dy, dz, idyz, idy, idz;
    double dbx = 1.0/nbx, dby = 1.0/nby, dbz;
    double dgx = 1.0/(gx-1.0), dgy = 1.0/(gy-1.0), dgz;
    if (bz > 1) {
        dbz = 1.0/nbz;
	dgz = 1.0/(gz-1.0);
    }

    // pass 1: find sparsity structure of map matrix
    int *rowlen = new int[glen];
    for (i = 0; i < glen; i++) rowlen[i] = 0;

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz ? (iz-1)*dbz:0.0);
	    z1 = (iz < bz-1 ? (iz+1)*dbz:1.0);
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy ? (iy-1)*dby:0.0);
	    y1 = (iy < by-1 ? (iy+1)*dby:1.0);
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x  = ix*dbx;
		x0 = (ix ? (ix-1)*dbx:0.0);
		x1 = (ix < bx-1 ? (ix+1)*dbx:1.0);
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));

		if (bz == 1) { // 2D case

		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 1.0) continue;
			idy = 1.0-dy;
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 1.0) continue;
			    rowlen[gidx]++;
			}
		    }

		} else { // 3D case

		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 1.0) continue;
			idz = 1.0-dz;
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 1.0) continue;
			    idyz = (1.0-dy)*idz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
				if (elref && elref[gidx] < 0) continue;
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 1.0) continue;
				rowlen[gidx]++;
			    }
			}
		    }
		}
	    }
	}
    }
    int *rowptr = new int[glen+1];
    rowptr[0] = 0;
    for (i = 0; i < glen; i++)
        rowptr[i+1] = rowptr[i] + rowlen[i];
    nzero = rowptr[glen];
    cerr << "LinPix2GridMatrix: found " << nzero << " nonzeros" << endl;

    // pass 2: construct map matrix
    int *colidx = new int[nzero];
    double *val = new double[nzero];
    for (i = 0; i < glen; i++) rowlen[i] = 0;

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz ? (iz-1)*dbz:0.0);
	    z1 = (iz < bz-1 ? (iz+1)*dbz:1.0);
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy ? (iy-1)*dby:0.0);
	    y1 = (iy < by-1 ? (iy+1)*dby:1.0);
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x  = ix*dbx;
		x0 = (ix ? (ix-1)*dbx:0.0);
		x1 = (ix < bx-1 ? (ix+1)*dbx:1.0);
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));

		if (bz == 1) { // 2D case

		    c = ix+iy*bx;
		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 1.0) continue;
			idy = 1.0-dy;
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 1.0) continue;
			    idx = rowptr[gidx]+rowlen[gidx]++;
			    colidx[idx] = c;
			    val[idx] = (1.0-dx)*idy;
			}
		    }

		} else { // 3D case

		    c = ix+(iy+iz*by)*bx;
		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 1.0) continue;
			idz = 1.0-dz;
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 1.0) continue;
			    idyz = (1.0-dy)*idz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
				if (elref && elref[gidx] < 0) continue;
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 1.0) continue;
				idx = rowptr[gidx]+rowlen[gidx]++;
				colidx[idx] = c;
				val[idx] = (1.0-dx)*idyz;
			    }
			}
		    }
		}
	    }
	}
    }

    RCompRowMatrix *IB2 = new RCompRowMatrix (glen, blen, rowptr, colidx, val);

    delete []rowlen;
    delete []rowptr;
    delete []colidx;
    delete []val;
    return IB2;
}

void SubsampleCubPixel (const RVector &gf, RVector &bf, const IVector &gdim,
    const IVector &bdim, const int *elref)
{
    int ix, iy, iz, gidx;
    int gx0, gx1, gy0, gy1, gz0=0, gz1=0, xx, yy, zz;
    int gx = gdim[0], gy = gdim[1], gz = (gdim.Dim() > 2 ? gdim[2]:1);
    int bx = bdim[0], by = bdim[1], bz = (bdim.Dim() > 2 ? bdim[2]:1);
    double nbx = bx-1.0,  nby = by-1.0,  nbz = bz-1.0;
    double dbx = 1.0/nbx, dby = 1.0/nby, dbz;
    double dgx = 1.0/(gx-1.0), dgy = 1.0/(gy-1.0), dgz;
    double x, y, z, x0, x1, y0, y1, z0, z1, dx, dy, dz, sx, sy, sz, syz;
    double t, sum, wgt, *sx_store;
    if (bz > 1) {
        dbz = 1.0/nbz;
	dgz = 1.0/(gz-1.0);
	sx_store = new double[(int)(dbx*4.0*gx)+4];
    }

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz-2)*dbz; if (z0 < 0.0) z0 = 0.0;
	    z1 = (iz+2)*dbz; if (z1 > 1.0) z1 = 1.0;
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy-2)*dby; if (y0 < 0.0) y0 = 0.0;
	    y1 = (iy+2)*dby; if (y1 > 1.0) y1 = 1.0;
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x = ix*dbx;
		x0 = (ix-2)*dbx; if (x0 < 0.0) x0 = 0.0;
		x1 = (ix+2)*dbx; if (x1 > 1.0) x1 = 1.0;
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));
		sum = wgt = 0.0;

		if (bz == 1) { // 2D case

		    // loop over fine grid support of pixel (ix,iy,iz)
		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 2.0) continue;
			else if (dy <= 1.0)
			    sy = 1.0 + dy*dy*(-2.5 + 1.5*dy);
			else
			    sy = 2.0 + dy*(-4.0 + dy*(2.5 - 0.5*dy));
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    // mesh support check
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 2.0) continue;
			    else if (dx <= 1.0)
			        sx = 1.0 + dx*dx*(-2.5 + 1.5*dx);
			    else
			        sx = 2.0 + dx*(-4.0 + dx*(2.5 - 0.5*dx));
			    t = sx*sy;
			    wgt += t;
			    sum += t*gf[gidx];
			}
		    }
		    if (wgt) bf[ix+iy*bx] = sum/wgt;
		    else bf[ix+iy*bx] = 0.0; // rather, pixel should be skipped

		} else { // 3D case

		    // pre-calculate the weights for the inner (x) loop
		    // for efficiency

		    for (xx = gx0; xx <= gx1; xx++) {
		        dx = fabs (xx*dgx - x)*nbx;
			if (dx >= 2.0)
			    sx = 0.0;
			else if (dx <= 1.0)
			    sx = 1.0 + dx*dx*(-2.5 + 1.5*dx);
			else
			    sx = 2.0 + dx*(-4.0 + dx*(2.5 - 0.5*dx));
			sx_store[xx-gx0] = sx;
		    }

		    // loop over fine grid support of pixel (ix,iy,iz)
		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 2.0) continue;
			else if (dz < 1.0)
			    sz = 1.0 + dz*dz*(-2.5 + 1.5*dz);
			else
			    sz = 2.0 + dz*(-4.0 + dz*(2.5 - 0.5*dz));
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 2.0) continue;
			    else if (dy <= 1.0)
			        sy = 1.0 + dy*dy*(-2.5 + 1.5*dy);
			    else
			        sy = 2.0 + dy*(-4.0 + dy*(2.5 - 0.5*dy));
			    syz = sy*sz;
			    for (xx = gx0, gidx = gx0 + (yy+zz*gy)*gx;
				 xx <= gx1; xx++, gidx++) {
			        //gidx = xx + (yy+zz*gy)*gx;
				if (elref && elref[gidx] < 0) continue;
				// mesh support check
				//dx = fabs (xx*dgx - x)*nbx;
				//if (dx >= 2.0) continue;
				//else if (dx <= 1.0)
				//    sx = 1.0 + dx*dx*(-2.5 + 1.5*dx);
				//else
				//    sx = 2.0 + dx*(-4.0 + dx*(2.5 - 0.5*dx));
				//t = sx*syz;
				t = sx_store[xx-gx0]*syz;
				wgt += t;
				sum += t*gf[gidx];
			    }
			}
		    }
		    if (wgt) bf[ix+(iy+iz*by)*bx] = sum/wgt;
		    else bf[ix+(iy+iz*by)*bx] = 0.0;
		}
	    }
	}
    }
    if (bz > 1) delete []sx_store;
}

void RasterCubPixel (RVector &gf, const RVector &bf, const IVector &gdim,
    const IVector &bdim, const int *elref)
{
    int ix, iy, iz, gidx;
    int gx0, gx1, gy0, gy1, gz0=0, gz1=0, xx, yy, zz;
    int gx = gdim[0], gy = gdim[1], gz = (gdim.Dim() > 2 ? gdim[2]:1);
    int bx = bdim[0], by = bdim[1], bz = (bdim.Dim() > 2 ? bdim[2]:1);
    double nbx = bx-1.0,  nby = by-1.0,  nbz = bz-1.0;
    double dbx = 1.0/nbx, dby = 1.0/nby, dbz;
    double dgx = 1.0/(gx-1.0), dgy = 1.0/(gy-1.0), dgz;
    double x, y, z, x0, x1, y0, y1, z0, z1, dx, dy, dz, sx, sy, sz;
    double t, val;
    if (bz > 1) {
        dbz = 1.0/nbz;
	dgz = 1.0/(gz-1.0);
    }
    gf.Clear(); // zero fine pixel grid

    // loop over coarse pixel grid
    for (iz = 0; iz < bz; iz++) {
        if (bz > 1) {
	    z = (double)iz*dbz;
	    z0 = (iz-2)*dbz; if (z0 < 0.0) z0 = 0.0;
	    z1 = (iz+2)*dbz; if (z1 > 1.0) z1 = 1.0;
	    gz0 = (int)floor (z0*(gz-1));
	    gz1 = (int)ceil  (z1*(gz-1));
	}
	for (iy = 0; iy < by; iy++) {
	    y = (double)iy*dby;
	    y0 = (iy-2)*dby; if (y0 < 0.0) y0 = 0.0;
	    y1 = (iy+2)*dby; if (y1 > 1.0) y1 = 1.0;
	    gy0 = (int)floor (y0*(gy-1));
	    gy1 = (int)ceil  (y1*(gy-1));
	    for (ix = 0; ix < bx; ix++) {
	        x = ix*dbx;
		x0 = (ix-2)*dbx; if (x0 < 0.0) x0 = 0.0;
		x1 = (ix+2)*dbx; if (x1 > 1.0) x1 = 1.0;
		gx0 = (int)floor (x0*(gx-1));
		gx1 = (int)ceil  (x1*(gx-1));

		if (bz == 1) { // 2D case

		    val = bf[ix+iy*bx];
		    // loop over fine grid support of pixel (ix,iy,iz)
		    for (yy = gy0; yy <= gy1; yy++) {
		        dy = fabs (yy*dgy - y)*nby;
			if (dy >= 2.0) continue;
			else if (dy <= 1.0)
			    sy = 1.0 + dy*dy*(-2.5 + 1.5*dy);
			else
			    sy = 2.0 + dy*(-4.0 + dy*(2.5 - 0.5*dy));
			for (xx = gx0; xx <= gx1; xx++) {
			    gidx = xx + yy*gx;
			    if (elref && elref[gidx] < 0) continue;
			    dx = fabs (xx*dgx - x)*nbx;
			    if (dx >= 2.0) continue;
			    else if (dx <= 1.0)
			        sx = 1.0 + dx*dx*(-2.5 + 1.5*dx);
			    else
			        sx = 2.0 + dx*(-4.0 + dx*(2.5 - 0.5*dx));
			    t = sx*sy;
			    gf[gidx] += t*val;
			}
		    }

		} else { // 3D case

		    val = bf[ix+(iy+iz*by)*bx];
		    // loop over fine grid support of pixel (ix,iy,iz)
		    for (zz = gz0; zz <= gz1; zz++) {
		        dz = fabs (zz*dgz - z)*nbz;
			if (dz >= 2.0) continue;
			else if (dz < 1.0)
			    sz = 1.0 + dz*dz*(-2.5 + 1.5*dz);
			else
			    sz = 2.0 + dz*(-4.0 + dz*(2.5 - 0.5*dz));
			for (yy = gy0; yy <= gy1; yy++) {
			    dy = fabs (yy*dgy - y)*nby;
			    if (dy >= 2.0) continue;
			    else if (dy <= 1.0)
			        sy = 1.0 + dy*dy*(-2.5 + 1.5*dy);
			    else
			        sy = 2.0 + dy*(-4.0 + dy*(2.5 - 0.5*dy));
			    for (xx = gx0; xx <= gx1; xx++) {
			        gidx = xx + (yy+zz*gy)*gx;
				if (elref && elref[gidx] < 0) continue;
				dx = fabs (xx*dgx - x)*nbx;
				if (dx >= 2.0) continue;
				else if (dx <= 1.0)
				    sx = 1.0 + dx*dx*(-2.5 + 1.5*dx);
				else
				    sx = 2.0 + dx*(-4.0 + dx*(2.5 - 0.5*dx));
				t = sx*sy*sz;
				gf[gidx] += t*val;
			    }
			}
		    }
		}
	    }
	}
    }
}

void CreateVoxelMesh (int ex, int ey, int ez, bool *egrid,
		      double dx, double dy, double dz, Mesh &mesh)
{
    // generate node grid
    int i, j, k, n, nel, nnd, idx, nd[8];
    int nx = ex+1, ny = ey+1, nz = ez+1; // why is this true?
    int *ngrid = new int[nx*ny*nz];
    for (i = 0; i < nx*ny*nz; i++) ngrid[i] = -1;

    // find number of elements and generate element list
    for (nel = i = 0; i < ex*ey*ez; i++)
        if (egrid[i]) nel++;
    cerr << "number of elements: " << nel << endl;
    Element **elist = new Element*[nel];
    for (i = 0; i < nel; i++)
        elist[i] = new Voxel8;
    mesh.elist.SetList (nel, elist);

    for (k = nel = nnd = 0; k < ez; k++) {
        for (j = 0; j < ey; j++) {
	    for (i = 0; i < ex; i++) {
	        if (egrid[i + j*ex + k*ex*ey]) {
		    nd[0] = i + j*nx + k*nx*ny;
		    nd[1] = nd[0]+1;
		    nd[2] = nd[0]+nx;
		    nd[3] = nd[2]+1;
		    nd[4] = nd[0]+nx*ny;
		    nd[5] = nd[4]+1;
		    nd[6] = nd[4]+nx;
		    nd[7] = nd[6]+1;
		    for (n = 0; n < 8; n++) {
		        if (ngrid[nd[n]] < 0) ngrid[nd[n]] = nnd++;
			elist[nel]->Node[n] = ngrid[nd[n]];
		    }
		    nel++;
		}
	    }
	}
    }

    // generate node list
    mesh.nlist.New (nnd);
    for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++) {
	        if ((idx = ngrid[i + j*nx + k*nx*ny]) >= 0) {
		    if (idx >= nnd) cerr << "Problems here!" << endl;
		    mesh.nlist[idx].New(3);
		    mesh.nlist[idx][0] = i*dx;
		    mesh.nlist[idx][1] = j*dy;
		    mesh.nlist[idx][2] = k*dz;
		}
	    }
	}
    }

    // generate parameter list
    mesh.plist.New (nnd);
    
    // label boundary nodes
    mesh.MarkBoundary ();

    // mesh setup
    mesh.Setup ();

    // cleanup
    delete []ngrid;
}
