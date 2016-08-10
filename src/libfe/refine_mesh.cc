#define FELIB_IMPLEMENTATION

#include "felib.h"
#include "toastdef.h"
#include "refine_mesh.h"

struct EDGE {   // structure for the edge of a 2-D element
    int nd[2];      // node indices
    int el0, el1;   // attached element indices (-1 for none)
    int sd0, sd1;   // corresponding local element side indices
    bool refine;    // flag for edge refinement
    int mnd;        // node index for midpoint node after edge refinement
};

int edgesort_comp (const void *arg1, const void *arg2)
{
    int k;
    const EDGE *e1 = (EDGE*)arg1;
    const EDGE *e2 = (EDGE*)arg2;
    for (k = 0; k < 2; k++) {
	if (e1->nd[k] < e2->nd[k]) return -1;
	else if (e1->nd[k] > e2->nd[k]) return 1;
    }
    return 0;
}

#ifdef UNDEF
void RefineTriangle3Mesh (Mesh *mesh, bool *elrefine)
{
    // Create global edge list
    int i, j, k, l, idx, tmp;
    int nel = mesh->elen();
    int nnd = mesh->nlen();
    int nedge_full = nel*3;
    EDGE *edge_full = new EDGE[nedge_full];

    for (i = idx = 0; i < nel; i++) {
	Element *pel = mesh->elist[i];
	int nsd = pel->nSide();
	for (j = 0; j < nsd; j++) {
	    int nsdnd = pel->nSideNode(j);
	    for (k = 0; k < nsdnd; k++)
		edge_full[idx].nd[k] = pel->Node[pel->SideNode(j,k)];
	    for (k = 0; k < nsdnd-1; k++) {
		for (l = k+1; l < nsdnd; l++) {
		    if (edge_full[idx].nd[k] > edge_full[idx].nd[l]) {
			tmp = edge_full[idx].nd[k];
			edge_full[idx].nd[k] = edge_full[idx].nd[l];
			edge_full[idx].nd[l] = tmp;
		    }
		}
	    }
	    edge_full[idx].el0 = i;
	    edge_full[idx].sd0 = j;
	    idx++;
	}
    }

    // order the edge list according to node indices
    qsort (edge_full, nedge_full, sizeof(EDGE), edgesort_comp);

    // now all non-boundary sides should be listed as consecutive pairs,
    // so condense all matching pairs into single entries
    int nedge = nedge_full;
    for (i = 0; i < nedge_full-1; i++) {
	if (edge_full[i].nd[0] == edge_full[i+1].nd[0] && 
	    edge_full[i].nd[1] == edge_full[i+1].nd[1])
	    nedge--;
    }

    EDGE *edge = new EDGE[nedge];
    for (i = idx = 0; i < nedge_full; i++) {
	edge[idx].nd[0] = edge_full[i].nd[0];
	edge[idx].nd[1] = edge_full[i].nd[1];
	edge[idx].el0 = edge_full[i].el0;
	edge[idx].sd0 = edge_full[i].sd0;
	if (i < nedge_full-1 &&
	  edge_full[i].nd[0] == edge_full[i+1].nd[0] && 
	  edge_full[i].nd[1] == edge_full[i+1].nd[1]) {
	    edge[idx].el1 = edge_full[i+1].el0;
	    edge[idx].sd1 = edge_full[i+1].sd0;
	    i++;
	} else {
	    edge[idx].el1 = -1;
	    edge[idx].sd1 = -1;
	}
	edge[idx].refine = false;
	idx++;
    }
    delete []edge_full;

    // create element->edge reference list
    typedef int TRI_EDGE[3];
    TRI_EDGE *el_edge = new TRI_EDGE[nel];

    for (i = 0; i < nedge; i++) {
	el_edge[edge[i].el0][edge[i].sd0] = i;
	if (edge[i].el1 >= 0)
	    el_edge[edge[i].el1][edge[i].sd1] = i;
    }

    // check which edges should be refined
    for (i = 0; i < nel; i++) {
	Element *pel = mesh->elist[i];
	int nsd = pel->nSide();
	if (elrefine[i]) {
	    for (j = 0; j < nsd; j++) {
		edge[el_edge[i][j]].refine = true;
	    }
	}
    }
    int nrefine = 0;
    for (i = 0; i < nedge; i++) {
	if (edge[i].refine) nrefine++;
    }

    // go through the list and check if any triangles try to refine
    // two of their sides, or if any single-side refinement is not refining
    // the longest side
    bool need_subdiv = false;
    do {
	need_subdiv = false;
	for (i = 0; i < nel; i++) {
	    Element *pel = mesh->elist[i];
	    int nsd = pel->nSide();
	    int nref = 0;
	    for (j = 0; j < nsd; j++)
		if (edge[el_edge[i][j]].refine) nref++;
	    if (nref == 2) {
		for (j = 0; j < nsd; j++)
		    edge[el_edge[i][j]].refine = true;
		nrefine++;
		need_subdiv = true;
	    } else if (nref == 1) {
		int subsd, maxsd;
		double len, lenmax = 0.0;
		for (j = 0; j < nsd; j++) {
		    if (edge[el_edge[i][j]].refine) subsd = j;
		    len = mesh->nlist[edge[el_edge[i][j]].nd[0]].Dist
			(mesh->nlist[edge[el_edge[i][j]].nd[1]]);
		    if (len > lenmax) lenmax = len, maxsd = j;
		}
		if (subsd != maxsd) {
		    for (j = 0; j < nsd; j++)
			edge[el_edge[i][j]].refine = true;
		    nrefine += 2;
		    need_subdiv = true;
		}
	    }
	}
    } while (need_subdiv);
	    
    // create the new nodes
    mesh->nlist.Append (nrefine);
    mesh->plist.Append (nrefine);
    for (i = 0; i < nedge; i++) {
	if (edge[i].refine) {
	    mesh->nlist[nnd].New(2);
	    for (j = 0; j < 2; j++)
		mesh->nlist[nnd][j] = 0.5 *
		    (mesh->nlist[edge[i].nd[0]][j] +
		     mesh->nlist[edge[i].nd[1]][j]);
	    edge[i].mnd = nnd++;
	}
    }

    // subdivide the elements
    for (i = 0, idx = nel; i < nel; i++) {
	Element *pel = mesh->elist[i];
	int nsd = pel->nSide();
	int nsub = 0;
	for (j = 0; j < nsd; j++)
	    if (edge[el_edge[i][j]].refine) nsub++;

	if (nsub) {
	    if (nsub == 3) { // split into 4 sub-triangles
		Triangle3 *tri;
		mesh->elist.Append (tri = new Triangle3);
		tri->Node[0] = edge[el_edge[i][0]].mnd;
		tri->Node[1] = pel->Node[1];
		tri->Node[2] = edge[el_edge[i][1]].mnd;
		mesh->elist.Append (tri = new Triangle3);
		tri->Node[0] = edge[el_edge[i][2]].mnd;
		tri->Node[1] = edge[el_edge[i][1]].mnd;
		tri->Node[2] = pel->Node[2];
		mesh->elist.Append (tri = new Triangle3);
		tri->Node[0] = edge[el_edge[i][0]].mnd;
		tri->Node[1] = edge[el_edge[i][1]].mnd;
		tri->Node[2] = edge[el_edge[i][2]].mnd;
		// now modify original triangle
		pel->Node[1] = edge[el_edge[i][0]].mnd;
		pel->Node[2] = edge[el_edge[i][2]].mnd;
	    } else if (nsub == 1) { // split into 2 sub-triangles
		Triangle3 *tri;
		mesh->elist.Append (tri = new Triangle3);
		for (j = 0; j < nsd; j++)
		    if (edge[el_edge[i][j]].refine) break;
		tri->Node[0] = pel->Node[(j+2)%3];
		tri->Node[1] = edge[el_edge[i][j]].mnd;
		tri->Node[2] = pel->Node[(j+1)%3];
		pel->Node[(j+1)%3] = edge[el_edge[i][j]].mnd;
	    }
	}
    }

    /*
    // jiggle nodes
    Point *pt = new Point[mesh->nlen()];
    for (i = 0; i < mesh->nlen(); i++)
	pt[i] = mesh->NeighbourBarycentre(i);
    for (i = 0; i < mesh->nlen(); i++)
	for (j = 0; j < 2; j++)
	    mesh->nlist[i][j] += (pt[i][j]-mesh->nlist[i][j]) * 0.5;
    delete []pt;
    */
    JiggleMesh (mesh, 0.1, 5);
    mesh->Setup();
}
#endif

// Generate a list of edges for a triangle mesh
void Triangle3_MakeEdgeList (Mesh *mesh, EDGE **triedge, int *ntriedge)
{
    // Create global edge list
    int i, j, k, l, idx, tmp;
    int nel = mesh->elen();
    int nedge_full = nel*3;
    EDGE *edge_full = new EDGE[nedge_full];

    // assemble edges for all elements (including duplicates)
    for (i = idx = 0; i < nel; i++) {
	Element *pel = mesh->elist[i];
	int nsd = pel->nSide();
	for (j = 0; j < nsd; j++) {
	    int nsdnd = pel->nSideNode(j);
	    for (k = 0; k < nsdnd; k++)
		edge_full[idx].nd[k] = pel->Node[pel->SideNode(j,k)];
	    for (k = 0; k < nsdnd-1; k++) {
		for (l = k+1; l < nsdnd; l++) {
		    if (edge_full[idx].nd[k] > edge_full[idx].nd[l]) {
			tmp = edge_full[idx].nd[k];
			edge_full[idx].nd[k] = edge_full[idx].nd[l];
			edge_full[idx].nd[l] = tmp;
		    }
		}
	    }
	    edge_full[idx].el0 = i;
	    edge_full[idx].sd0 = j;
	    idx++;
	}
    }

    // order the edge list according to node indices
    qsort (edge_full, nedge_full, sizeof(EDGE), edgesort_comp);

    // now all non-boundary sides should be listed as consecutive pairs,
    // so condense all matching pairs into single entries
    int nedge = nedge_full;
    for (i = 0; i < nedge_full-1; i++) {
	if (edge_full[i].nd[0] == edge_full[i+1].nd[0] && 
	    edge_full[i].nd[1] == edge_full[i+1].nd[1])
	    nedge--;
    }

    EDGE *edge = new EDGE[nedge];
    for (i = idx = 0; i < nedge_full; i++) {
	edge[idx].nd[0] = edge_full[i].nd[0];
	edge[idx].nd[1] = edge_full[i].nd[1];
	edge[idx].el0 = edge_full[i].el0;
	edge[idx].sd0 = edge_full[i].sd0;
	if (i < nedge_full-1 &&
	    edge_full[i].nd[0] == edge_full[i+1].nd[0] && 
	    edge_full[i].nd[1] == edge_full[i+1].nd[1]) {
	    edge[idx].el1 = edge_full[i+1].el0;
	    edge[idx].sd1 = edge_full[i+1].sd0;
	    i++;
	} else {
	    edge[idx].el1 = -1;
	    edge[idx].sd1 = -1;
	}
	edge[idx].refine = false;
	idx++;
    }
    delete []edge_full;
    *triedge = edge;
    *ntriedge = nedge;
}

Element *Triangle3_FindNeighbour (Mesh *mesh, Element *e, int side)
{
    int i, j, sdnd[2];
    int nel = mesh->elen();

    for (i = 0; i < 2; i++)
	sdnd[i] = e->Node[e->SideNode(side,i)];

    for (i = 0; i < nel; i++) {
	Element *pel = mesh->elist[i];
	if (pel == e) continue;
	int nfound = 0;
	for (j = 0; j < 3; j++) {
	    if (pel->Node[j] == sdnd[0] || pel->Node[j] == sdnd[1])
		nfound++;
	}
	if (nfound == 2)
	    return pel;
    }
    return NULL;
}

// Subdivide element 'e' using the midpoint node 'sidend' on side 'side'
// Note: The rest of the mesh (except for element 'e') is assumed to be valid.
// In particular, the neighbour element sharing 'side' (if any) is assumed to
// be subdivided already, i.e. the midpoint node 'sidend' has already been
// accommodated there.

void Triangle3_SplitElement (Mesh *mesh, Element *e, int sidend)
{
    int i, j, side;
    int gidx[3]; // global node indices
    double dst, maxdst;
    for (i = 0; i < 3; i++)
	gidx[i] = e->Node[i];
    
    // sidend is supposed to point to a side midpoint. So let's first identify
    // the side
    for (i = 0; i < 3; i++) {
	Point midpoint(2);
	for (j = 0; j < 2; j++)
	    midpoint[j] = ((mesh->nlist[gidx[e->SideNode(i,0)]])[j] +
			   (mesh->nlist[gidx[e->SideNode(i,1)]])[j]) * 0.5;
	dst = midpoint.Dist (mesh->nlist[sidend]);
	if (!i || dst < maxdst) {
	    maxdst = dst;
	    side = i;
	}
    }

    // the other triangle sides
    int offside[2];
    for (i = j = 0; i < 3; i++)
	if (i != side) offside[j++] = i;

    // opposing node
    int offnode;
    for (i = 0; i < 3; i++)
	if (i != e->SideNode(side,0) && i != e->SideNode(side,1)) {
	    offnode = i;
	    break;
	}

    // side lengths
    double sidelen[3];
    for (i = 0; i < 3; i++)
	sidelen[i] = mesh->nlist[gidx[e->SideNode(i,0)]].Dist
	    (mesh->nlist[gidx[e->SideNode(i,1)]]);

    if (sidelen[side] >= sidelen[offside[0]]) {
	if (sidelen[side] >= sidelen[offside[1]]) {
	    // subdividing the longest side, so no other sides need splitting
	    Element *ne = new Triangle3;
	    ne->Node[0] = e->Node[offnode];
	    ne->Node[1] = sidend;
	    int flipnode = (offnode+2) % 3;
	    ne->Node[2] = e->Node[flipnode];
	    e->Node[flipnode] = sidend;
	    mesh->elist.Append (ne);
	} else {
	    Element *nbor = Triangle3_FindNeighbour (mesh, e, offside[1]);
	    // also need to subdivide offside[1]
	    int nidx = mesh->nlen();
	    mesh->nlist.Append(1);
	    mesh->nlist[nidx].New(2);
	    for (i = 0; i < 2; i++)
		mesh->nlist[nidx][i] = 0.5*(
		      mesh->nlist[gidx[e->SideNode(offside[1],0)]][i] +
		      mesh->nlist[gidx[e->SideNode(offside[1],1)]][i]);
	    int vtxnode[2];
	    vtxnode[0] = e->SideNode(side,0);
	    vtxnode[1] = e->SideNode(side,1);
	    for (i = 0; i < 2; i++)
		if (e->SideNode(offside[0],i) == vtxnode[1]) {
		    vtxnode[0] = e->SideNode(side,1);
		    vtxnode[1] = e->SideNode(side,0);
		    break;
		}
	    double dst1 = mesh->nlist[sidend].Dist (mesh->nlist[gidx[offnode]]);
	    double dst2 = mesh->nlist[nidx].Dist
		(mesh->nlist[gidx[vtxnode[0]]]);
	    Element *ne1 = new Triangle3;
	    Element *ne2 = new Triangle3;
	    if (dst1 < dst2) {
		ne1->Node[0] = sidend;
		ne1->Node[1] = nidx;
		ne1->Node[2] = gidx[offnode];
		ne2->Node[0] = sidend;
		ne2->Node[1] = nidx;
		ne2->Node[2] = gidx[vtxnode[0]];
	    } else {
		ne1->Node[0] = sidend;
		ne1->Node[1] = nidx;
		ne1->Node[2] = gidx[vtxnode[0]];
		ne2->Node[0] = nidx;
		ne2->Node[1] = gidx[offnode];
		ne2->Node[2] = gidx[vtxnode[0]];
	    }
	    e->Node[0] = sidend;
	    e->Node[1] = gidx[vtxnode[1]];
	    e->Node[2] = nidx;
	    mesh->elist.Append (ne1);
	    mesh->elist.Append (ne2);

	    // now we need to recursively split the corresponding
	    // neighbour element
	    if (nbor)
		Triangle3_SplitElement (mesh, nbor, nidx);
	}
    } else {  // sidelen[side] < sidelen[offside[0]]
	if (sidelen[side] >= sidelen[offside[1]]) {
	    // also need to subdivide offside[0]
	    Element *nbor = Triangle3_FindNeighbour (mesh, e, offside[0]);
	    int nidx = mesh->nlen();
	    mesh->nlist.Append(1);
	    mesh->nlist[nidx].New(2);
	    for (i = 0; i < 2; i++)
		mesh->nlist[nidx][i] = 0.5*(
		      mesh->nlist[gidx[e->SideNode(offside[0],0)]][i] +
		      mesh->nlist[gidx[e->SideNode(offside[0],1)]][i]);
	    int vtxnode[2];
	    vtxnode[0] = e->SideNode(side,0);
	    vtxnode[1] = e->SideNode(side,1);
	    for (i = 0; i < 2; i++)
		if (e->SideNode(offside[0],i) == vtxnode[1]) {
		    vtxnode[0] = e->SideNode(side,1);
		    vtxnode[1] = e->SideNode(side,0);
		    break;
		}
	    double dst1 = mesh->nlist[sidend].Dist (mesh->nlist[gidx[offnode]]);
	    double dst2 = mesh->nlist[nidx].Dist
		(mesh->nlist[gidx[vtxnode[1]]]);
	    Element *ne1 = new Triangle3;
	    Element *ne2 = new Triangle3;
	    if (dst1 < dst2) {
		ne1->Node[0] = sidend;
		ne1->Node[1] = gidx[offnode];
		ne1->Node[2] = nidx;
		ne2->Node[0] = sidend;
		ne2->Node[1] = gidx[vtxnode[1]];
		ne2->Node[2] = gidx[offnode];
	    } else {
		ne1->Node[0] = sidend;
		ne1->Node[1] = gidx[vtxnode[1]];
		ne1->Node[2] = nidx;
		ne2->Node[0] = nidx;
		ne2->Node[1] = gidx[vtxnode[1]];
		ne2->Node[2] = gidx[offnode];
	    }
	    e->Node[0] = sidend;
	    e->Node[1] = nidx;
	    e->Node[2] = gidx[vtxnode[0]];
	    mesh->elist.Append (ne1);
	    mesh->elist.Append (ne2);

	    // now we need to recursively split the corresponding
	    // neighbour element
	    if (nbor)
	    	Triangle3_SplitElement (mesh, nbor, nidx);
	} else {
	    // need to split both sides
	    Element *nbor1 = Triangle3_FindNeighbour (mesh, e, offside[0]);
	    Element *nbor2 = Triangle3_FindNeighbour (mesh, e, offside[1]);

	    int nidx1 = mesh->nlen();
	    mesh->nlist.Append(1);
	    mesh->nlist[nidx1].New(2);
	    for (i = 0; i < 2; i++)
		mesh->nlist[nidx1][i] = 0.5*(
		      mesh->nlist[gidx[e->SideNode(offside[0],0)]][i] +
		      mesh->nlist[gidx[e->SideNode(offside[0],1)]][i]);
	    int nidx2 = mesh->nlen();
	    mesh->nlist.Append(1);
	    mesh->nlist[nidx2].New(2);
	    for (i = 0; i < 2; i++)
		mesh->nlist[nidx2][i] = 0.5*(
		      mesh->nlist[gidx[e->SideNode(offside[1],0)]][i] +
		      mesh->nlist[gidx[e->SideNode(offside[1],1)]][i]);
  
	    int vtxnode[2];
	    vtxnode[0] = e->SideNode(side,0);
	    vtxnode[1] = e->SideNode(side,1);
	    for (i = 0; i < 2; i++)
		if (e->SideNode(offside[0],i) == vtxnode[1]) {
		    vtxnode[0] = e->SideNode(side,1);
		    vtxnode[1] = e->SideNode(side,0);
		    break;
		}
	    Element *ne1 = new Triangle3;
	    ne1->Node[0] = sidend;
	    ne1->Node[1] = nidx1;
	    ne1->Node[2] = gidx[vtxnode[0]];
	    mesh->elist.Append (ne1);
	    Element *ne2 = new Triangle3;
	    ne2->Node[0] = sidend;
	    ne2->Node[1] = gidx[vtxnode[1]];
	    ne2->Node[2] = nidx2;
	    mesh->elist.Append (ne2);
	    Element *ne3 = new Triangle3;
	    ne3->Node[0] = gidx[offnode];
	    ne3->Node[1] = nidx1;
	    ne3->Node[2] = nidx2;
	    mesh->elist.Append (ne3);
	    e->Node[0] = sidend;
	    e->Node[1] = nidx2;
	    e->Node[2] = nidx1;
	    if (nbor1) Triangle3_SplitElement (mesh, nbor1, nidx1);
	    if (nbor2) Triangle3_SplitElement (mesh, nbor2, nidx2);
	}
    }

    // check mesh integrity
    mesh->Setup(false);
    for (i = 0; i < mesh->elen(); i++) {
	if (mesh->elist[i]->Size() < 0.0) {
	    cerr << "flipping element " << i << endl;
	    int t = mesh->elist[i]->Node[1];
	    mesh->elist[i]->Node[1] = mesh->elist[i]->Node[2];
	    mesh->elist[i]->Node[2] = t;
	}
    }
}

// Subdivide element 'el' in 'mesh'. If necessary, recursively refine the
// neighbours as well.
FELIB void Mesh_SplitTriangle3 (Mesh *mesh, int el)
{
    Element *pel = mesh->elist[el];
    if (pel->Type() != ELID_TRI3)
	return; // error

    int i, s, smax;
    double lmax = 0.0;
    int gidx[3];  // global node indices of the triangle
    for (i = 0; i < 3; i++)
	gidx[i] = pel->Node[i];

    // find the longest side
    for (s = 0; s < 3; s++) {
	int n0 = gidx[pel->SideNode(s,0)];
	int n1 = gidx[pel->SideNode(s,1)];
	double len = mesh->nlist[n0].Dist (mesh->nlist[n1]);
	if (len >= lmax) {
	    lmax = len;
	    smax = s;
	}
    }
    // Find neighbour element
    Element *nbor = Triangle3_FindNeighbour (mesh, pel, smax);

    // create new midpoint node
    int midx = mesh->nlen(); // midpoint node index
    mesh->nlist.Append(1);
    mesh->nlist[midx].New(2);
    for (i = 0; i < 2; i++)
	mesh->nlist[midx][i] =
    	    ((mesh->nlist[gidx[pel->SideNode(smax,0)]])[i] +
    	     (mesh->nlist[gidx[pel->SideNode(smax,1)]])[i]) * 0.5;

    // create a new element, and adjust the indices
    Element *ne = new Triangle3;
    mesh->elist.Append(ne);
    ne->Node[0] = midx;
    ne->Node[1] = gidx[(smax+1)%3];
    ne->Node[2] = gidx[(smax+2)%3];

    pel->Node[0] = midx;
    pel->Node[1] = gidx[(smax+2)%3];
    pel->Node[2] = gidx[smax];

    // Now we need to subdivide the neighbour
    if (nbor)
	Triangle3_SplitElement (mesh, nbor, midx); 

    // Refresh mesh to register changes
    mesh->Setup();
}

FELIB void RefineTriangle3Mesh (Mesh *mesh, bool *elrefine)
{
    int i, j, idx;
    int nel = mesh->elen();
    int nnd = mesh->nlen();

    // create the edge list for the mesh
    EDGE *edge;
    int nedge;
    Triangle3_MakeEdgeList (mesh, &edge, &nedge);

    // create element->edge reference list
    typedef int TRI_EDGE[3];
    TRI_EDGE *el_edge = new TRI_EDGE[nel];

    for (i = 0; i < nedge; i++) {
	el_edge[edge[i].el0][edge[i].sd0] = i;
	if (edge[i].el1 >= 0)
	    el_edge[edge[i].el1][edge[i].sd1] = i;
    }

    // check which edges should be refined
    for (i = 0; i < nel; i++) {
	Element *pel = mesh->elist[i];
	int nsd = pel->nSide();
	if (elrefine[i]) {
	    for (j = 0; j < nsd; j++) {
		edge[el_edge[i][j]].refine = true;
	    }
	}
    }
    int nrefine = 0;
    for (i = 0; i < nedge; i++) {
	if (edge[i].refine) nrefine++;
    }

    // go through the list and check if any triangles try to refine
    // two of their sides, or if any single-side refinement is not refining
    // the longest side
    bool need_subdiv = false;
    do {
	need_subdiv = false;
	for (i = 0; i < nel; i++) {
	    Element *pel = mesh->elist[i];
	    int nsd = pel->nSide();
	    int nref = 0;
	    for (j = 0; j < nsd; j++)
		if (edge[el_edge[i][j]].refine) nref++;
	    if (nref == 2) {
		for (j = 0; j < nsd; j++)
		    edge[el_edge[i][j]].refine = true;
		nrefine++;
		need_subdiv = true;
	    } else if (nref == 1) {
		int subsd, maxsd;
		double len, lenmax = 0.0;
		for (j = 0; j < nsd; j++) {
		    if (edge[el_edge[i][j]].refine) subsd = j;
		    len = mesh->nlist[edge[el_edge[i][j]].nd[0]].Dist
			(mesh->nlist[edge[el_edge[i][j]].nd[1]]);
		    if (len > lenmax) lenmax = len, maxsd = j;
		}
		if (subsd != maxsd) {
		    for (j = 0; j < nsd; j++)
			edge[el_edge[i][j]].refine = true;
		    nrefine += 2;
		    need_subdiv = true;
		}
	    }
	}
    } while (need_subdiv);
	    
    // create the new nodes
    mesh->nlist.Append (nrefine);
    for (i = 0; i < nedge; i++) {
	if (edge[i].refine) {
	    mesh->nlist[nnd].New(2);
	    for (j = 0; j < 2; j++)
		mesh->nlist[nnd][j] = 0.5 *
		    (mesh->nlist[edge[i].nd[0]][j] +
		     mesh->nlist[edge[i].nd[1]][j]);
	    edge[i].mnd = nnd++;
	}
    }

    // subdivide the elements
    for (i = 0, idx = nel; i < nel; i++) {
	Element *pel = mesh->elist[i];
	int nsd = pel->nSide();
	int nsub = 0;
	for (j = 0; j < nsd; j++)
	    if (edge[el_edge[i][j]].refine) nsub++;

	if (nsub) {
	    if (nsub == 3) { // split into 4 sub-triangles
		Triangle3 *tri;
		mesh->elist.Append (tri = new Triangle3);
		tri->Node[0] = edge[el_edge[i][0]].mnd;
		tri->Node[1] = pel->Node[1];
		tri->Node[2] = edge[el_edge[i][1]].mnd;
		mesh->elist.Append (tri = new Triangle3);
		tri->Node[0] = edge[el_edge[i][2]].mnd;
		tri->Node[1] = edge[el_edge[i][1]].mnd;
		tri->Node[2] = pel->Node[2];
		mesh->elist.Append (tri = new Triangle3);
		tri->Node[0] = edge[el_edge[i][0]].mnd;
		tri->Node[1] = edge[el_edge[i][1]].mnd;
		tri->Node[2] = edge[el_edge[i][2]].mnd;
		// now modify original triangle
		pel->Node[1] = edge[el_edge[i][0]].mnd;
		pel->Node[2] = edge[el_edge[i][2]].mnd;
	    } else if (nsub == 1) { // split into 2 sub-triangles
		Triangle3 *tri;
		mesh->elist.Append (tri = new Triangle3);
		for (j = 0; j < nsd; j++)
		    if (edge[el_edge[i][j]].refine) break;
		tri->Node[0] = pel->Node[(j+2)%3];
		tri->Node[1] = edge[el_edge[i][j]].mnd;
		tri->Node[2] = pel->Node[(j+1)%3];
		pel->Node[(j+1)%3] = edge[el_edge[i][j]].mnd;
	    }
	}
    }

    /*
    // jiggle nodes
    Point *pt = new Point[mesh->nlen()];
    for (i = 0; i < mesh->nlen(); i++)
	pt[i] = mesh->NeighbourBarycentre(i);
    for (i = 0; i < mesh->nlen(); i++)
	for (j = 0; j < 2; j++)
	    mesh->nlist[i][j] += (pt[i][j]-mesh->nlist[i][j]) * 0.5;
    delete []pt;
    */
    JiggleMesh (mesh, 0.1, 5);
    mesh->Setup();
}


FELIB void RefineTetrahedron4Mesh (Mesh *mesh, bool *elrefine)
{
}

void JiggleMesh (Mesh *mesh, double scale, int iterations)
{
    int it, i, j;

    for (it = 0; it < iterations; it++) {
	Point *pt = new Point[mesh->nlen()];
	for (i = 0; i < mesh->nlen(); i++)
	    pt[i] = mesh->NeighbourBarycentre(i);
	for (i = 0; i < mesh->nlen(); i++)
	    for (j = 0; j < 2; j++)
		mesh->nlist[i][j] += (pt[i][j]-mesh->nlist[i][j]) * scale;
	mesh->Setup();
	
	bool ok = true;
	for (i = 0; i < mesh->elen(); i++) {
	    if (mesh->ElSize(i) <= 0.0) {
		ok = false;
		break;
	    }
	}
	if (!ok) {
	    std::cerr << "Jiggle: Negative element size detected. Reverting"
		      << std::endl;
	    for (i = 0; i < mesh->nlen(); i++)
		for (j = 0; j < 2; j++)
		    mesh->nlist[i][j] -= (pt[i][j]-mesh->nlist[i][j]) * scale;
	    mesh->Setup();
	    
	}
	delete []pt;
	if (!ok) break;
    }
}
