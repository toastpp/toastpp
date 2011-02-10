#include "felib.h"

typedef struct EDGE {
    int nd[2];    // node indices
    int el0, el1; // attached element indices (-1 for none)
    int sd0, sd1; // corresponding local element side indices
    bool refine;
    int mnd;      // node index for midpoint node (only if refine==true)
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

    // jiggle nodes
    Point *pt = new Point[mesh->nlen()];
    for (i = 0; i < mesh->nlen(); i++)
	pt[i] = mesh->NeighbourBarycentre(i);
    for (i = 0; i < mesh->nlen(); i++)
	for (j = 0; j < 2; j++)
	    mesh->nlist[i][j] += (pt[i][j]-mesh->nlist[i][j]) * 0.5;
    delete []pt;

    mesh->Setup();
}
