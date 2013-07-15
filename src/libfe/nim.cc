// -*-C++-*-
// Utilities for nodal images

#include "nim.h"

// ===========================================================================
// Local prototypes

int ProcessNeighbours (const Mesh *mesh, int **nbrs, int *nnbrs,
    RVector &phase, TVector<bool> &processed, TVector<int> &front, int &nfront);


// ===========================================================================
// Unwrap nodal phase vector 'phase', starting from point 'seed'
// Return the number of unwrapped nodes

int NimPhaseUnwrap (const Mesh *mesh, RVector &phase, Point seed)
{
    // Get the node neighbour list
    int *nnbrs, **nbrs;
    int i, nlen = mesh->nlen();
    mesh->NodeNeighbourList (&nnbrs, &nbrs);

    // Find node closest to seed
    double dstmin = seed.Dist (mesh->nlist[0]);
    int seednd = 0;
    for (i = 1; i < nlen; i++) {
	double dst = seed.Dist (mesh->nlist[i]);
	if (dst < dstmin) {
	    dstmin = dst;
	    seednd = i;
	}
    }

    TVector<bool> processed(nlen);
    TVector<int> frontline(nlen);
    frontline[0] = seednd;
    int nfront = 1;
    int nunwrap = 0;

    for (i = 0; i < nlen; i++)
	processed[i] = false;
    processed[seednd] = true;

    while (nfront) {
	nunwrap += ProcessNeighbours (mesh, nbrs, nnbrs, phase, processed,
				      frontline, nfront);
    }

    // cleanup
    for (i = 0; i < nlen; i++)
	delete []nbrs[i];
    delete []nbrs;
    delete []nnbrs;

    return nunwrap;
}

int ProcessNeighbours (const Mesh *mesh, int **nbrs, int *nnbrs,
    RVector &phase, TVector<bool> &processed, TVector<int> &front, int &nfront)
{
    int nunwrap = 0;
    int i, j, noldfront = nfront;
    TVector<int> oldfront(noldfront);
    for (i = 0; i < noldfront; i++)
	oldfront[i] = front[i];

    nfront = 0;
    for (i = 0; i < noldfront; i++) {
	int s = oldfront[i];
	for (j = 0; j < nnbrs[s]; j++) {
	    int nb = nbrs[s][j];
	    if (!processed[nb]) {
	        if (phase[nb] >= phase[s]+Pi) {
		    while (phase[nb] >= phase[s]+Pi)
		        phase[nb] -= Pi2;
		    nunwrap++;
		} else if (phase[nb] < phase[s]-Pi) {
		    while (phase[nb] < phase[s]-Pi)
		        phase[nb] += Pi2;
		    nunwrap++;
		}
		processed[nb] = true;
		front[nfront++] = nb;
	    }
	}
    }
    return nunwrap;
}
