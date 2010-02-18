#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"
#include "toastmex.h"

using namespace std;

void Reorder (Mesh &mesh, int *perm)
{
    int i, j, ii, ij, nds = mesh.nlen(), els = mesh.elen();
    int *iperm = new int[nds];

    for (i = 0; i < nds; i++) iperm[perm[i]] = i;

    // reorder node indices in element list

    for (i = 0; i < els; i++)
        for (j = 0; j < mesh.elist[i]->nNode(); j++)
	    mesh.elist[i]->Node[j] = iperm[mesh.elist[i]->Node[j]];

    // reorder node and parameter lists

    for (i = 0; i < nds; i++) {
        j = perm[i];
	ii = iperm[i], ij = iperm[j];
        if (i == j) continue;
	mesh.nlist.Swap (i, j);
	mesh.plist.Swap (i, j);
	perm[ii] = j, perm[ij] = i;
	iperm[i] = ij, iperm[j] = ii;
    }

    delete []iperm;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    RVector fperm;
    CopyVector (fperm, prhs[1]);
    int i, len = fperm.Dim();
    int *perm = new int[len];

    for (i = 0; i < len; i++)
	perm[i] = (int)(fperm[i]+0.5);

    Reorder (*mesh, perm);

    delete []perm;
}
