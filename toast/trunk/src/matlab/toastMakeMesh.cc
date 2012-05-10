// =========================================================================
// toastMakeMesh
// Creates a TOAST mesh from vertex coordinates and element connectivity
// data.
//
// RH parameters:
//     1: vertex array (double, n x 2 or n x 3, where n is number of nodes)
//     2: element index list (integer, m x v, where m is number of
//        elements, and v is max. number of nodes in an element). The list is
//        1-based.
//     3: element type list (integer, m x 1). Supported types:
//           3: 4-noded tetrahedra
//     4: array of nodal parameter values (double, n x (1-3)), optional
// LH parameters:
//     1: mesh handle (pointer)
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, k;
    int nvtx = (int)mxGetM(prhs[0]);
    int nel  = (int)mxGetM(prhs[1]);
    int dim  = (int)mxGetN(prhs[0]);
    int nnd0 = (int)mxGetN(prhs[1]);
    double *vtx = mxGetPr (prhs[0]);
    double *idx = mxGetPr (prhs[1]);
    double *etp = mxGetPr (prhs[2]);

    Mesh *mesh = new QMMesh;;

    // create node list
    mesh->nlist.New (nvtx);
    for (i = 0; i < nvtx; i++) {
	mesh->nlist[i].New(dim);
	mesh->nlist[i].SetBndTp (BND_NONE); // don't know
    }
    for (j = k = 0; j < dim; j++) {
	for (i = 0; i < nvtx; i++) {
	    mesh->nlist[i][j] = vtx[k++];
	}
    }

    // create element list
    Element *el, **list = new Element*[nel];
    for (i = 0; i < nel; i++) {
	int eltp = (int)(etp[i]+0.5);
	switch (eltp) {
	case ELID_TRI3OLD:
	    list[i] = new Triangle3old;
	    break;
	case ELID_TET4:
	    list[i] = new Tetrahedron4;
	    break;
	case ELID_WDG6:
	    list[i] = new Wedge6;
	    break;
	case ELID_VOX8:
	    list[i] = new Voxel8;
	    break;
	case ELID_TRI6:
	    list[i] = new Triangle6;
	    break;
	case ELID_TET10:
	    list[i] = new Tetrahedron10;
	    break;
	case ELID_TRI6_IP:
	    list[i] = new Triangle6_ip;
	    break;
	case ELID_TRI10:
	    list[i] = new Triangle10;
	    break;
	case ELID_TRI10_IP:
	    list[i] = new Triangle10_ip;
	    break;
	case ELID_TET10_IP:
	    list[i] = new Tetrahedron10_ip;
	    break;
	case ELID_PIX4:
	    list[i] = new Pixel4;
	    break;
	case ELID_TRI3:
	    list[i] = new Triangle3;
	    break;
	case ELID_TRI3D3:
	    list[i] = new Triangle3D3;
	    break;
	case ELID_TRI3D6:
	    list[i] = new Triangle3D6;
	    break;
	default:
	  mexPrintf ("Element type %d not supported!\n", eltp);
	    list[i] = 0;
	    break;
	}
    }
    mesh->elist.SetList (nel, list);
    delete []list;

    for (j = k = 0; j < nnd0; j++) {
	for (i = 0; i < nel; i++) {
	    if (el = mesh->elist[i]) {
		if (j < el->nNode())
		    el->Node[j] = (int)(idx[k]-0.5);
	    }
	    k++;
	   }
    }

    // create dummy parameter list
    mesh->plist.New (nvtx);
    mesh->plist.SetMua (0.01);
    mesh->plist.SetMus (1);
    mesh->plist.SetN (1);


    // copy user-provided list if available
    if (nrhs >= 4) {
	RVector prm(nvtx);
	double *pprm = mxGetPr (prhs[3]);
	int nprm = mxGetN (prhs[3]);
	if (nprm >= 1) {
	    for (i = 0; i < nvtx; i++) prm[i] = *pprm++;
	    mesh->plist.SetMua (prm);
	}
	if (nprm >= 2) {
	    for (i = 0; i < nvtx; i++) prm[i] = *pprm++;
	    mesh->plist.SetMus (prm);
	}
	if (nprm >= 3) {
	    for (i = 0; i < nvtx; i++) prm[i] = *pprm++;
	    mesh->plist.SetN (prm);
	}
    }
    
    // set up mesh
    mesh->Setup();

    plhs[0] = mxCreateDoubleScalar (Ptr2Handle (mesh));

    mexPrintf ("Mesh: %d nodes, %d elements, dimension %d\n", nvtx, nel, dim);
    mexPrintf ("      %d boundary nodes\n", mesh->nbnd());
}
