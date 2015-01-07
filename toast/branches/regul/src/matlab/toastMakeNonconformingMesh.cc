// =========================================================================
// toastMakeNonconformingMesh
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
#include "nonconformingMesh.h"
#include "util.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, k;
    size_t nvtx = mxGetM(prhs[0]);
    size_t nel  = mxGetM(prhs[1]);
    size_t dim  = mxGetN(prhs[0]);
    size_t nnd0 = mxGetN(prhs[1]);
    double *vtx = mxGetPr (prhs[0]);
    double *idx = mxGetPr (prhs[1]);
    double *etp = mxGetPr (prhs[2]);
   

    NonconformingMesh *mesh = new NonconformingMesh();

    // create node list
    mesh->nlist.New ((int)nvtx);
    for (i = 0; i < (int)nvtx; i++) {
	mesh->nlist[i].New((int)dim);
	mesh->nlist[i].SetBndTp (BND_NONE); // don't know
    }
    for (j = k = 0; j < (int)dim; j++) {
	for (i = 0; i < (int)nvtx; i++) {
	    mesh->nlist[i][j] = vtx[k++];
	}
    }

    // create element list
    Element *el, **list = new Element*[nel];
    for (i = 0; i < (int)nel; i++) {
	int eltp = (int)(etp[i]+0.5);
	switch (eltp) {
	case ELID_TET4:
	    list[i] = new Tetrahedron4;
	    break;
	case ELID_TRI3:
	    list[i] = new Triangle3;
	    break;
	default:
	    mexPrintf("Nonconforming mesh class is only supported for Tri3 and Tet4 elements.\n");
	    return;
	}
    }
    mesh->elist.Clear();
    mesh->elist.AppendList ((int)nel, list);
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

    std::vector<int> vec1(2), vec2(3);
    
    // set up mesh
    if(nrhs == 3)
    {
	//mesh->MarkBoundary();
    	mesh->SetupEdgeTables();
     }
     else
     {
	double *ielist = mxGetPr(prhs[3]);
    	double *inlist = mxGetPr(prhs[4]);
    	double *istate = mxGetPr(prhs[5]);
    	double *belist = mxGetPr(prhs[6]);
    	double *bnlist = mxGetPr(prhs[7]);
    	double *bstate = mxGetPr(prhs[8]);
 	size_t nie = mxGetM(prhs[3]);
    	size_t nbe = mxGetM(prhs[6]);


	mexPrintf("SetupEdgeTables function skipped ...\n");
	j=0; k=0; 
        int l=0;
	for(i=0; i <nie; i++)
	{
	  vec1.clear();vec2.clear();
	  vec1.push_back((int)(ielist[j]-0.5)); vec1.push_back((int)(ielist[j+nie]-0.5)); 
	  j++; 
	  vec2.push_back((int)(inlist[k]-0.5)); vec2.push_back((int)(inlist[k+nie]-0.5)); vec2.push_back((int)(inlist[k+2*nie]-0.5)); 
	  k++;

	  mesh->iedge_elist.insert(mesh->iedge_elist.end(), vec1);
	  mesh->iedge_nlist.insert(mesh->iedge_nlist.end(), vec2);
	  mesh->iedge_state.push_back((short)istate[l++]);
	}

	j=0; k=0; l=0;
	for(i=0; i <nbe; i++)
	{
	  vec2.clear();
	  vec2.push_back((int)(bnlist[k]-0.5)); vec2.push_back((int)(bnlist[k+nbe]-0.5)); vec2.push_back((int)(bnlist[k+2*nbe]-0.5)); 
	  k++;
	  mesh->bedge_elist.push_back((int)(belist[j++]-0.5));
	  mesh->bedge_nlist.insert(mesh->bedge_nlist.end(), vec2);
	  mesh->bedge_state.push_back((short)bstate[l++]);	
	}
	  //mesh->MarkBoundary();
	  mesh->Setup();
		
	}
  
    plhs[0] = mxCreateDoubleScalar (Ptr2Handle (mesh));

    mexPrintf ("Mesh nodes:%d elements:%d dimension:%d number of interior edges:%d number of boundary edges:%d\n", nvtx, nel, dim, mesh->iedge_elist.size(), mesh->bedge_elist.size());
    /*for(int i =0 ; i<10; i++)
    { 
	std::set<elems>::iterator it1;
	std::set<nodes>::iterator it2;
	it1 = mesh->iedgelist.elems.begin();
	it2 = mesh->iedgelist.nodes.begin();
	std::vector<int> vec1 = *it1;
	std::vector<int> vec2 = *it2;
        int e = vec1.back();
	vec1.pop_back();
	int el = vec1.back();
	vec1.pop_back();

	int n1 = vec2.back(); vec2.pop_back();
	int n2 = vec2.back(); vec2.pop_back();

	mexPrintf("%d %d %d %d\n",e, el, n1, n2);
	it1++;
	it2++;
    }*/
}
