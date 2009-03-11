// =======================================================================
// mex version
// takes a mesh handle as input and 
// Converts TRI3 -> TRI6
// and      TET4 -> TET10
// Other element types are not supported
// the same handle points to the new mesh
// =======================================================================


#include "mex.h"
#include <mathlib.h>
#include <felib.h>
#include "util.h"


using namespace std;

void Lin2Quad (Mesh &mesh);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
	cerr<<" Lin2Quad-----------------------------------------------------------";
   int elen = mesh->elen();
   int nlen = mesh->nlen();

    // make sure you have the correct element types
  
 //       switch (mesh->elist[0]->Type()) {
//	case ELID_TRI3OLD:
//	    Lin2Quad (*mesh);
//	case ELID_TRI3:
//	    Lin2Quad (*mesh);
//	case ELID_TET4:
//	    Lin2Quad (*mesh);
//	default:
//	    cerr << "lin2quadmesh: Unsupported element type detected.\n";
//	    cerr << "lin2quadmesh: Giving up.\n";
//	    exit (1);
//	}
    


    Lin2Quad (*mesh);
	cerr<<" done-----------------------------------------------------------";
    //mesh.Setup();
    // set up mesh
    //mesh->MarkBoundary();
   // mesh->Setup();

   // plhs[0] = mxCreateScalarDouble (Ptr2Handle (mesh));

   // mexPrintf ("Mesh: %d nodes, %d elements, dimension %d\n", nvtx, nel, dim);
}

//int main (void)
//{
//    CHECK_EXPIRED();

//    Mesh mesh;

 //   cin >> mesh;
//    mesh.Setup();
//    Lin2Quad (mesh);
//    mesh.Setup();
//    cout << mesh;
//    return 0;
//}

void Lin2Quad (Mesh &mesh)
{
    const int chunksize = 64;
    int el, elen, nlen, i, j, k;
    int nn, nbj, *nb, *nd;
    int nnlen = 0, maxnew = 0;
    int **elref, *nelref, *bufsize;
    int dim = mesh.Dimension();

    elen = mesh.elen();
    nlen = mesh.nlen();

    // find max number of new nodes
    for (el = 0; el < elen; el++) {
        switch (mesh.elist[el]->Type()) {
	case ELID_TRI3OLD:
	    maxnew += 3;
	    break;
	case ELID_TRI3:
	    maxnew += 3;
	    break;
	case ELID_TET4:
	    maxnew += 6;
	    break;
	default:
	    cerr << "lin2quadmesh: Unsupported element type detected.\n";
	    cerr << "lin2quadmesh: Giving up.\n";
	    exit (1);
	}
    }

    // allocate lists for new nodes
    NodeList nnlist (maxnew);
    ParameterList nplist (maxnew);

    // find neighbour nodes for each node in the mesh
    cout << "Setting up node neighbour list\n";
    int *nd_nnbhr, **nd_nbhr;
    mesh.NodeNeighbourList (&nd_nnbhr, &nd_nbhr);

    // build node->element reference list
    cout << "Building node->element reference list\n";
    elref   = new int*[nlen];
    nelref  = new int[nlen];
    bufsize = new int[nlen];
    for (i = 0; i < nlen; i++) {
        elref[i] = new int[bufsize[i]=chunksize];
	nelref[i] = 0;
    }
    for (el = 0; el < elen; el++) {
        nn = mesh.elist[el]->nNode();
        nd = mesh.elist[el]->Node;
	for (i = 0; i < nn; i++) {
	    int n = nd[i];
	    if (nelref[n] == bufsize[n]) { // reallocate
	        int *tmp = new int[bufsize[n]+chunksize];
		memcpy (tmp, elref[n], bufsize[n]*sizeof(int));
		delete []elref[n];
		elref[n] = tmp;
		bufsize[n] += chunksize;
		cout<< "Increasing elref buffer to " << bufsize[n] << endl;
	    }
	    elref[n][nelref[n]++] = el;
	}
    }
    delete []bufsize;
    
    // copy element list
    cout << "Copying element list\n";
    for (el = 0; el < elen; el++) {
        Element *nel;
	switch (mesh.elist[el]->Type()) {
	case ELID_TRI3OLD:
	    nel = new Triangle6;
	    // node order in TRI3OLD is inconsistent
	    // the following hack accounts for this
	    i = mesh.elist[el]->Node[1];
	    mesh.elist[el]->Node[1] = mesh.elist[el]->Node[2];
	    mesh.elist[el]->Node[2] = i;
	    break;
	case ELID_TRI3:
	    nel = new Triangle6;
	    break;
	case ELID_TET4:
	    nel = new Tetrahedron10;
	    break;
	}
	nn = mesh.elist[el]->nNode();
	for (i = 0; i < nn; i++) nel->Node[i] = mesh.elist[el]->Node[i];
	delete mesh.elist[el];
	mesh.elist[el] = nel;
    }

    for (i = 0; i < nlen; i++) {
        if (!(i%1000)) cout << "Processing node " << i
			    << " of " << nlen << endl;
        nn = nd_nnbhr[i];
	nb = nd_nbhr[i];
	Node &ni = mesh.nlist[i];
	Parameter &pi = mesh.plist[i];
	for (j = 0; j < nn; j++) {
	    if ((nbj = nb[j]) <= i) continue; // to avoid counting twice
	    Node &nj = mesh.nlist[nbj];
	    Parameter &pj = mesh.plist[nbj];
	    // create a new node between i and its jth neighbour
	    nnlist[nnlen].New(dim);
	    for (k = 0; k < dim; k++)
	        nnlist[nnlen][k] = 0.5*(ni[k]+nj[k]);
	    if (ni.BndTp() == nj.BndTp()) nnlist[nnlen].SetBndTp (ni.BndTp());
	    else                          nnlist[nnlen].SetBndTp (BND_NONE);

	    nplist[nnlen].SetMua ((pi.Mua()+pj.Mua())*0.5);
	    nplist[nnlen].SetKappa ((pi.Kappa()+pj.Kappa())*0.5);
	    nplist[nnlen].SetN ((pi.N()+pj.N())*0.5);
	    nplist[nnlen].SetA ((pi.A()+pj.A())*0.5);

	    // now find all elements sharing the new node
	    for (k = 0; k < nelref[i]; k++) {
	        el = elref[i][k];
	        int *node = mesh.elist[el]->Node;
		switch (mesh.elist[el]->Type()) {
		case ELID_TRI6:
		    if (i == node[0]) {
		        if      (nbj == node[1]) node[3] = nlen+nnlen;
			else if (nbj == node[2]) node[5] = nlen+nnlen;
		    } else if (i == node[1]) {
		        if      (nbj == node[0]) node[3] = nlen+nnlen;
			else if (nbj == node[2]) node[4] = nlen+nnlen;
		    } else if (i == node[2]) {
		        if      (nbj == node[0]) node[5] = nlen+nnlen;
			else if (nbj == node[1]) node[4] = nlen+nnlen;
		    } else
		        cerr << "Panic!\n";
		    break;
		case ELID_TET10:
		    if (i == node[0]) {
		        if      (nbj == node[1]) node[4] = nlen+nnlen;
			else if (nbj == node[2]) node[5] = nlen+nnlen;
			else if (nbj == node[3]) node[6] = nlen+nnlen;
		    } else if (i == node[1]) {
		        if      (nbj == node[0]) node[4] = nlen+nnlen;
			else if (nbj == node[2]) node[7] = nlen+nnlen;
			else if (nbj == node[3]) node[8] = nlen+nnlen;
		    } else if (i == node[2]) {
		        if      (nbj == node[0]) node[5] = nlen+nnlen;
			else if (nbj == node[1]) node[7] = nlen+nnlen;
			else if (nbj == node[3]) node[9] = nlen+nnlen;
		    } else if (i == node[3]) {
		        if      (nbj == node[0]) node[6] = nlen+nnlen;
			else if (nbj == node[1]) node[8] = nlen+nnlen;
			else if (nbj == node[2]) node[9] = nlen+nnlen;
		    } else
		        cerr << "Panic!\n";
		    break;
		}
	    }
	    nnlen++;
	}
    }

    cout << "Finalising mesh\n";
    if (nnlen) {
        mesh.nlist.Append (nnlen);
	mesh.plist.Append (nnlen);
	for (i = 0; i < nnlen; i++) {
	    mesh.nlist[nlen+i].Copy (nnlist[i]);
	    mesh.plist[nlen+i] = nplist[i];
	}
    }

    // cleanup
    for (i = 0; i < nlen; i++) delete elref[i];
    delete []elref;
    delete []nelref;
}
