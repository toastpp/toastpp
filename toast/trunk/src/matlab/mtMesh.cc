// ========================================================================
// Implementation of class MatlabToast
// mesh-related methods
// NOTE: The core mesh optimisation methods should be moved from here into
// the Mesh class.
// ========================================================================

#include "matlabtoast.h"

using namespace std;
using namespace toast;

// =========================================================================

void MatlabToast::ReadMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char meshname[256];

    if (mxIsChar (prhs[0]))
	mxGetString (prhs[0], meshname, 256);
    else
	mexErrMsgTxt ("ReadMesh: Argument 1: file name expected.");
	
    if (fileExists (meshname) != 1)
	mexErrMsgTxt ("ReadMesh: Mesh file not found.");

    QMMesh *mesh = new QMMesh;
    ifstream ifs(meshname);
    ifs >> *mesh;
    if (!ifs.good())
	mexErrMsgTxt ("ReadMesh: Mesh file invalid format.");
    mesh->Setup();

    Mesh **tmp = new Mesh*[nmesh+1];
    if (nmesh) {
	memcpy (tmp, meshlist, nmesh*sizeof(Mesh*));
	delete []meshlist;
    }
    meshlist = tmp;
    meshlist[nmesh++] = mesh;

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT32_CLASS, mxREAL);
    unsigned int *ptr = (unsigned int*)mxGetData (plhs[0]);
    *ptr = nmesh;

    if (verbosity >= 1) {
	mexPrintf("Mesh: %d nodes, %d elements\n", mesh->nlen(), mesh->elen());	
    }
}

// =========================================================================

void MatlabToast::MakeMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j, k;
    size_t nvtx = mxGetM(prhs[0]);
    size_t nel  = mxGetM(prhs[1]);
    size_t dim  = mxGetN(prhs[0]);
    size_t nnd0 = mxGetN(prhs[1]);
    double *vtx = mxGetPr (prhs[0]);
    double *idx = mxGetPr (prhs[1]);
    double *etp = mxGetPr (prhs[2]);

    Mesh *mesh = new QMMesh;;

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
	    mexErrMsgTxt ("Element type not supported!\n");
	    list[i] = 0;
	    break;
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

    
    // create dummy parameter list
    mesh->plist.New ((int)nvtx);
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
    mesh->MarkBoundary();
    mesh->Setup();

    // insert mesh into mesh list
    Mesh **tmp = new Mesh*[nmesh+1];
    if (nmesh) {
	memcpy (tmp, meshlist, nmesh*sizeof(Mesh*));
	delete []meshlist;
    }
    meshlist = tmp;
    meshlist[nmesh++] = mesh;

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT32_CLASS, mxREAL);
    unsigned int *ptr = (unsigned int*)mxGetData (plhs[0]);
    *ptr = nmesh;

    if (verbosity >= 1)
	mexPrintf ("Mesh: %d nodes, %d elements, dimension %d\n", nvtx, nel, dim);
}

// =========================================================================

void MatlabToast::WriteMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    char meshname[256];
    
    if (mxIsChar (prhs[1]))
	mxGetString (prhs[1], meshname, 256);
    else
	mexErrMsgTxt ("WriteMesh: Argument 2: file name expected.");
	
    ofstream ofs (meshname);
    ofs << *mesh;

    if (verbosity >= 1)
	mexPrintf("Mesh: written to %s\n", meshname);
}

// =========================================================================

void MatlabToast::MeshOpt (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, res, len, *perm;
    char optmode[256] = "\0";

	Mesh *mesh = GETMESH_SAFE(0);
    ASSERTARG(mxGetString (prhs[1], optmode, 256), 2, "expected string");

    len = mesh->nlen();
    perm = new int[len];
    for (i = 0; i < len; i++) perm[i] = i;

    if (!strcasecmp(optmode, "mmd")) {
	res = Optimise_MMD (*mesh, perm, 0, len);
    } else if (!strcasecmp(optmode, "bandw")) {
	res = Optimise_MinBandwidth (*mesh, perm, 0, len);
    } else if (!strcasecmp (optmode, "tinney2")) {
	res = Optimise_Tinney2 (*mesh, perm, 0, len);
    } else {
	res = 1;
    }

    ASSERTARG(res == 0, 2, "optimisation reports failure");

    plhs[0] = mxCreateDoubleMatrix (len,1,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < len; i++)
	*pr++ = perm[i]+1; // switch to 1-based

    delete []perm;
}

// =========================================================================

void MatlabToast::MeshNodeCount (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
	Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateScalarDouble (mesh->nlen());
}

// =========================================================================

void MatlabToast::MeshElementCount (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
	Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateScalarDouble (mesh->elen());
}

// =========================================================================

void MatlabToast::ClearMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
	Mesh *mesh = GETMESH_SAFE(0);
	unsigned int meshid = *(unsigned int*)mxGetData(prhs[0])-1;
	if (meshlist[meshid]) {
		delete meshlist[meshid];
		meshlist[meshid] = 0;
		if (verbosity >= 1)
			mexPrintf("%s: Mesh %d cleared.\n", __FUNCTION__, meshid+1);
	}
}

// =========================================================================

void MatlabToast::MeshData (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    int i, j;
    double *pr;

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // vertex coordinate list
    mxArray *vtx = mxCreateDoubleMatrix (nlen, dim, mxREAL);
    pr = mxGetPr (vtx);

    for (i = 0; i < dim; i++)
	for (j = 0; j < nlen; j++)
	    *pr++ = mesh->nlist[j][i];

    // max number of nodes per element
    int nnd = mesh->elist[0]->nNode();
    for (i = 0; i < elen; i++)
	nnd = ::max (nnd, mesh->elist[i]->nNode());

    // element index list
    // (1-based; value 0 indicates unused matrix entry)
    mxArray *idx = mxCreateDoubleMatrix (elen, nnd, mxREAL);
    pr = mxGetPr (idx);
    
    for (i = 0; i < nnd; i++)
	for (j = 0; j < elen; j++)
	    if (i < mesh->elist[j]->nNode())
		*pr++ = mesh->elist[j]->Node[i]+1;
	    else
		*pr++ = 0;

    plhs[0] = vtx;
    plhs[1] = idx;

    // element type list
    if (nlhs > 2) {
	mxArray *eltp = mxCreateDoubleMatrix (elen, 1, mxREAL);
	pr = mxGetPr (eltp);

	for (i = 0; i < elen; i++)
	    *pr++ = (double)mesh->elist[i]->Type();

	plhs[2] = eltp;
    }
}

// =========================================================================

void MatlabToast::SurfData (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    int i, j, k;
    double *pr;

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();
    int nbnd = mesh->nlist.NumberOf (BND_ANY);

    // vertex coordinate list
    mxArray *vtx = mxCreateDoubleMatrix (nbnd, dim, mxREAL);
    pr = mxGetPr (vtx);

    for (i = 0; i < dim; i++)
	for (j = 0; j < nlen; j++)
	    if (mesh->nlist[j].isBnd())
		*pr++ = mesh->nlist[j][i];

	if (nlhs >= 2) {
	    int *bndidx = new int[nlen];

		for (j = k = 0; j < nlen; j++)
			bndidx[j] = (mesh->nlist[j].isBnd() ? k++ : -1);
    
		// boundary element index list
		// note: this currently assumes that all elements contain the
		// same number of vertices!

		int nnd = 0, nface, sd, nn, nd, bn, *bndellist, *bndsdlist;
		nface = mesh->BoundaryList (&bndellist, &bndsdlist);
		for (j = 0; j < nface; j++)
			nnd = ::max (nnd, mesh->elist[bndellist[j]]->nSideNode(bndsdlist[j]));
		mxArray *idx = mxCreateDoubleMatrix (nface, nnd, mxREAL);
		pr = mxGetPr (idx);
    
		for (i = 0; i < nnd; i++)
			for (j = 0; j < nface; j++) {
				Element *pel = mesh->elist[bndellist[j]];
				sd = bndsdlist[j];
				nn = pel->nSideNode (sd);
				if (i < nn) {
				nd = pel->Node[pel->SideNode (sd, i)];
				bn = bndidx[nd]+1;
			} else bn = 0;
			*pr++ = bn;
		}

		plhs[0] = vtx;
		plhs[1] = idx;

		if (nlhs >= 3) {
			// generate nodal permutation index list
			mxArray *perm = mxCreateDoubleMatrix (nbnd, 1, mxREAL);
			pr = mxGetPr (perm);
			for (i = 0; i < nlen; i++) {
				if (bndidx[i] >= 0)
				*pr++ = i+1;
			}
			plhs[2] = perm;
		}

		// cleanup
		delete []bndidx;
	}
}

// =========================================================================

void MatlabToast::MarkMeshBoundary (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    const char id[3] = {0,3,2};
    const double rid[4] = {0,-1,2,1};

    Mesh *mesh = GETMESH_SAFE(0);

    int i, n = mesh->nlen();
    RVector bnd(n);
    if (nrhs > 1) {
	CopyVector (bnd,prhs[1]);
	for (i = 0; i < n; i++)
	    mesh->nlist[i].SetBndTp (id[(int)(bnd[i]+0.5)]);
    } else {
	mesh->MarkBoundary();
    }
    if (nlhs > 0) {
	for (i = 0; i < n; i++)
	    bnd[i] = rid[mesh->nlist[i].BndTp()];
	CopyVector (&plhs[0], bnd);
    }
}

// =========================================================================

void MatlabToast::MeshBB (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i;

    Mesh *mesh = GETMESH_SAFE(0);

    int dim = mesh->Dimension();
    Point pmin(dim), pmax(dim);
    mesh->BoundingBox (pmin, pmax);

    plhs[0] = mxCreateDoubleMatrix (2, dim, mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < dim; i++) {
	*pr++ = pmin[i];
	*pr++ = pmax[i];
    }
}

// =========================================================================

void MatlabToast::MeshSize (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateScalarDouble (mesh->FullSize());
}

// =========================================================================

void MatlabToast::MeshDimension (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateScalarDouble (mesh->Dimension()); 
}

// =========================================================================

void MatlabToast::ElementSize (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    int n = mesh->elen();
    plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    double *pr = mxGetPr (plhs[0]);
    for (int i = 0; i < n; i++)
	pr[i] = mesh->ElSize(i);
}

// =========================================================================

void MatlabToast::ElementData (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;
    double *pr;

    Mesh *mesh = GETMESH_SAFE(0);

    // element index
    int el = (int)mxGetScalar (prhs[1]) - 1;
    if (el < 0 || el >= mesh->elen()) { // invalid index
	plhs[0] = mxCreateDoubleMatrix (0, 0, mxREAL);
	return;
    }

    Element *pel = mesh->elist[el];
    int dim = pel->Dimension();
    int nnd = pel->nNode();
    int nsd = pel->nSide();
    int nsn = pel->nSideNode(0);
    for (i = 1; i < nsd; i++)
	nsn = ::max (nsn, pel->nSideNode(i));

    mxArray *vtx = mxCreateDoubleMatrix (nnd, dim, mxREAL);
    pr = mxGetPr (vtx);
    for (i = 0; i < dim; i++)
	for (j = 0; j < nnd; j++)
	    *pr++ = mesh->nlist[pel->Node[j]][i];

    mxArray *idx = mxCreateDoubleMatrix (nsd, nsn, mxREAL);
    pr = mxGetPr (idx);
    for (i = 0; i < nsn; i++)
	for (j = 0; j < nsd; j++)
	    if (i < pel->nSideNode(j))
		*pr++ = pel->Node[pel->SideNode(j,i)]+1;
	    else
		*pr++ = 0;
    
    plhs[0] = vtx;
    plhs[1] = idx;
}
