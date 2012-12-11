// ========================================================================
// Implementation of class MatlabToast
// mesh-related methods
// NOTE: The core mesh optimisation methods should be moved from here into
// the Mesh class.
// ========================================================================

#include "matlabtoast.h"
#include "toastmex.h"

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

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(mesh);
}

// =========================================================================

void MatlabToast::MakeMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
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
	    mexErrMsgTxt ("Element type not supported!\n");
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

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(mesh);
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
    plhs[0] = mxCreateDoubleScalar (mesh->nlen());
}

// =========================================================================

void MatlabToast::MeshElementCount (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
	Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateDoubleScalar (mesh->elen());
}

// =========================================================================

void MatlabToast::ClearMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    delete mesh;
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
    plhs[0] = mxCreateDoubleScalar (mesh->FullSize());
}

// =========================================================================

void MatlabToast::MeshDimension (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateDoubleScalar (mesh->Dimension()); 
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
    if (nlhs > 2)
        plhs[2] = mxCreateDoubleScalar(pel->Type());
}

// =========================================================================

void MatlabToast::FindElement (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    int i, j, n, dim = mesh->Dimension();
    RDenseMatrix pt;
    CopyMatrix (pt, prhs[1]);
    ASSERTARG(pt.nCols() == dim, 2, "Invalid dimension");
    n = pt.nRows();
    plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    double *pr = mxGetPr(plhs[0]);
    Point p(dim);
    for (i = 0; i < n; i++) {
	for (j = 0; j < dim; j++) p[j] = pt(i,j);
	*pr++ = mesh->ElFind(p)+1;
    }
}

// =========================================================================

void MatlabToast::Elmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];
    mxArray *elmat;
    double *pr;

    Mesh *mesh = GETMESH_SAFE(0);

    int idx = (int)mxGetScalar(prhs[1]) - 1; // convert to zero-based
    Element *pel = mesh->elist[idx];
    int i, j, k, l, nnd = pel->nNode();
    int dim = mesh->Dimension();

    mxGetString (prhs[2], cbuf, 256);

    if (!strcmp(cbuf, "F")) {
	elmat = mxCreateDoubleMatrix (nnd, 1, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    pr[i] = pel->IntF(i);
    } else if (!strcmp(cbuf, "FF")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = pel->IntFF(i,i);
	    for (j = 0; j < i; j++)
		pr[i*nnd+j] = pr[j*nnd+i] = pel->IntFF(i,j);
	}
    } else if (!strcmp(cbuf, "FFF")) {
	mwSize dims[3] = {nnd,nnd,nnd};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    for (j = 0; j < nnd; j++) {
		for (k = 0; k < nnd; k++) {
		    pr[(i*nnd+j)*nnd+k] = pel->IntFFF(i,j,k);
		}
	    }
	}
    } else if (!strcmp(cbuf, "DD")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = pel->IntDD(i,i);
	    for (j = 0; j < i; j++)
		pr[i*nnd+j] = pr[j*nnd+i] = pel->IntDD(i,j);
	}
    } else if (!strcmp(cbuf, "FD")) {
	mwSize dims[3] = {nnd, nnd, dim};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    for (j = 0; j < nnd; j++) {
		RVector fd = pel->IntFD(i,j);
		if (fd.Dim() == 0) {
		    plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		    return;
		}
		for (k = 0; k < dim; k++)
		    pr[i + nnd*(j + k*nnd)] = fd[k];
	    }
    } else if (!strcmp(cbuf, "FDD")) {
	mwSize dims[3] = {nnd, nnd, nnd};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    for (j = 0; j < nnd; j++)
		for (k = 0; k < nnd; k++)
		    pr[i+nnd*(j+nnd*k)] = pel->IntFDD(i,j,k);
    } else if (!strcmp(cbuf, "dd")) {
	mwSize dims[4] = {nnd, dim, nnd, dim};
	elmat = mxCreateNumericArray (4, dims, mxDOUBLE_CLASS, mxREAL);
	RSymMatrix intdd = pel->Intdd();
	pr = mxGetPr(elmat);
	for (l = 0; l < dim; l++)
	    for (k = 0; k < nnd; k++)
		for (j = 0; j < dim; j++)
		    for (i = 0; i < nnd; i++)
			*pr++ = intdd(i*dim+j,k*dim+l);

    } else if (!strcmp(cbuf, "PFF")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	RVector prm;
	CopyVector (prm, prhs[3]);
	RSymMatrix intPFF = pel->IntPFF(prm);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = intPFF(i,i);
	    for (j = 0; j < i; j++)
	        pr[i*nnd+j] = pr[j*nnd+i] = intPFF(i,j);
	}
    } else if (!strcmp(cbuf, "PDD")) {
        elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	RVector prm;
	CopyVector (prm, prhs[3]);
	RSymMatrix intPDD = pel->IntPDD(prm);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = intPDD(i,i);
	    for (j = 0; j < i; j++)
	        pr[i*nnd+j] = pr[j*nnd+i] = intPDD(i,j);
	}

    } else if (!strcmp(cbuf, "BndF")) {

	// for now, only integrals over all boundary sides are supported
	elmat = mxCreateDoubleMatrix (nnd, 1, mxREAL);
	pr = mxGetPr(elmat);
	RVector bndintf = pel->BndIntF();
	for (i = 0; i < pel->nNode(); i++)
	    pr[i] = bndintf[i];

    } else if (!strcmp(cbuf, "BndFF")) {
	int ii, jj, sd;
	if (nrhs > 3) sd = (int)mxGetScalar(prhs[3]) - 1; // side index
	else sd = -1;
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);

	if (sd >= 0) { // integral over a single side
	    for (ii = 0; ii < pel->nSideNode(sd); ii++) {
		i = pel->SideNode(sd,ii);
		pr[i*nnd+i] = pel->BndIntFFSide(i,i,sd);
		for (jj = 0; jj < ii; jj++) {
		    j = pel->SideNode(sd,jj);
		    pr[i*nnd+j] = pr[j*nnd+i] = pel->BndIntFFSide(i,j,sd);
		}
	    }
	} else { // integral over all boundary sides
	    for (sd = 0; sd < pel->nSide(); sd++) {
		if (!pel->IsBoundarySide (sd)) continue;
		for (ii = 0; ii < pel->nSideNode(sd); ii++) {
		    i = pel->SideNode(sd,ii);
		    pr[i*nnd+i] += pel->BndIntFFSide(i,i,sd);
		    for (jj = 0; jj < ii; jj++) {
			j = pel->SideNode(sd,jj);
			pr[i*nnd+j] = pr[j*nnd+i] += pel->BndIntFFSide(i,j,sd);
		    }
		}
	    }
	}
    } else if (!strcmp(cbuf, "BndPFF")) {
        int sd;
        if (nrhs > 4) sd =  (int)mxGetScalar(prhs[4]) - 1; // side index
	else sd = -1;
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	RVector prm;
	CopyVector (prm, prhs[3]);

	if (sd >= 0) { // integral over a single side
	    xERROR("Not implemented yet!");
	} else { // integral over all boundary sides
	    for (i = 0; i < nnd; i++) {
	        for (j = 0; j < nnd; j++)
		    pr[i*nnd+j] = pel->BndIntPFF(i,j,prm);
	    }
	}
    } else {
	mexErrMsgTxt ("Elmat: Integral type string not recognised");
    }
    plhs[0] = elmat;
}

// =========================================================================

void MatlabToast::ReadQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char qmname[256];

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    if (mxIsChar (prhs[1]))
	mxGetString (prhs[1], qmname, 256);
    else
	mexErrMsgTxt ("ReadQM: Argument 1: file name expected.");
	
    if (fileExists (qmname) != 1)
	mexErrMsgTxt ("ReadQM: QM file not found.");

    ifstream ifs(qmname);
    mesh->LoadQM (ifs);

    if (verbosity >= 1)
	mexPrintf ("QM: %d sources, %d detectors, %d measurements\n",
	       mesh->nQ, mesh->nM, mesh->nQM);

}

// =========================================================================

void MatlabToast::SetQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    size_t i, j, d;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    size_t dim = (size_t)mesh->Dimension();

    // Copy source positions
    size_t nq = mxGetM(prhs[1]);
    size_t dimq = mxGetN(prhs[1]);
    d = std::min(dim,dimq);
    if (dim != dimq) {
	cerr << "Warning: toastSetQM: param 2:" << endl;
	cerr << "source array size(" << nq << "," << dimq
	     << ") does not correspond\nto mesh dimension ("
	     << dim << ")" << endl;
    }
    RDenseMatrix Q(nq,dimq);
    CopyMatrix (Q, prhs[1]);
    Point *pQ = new Point[nq];
    for (i = 0; i < nq; i++) {
	pQ[i].New(dim);
	for (j = 0; j < d; j++) pQ[i][j] = Q(i,j);
    }

    // Copy detector positions
    size_t nm = mxGetM(prhs[2]);
    size_t dimm = mxGetN(prhs[2]);
    d = std::min(dim,dimm);
    if (dim != dimm) {
	cerr << "Warning: toastSetQM: param 3:" << endl;
	cerr << "detector array size(" << nm << "," << dimm
	     << ") does not correspond\nto mesh dimension ("
	     << dim << ")" << endl;
    }
    RDenseMatrix M(nm,dimm);
    CopyMatrix (M, prhs[2]);
    Point *pM = new Point[nm];
    for (i = 0; i < nm; i++) {
	pM[i].New(dim);
	for (j = 0; j < d; j++) pM[i][j] = M(i,j);
    }

    // Assign Q and M arrays to mesh
    mesh->SetupQM (pQ, nq, pM, nm);

    if (verbosity >= 1) {
	mexPrintf ("QM: %d sources, %d detectors, %d measurements\n",
		   mesh->nQ, mesh->nM, mesh->nQM);
    }
}

// ============================================================================

void Integrate_Lin_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_u1)
{
    double arg1 = 2.0*a / (Pi*Pi*(x0-x1));

    int_cos_u0 = arg1 * (-2*a * cos(Pi*(d-x0)/(2*a)) +
			 2*a * cos(Pi*(d-x1)/(2*a)) +
			 Pi * (x0-x1) * sin(Pi*(d-x0)/(2*a)));
    int_cos_u1 = arg1 * (2*a * cos(Pi*(d-x0)/(2*a)) -
			 2*a * cos(Pi*(d-x1)/(2*a)) -
			 Pi * (x0-x1) * sin(Pi*(d-x1)/(2*a)));
}

// ----------------------------------------------------------------------------

CVector CompleteTrigSourceVector (const Mesh &mesh, int order)
{
    // currently only works with 2D circular mesh centered at origin
    int el, sd, nnode, *node;
    int n = mesh.nlen();
    double phi0, phi1, rad, a, f0, f1;
    Element *pel;
    CVector qvec (n);

    for (el = 0; el < mesh.elen(); el++) {
	pel = mesh.elist[el];
	xASSERT (pel->Type() == ELID_TRI3, "Element type not supported");
	nnode = pel->nNode();
	node  = pel->Node;
	for (sd = 0; sd < pel->nSide(); sd++) {
	    if (!pel->IsBoundarySide (sd)) continue;
	    Node &nd0 = mesh.nlist[node[pel->SideNode (sd, 0)]];
	    Node &nd1 = mesh.nlist[node[pel->SideNode (sd, 1)]];
	    phi0 = atan2 (nd0[1], nd0[0]);
	    phi1 = atan2 (nd1[1], nd1[0]);

	    if (fabs (phi0-phi1) > Pi) {
		if (phi1 > phi0) phi0 += 2.0*Pi;
		else             phi1 += 2.0*Pi;
	    }
	    if (order) {
		a    = 2.0*Pi/4.0/order;
		Integrate_Lin_Cosine (0, a, phi0, phi1, f0, f1);
	    } else {
		f0 = f1 = 0.0;
	    }
	    f0 += fabs (phi1-phi0);
	    f1 += fabs (phi1-phi0);
	    qvec[node[pel->SideNode(sd,0)]] += f0;
	    qvec[node[pel->SideNode(sd,1)]] += f1;
	}
    }
    return qvec;
}

// ----------------------------------------------------------------------------

void MatlabToast::Qvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char typestr[256] = "";
    char profstr[256] = "";
    double w = 0.0;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    int i, j, idx, n = mesh->nlen(), nQ = mesh->nQ;

    // read source parameters from function parameters
    if (mxIsChar(prhs[1])) mxGetString (prhs[1], typestr, 256);
    if (mxIsChar(prhs[2])) mxGetString (prhs[2], profstr, 256);
    if (nrhs >= 4) {
	idx = 3;
	if (mxIsNumeric(prhs[idx]) && mxGetNumberOfElements(prhs[idx])==1){
	    w = mxGetScalar (prhs[idx]);
	    idx++;
	}
	// additional optional parameters to go here
    }

    SourceMode qtype;
    if      (!strcasecmp (typestr, "Neumann"))   qtype = SRCMODE_NEUMANN;
    else if (!strcasecmp (typestr, "Isotropic")) qtype = SRCMODE_ISOTROPIC;
    else    mexErrMsgTxt ("toastQvec: Invalid source type");

    SRC_PROFILE qprof;
    if      (!strcasecmp (profstr, "Point"))     qprof = PROF_POINT;
    else if (!strcasecmp (profstr, "Gaussian"))  qprof = PROF_GAUSSIAN;
    else if (!strcasecmp (profstr, "Cosine"))    qprof = PROF_COSINE;
    else if (!strcasecmp (profstr, "TrigBasis")) qprof = PROF_COMPLETETRIG;
    else    mexErrMsgTxt ("toastQvec: Invalid source profile");

    double qwidth;
    if (qprof != PROF_POINT) {
	if   (w > 0) qwidth = w;
	else mexErrMsgTxt ("toastQvec: Invalid source width");
    }

    CCompRowMatrix qvec;

    // build the source vectors
    qvec.New (nQ, n);

    for (i = 0; i < nQ; i++) {
	CVector q(n);
	switch (qprof) {
	case PROF_POINT:
	    SetReal (q, QVec_Point (*mesh, mesh->Q[i], qtype));
	    break;
	case PROF_GAUSSIAN:
	    SetReal (q, QVec_Gaussian (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COSINE:
	    SetReal (q, QVec_Cosine (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COMPLETETRIG:
	    q = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	qvec.SetRow (i, q);
    }

    // return source vectors as matrix columns to matlab
    CopyTMatrix (&plhs[0], qvec);
}

// =========================================================================

void MatlabToast::Mvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];

    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    if (!mesh) mexErrMsgTxt ("DataLinkList: Mesh not found.");

    SRC_PROFILE mprof;
    mxGetString (prhs[1], cbuf, 256);
    if      (!strcasecmp (cbuf, "Gaussian"))  mprof = PROF_GAUSSIAN;
    else if (!strcasecmp (cbuf, "Cosine"))    mprof = PROF_COSINE;
    else if (!strcasecmp (cbuf, "TrigBasis")) mprof = PROF_COMPLETETRIG;
    else mexErrMsgTxt ("Invalid measurement profile");

    double mwidth = mxGetScalar (prhs[2]);

    int n, nM;
    int i, j, k, idx, ofs, resettp, cmd;

    n = mesh->nlen();
    nM = mesh->nM;
    CCompRowMatrix mvec;

    // build the measurement vectors
    mvec.New (nM, n);
    for (i = 0; i < nM; i++) {
	CVector m(n);
	switch (mprof) {
	case PROF_GAUSSIAN:
	    SetReal (m, QVec_Gaussian (*mesh, mesh->M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case PROF_COSINE:
	    SetReal (m, QVec_Cosine (*mesh, mesh->M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	case PROF_COMPLETETRIG:
	    m = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	for (j = 0; j < n; j++) m[j] *= mesh->plist[j].C2A();
	mvec.SetRow (i, m);
    }

    // return source vectors as matrix columns to matlab
    CopyTMatrix (&plhs[0], mvec);
}

// =========================================================================

void MatlabToast::QPos (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    
    int nq = mesh->nQ;
    int dim = mesh->Dimension();
    
    mxArray *qpos = mxCreateDoubleMatrix (nq, dim, mxREAL);
    double *pr = mxGetPr (qpos);

    for (j = 0; j < dim; j++) {
	for (i = 0; i < nq; i++) {
	    *pr++ = mesh->Q[i][j];
	}
    }

    plhs[0] = qpos;
}

// =========================================================================

void MatlabToast::MPos (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    
    int nm = mesh->nM;
    int dim = mesh->Dimension();
    
    mxArray *mpos = mxCreateDoubleMatrix (nm, dim, mxREAL);
    double *pr = mxGetPr (mpos);

    for (j = 0; j < dim; j++) {
	for (i = 0; i < nm; i++) {
	    *pr++ = mesh->M[i][j];
	}
    }

    plhs[0] = mpos;
}

