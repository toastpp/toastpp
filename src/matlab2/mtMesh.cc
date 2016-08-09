
// ========================================================================
// Implementation of class MatlabToast
// mesh-related methods
// NOTE: The core mesh optimisation methods should be moved from here into
// the Mesh class.
// ========================================================================

#include "matlabtoast.h"
#include "mex.h"
#include "mexutil.h"

using namespace std;

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
    mexLock(); // prevent mex file unloading while mesh is allocated

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

    // check mesh consistency
    if (mesh->Shrink()) {
	mexWarnMsgTxt ("toastMesh: removed unused nodes");
	nvtx = mesh->nlen();
    }
    
    // set up mesh
    mesh->Setup();

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(mesh);
}

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
    mexLock(); // prevent mex file unloading while mesh is allocated

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

void MatlabToast::WriteMeshVtk (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    char meshname[256];
    
    if (mxIsChar (prhs[1]))
	mxGetString (prhs[1], meshname, 256);
    else
	mexErrMsgTxt ("WriteMesh: Argument 2: file name expected.");
	
    RVector nim(mesh->nlen());
    if (nrhs > 2)
	CopyVector(nim, prhs[2]);

    ofstream ofs (meshname);
    mesh->WriteVtk(ofs, nim);

    if (verbosity >= 1)
	mexPrintf("Mesh: Vtk format written to %s\n", meshname);
}

// =========================================================================

void MatlabToast::MeshOpt (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, res, len, *perm;
    char optmode[256] = "\0";

    Mesh *mesh = GETMESH_SAFE(0);
    ASSERTARG_CHAR(1);
    mxGetString (prhs[1], optmode, 256);

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

void MatlabToast::MeshReorder (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    double *pr = mxGetPr(prhs[1]);

    int nlen = mesh->nlen();
    IVector perm(nlen);
    for (int i = 0; i < nlen; i++)
	perm[i] = (int)(pr[i]-0.5); // switch to 0-based
    
    mesh->Reorder(perm);
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

void MatlabToast::MeshDimension (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    plhs[0] = mxCreateDoubleScalar (mesh->Dimension()); 
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

void MatlabToast::ElementNeighbours (int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    int ne = mesh->elen();
    mesh->PopulateNeighbourLists();
    int i, j;
    int nnbmax = 0;
    for (i = 0; i < ne; i++) {
	if (mesh->elist[i]->nSide() > nnbmax)
	    nnbmax = mesh->elist[i]->nSide();
    }
    RDenseMatrix nb(ne, nnbmax);
    for (i = 0; i < ne; i++) {
	Element *pel = mesh->elist[i];
	int ns = pel->nSide();
	for (j = 0; j < ns; j++) {
	    nb(i, j) = pel->sdnbhridx[j]+1; // one-based
	}
    }
    CopyMatrix(&plhs[0], nb);
}

// =========================================================================

void MatlabToast::SysmatSparsityStructure (int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    int n = mesh->nlen();
    int *rowptr, *colidx, nzero;
    double *ndata;
    
    mesh->SparseRowStructure (rowptr, colidx, nzero);
    
    ndata = new double[nzero];
    for(int i=0; i<nzero; i++) ndata[i] = 1.0;
    
    RCompRowMatrix F(n, n, rowptr, colidx, ndata);
    
    delete []rowptr;
    delete []colidx;
    delete []ndata;
    
    CopyMatrix ( &plhs[0], F);
}


// =========================================================================

void CalcSysmatComponent (QMMesh *mesh, RVector &prm, char *integral_type,
    bool elbasis, bool prmint, mxArray **res, RCompRowMatrix &F)
{
    int n = mesh->nlen();

    // If the supplied matrix has the correct dimension, assume it is
    // of the correct sparsity structure, else compute and reinitialise.
    if((F.nRows() == n) && (F.nCols() == n)) {
        
        F.Zero();
        
    }
    else {
        
        int *rowptr, *colidx, nzero;
        
        mesh->SparseRowStructure (rowptr, colidx, nzero);
        F.New(n, n);
        F.Initialise(rowptr, colidx);
        delete []rowptr;
        delete []colidx;
        
    }
    
    if(!prmint) {
        
        if (!strcasecmp (integral_type, "FF")) {
            AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_FF);
        } else if (!strcasecmp (integral_type, "DD")) {
            AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_DD);
        } else {
            mexErrMsgTxt ("Unknown integral string for non-parameter integration");
        }
        
    } else {
        
        if (!strcasecmp (integral_type, "PFF")) {
            AddToSysMatrix (*mesh, F, &prm,
                            elbasis ? ASSEMBLE_PFF_EL : ASSEMBLE_PFF);
        } else if (!strcasecmp (integral_type, "PDD")) {
            AddToSysMatrix (*mesh, F, &prm,
                            elbasis ? ASSEMBLE_PDD_EL : ASSEMBLE_PDD);
        } else if (!strcasecmp (integral_type, "BNDPFF")) {
            AddToSysMatrix (*mesh, F, &prm,
                            elbasis ? ASSEMBLE_BNDPFF_EL : ASSEMBLE_BNDPFF);
        } else {
            mexErrMsgTxt ("Unknown integral string for parameter integration");
        }
        
    }
    
    // Return system matrix to MATLAB
    CopyMatrix (res, F);
}

void MatlabToast::SysmatComponent (int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    RVector prm;
    RCompRowMatrix F;
    bool prmint = false;            // Is this a parameter integral
    bool elbasis = false;
    char integral_type[32] = "";

    // First argument is the required integral type
    if (nrhs >= 2 && mxIsChar(prhs[1])) {
	mxGetString (prhs[1], integral_type, 32);
    } else {
	mexErrMsgTxt ("Parameter 2: string expected");
    }
    
    // Second argument is the integration parameter, or empty
    if(nrhs >=3) {
        
        prmint = true;
        CopyVector (prm, prhs[2]);

        if(prm.Dim() == mesh->nlen())  elbasis = false;
        else if (prm.Dim() == mesh->elen()) elbasis = true;
        else prmint = false;
   
    }
    
    // Thurd argument is the system matrix sparsity pattern, or empty
    if (nrhs == 4 && mxIsSparse(prhs[3])) CopyTMatrix(F, prhs[3]);
    
    CalcSysmatComponent (mesh, prm, integral_type, elbasis, prmint, &plhs[0], F);
}

// =========================================================================

void MatlabToast::Massmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);
    RCompRowMatrix *B = mesh->MassMatrix();
    CopyMatrix (&plhs[0], *B);
    delete B;
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
}

// =========================================================================

void MatlabToast::WriteQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    RDenseMatrix qpos, mpos;
    RCompRowMatrix lnk;
    char qmname[256];
    int i, j;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    mxGetString (prhs[1], qmname, 256);

    CopyMatrix (qpos, prhs[2]);
    if (qpos.nRows() == 0 && qpos.nCols() == 0) {
	qpos.New(mesh->nQ,mesh->Dimension());
	for (i = 0; i < mesh->nQ; i++)
	    for (j = 0; j < mesh->Dimension(); j++)
		qpos(i,j) = mesh->Q[i][j];
    }

    CopyMatrix (mpos, prhs[3]);
    if (mpos.nRows() == 0 && mpos.nCols() == 0) {
	mpos.New(mesh->nM,mesh->Dimension());
	for (i = 0; i < mesh->nM; i++)
	    for (j = 0; j < mesh->Dimension(); j++)
		mpos(i,j) = mesh->M[i][j];
    }

    int dim = qpos.nCols();
    int nq = qpos.nRows();
    int nm = mpos.nRows();
    
    if (mxGetM(prhs[4]) && mxGetN(prhs[4])) {
	ASSERTARG(mxIsSparse(prhs[4]), 4, "expected sparse matrix.");
	CopyTMatrix (lnk, prhs[4]);

	if (nq != lnk.nRows() || nm != lnk.nCols()) {
	    char cbuf[256];
	    sprintf (cbuf, "Invalid dimension (was: %d x %d, expected %d x %d)",
		     lnk.nCols(), lnk.nRows(), nm, nq);
	    ASSERTARG(0, 4, cbuf);
	}
    } else {
	lnk.New(nq, nm);
	for (i = 0; i < nq; i++) {
	    RVector qlnk(nm);
	    for (j = 0; j < nm; j++)
		if (mesh->Connected (i,j))
		    qlnk[j] = 1.0;
	    lnk.SetRow (i, qlnk);
	}
    }

    ofstream ofs(qmname);
    ofs << "QM file" << endl;
    ofs << "Dimension " << dim << endl << endl;
    ofs << "SourceList " << nq << endl;
    for (i = 0; i < nq; i++)
	for (j = 0; j < dim; j++)
	    ofs << qpos(i,j) << (j == dim-1 ? '\n':' ');
    ofs << endl;
    ofs << "MeasurementList " << nm << endl;
    for (i = 0; i < nm; i++)
	for (j = 0; j < dim; j++)
	    ofs << mpos(i,j) << (j == dim-1 ? '\n':' ');
    ofs << endl;

    ofs << "LinkList" << endl;
    int *ci = new int[nm];
    double *idx = new double[nm];

    for (i = 0; i < nq; i++) {
	int nz = lnk.SparseRow (i, ci, idx);
	ofs << nz << ':';
	for (j = 0; j < nz; j++)
	    ofs << ci[j] << (j == nz-1 ? '\n':' ');
    }
}

// =========================================================================

void MatlabToast::GetQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);
    int q, m, n, idx;
    int nq = mesh->nQ;
    int nm = mesh->nM;
    int nqm = mesh->nQM;

    mxArray *lnk = mxCreateSparse (nm, nq, nqm, mxREAL);
    double  *pr = mxGetPr (lnk);
    mwIndex *ir = mxGetIr (lnk);
    mwIndex *jc = mxGetJc (lnk);

    *jc++ = 0;
    for (q = idx = 0; q < nq; q++) {
	for (m = n = 0; m < nm; m++) {
	    if (!mesh->Connected (q,m)) continue;
	    *pr++ = ++idx;
	    *ir++ = m;
	    n++;
	}
	*jc = *(jc-1) + n;
	jc++;
    }

    plhs[0] = lnk;
}

// =========================================================================

void MatlabToast::DataLinkList (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    int q, m, idx;
    plhs[0] = mxCreateDoubleMatrix (1, mesh->nQM, mxREAL);
    double *pr = mxGetPr (plhs[0]);
    for (q = idx = 0; q < mesh->nQ; q++)
	for (m = 0; m < mesh->nM; m++)
	    if (mesh->Connected (q, m))
		pr[idx++] = q*mesh->nM + m + 1;    
}

// =========================================================================

void MatlabToast::ClearMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);
    delete mesh;
    mexUnlock();

    if (verbosity >= 1)
        mexPrintf ("<Mesh object deleted>\n");
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

void MatlabToast::SetQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    size_t i, j, d;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    size_t dim = (size_t)mesh->Dimension();

    // Copy source positions
    mwSize nq = mxGetM(prhs[1]);
    mwSize dimq = mxGetN(prhs[1]);
    d = std::min(dim, dimq);
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

// =========================================================================

void MatlabToast::QPos (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);
    
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

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);
    
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

// =========================================================================

void MatlabToast::MeshRefine (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    void(*RefineFunc)(Mesh*,bool*);
    Mesh *mesh = (Mesh*)GETMESH_SAFE(0);
    int i, nv, nlen = mesh->nlen(), elen = mesh->elen();
    BYTE eltp = mesh->elist[0]->Type();
    switch (eltp) {
    case ELID_TRI3:
        RefineFunc = RefineTriangle3Mesh;
	break;
    case ELID_TET4:
        RefineFunc = RefineTetrahedron4Mesh;
	break;
    default:
        mexErrMsgTxt("MeshRefine: Element type not supported");
	return;
    }
    for (i = 1; i < nlen; i++)
        if (mesh->elist[i]->Type() != eltp) {
	    mexErrMsgTxt("MeshRefine: Mesh must consist of elements of only one type");
	}

    bool *elrefine = new bool[elen];
    for (i = 0; i < elen; i++)
        elrefine[i] = false;

    if (nrhs > 1) {
        RVector v;
	CopyVector (v, prhs[1]);
	nv = v.Dim();
	for (i = 0; i < nv; i++) {
	    int el = (int)(v[i]-0.5);
	    elrefine[el] = true;
	}
    } else {
        for (i = 0; i < elen; i++)
	    elrefine[i] = true;
    }

    RefineFunc (mesh, elrefine);
}

// =========================================================================

void MatlabToast::SplitElement (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // currently only works for Triangle3 meshes.

    Mesh *mesh = (Mesh*)GETMESH_SAFE(0);
    int el = (int)(mxGetScalar(prhs[1])-0.5); // switch to 0-based

    Mesh_SplitTriangle3 (mesh, el);
}

// =========================================================================

void MatlabToast::NodeNeighbour (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int *nnbrs, **nbrs;
    Mesh *mesh = (Mesh*)GETMESH_SAFE(0);
    int i, j, nlen = mesh->nlen();
    mesh->NodeNeighbourList (&nnbrs, &nbrs);

    int maxnbr = 0;
    for (i = 0; i < nlen; i++)
	if (nnbrs[i] > maxnbr)
	    maxnbr = nnbrs[i];

    if (nlhs > 0) {
	plhs[0] = mxCreateNumericMatrix (nlen, maxnbr, mxINT32_CLASS, mxREAL);
	int *ptr = (int*)mxGetData (plhs[0]);
    
	for (i = 0; i < maxnbr; i++) {
	    for (j = 0; j < nlen; j++) {
		*ptr++ = (i < nnbrs[j] ? nbrs[j][i]+1 : 0);
	    }
	}
    }

    for (i = 0; i < nlen; i++)
	delete []nbrs[i];
    delete []nbrs;
    delete []nnbrs;
}

// =========================================================================

void MatlabToast::UnwrapPhase (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    if (nlhs < 1) return; // sanity check

    Mesh *mesh = (Mesh*)GETMESH_SAFE(0);
    int i;
    int nlen = mesh->nlen();
    int dim = mesh->Dimension();
    RVector phase;
    RVector seed;
    CopyVector (phase, prhs[1]);
    CopyVector (seed, prhs[2]);

    if (seed.Dim() != dim)
	mexErrMsgTxt ("Invalid dimension for seed point");

    Point seedpt(dim);
    for (i = 0; i < dim; i++)
	seedpt[i] = seed[i];

    NimPhaseUnwrap (mesh, phase, seedpt);
    plhs[0] = mxCreateDoubleMatrix(nlen,1,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < nlen; i++)
	*pr++ = phase[i];
}

// =========================================================================

template<class T>
TVector<T> IntFG (const Mesh &mesh, const TVector<T> &f, const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), "Wrong vector size");
    dASSERT(g.Dim() == mesh.nlen(), "Wrong vector size");

    int el, nnode, *node, i, j, k, nj, nk, bs;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.nlen());

    for (el = 0; el < mesh.elen(); el++) {
        pel   = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    bs = node[i];
	    for (j = 0; j < nnode; j++) {
	        nj = node[j];
		sum = (f[nj] * g[nj]) * pel->IntFFF(i,j,j);
		for (k = 0; k < j; k++) {
		    nk = node[k];
		    sum += (f[nj]*g[nk] + f[nk]*g[nj]) * pel->IntFFF(i,j,k);
		}
		tmp[bs] += sum;
	    }
	}
    }
    return tmp;
}

// -------------------------------------------------------------------------

void MatlabToast::IntFG (int nlhs, mxArray *plhs[], int nrhs,
                        const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    bool isCplx = mxIsComplex (prhs[1]);
    if (isCplx != mxIsComplex (prhs[2]))
	mexErrMsgTxt
	   ("toastIntGradFGradG: Arguments must be both real or both complex");

    if (isCplx) {
	CVector fvec, gvec, rvec;
	CopyVector (fvec, prhs[1]);
	CopyVector (gvec, prhs[2]);
	rvec = ::IntFG (*mesh, fvec, gvec);
	CopyVector (&plhs[0], rvec);
    } else {
	RVector fvec, gvec, rvec;
	CopyVector (fvec, prhs[1]);
	CopyVector (gvec, prhs[2]);
	rvec= ::IntFG (*mesh, fvec, gvec);
	CopyVector (&plhs[0],rvec);
    }
}

// =========================================================================

template<class T>
TVector<T> IntGradFGradG (const Mesh &mesh, const TVector<T> &f,
    const TVector<T> &g)
{
    dASSERT(f.Dim() == mesh.nlen(), "Wrong vector size");
    dASSERT(g.Dim() == mesh.nlen(), "Wrong vector size");

    int el, nnode, *node, i, j, k, nj, nk, bs;
    T sum;
    Element *pel;
    TVector<T> tmp(mesh.nlen());

    for (el = 0; el < mesh.elen(); el++) {
	pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    bs = node[i];
	    for (j = 0; j < nnode; j++) {
		nj = node[j];
		sum = (f[nj] * g[nj]) * pel->IntFDD (i,j,j);
		for (k = 0; k < j; k++) {
		    nk = node[k];
		    sum += (f[nj]*g[nk] + f[nk]*g[nj]) * pel->IntFDD (i,j,k);
		}
		// we exploit the fact that IntFDD(i,j,k) is symmetric in
		// j and k: IntFDD(i,j,k) = IntFDD(i,k,j), so that two terms
		// can be combined in each pass of the inner (k) loop
		tmp[bs] += sum;
	    }
	}
    }
    return tmp;
}

// -------------------------------------------------------------------------

void MatlabToast::IntGradFGradG (int nlhs, mxArray *plhs[], int nrhs,
                        const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    bool isCplx = mxIsComplex (prhs[1]);
    if (isCplx != mxIsComplex (prhs[2]))
	mexErrMsgTxt
	   ("toastIntGradFGradG: Arguments must be both real or both complex");

    if (isCplx) {
	CVector fvec, gvec, rvec;
	CopyVector (fvec, prhs[1]);
	CopyVector (gvec, prhs[2]);
	rvec = ::IntGradFGradG (*mesh, fvec, gvec);
	CopyVector (&plhs[0], rvec);
    } else {
	RVector fvec, gvec, rvec;
	CopyVector (fvec, prhs[1]);
	CopyVector (gvec, prhs[2]);
	rvec = ::IntGradFGradG (*mesh, fvec, gvec);
	CopyVector (&plhs[0], rvec);
    }
}
