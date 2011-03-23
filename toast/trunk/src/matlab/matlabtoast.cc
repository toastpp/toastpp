#include "mex.h"
#include "stoastlib.h"
#include "util.h"
#include "toastmex.h"
#include "matlabtoast.h"
#include "timing.h"

using namespace std;
using namespace toast;

// =========================================================================

void AssertArg (bool cond, const char *func, int argno, const char *errmsg)
{
    if (!cond) {
	char str[1024];
	if (argno)
	    sprintf (str, "%s: argument %d: %s.", func, argno, errmsg);
	else
	    sprintf (str, "%s: %s.", func, errmsg);
	mexErrMsgTxt (str);
    }
}

// =========================================================================

bool fileExists(const std::string& fileName)
{
  std::fstream fin;
  fin.open(fileName.c_str(),std::ios::in);
  if( fin.is_open() )
  {
    fin.close();
    return true;
  }
  fin.close();
  return false;
}

// =========================================================================

bool ReadNim (char *name, int idx, RVector &img)
{
    char cbuf[256];
    int i, j = 0, imgsize = 0;

    ifstream ifs (name);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM") && strcmp (cbuf, "RIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    for (;;) {
	do {
	    ifs.getline (cbuf, 256);
	} while (ifs.good() && strncasecmp (cbuf, "Image", 5));
	if (!ifs.good()) break;
	for (i = 0; i < imgsize; i++)
	    ifs >> img[i];
	if (++j == idx) break;
    }
    return true;
}

// =========================================================================
// =========================================================================
// MatlabToast class implementation

MatlabToast::MatlabToast ()
{
    SetErrorhandler (ErrorHandler);
    nmesh = 0;
    nbasis = 0;
    nreg = 0;
    verbosity = 0;
}

// =========================================================================

MatlabToast::~MatlabToast ()
{
    if (nmesh) delete []meshlist;
    if (nbasis) delete []basislist;
    if (nreg) delete []reglist;
}

// =========================================================================

void MatlabToast::ErrorHandler (char *msg)
{
    mexErrMsgTxt (msg);
}

// =========================================================================

Mesh *MatlabToast::GetMesh (const mxArray *idx, int *errid)
{
    unsigned int meshid;
	if (errid) *errid = 0;

    if (mxIsUint32(idx)) {
		meshid = *(unsigned int*)mxGetData(idx) - 1;
		if (meshid >= nmesh) {
			if (errid) {
				*errid = MESH_INDEXOUTOFRANGE;
			} else if (nmesh) {
				char cbuf[256];
				sprintf (cbuf, "GetMesh: mesh index out of range (expected 1..%d).\n", nmesh);
				mexErrMsgTxt (cbuf);
			} else {
				mexErrMsgTxt ("GetMesh: mesh index out of range (no meshes registered).\n");
			}
			return 0;
		}
		if (!meshlist[meshid]) {
			if (errid) *errid = MESH_INDEXCLEARED;
			else mexErrMsgTxt ("GetMesh: mesh index points to cleared mesh.\n");
		}
		return meshlist[meshid];
    } else {
		if (errid) *errid = MESH_INDEXFORMAT;
		else mexErrMsgTxt ("GetMesh: Invalid mesh index format (expected uint32).\n");
		return 0;
    }
}

// =========================================================================

Mesh *MatlabToast::GetMesh_Safe (const mxArray *arr, const char *func, int argno)
{
	int err;
	Mesh *meshptr = GetMesh(arr, &err);
	switch(err) {
		case MESH_INDEXOUTOFRANGE:
			AssertArg (0, func, argno, "Mesh index out of range");
			break;
		case MESH_INDEXCLEARED:
			AssertArg (0, func, argno, "Index refers to cleared mesh");
			break;
		case MESH_INDEXFORMAT:
			AssertArg (0, func, argno, "Wrong index format (expected uint32)");
			break;
	}
	return meshptr;
}

// =========================================================================

Raster *MatlabToast::GetBasis (const mxArray *idx)
{
    unsigned int basisid;

    if (mxIsUint32(idx)) {
	basisid = *(unsigned int*)mxGetData(idx) - 1;
	if (basisid >= nbasis) {
	    if (nbasis)
		mexPrintf ("GetBasis: index out of range (expected 1..%d).\n", nbasis);
	    else
		mexPrintf ("GetBasis: index out of range (no basis instances registered).\n");
	    return 0;
	}
	if (!basislist[basisid]) {
	    mexPrintf ("GetBasis: index points to cleared basis instance.\n");
	}
	return basislist[basisid];
    } else {
	mexPrintf ("GetBasis: Invalid index format (expected uint32).\n");
	return 0;
    }
}

// =========================================================================

Regularisation *MatlabToast::GetRegul (const mxArray *idx)
{
    unsigned int regid;

    if (mxIsUint32(idx)) {
	regid = *(unsigned int*)mxGetData(idx);
	if (!regid) return 0; // "no regularisation"
	regid--;
	if (regid >= nreg) {
	    if (nreg)
		mexPrintf ("GeRegul: index out of range (expected 1..%d).\n", nbasis);
	    else
		mexPrintf ("GetRegul: index out of range (no regularisation instances registered).\n");
	    return 0;
	}
	if (!reglist[regid]) {
	    mexPrintf ("GetRegul: index points to cleared regularisation instance.\n");
	}
	return reglist[regid];
    } else {
	mexPrintf ("GetRegul: Invalid index format (expected uint32).\n");
	return 0;
    }
}

// =========================================================================

void MatlabToast::SetVerbosity (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    const unsigned int max_verbosity = 1;

    if (mxIsUint32 (prhs[0])) {
	unsigned int verb = *(unsigned int*)mxGetData (prhs[0]);
	if (verb == verbosity) return; // nothing to do
	verbosity = std::min (verb, max_verbosity);
    } else {
	mexErrMsgTxt ("SetVerbosity: Argument 1: invalid format (uint32 expected).");
    }

    if (verbosity >= 1) {
	mexPrintf("Verbosity level set to %d.\n", verbosity);
    }
}

// =========================================================================

void MatlabToast::MeshLin2Quad (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG (mesh, 1, "Mesh not found");

    Mesh *quadmesh = Lin2Quad (*mesh);

    // insert mesh into mesh list
    Mesh **tmp = new Mesh*[nmesh+1];
    if (nmesh) {
	memcpy (tmp, meshlist, nmesh*sizeof(Mesh*));
	delete []meshlist;
    }
    meshlist = tmp;
    meshlist[nmesh++] = quadmesh;

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT32_CLASS, mxREAL);
    unsigned int *ptr = (unsigned int*)mxGetData (plhs[0]);
    *ptr = nmesh;

    if (verbosity >= 1)
	mexPrintf ("Mesh: %d nodes, %d elements, dimension %d\n",
		   quadmesh->nlen(), quadmesh->elen(), quadmesh->Dimension());
}

// =========================================================================

void MatlabToast::ReadQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char qmname[256];

    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

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

    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");
    
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

// =========================================================================

void MatlabToast::GetQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

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

void MatlabToast::WriteQM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    RDenseMatrix qvec, mvec;
    RCompRowMatrix lnk;
    char qmname[256];
    int i, j;

    CopyMatrix (qvec, prhs[0]);
    CopyMatrix (mvec, prhs[1]);

    ASSERTARG(mxIsSparse(prhs[2]), 3, "expected sparse matrix.");
    CopyTMatrix (lnk, prhs[2]);

    int dim = qvec.nCols();
    int nq = qvec.nRows();
    int nm = mvec.nRows();
    
    if (nq != lnk.nRows() || nm != lnk.nCols()) {
	char cbuf[256];
	sprintf (cbuf, "Invalid dimension (was: %d x %d, expected %d x %d)",
		 lnk.nCols(), lnk.nRows(), nm, nq);
	ASSERTARG(0, 3, cbuf);
    }

    mxGetString (prhs[3], qmname, 256);

    ofstream ofs(qmname);
    ofs << "QM file" << endl;
    ofs << "Dimension " << dim << endl << endl;
    ofs << "SourceList " << nq << endl;
    for (i = 0; i < nq; i++)
	for (j = 0; j < dim; j++)
	    ofs << qvec(i,j) << (j == dim-1 ? '\n':' ');
    ofs << endl;
    ofs << "MeasurementList " << nm << endl;
    for (i = 0; i < nm; i++)
	for (j = 0; j < dim; j++)
	    ofs << mvec(i,j) << (j == dim-1 ? '\n':' ');
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

void MatlabToast::DataLinkList (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found.");

    int q, m, idx;
    plhs[0] = mxCreateDoubleMatrix (1, mesh->nQM, mxREAL);
    double *pr = mxGetPr (plhs[0]);
    for (q = idx = 0; q < mesh->nQ; q++)
	for (m = 0; m < mesh->nM; m++)
	    if (mesh->Connected (q, m))
		pr[idx++] = q*mesh->nM + m + 1;    
}

// =========================================================================

void MatlabToast::FindElement (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found.");

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

void MatlabToast::ShapeFunc (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;
    double *pr;

    // mesh handle
    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found.");

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // element index
    int idx = (int)(floor(mxGetScalar (prhs[1])-0.5));
    ASSERTARG(idx >= 0 && idx < elen, 2, "Invalid element index");

    // global point
    Point glob(dim);
    const mwSize *gdim = mxGetDimensions (prhs[2]);
    ASSERTARG(gdim[0]*gdim[1] == dim, 3, "Invalid point dimension");
    pr = mxGetPr (prhs[2]);
    for (i = 0; i < dim; i++) glob[i] = pr[i];

    // calculate shape functions
    RVector fun = mesh->elist[idx]->GlobalShapeF (mesh->nlist, glob);
    int nn = fun.Dim();

    // vertex coordinate list
    mxArray *sf = mxCreateDoubleMatrix (1, nn, mxREAL);
    pr = mxGetPr (sf);
    for (i = 0; i < nn; i++) pr[i] = fun[i];

    plhs[0] = sf;    
}

// =========================================================================

void MatlabToast::ShapeGrad (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;
    double *pr;

    // mesh handle
    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found.");

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // element index
    int idx = (int)(floor(mxGetScalar (prhs[1])-0.5));
    ASSERTARG(idx >= 0 && idx < elen, 2, "Invalid element index");

    // global point
    Point glob(dim);
    const mwSize *gdim = mxGetDimensions (prhs[2]);
    ASSERTARG(gdim[0]*gdim[1] == dim, 3, "Invalid point dimension");
    pr = mxGetPr (prhs[2]);
    for (i = 0; i < dim; i++) glob[i] = pr[i];

    // calculate shape function derivatives
    RDenseMatrix fgrad = mesh->elist[idx]->GlobalShapeD (mesh->nlist, glob);
    int nn = fgrad.nCols();

    // vertex coordinate list
    mxArray *sf = mxCreateDoubleMatrix (dim, nn, mxREAL);
    CopyMatrix (&sf, fgrad);

    plhs[0] = sf;
}

// =========================================================================

void MatlabToast::ReadNIM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char nimname[256];
    int idx;

    if (mxIsChar (prhs[0]))
	mxGetString (prhs[0], nimname, 256);
    else
	mexErrMsgTxt ("ReadNIM: Argument 1: file name expected.");

    if (fileExists (nimname) != 1)
	mexErrMsgTxt ("ReadNIM: Image file not found.");

    if (nrhs < 2) idx = 1;
    else          idx = (int)mxGetScalar (prhs[1]);

    RVector img;
    ReadNim (nimname, idx, img);

    plhs[0] = mxCreateDoubleMatrix (img.Dim(), 1, mxREAL);
    memcpy (mxGetPr (plhs[0]), img.data_buffer(), img.Dim()*sizeof(double));

    if (verbosity >= 1)
	mexPrintf ("NIM: size = %d\n", img.Dim());
}

// ============================================================================

void MatlabToast::WriteNIM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char dname[256];
    mxGetString (prhs[0], dname, 256);
    ofstream ofs (dname);
    ofs << "NIM" << endl;

    mxGetString (prhs[1], dname, 256);
    ofs << "Mesh = " << dname << endl;

    ofs << "SolutionType = N/A" << endl;
    
    int n = mxGetN(prhs[2]) * mxGetM(prhs[2]);
    ofs << "ImageSize = " << n << endl;

    ofs << "EndHeader" << endl;

    ofs << "Image 0" << endl;

    double *val = mxGetPr (prhs[2]);
    ofs.precision(12);
    ofs.setf (ios::scientific);
    for (int i = 0; i < n; i++)
	ofs << val[i] << ' ';
    ofs << endl;
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
	xASSERT (pel->Type() == ELID_TRI3, Element type not supported);
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

// =========================================================================

void MatlabToast::Qvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char typestr[256] = "";
    char profstr[256] = "";
    double w = 0.0;

    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    if (!mesh) mexErrMsgTxt ("DataLinkList: Mesh not found.");
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

// =========================================================================

void CalcSysmat (QMMesh *mesh, RVector &mua, RVector &mus, RVector &ref,
		 double freq, bool elbasis, mxArray **res)
{
    int n = (elbasis ? mesh->elen() : mesh->nlen());

    // Set optical coefficients
    Solution sol (OT_NPARAM, n);
    sol.SetParam (OT_CMUA, mua*c0/ref);
    sol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (int i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    sol.SetParam (OT_C2A, c2a);

    // Create forward solver to initialise system matrix
    CFwdSolver FWS (LSOLVER_DIRECT, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    double omega = freq * 2.0*Pi*1e-6;
    
    FWS.Allocate (*mesh);
    FWS.AssembleSystemMatrix (sol, omega, elbasis);

    // Return system matrix to MATLAB
    CopyMatrix (res, *FWS.F);
}

void CalcBndSysmat (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    CFwdSolver FWS (LSOLVER_DIRECT, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    FWS.Allocate (*mesh);
    RVector prm(n);
    for (int i = 0; i < n; i++) prm[i] = c0/ref[i];
    AddToSysMatrix (*mesh, *FWS.F, &prm, ASSEMBLE_BNDPFF);
    CopyMatrix (res, *FWS.F);
}

void CalcBndFactors (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int n = mesh->nlen();
    RVector prm(n);
    for (int i = 0; i < n; i++)
	prm[i] = A_Keijzer (ref[i]);
    CopyVector (res, prm);
}

void MatlabToast::Sysmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    if (!mesh) mexErrMsgTxt ("Sysmat: Mesh not found.");

    RVector mua, mus, ref;
    double freq;
    int n;
    bool elbasis = false;

    CopyVector (mua, prhs[1]);
    CopyVector (mus, prhs[2]);
    CopyVector (ref, prhs[3]);
    freq = mxGetScalar (prhs[4]);

    if (nrhs >= 6 && mxIsChar(prhs[5])) {
	char cbuf[32];
	mxGetString (prhs[5], cbuf, 32);
	elbasis = (strcasecmp (cbuf, "EL") == 0);
    }

    if (elbasis) n = mesh->elen();
    else         n = mesh->nlen();
    ASSERTARG (n == mua.Dim(), 2, "wrong size");
    ASSERTARG (n == mus.Dim(), 3, "wrong size");
    ASSERTARG (n == ref.Dim(), 4, "wrong size");

    CalcSysmat (mesh, mua, mus, ref, freq, elbasis, &plhs[0]);

    if (nlhs > 1) { // separately provide boundary integral matrix
	CalcBndSysmat (mesh, ref, &plhs[1]);

	if (nlhs > 2) { // boundary matrix prefactors
	    CalcBndFactors (mesh, ref, &plhs[2]);
	}
    }
}

// =========================================================================

void MatlabToast::Massmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n = mesh->nlen();

    // Create forward solver to initialise system matrix
    RFwdSolver FWS (LSOLVER_DIRECT, 1e-10);
    FWS.AssembleMassMatrix (mesh);

    // Return system matrix to MATLAB
    CopyMatrix (&plhs[0], *FWS.B);
}

// =========================================================================

void CalcVolmat (Mesh *mesh, const char *intstr, RVector &prm, bool elbasis,
    mxArray **res)
{
    // sysmatrix structure
    int n = mesh->nlen();
    int *rowptr, *colidx, nzero;
    mesh->SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    if (!strcasecmp (intstr, "FF")) {
	AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_FF);
    } else if (!strcasecmp (intstr, "DD")) {
	AddToSysMatrix (*mesh, F, &prm, ASSEMBLE_DD);
    } else if (!strcasecmp (intstr, "PFF")) {
	AddToSysMatrix (*mesh, F, &prm,
			elbasis ? ASSEMBLE_PFF_EL : ASSEMBLE_PFF);
    } else if (!strcasecmp (intstr, "PDD")) {
	AddToSysMatrix (*mesh, F, &prm,
			elbasis ? ASSEMBLE_PDD_EL : ASSEMBLE_PDD);
    } else if (!strcasecmp (intstr, "BndPFF")) {
	AddToSysMatrix (*mesh, F, &prm,
			elbasis ? ASSEMBLE_BNDPFF_EL : ASSEMBLE_BNDPFF);
    } else if (!strcasecmp (intstr, "BndFF")) {
	RVector dummy(n);
	dummy = 1.0;
	AddToSysMatrix (*mesh, F, &dummy, ASSEMBLE_BNDPFF);
    }
    CopyMatrix (res, F);
}

void MatlabToast::Volmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    char intstr[256];
    RVector prm;
    bool elbasis = false;

    mxGetString (prhs[1], intstr, 256);
    
    if (nrhs > 2) {
	CopyVector (prm, prhs[2]);
	if (nrhs > 3 && mxIsChar (prhs[3])) {
	    char cbuf[32];
	    mxGetString (prhs[3], cbuf, 32);
	    elbasis = (strcasecmp (cbuf, "EL") == 0);
	}
	ASSERTARG((elbasis && prm.Dim()==mesh->elen()) ||
		  (!elbasis && prm.Dim()==mesh->nlen()),3,"Invalid length");
    }

    CalcVolmat (mesh, intstr, prm, elbasis, &plhs[0]);
}

// =========================================================================

#define EQ 0
#define GT 1
#define GE 2
#define LT 3
#define LE 4

bool Inside (Mesh *mesh, int nd, int dim, int rel, double v)
{
    const double EPS = 1e-8;
    double nv;
    if (dim < 3)
	nv = mesh->nlist[nd][dim];
    else
	nv = length(mesh->nlist[nd]);

    switch (rel) {
    case EQ:
	return (fabs(nv-v) < EPS);
    case GT:
	return (nv > v);
    case GE:
	return (nv > v-EPS);
    case LT:
	return (nv < v);
    case LE:
	return (nv < v+EPS);
    default:
	return false;
    }
}

void CalcBndmat (Mesh *mesh, char *intstr, int dim, int rel, double v,
    mxArray **res)
{
    // sysmatrix structure
    int el, i, j, is, js, nside, nnode;
    int n = mesh->nlen();
    int nel = mesh->elen();
    int *rowptr, *colidx, nzero;
    double val;
    bool subreg = (dim >= 0);
    mesh->SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    if (!strcasecmp (intstr, "FF")) {
	for (el = 0; el < nel; el++) {
	    Element *pel = mesh->elist[el];
	    nside = pel->nSide();
	    nnode = pel->nNode();
	    for (i = 0; i < nnode; i++) {
		is = pel->Node[i];
		if (subreg && !Inside (mesh, is, dim, rel, v)) continue;
		for (j = 0; j < nnode; j++) {
		    js = pel->Node[j];
		    if (subreg && !Inside (mesh, js, dim, rel, v)) continue;
		    val = pel->BndIntFF (i, j);
		    F.Add (is, js, val);
		}
	    }
	}
    }
    CopyMatrix (res, F);
}

void MatlabToast::Bndmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    static const char *Relstr[5] = {"EQ", "GT", "GE", "LT", "LE"};
    static const char *Dimstr[4] = {"X", "Y", "Z", "R"};

    char intstr[256];
    RVector prm;
    int n, i;
    bool elbasis = false;
    int rel = -1;
    int dim = -1;
    double val;

    mxGetString (prhs[1], intstr, 256);
    
    if (nrhs > 2) {
	char regionstr[256];
	char dimstr[32], relstr[32];
	mxGetString (prhs[2], regionstr, 256);
	sscanf (regionstr, "%s %s %lf", dimstr, relstr, &val);
	//cerr << relstr << endl;
	switch (toupper(dimstr[0])) {
	case 'X': dim = 0; break;
	case 'Y': dim = 1; break;
	case 'Z': dim = 2; break;
	case 'R': dim = 3; break;
	}
	for (i = 0; i <= 4; i++) {
	    if (!strcasecmp (relstr, Relstr[i])) {
		rel = i; break;
	    }
	}
    }

    if (dim >= 0)
	cerr << "toastBndmat: region limit: " << Dimstr[dim] << ' '
	     << Relstr[rel] << ' ' << val << endl;

    CalcBndmat (mesh, intstr, dim, rel, val, &plhs[0]);
}

// =========================================================================

void MatlabToast::Elmat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];
    mxArray *elmat;
    double *pr;

    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found.");

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
    } else {
	mexErrMsgTxt ("Elmat: Integral type string not recognised");
    }
    plhs[0] = elmat;
}

// =========================================================================

void MatlabToast::SampleField (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j, nsample, el;
    double *pr;

    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found.");

    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();

    // field coefficients
    RVector phi;
    CopyVector (phi, prhs[1]);
    ASSERTARG(phi.Dim() == nlen, 2, "Invalid array dimension");

    // list of sampling points
    RDenseMatrix pt;
    CopyMatrix (pt, prhs[2]);
    nsample = pt.nRows();
    ASSERTARG(pt.nCols() == dim, 3, "Invalid array dimension");

    // create the output array
    mxArray *sample = mxCreateDoubleMatrix (nsample, 1, mxREAL);
    pr = mxGetPr (sample);

    for (i = 0; i < nsample; i++) {
	Point p(dim);
	for (j = 0; j < dim; j++) p[j] = pt(i,j);
	el = mesh->ElFind (p);
	if (el < 0) { // didn't find an element
	    double dst, dstmin = 1e10;
	    for (j = 0; j < elen; j++) {
		dst = p.Dist (mesh->ElCentre(j));
		if (dst < dstmin) dstmin = dst, el = j;
	    }
	}
	Element *pel = mesh->elist[el];
	RVector fun = pel->GlobalShapeF (mesh->nlist, p);

	double val = 0.0;
	for (j = 0; j < fun.Dim(); j++)
	    val += fun[j] * phi[pel->Node[j]];
	*pr++ = val;
    }

    plhs[0] = sample;
}

// =========================================================================

void MatlabToast::BndReflectionTerm (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];
    double refind = mxGetScalar (prhs[0]);
    mxGetString (prhs[1], cbuf, 256);

    if (!strcasecmp (cbuf, "Keijzer")) {
	plhs[0] = mxCreateDoubleScalar(A_Keijzer(refind));
    } else if (!strcasecmp (cbuf, "Contini")) {
	plhs[0] = mxCreateDoubleScalar(A_Contini(refind));
    } else {
	plhs[0] = mxCreateDoubleScalar(1.0);
    }
}

// =========================================================================

bool ReadRVector (char *name, RVector &data, int idx)
{
    char c;

    ifstream ifs(name);

    if (!idx) { // read last vector
	RVector tmp;
	ifs >> tmp;
	while (ifs.good ()) {
	    data.Copy (tmp);
	    ifs >> tmp;
	}
    } else {
	for (int i = 0; i < idx; i++) {
	    ifs >> data;
	    if (!ifs.good()) {
		data.New(0);
		return false;
	    }
	}
    }
    return true;
}    

bool ReadCVector (char *name, CVector &data, int idx)
{
    char c;

    ifstream ifs(name);

    if (!idx) { // read last vector
	CVector tmp;
	ifs >> tmp;
	while (ifs.good ()) {
	    data.Copy (tmp);
	    ifs >> tmp;
	}
    } else {
	for (int i = 0; i < idx; i++) {
	    ifs >> data;
	    if (!ifs.good()) {
		data.New(0);
		return false;
	    }
	}
    }
    return true;
}    

void MatlabToast::ReadVector (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int idx = 1;
    char dname[256];

    ASSERTARG(mxIsChar(prhs[0]), 1, "File name expected.");
    mxGetString (prhs[0], dname, 256);

    if (nrhs > 1) {
	double v = mxGetScalar (prhs[1]);
	idx = (int)(v+0.5);
    }

    RVector rdata;
    CVector cdata;
    bool rok = false, cok = false;
    int dim;
    if (idx >= 0) {
	rok = ReadRVector (dname, rdata, idx);
	if (!rok || (dim = rdata.Dim()) == 0) {
	    cok = ReadCVector (dname, cdata, idx);
	    dim = cdata.Dim();
	}
	ASSERTARG((rok || cok)  && dim > 0, 1, "Error reading vector.");
    }

    if (cok) {
	plhs[0] = mxCreateDoubleMatrix (cdata.Dim(), 1, mxCOMPLEX);
	double *pr = mxGetPr (plhs[0]);
	double *pi = mxGetPi (plhs[0]);
	toast::complex *buf = cdata.data_buffer();
	for (int i = 0; i < cdata.Dim(); i++) {
	    pr[i] = buf[i].re;
	    pi[i] = buf[i].im;
	}
	    
    } else if (rok) {
	plhs[0] = mxCreateDoubleMatrix (rdata.Dim(), 1, mxREAL);
	double *pr  = mxGetPr (plhs[0]);
	double *buf = rdata.data_buffer();
	for (int i = 0; i < rdata.Dim(); i++)
	    pr[i] = buf[i];
    }
}

// =========================================================================

void MatlabToast::WriteVector (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char dname[256];
    ASSERTARG(mxGetString (prhs[0], dname, 256) == 0, 1, "String expected");
    ASSERTARG(mxIsNumeric (prhs[1]), 2, "Expected numeric array");

    ofstream ofs(dname);
    ofs.precision(12);
    ofs.setf (ios::scientific);

    if (mxIsComplex(prhs[1])) {
	CVector data;
	CopyVector (data, prhs[1]);
	ofs << data;
    } else {
	RVector data;
	CopyVector (data, prhs[1]);
	ofs << data;
    }
}

// ==========================================================================
// This data structure defines the Hessian implicitly

struct HESS_DATA {
    const RMatrix *J;                // Jacobian
    const RVector *M;                // normalisation diagonal matrix
    const double *lambda;            // diagonal scaling factor
    const Regularisation *reg;       // regularisation
    const RCompRowMatrix *RHess;     // Hessian of regularisation operator
};

// ==========================================================================

static RVector JTJx_clbk (const RVector &x, void *context)
{
    // Computes (J^T J + M P'' M + lambda I) x
    // where J is the column-normalised Jacobian,
    // P'' is the second derivative of the prior term,
    // M is the diagonal scaling matrix, and
    // lambda is a scalar

    // unpack the data
    HESS_DATA *data             = (HESS_DATA*)context;
    const RMatrix *J            =  data->J;
    const RCompRowMatrix *RHess =  data->RHess;
    const RVector &M            = *data->M;
    const double lambda         = *data->lambda;
    int m = J->nRows(), n = J->nCols();

    RVector Px(n);

    // add prior to Hessian
    if (RHess) {
	int i, j, k, nz, *colidx = new int[n];
	double *val = new double[n];
	for (i = 0; i < n; i++) {
	    nz = RHess->SparseRow (i, colidx, val);
	    for (k = 0; k < nz; k++) {
		j = colidx[k];
		Px[i] += val[k] * x[j] * M[i]*M[j];
	    }
	}
	delete []colidx;
	delete []val;
    }

    return ATx (*J, Ax (*J, x)) + Px + lambda*x;
}

// ==========================================================================
// Build the Hessian of the prior from individual parameter contributions

RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x,
    int nprm)
{
    int n = x.Dim();
    int n0 = n/nprm;
    int i, j;

    RCompRowMatrix H, Hi, Hij;
    for (i = 0; i < nprm; i++) {
	for (j = 0; j < nprm; j++) {
	    if (!j) {
		Hi.New(n0,n0);
		if (j==i) reg->SetHess1 (Hi, x, j);
	    } else {
		Hij.New(n0,n0);
		if (j==i) reg->SetHess1 (Hij, x, j);
		Hi = cath (Hi, Hij);
	    }
	}
	if (!i) H = Hi;
	else    H = catv (H, Hi);
    }
    return H;
}

void MatlabToast::Krylov (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int m = mxGetM (prhs[1]);
    int n = mxGetN (prhs[1]);
    int nprm = 2; // for now

    // copy current solution
    RVector x(n);
    CopyVector (x, prhs[0]);

    // copy Jacobian
    RDenseMatrix J(m,n);
    CopyMatrix (J, prhs[1]);

    // copy Gradient
    RVector g(n);
    CopyVector (g, prhs[2]);

    // copy M
    RVector M(n);
    CopyVector (M, prhs[3]);

    // copy lambda
    double lambda = mxGetScalar (prhs[4]);

    // regularisation handle
    Regularisation *reg = GetRegul(prhs[5]);
    ASSERTARG(reg, 6, "Regularisation instance not found");
    RCompRowMatrix *RHess = 0;
    if (reg)
	RHess = new RCompRowMatrix (BuildRHessian (reg, x, nprm));

    // copy tolerance
    double tol = mxGetScalar (prhs[6]);

    RVector z(n);
    
    HESS_DATA hdata = {&J, &M, &lambda, reg, RHess};
    static RPrecon_Diag precon;
    precon.ResetFromDiagonal (M);
    int iter;
    double res, time;
    tic();
    GMRES (JTJx_clbk, &hdata, g, z, tol, &precon, 30, &iter, &res);
    //GMRES (JTJx_clbk, &hdata, g, z, tol, (RPreconditioner*)0, 30, &iter, &res);
    time = toc();

    if (RHess) delete RHess;

    // copy result back to MATLAB
    CopyVector (&plhs[0], z);

    if (nlhs > 1) {
	static const char *fnames[3] = {"iter","res","time"};
	plhs[1] = mxCreateStructMatrix (1, 1, 3, fnames);
	mxSetField (plhs[1], 0, "iter", mxCreateDoubleScalar(iter));
	mxSetField (plhs[1], 0, "res", mxCreateDoubleScalar(res));
	mxSetField (plhs[1], 0, "time", mxCreateDoubleScalar(time));
    }
}

