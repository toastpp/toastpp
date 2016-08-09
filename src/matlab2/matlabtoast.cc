#include "matlabtoast.h"
#include "mexutil.h"
#include "util.h"
#include "timing.h"

using namespace std;

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

// -------------------------------------------------------------------------

void AssertArg_Char (const mxArray *arg, const char *func, int argno)
{
    AssertArg (mxIsChar(arg), func, argno, "Expected character array");
}

// -------------------------------------------------------------------------

void GetString_Safe (const mxArray *arg, char *pc, int len,
    const char *func, int argno)
{
    AssertArg_Char (arg, func, argno);
    mxGetString (arg, pc, len);
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
// Convert a pointer to a 64-bit handle on both 32 and 64 bit systems

uint64_T Ptr2Handle (void *ptr)
{
    switch (sizeof(ptr)) {
    case 4: { // 32-bit system
        union {
	    void *p[2];
	    uint64_T handle;
	} h;
	h.p[0] = ptr;
	h.p[1] = NULL;
	return h.handle;
        }
    case 8: { // 64-bit system
        union {
	    void *p;
	    uint64_T handle;
	} h;
	h.p = ptr;
	return h.handle;
        }
    default:
        xERROR ("Unsupported pointer size");
	return 0;
    }
}

// =========================================================================
// Convert a 64-bit handle to a pointer on both 32 and 64 bit systems

void *Handle2Ptr (uint64_T handle)
{
    switch (sizeof(void*)) {
    case 4: { // 32-bit system
        union {
	    void *p[2];
	    uint64_T handle;
	} h;
	h.handle = handle;
	return h.p[0];
        }
    case 8: { // 64-bit system
        union {
	    void *p;
	    uint64_T handle;
	} h;
	h.handle = handle;
	return h.p;
        }
    default:
        xERROR ("Unsupported pointer size");
	return NULL;
    }
}

// =========================================================================
// =========================================================================
// MatlabToast class implementation

MatlabToast::MatlabToast ()
{
    SetErrorhandler (ErrorHandler);
    nbasis = 0;
    verbosity = 0;
}

// =========================================================================

MatlabToast::~MatlabToast ()
{
    if (nbasis) delete []basislist;
}

// =========================================================================

void MatlabToast::ErrorHandler (char *msg)
{
    mexErrMsgTxt (msg);
}

// =========================================================================

Mesh *MatlabToast::GetMesh (const mxArray *idx, int *errid)
{
    if (errid) *errid = 0;

    if (!mxIsUint64(idx))
        mexErrMsgTxt("GetMesh: Invalid handle format (expected uint64)");

    uint64_T meshidx = *(uint64_T*)mxGetData(idx);
    Mesh *mesh = (Mesh*)Handle2Ptr(meshidx);
    if (!mesh)
        mexErrMsgTxt("GetMesh: Index out of range");

    return mesh;
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

Raster *MatlabToast::GetBasis (const mxArray *idx, int *errid, bool allownull)
{
    if (errid) *errid = 0;

    if (allownull) {
        if (mxIsNumeric(idx)) {
	    double v = mxGetScalar(idx);
	    if (!v) return NULL;
	}
    }

    if (!mxIsUint64(idx)) {
        if (errid) {
	    *errid = 1;
	    return 0;
	} else
	    mexErrMsgTxt("GetBasis: Invalid handle format (expected uint64)");
    }

    uint64_T basisidx = *(uint64_T*)mxGetData(idx);
    Raster *raster = (Raster*)Handle2Ptr(basisidx);
    if (!raster) {
        if (errid) {
	    *errid = 2;
	    return 0;
	} else
	    mexErrMsgTxt("GetBasis: Index out of range");
    }

    return raster;
}

// =========================================================================

Raster *MatlabToast::GetBasis_Safe (const mxArray *arr, const char *func,
    int argno)
{
    int err;
    Raster *rasterptr = GetBasis(arr, &err, true);
    switch(err) {
    case 1:
	AssertArg (0, func, argno, "Invalid handle format (expected uint64)");
	break;
    case 2:
	AssertArg (0, func, argno, "Index out of range");
	break;
    }
    return rasterptr;
}

// =========================================================================

Regularisation *MatlabToast::GetRegul (const mxArray *idx)
{
    uint64_T handle = *(uint64_T*)mxGetData(idx);
    Regularisation *reg = (Regularisation*)Handle2Ptr (handle);
    return reg;
    // put some more safeguards into this
}

// =========================================================================

void MatlabToast::SetVerbosity (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    const unsigned int max_verbosity = 1;

    if (mxIsUint32 (prhs[0])) {
	unsigned int verb = *(unsigned int*)mxGetData (prhs[0]);
	if (verb == verbosity) return; // nothing to do
	toastVerbosity = verbosity = std::min (verb, max_verbosity);
    } else {
	mexErrMsgTxt ("SetVerbosity: Argument 1: invalid format (uint32 expected).");
    }
}

// =========================================================================

void MatlabToast::ThreadCount (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int nthread_get = 1;
    int ncore = 1;
#ifdef TOAST_THREAD
    nthread_get = Task::GetThreadCount();
    ncore = Task::nProcessor();
    if (nrhs > 0) {
	int nthread_set = (int)(mxGetScalar(prhs[0])+0.5);
	Task::SetThreadCount (nthread_set);
    }
#endif
    if (nlhs > 0) {
	plhs[0] = mxCreateDoubleScalar(nthread_get);
	if (nlhs > 1) {
	    plhs[1] = mxCreateDoubleScalar(ncore);
	}
    }
}

// =========================================================================

void MatlabToast::MeshLin2Quad (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GetMesh(prhs[0]);
    ASSERTARG (mesh, 1, "Mesh not found");

    Mesh *quadmesh = Lin2Quad (*mesh);

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(quadmesh);

    if (verbosity >= 1)
	mexPrintf ("Mesh: %d nodes, %d elements, dimension %d\n",
		   quadmesh->nlen(), quadmesh->elen(), quadmesh->Dimension());
}

// =========================================================================

void MatlabToast::ReadNIM (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char nimname[256], meshname[256];
    int idx;
	bool read_all = false;

    if (mxIsChar (prhs[0]))
	mxGetString (prhs[0], nimname, 256);
    else
	mexErrMsgTxt ("ReadNIM: Argument 1: file name expected.");

    if (fileExists (nimname) != 1)
	mexErrMsgTxt ("ReadNIM: Image file not found.");

	if (nrhs < 2) {
		idx = 1;
	} else if (mxIsChar(prhs[1])) {
		char cbuf[256];
		mxGetString (prhs[1], cbuf, 256);
		if (!strcasecmp (cbuf, "all"))
			read_all = true;
		else
			mexErrMsgTxt ("Argument 2 must be index or 'all'");
	} else {
		idx = (int)mxGetScalar (prhs[1]);
	}

	if (read_all) {
		RDenseMatrix img;
		ReadNimAll (nimname, img);
		CopyMatrix (&plhs[0], img);

		if (verbosity >= 1) {
			mexPrintf ("Nodal image:\n");
			mexPrintf ("--> Size............%d x %d\n", img.nRows(), img.nCols());
		}
	} else {
	    RVector img;
		ReadNim (nimname, idx, img, meshname);

		plhs[0] = mxCreateDoubleMatrix (img.Dim(), 1, mxREAL);
		memcpy (mxGetPr (plhs[0]), img.data_buffer(), img.Dim()*sizeof(double));

		if (verbosity >= 1) {
			mexPrintf ("Nodal image:\n");
			mexPrintf ("--> Size............%d\n", img.Dim());
		}
    }

    if (nlhs > 1) {
        plhs[1] = mxCreateString(meshname);
    }
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

// =========================================================================

void CalcSysmat (QMMesh *mesh, RVector &mua, RVector &mus, RVector &ref,
		 double freq, bool elbasis, mxArray **res)
{
    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int nz, n = (elbasis ? elen : nlen);
    idxtype *rowptr, *colidx;

    mesh->SparseRowStructure (rowptr, colidx, nz);

    // Set optical coefficients
    RVector cmua = mua*c0/ref;
    RVector ckappa = c0/(3.0*ref*(mua+mus));
    RVector c2a(n);
    for (int i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));

    if (freq) { // complex-valued system matrix
	double omega = freq * 2.0*Pi*1e-6;
	CCompRowMatrix F(nlen, nlen, rowptr, colidx);
	AddToSysMatrix (*mesh, F, &cmua,
			elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
	AddToSysMatrix (*mesh, F, &ckappa,
			elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
	AddToSysMatrix (*mesh, F, &c2a,
			elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
	AddToSysMatrix (*mesh, F, omega, ASSEMBLE_iCFF);
	CopyMatrix (res, F);
    } else { // real-valued problem
	RCompRowMatrix F(nlen, nlen, rowptr, colidx);
	AddToSysMatrix (*mesh, F, &cmua,
			elbasis ? ASSEMBLE_PFF_EL:ASSEMBLE_PFF);
	AddToSysMatrix (*mesh, F, &ckappa,
			elbasis ? ASSEMBLE_PDD_EL:ASSEMBLE_PDD);
	AddToSysMatrix (*mesh, F, &c2a,
			elbasis ? ASSEMBLE_BNDPFF_EL:ASSEMBLE_BNDPFF);
	CopyMatrix (res, F);
    }

    delete []rowptr;
    delete []colidx;
}

void CalcBndSysmat (QMMesh *mesh, RVector &ref, mxArray **res)
{
    int nz, nlen = mesh->nlen();
    idxtype *rowptr, *colidx;

    mesh->SparseRowStructure (rowptr, colidx, nz);
    RCompRowMatrix BF(nlen, nlen, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    RVector prm(nlen);
    for (int i = 0; i < nlen; i++)
	prm[i] = c0/ref[i];

    AddToSysMatrix (*mesh, BF, &prm, ASSEMBLE_BNDPFF);
    CopyMatrix (res, BF);
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
  QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

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
    int i;
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
    double refind = mxGetScalar (prhs[0]);
    bool match = false;
    if (nrhs > 1 && mxIsChar(prhs[1])) {
	char cbuf[256];
        mxGetString (prhs[1], cbuf, 256);
	if (!strcasecmp (cbuf, "Keijzer")) {
	    plhs[0] = mxCreateDoubleScalar(A_Keijzer(refind));
	    match = true;
	} else if (!strcasecmp (cbuf, "Contini")) {
	    plhs[0] = mxCreateDoubleScalar(A_Contini(refind));
	    match = true;
	}
    }
    if (!match)
	plhs[0] = mxCreateDoubleScalar(1.0);
}

// =========================================================================

bool ReadRVector (char *name, RVector &data, int idx)
{
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
	std::complex<double> *buf = cdata.data_buffer();
	for (int i = 0; i < cdata.Dim(); i++) {
	    pr[i] = buf[i].real();
	    pi[i] = buf[i].imag();
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
    int n = J->nCols();

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
    mwSize m = mxGetM (prhs[1]);
    mwSize n = mxGetN (prhs[1]);

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
    RCompRowMatrix *RHess = 0;
    if (reg) {
        int nprm = reg->GetNParam();
	RHess = new RCompRowMatrix (BuildRHessian (reg, x, nprm));
    }

    // copy tolerance
    double tol = mxGetScalar (prhs[6]);

    RVector z(n);
    
    HESS_DATA hdata = {&J, &M, &lambda, reg, RHess};
    static RPrecon_Diag precon;
    precon.ResetFromDiagonal (M);
    int iter;
    double res, time;
    tic();
    GMRES (JTJx_clbk, &hdata, g, z, tol, &precon, 30, 0, &iter, &res);
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
