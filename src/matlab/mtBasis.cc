// ========================================================================
// Implementation of class MatlabToast
// Basis-related methods
// ========================================================================

#include "matlabtoast.h"

using namespace std;
using namespace toast;

// =========================================================================
// Matlab interface
// =========================================================================

// =========================================================================

void MatlabToast::SetBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char basistype[256] = "LINEAR";
    int i, j, arg = 0, basistp;
    RDenseMatrix *bb = 0;

    // basis type
    ASSERTARG (nrhs > arg, 0, "Too few parameters provided");

    if (mxIsChar (prhs[arg])) {
	mxGetString (prhs[arg], basistype, 256);
	arg++;
    }
    if (!strcasecmp (basistype, "LINEAR")) {
	basistp = 0;
    } else if (!strcasecmp (basistype, "CUBIC")) {
	basistp = 1;
    } else {
	basistp = -1;
    }
    ASSERTARG(basistp >= 0, 1, "unexpected value");

    // mesh

    ASSERTARG(nrhs > arg, 0, "Too few parameters provided");
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[arg++]);
    ASSERTARG(mesh, 2, "Mesh not found");

    // solution basis dimensions
    ASSERTARG(nrhs > arg, 0, "Too few parameters provided");
    int dim = mxGetNumberOfElements (prhs[arg]);
    ASSERTARG(dim == mesh->Dimension(), 3, "Invalid basis dimensions");
    IVector bdim(dim);
    for (i = 0; i < dim; i++)
	bdim[i] = (int)mxGetPr (prhs[arg])[i];
    arg++;
    IVector gdim(bdim);

    if (nrhs > arg) {
	// sampling grid dimensions
	if (mxGetNumberOfElements (prhs[arg]) == dim) {
	    // argument is intermediate grid dimension
	    for (i = 0; i < dim; i++)
		gdim[i] = (int)mxGetPr (prhs[arg])[i];
	    arg++;
	}
    }
    if (nrhs > arg) {
	if (mxGetM (prhs[arg]) == dim && mxGetN (prhs[arg]) == 2) {
	    // argument is grid bounding box
	    bb = new RDenseMatrix (2,dim);
	    CopyTMatrix (*bb, prhs[arg]);
	    arg++;
	}
    }

    Raster *raster;
    switch (basistp) {
    case 0:
	raster = new Raster_Pixel (bdim, gdim, mesh, bb);
	break;
    case 1:
	raster = new Raster_CubicPixel (bdim, gdim, mesh, bb);
	break;
    }
	
    Raster **tmp = new Raster*[nbasis+1];
    if (nbasis) {
	memcpy (tmp, basislist, nbasis*sizeof(Raster*));
	delete []basislist;
    }
    basislist = tmp;
    basislist[nbasis++] = raster;

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT32_CLASS, mxREAL);
    unsigned int *ptr = (unsigned int*)mxGetData (plhs[0]);
    *ptr = nbasis;

    if (verbosity >= 1) {
	char cbuf[256];
	mexPrintf ("basis: Type:       %s-%s\n",
		   mesh->Dimension()==2 ? "BI":"TRI", basistype);
	sprintf (cbuf, "basis: Grid size:  %d [", raster->GLen());
	for (i = 0; i < dim; i++)
	    sprintf (cbuf+strlen(cbuf), "%d%c", gdim[i], i==dim-1 ? ']':'x');
	mexPrintf ("%s\n", cbuf);
	
	mexPrintf ("basis: Sol. size:  %d\n", raster->SLen());
	mexPrintf ("basis: Mesh size:  %d\n", mesh->nlen());
	if (bb) {
	    mexPrintf ("basis: Grid bounding box:\n");
	    for (i = 0; i < 2; i++) {
		for (j = 0; j < dim; j++)
		    mexPrintf ("  %12g", bb->Get(i,j));
		mexPrintf ("\n");
	    }
	}
    }

    if (bb) delete bb;
}

// =========================================================================

void MatlabToast::ClearBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    unsigned int basisid;

    if (mxIsUint32(prhs[0])) {
	basisid = *(unsigned int*)mxGetData(prhs[0]) - 1;
	if (basisid >= nbasis) {
	    if (nbasis)
		mexPrintf ("ClearBasis: basis index out of range (expected 1..%d).\n", nbasis);
	    else
		mexPrintf ("ClearBasis: basis index out of range (no basis instances registered).\n");
	    return;
	}
    } else {
	mexPrintf ("ClearBasis: Invalid basis index format (expected uint32).\n");
	return;
    }

    if (basislist[basisid]) {
	delete basislist[basisid];
	basislist[basisid] = 0;
    } else {
	mexPrintf ("ClearBasis: basis already cleared.\n");
    }
}

// =========================================================================

void MatlabToast::GetBasisSize (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int i, dim = raster->Dim();
    if (nlhs >= 1) {
	RVector tmp(dim);
	for (i = 0; i < dim; i++) tmp[i] = raster->BDim()[i];
	CopyVector (&plhs[0], tmp);
    }
    if (nlhs >= 2) {
	RVector tmp(dim);
	for (i = 0; i < dim; i++) tmp[i] = raster->GDim()[i];
	CopyVector (&plhs[1], tmp);
    }
}

// =========================================================================

void MatlabToast::MapBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // basis handle
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    // source and target basis representations
    ASSERTARG(mxIsChar (prhs[1]), 2, "Expected string.");

    char basisstr[256], srcid, tgtid;
    mxGetString (prhs[1], basisstr, 256);
    ASSERTARG(strlen(basisstr) == 4, 2, "Format not recognised.");
    if (!strncmp (basisstr+1, "->", 2)) {
	srcid = toupper(basisstr[0]);
	tgtid = toupper(basisstr[3]);
    } else if (!strncmp (basisstr+1, "<-", 2)) {
	srcid = toupper(basisstr[3]);
	tgtid = toupper(basisstr[0]);
    } else {
	ASSERTARG(0, 2, "Format not recognised.");
    }
    
    // source field
    mwIndex m = mxGetM (prhs[2]);
    mwIndex n = mxGetN (prhs[2]);
    int nsrc = (int)(m*n);

    if (mxIsComplex (prhs[2])) {

	double *pr = mxGetPr(prhs[2]);
	double *pi = mxGetPi(prhs[2]);
	CVector src(nsrc), tgt;
	toast::complex *vptr = src.data_buffer();
	for (int i = 0; i < nsrc; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	switch (srcid) {
	case 'M':
	    ASSERTARG(raster->mesh().nlen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_MeshToGrid (src, tgt);
		break;
	    case 'B':
		raster->Map_MeshToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_MeshToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'G':
	    ASSERTARG(raster->GLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'M':
		raster->Map_GridToMesh (src, tgt);
		break;
	    case 'B':
		raster->Map_GridToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_GridToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'B':
	    ASSERTARG(raster->BLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_BasisToGrid (src, tgt);
		break;
	    case 'S':
		raster->Map_BasisToSol (src, tgt);
		break;
	    case 'M':
		raster->Map_BasisToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'S':
	    ASSERTARG(raster->SLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'B':
		raster->Map_SolToBasis (src, tgt);
		break;
	    case 'G':
		raster->Map_SolToGrid (src, tgt);
		break;
	    case 'M':
		raster->Map_SolToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	default:
	    ASSERTARG(0, 2, "Source id not recognised.");
	    return;
	}

	CopyVector (&plhs[0], tgt);

    } else {

	RVector src (nsrc, mxGetPr (prhs[2]));
	RVector tgt;

	switch (srcid) {
	case 'M':
	    ASSERTARG(raster->mesh().nlen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_MeshToGrid (src, tgt);
		break;
	    case 'B':
		raster->Map_MeshToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_MeshToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'G':
	    ASSERTARG(raster->GLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'M':
		raster->Map_GridToMesh (src, tgt);
		break;
	    case 'B':
		raster->Map_GridToBasis (src, tgt);
		break;
	    case 'S':
		raster->Map_GridToSol (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'B':
	    ASSERTARG(raster->BLen() == nsrc, 3,
	        "Source vector unexpected length.");
	    switch (tgtid) {
	    case 'G':
		raster->Map_BasisToGrid (src, tgt);
		break;
	    case 'S':
		raster->Map_BasisToSol (src, tgt);
		break;
	    case 'M':
		raster->Map_BasisToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	case 'S':
	    ASSERTARG(raster->SLen() == nsrc, 3,
		"Source vector unexpected length.");
	    switch (tgtid) {
	    case 'B':
		raster->Map_SolToBasis (src, tgt);
		break;
	    case 'G':
		raster->Map_SolToGrid (src, tgt);
		break;
	    case 'M':
		raster->Map_SolToMesh (src, tgt);
		break;
	    default:
		ASSERTARG(0, 2, "Target id not recognised.");
		return;
	    }
	    break;
	default:
	    ASSERTARG(0, 2, "Source id not recognised.");
	    return;
	}

	CopyVector (&plhs[0], tgt);
    }    
}

// =========================================================================

void MatlabToast::MapMeshToBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int nn = raster->mesh().nlen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == nn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector nim (nn);
	CVector img (bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = nim.data_buffer();
	for (int i = 0; i < nn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_MeshToBasis (nim, img);
	CopyVector (&plhs[0], img);

    } else {

	RVector nim (nn, mxGetPr (prhs[1]));
	RVector img (bn);
	raster->Map_MeshToBasis (nim, img);
	CopyVector (&plhs[0], img);

    }
}

// =========================================================================

void MatlabToast::MapMeshToGrid (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int nn = raster->mesh().nlen();
    int gn = raster->GLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == nn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector nim (nn);
	CVector img (gn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = nim.data_buffer();
	for (int i = 0; i < nn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_MeshToGrid (nim, img);
	CopyVector (&plhs[0], img);

    } else {

	RVector nim (nn, mxGetPr (prhs[1]));
	RVector img (gn);
	raster->Map_MeshToGrid (nim, img);
	CopyVector (&plhs[0], img);

    }    
}

// =========================================================================

void MatlabToast::MapMeshToSol (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int nn = raster->mesh().nlen();
    int sn = raster->SLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == nn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector nim (nn);
	CVector img (sn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = nim.data_buffer();
	for (int i = 0; i < nn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_MeshToSol (nim, img);
	CopyVector (&plhs[0], img);

    } else {

	RVector nim (nn, mxGetPr (prhs[1]));
	RVector img (sn);
	raster->Map_MeshToSol (nim, img);
	CopyVector (&plhs[0], img);

    }
}

// =========================================================================

void MatlabToast::MapBasisToMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int nn = raster->mesh().nlen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == bn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (bn);
	CVector nim (nn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < bn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_BasisToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    } else {

	RVector img (bn, mxGetPr(prhs[1]));
	RVector nim (nn);
	raster->Map_BasisToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    }
}

// =========================================================================

void MatlabToast::MapSolToMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int nn = raster->mesh().nlen();
    int sn = raster->SLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == sn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (sn);
	CVector nim (nn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < sn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_SolToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    } else {

	RVector img (sn, mxGetPr(prhs[1]));
	RVector nim (nn);
	raster->Map_SolToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    }
}

// =========================================================================

void MatlabToast::MapSolToBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int sn = raster->SLen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == sn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (sn);
	CVector bimg(bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < sn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_SolToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    } else {

	RVector img (sn, mxGetPr(prhs[1]));
	RVector bimg(bn);
	raster->Map_SolToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    }
}

// =========================================================================

void MatlabToast::MapSolToGrid (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int sn = raster->SLen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == sn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (sn);
	CVector bimg(bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < sn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_SolToGrid (img, bimg);
	CopyVector (&plhs[0], bimg);

    } else {

	RVector img (sn, mxGetPr(prhs[1]));
	RVector bimg(bn);
	raster->Map_SolToGrid (img, bimg);
	CopyVector (&plhs[0], bimg);

    }
}

// =========================================================================

void MatlabToast::MapGridToMesh (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int nn = raster->mesh().nlen();
    int gn = raster->GLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == gn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (gn);
	CVector nim (nn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < gn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_GridToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    } else {

	RVector img (gn, mxGetPr(prhs[1]));
	RVector nim (nn);
	raster->Map_BasisToMesh (img, nim);
	CopyVector (&plhs[0], nim);

    }
}

// =========================================================================

void MatlabToast::MapGridToBasis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int gn = raster->GLen();
    int bn = raster->BLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == gn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (gn);
	CVector bimg(bn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < gn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_GridToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    } else {

	RVector img (gn, mxGetPr(prhs[1]));
	RVector bimg(bn);
	raster->Map_GridToBasis (img, bimg);
	CopyVector (&plhs[0], bimg);

    }
}

// =========================================================================

void MatlabToast::MapGridToSol (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int gn = raster->GLen();
    int sn = raster->SLen();

    // nodal image
    mwIndex m = mxGetM (prhs[1]);
    mwIndex n = mxGetN (prhs[1]);
    ASSERTARG((int)(m*n) == gn, 2, "Invalid image dimensions");

    if (mxIsComplex (prhs[1])) {

	CVector img (gn);
	CVector simg(sn);
	double *pr = mxGetPr(prhs[1]);
	double *pi = mxGetPi(prhs[1]);
	toast::complex *vptr = img.data_buffer();
	for (int i = 0; i < gn; i++) {
	    vptr[i].re = pr[i];
	    vptr[i].im = pi[i];
	}
	raster->Map_GridToSol (img, simg);
	CopyVector (&plhs[0], simg);

    } else {

	RVector img (gn, mxGetPr(prhs[1]));
	RVector simg(sn);
	raster->Map_GridToSol (img, simg);
	CopyVector (&plhs[0], simg);

    }
}

// =========================================================================

void MatlabToast::BasisToMeshMatrix (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    // copy mapping matrix
    RCompRowMatrix &map = 
	(RCompRowMatrix&)raster->Basis2MeshMatrix(); // dodgy cast

    CopyMatrix (&plhs[0], map);    
}

// =========================================================================

void MatlabToast::MeshToBasisMatrix (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    // copy mapping matrix
    RCompRowMatrix &map = 
	(RCompRowMatrix&)raster->Mesh2BasisMatrix(); // dodgy cast

    CopyMatrix (&plhs[0], map);
}

// =========================================================================

void MatlabToast::SolutionMask (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");
    int slen = raster->SLen();

    plhs[0] = mxCreateDoubleMatrix (1, slen, mxREAL);
    double *pr = mxGetPr (plhs[0]);

    for (int i = 0; i < slen; i++)
	pr[i] = raster->Sol2Basis (i) + 1;
        // "+1" to adjust to matlab's 1-based array indices
}

// =========================================================================

void MatlabToast::GridElref (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    int i, n = raster->GLen();
    int *elref = raster->Elref();

    plhs[0] = mxCreateDoubleMatrix (n,1,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (i = 0; i < n; i++) {
      pr[i] = elref[i] + 1;
        // make 1-based. Gridpoints outside mesh support are set to zero
    }    
}

// =========================================================================

#ifdef UNDEF
template<class T>
void ImageGradient (const IVector &gdim, const RVector &gsize,
    const TVector<T> &im, TVector<T> *grad, const int *mask)
{
    // this should be done by Fourier transform
    int x, y, z, idx, dim = gdim.Dim();
    int nx = gdim[0], ny = gdim[1], nz = (dim >= 3 ? gdim[2]:1);
    int n = nx*ny*nz;

    // x gradient
    double dx = gsize[0]/nx, ix = 1.0/dx, i2x = 0.5*ix;
    TVector<T> &gradx = grad[0];
    gradx.New (n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt   = (!mask || (mask[idx] >= 0));
		bool bleft  = (x > 0 && (!mask || (mask[idx-1] >= 0)));
		bool bright = (x < nx-1 && (!mask || (mask[idx+1] >= 0)));
		if (bleft && bright) {
		    gradx[idx] = (im[idx+1]-im[idx-1]) * i2x;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    gradx[idx] = 0.0;
		} else if (bleft) {
		    gradx[idx] = (im[idx]-im[idx-1]) * ix;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (bright) {
		    gradx[idx] = (im[idx+1]-im[idx]) * ix;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else {
		    gradx[idx] = 0.0;
		}
		idx++;
	    }
	}
    }

    // y gradient
    double dy = gsize[1]/ny, iy = 1.0/dy, i2y = 0.5*iy;
    TVector<T> &grady = grad[1];
    grady.New(n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt  = (!mask || (mask[idx] >= 0));
		bool bup   = (y > 0 && (!mask || (mask[idx-nx] >= 0)));
		bool bdown = (y < ny-1 && (!mask || (mask[idx+nx] >= 0)));
		if (bup && bdown) {
		    grady[idx] = (im[idx+nx]-im[idx-nx]) * i2y;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    grady[idx] = 0.0;
		} else if (bup) {
		    grady[idx] = (im[idx]-im[idx-nx]) * iy;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (bdown) {
		    grady[idx] = (im[idx+nx]-im[idx]) * iy;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else {
		    grady[idx] = 0.0;
		}
		idx++;
	    }
	}
    }
    if (dim < 3) return;

    // z gradient
    double dz = gsize[2]/nz, iz = 1.0/dz, i2z = 0.5*iz;
    int stridez = nx*ny;
    TVector<T> &gradz = grad[2];
    gradz.New(n);
    for (z = idx = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
	    for (x = 0; x < nx; x++) {
	        bool bcnt = (!mask || (mask[idx] >= 0));
		bool bfront = (z > 0 && !(mask || (mask[idx-stridez] >= 0)));
		bool bback  = (z < nz-1 && (!mask || (mask[idx+stridez] >= 0)));
	        if (bfront && bback) {
		    gradz[idx] = (im[idx+stridez]-im[idx-stridez]) * i2z;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (!bcnt) {
		    gradz[idx] = 0.0;
		} else if (bfront) {
		    gradz[idx] = (im[idx]-im[idx-stridez]) * iz;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else if (bback) {
		    gradz[idx] = (im[idx+stridez]-im[idx]) * iz;
		    INC_FLOPS_ADD(1);
		    INC_FLOPS_MUL(1);
		} else {
		    gradz[idx] = 0.0;
		}
		idx++;
	    }
	}
    }
}
#endif

void MatlabToast::ImageGradient (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // raster
    Raster *raster = GetBasis(prhs[0]);
    ASSERTARG(raster, 1, "Basis not found");

    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    int glen = raster->GLen();
    int d, dim = gdim.Dim();

    if (mxIsComplex (prhs[1])) {
	CVector img;
	CVector *grad = new CVector[dim];
	for (d = 0; d < dim; d++) grad[d].New (glen);
	CopyVector (img, prhs[1]);
	::ImageGradient (gdim, gsize, img, grad, elref);
	CDenseMatrix g (dim, glen);
	for (d = 0; d < dim; d++) g.SetRow (d, grad[d]);
	delete []grad;
	CopyMatrix (&plhs[0], g);
    } else {
	RVector img;
	RVector *grad = new RVector[dim];
	for (d = 0; d < dim; d++) grad[d].New (glen);
	CopyVector (img, prhs[1]);
	::ImageGradient (gdim, gsize, img, grad, elref);
	RDenseMatrix g (dim, glen);
	for (d = 0; d < dim; d++) g.SetRow (d, grad[d]);
	delete []grad;
	CopyMatrix (&plhs[0], g);
    }
}
