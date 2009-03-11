// =========================================================================
// toastGradient
// gradient of objective function (Frechet derivative) at a given set
// of parameters. This includes only the data part, not any priors
//
// RH parameters:
//     1: mesh handle (double)
//     2: basis mapper handle (double)
//     3: source vectors (columns in complex sparse matrix)
//     4: measurement vectors (columns in complex sparse matrix)
//     5: mua (real vector)
//     6: mus (real vector)
//     7: refractive index (real vector)
//     8: modulation frequency [MHz]
//     9: data vector (real)
//    10: sd vector (real)
//    11: linear solver (string)
//    12: iterative solver tolerance (double) (optional)
// LH parameters:
//     1: gradient of objective function (real vector)
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

using namespace toast;

// =========================================================================

void AddDataGradient (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad)
{
    int i, j, q, m, n, idx, ofs_mod, ofs_arg;
    double term;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    CVector wqa (mesh->nlen());
    RVector wqb (mesh->nlen());
    CVector dgrad (slen);
    ofs_mod = 0;
    ofs_arg = mesh->nQM;
    RVector grad_cmua (grad, 0, slen);
    RVector grad_ckappa (grad, slen, slen);

    for (q = 0; q < mesh->nQ; q++) {

	// expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh->nQMref[q];

	RVector y_mod (data, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector ypm_mod (proj, ofs_mod, n);
	RVector b_mod(n);
	b_mod = (y_mod-ypm_mod)/s_mod;

	RVector y_arg (data, ofs_arg, n);
	RVector s_arg (sd, ofs_arg, n);
	RVector ypm_arg (proj, ofs_arg, n);
	RVector b_arg(n);
	b_arg = (y_arg-ypm_arg)/s_arg;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	//Project_cplx (*mesh, q, dphi[q], cproj);
	cproj = ProjectSingle (mesh, q, mvec, dphi[q]);
	wqa = complex(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].re;
	    double ip = cproj[idx].im;
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * complex (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = -2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]);
	    wqa += qs * complex (-term*ip*dn, -term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh->nlen());
	FWS.CalcField (wqa, wphia);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster->Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster->Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster->Map_GridToSol (gk, dgrad);
	grad_ckappa -= Re(dgrad);

	delete []cdfield_grad;
	delete []cafield_grad;

	ofs_mod += n; // step to next source
	ofs_arg += n;
    }
}


void GetGradient (QMMesh *mesh, Raster *raster, CFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref, double freq,
    const RVector &data, const RVector &sd,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    RVector &grad)
{
    const double c0 = 0.3;
    int i, n = mesh->nlen();
    int nQ = mesh->nQ;
    CVector *dphi;
    RVector proj(data.Dim());

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);

    // Calculate fields
    FWS.Allocate (*mesh);
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

    AddDataGradient (mesh, raster, FWS, proj, data, sd, dphi, mvec, grad);

    delete []dphi;
}

// ============================================================================
// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int n = mesh->nlen();

    Raster *raster = (Raster*)Handle2Ptr (mxGetScalar (prhs[1]));

    // source vectors
    CCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    CCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));

    // modulation frequency
    double freq = mxGetScalar (prhs[7]);
    
    // data & sd vectors
    RVector data (mesh->nQM*2, mxGetPr (prhs[8]));
    RVector sd (mesh->nQM*2, mxGetPr (prhs[9]));
    
    // linear solver parameters
    char solver[128];
    double tol = 1e-10;
    mxGetString (prhs[10], solver, 128);
    if (nrhs >= 12) tol = mxGetScalar (prhs[11]);
    CFwdSolver FWS (solver, tol);
    FWS.SetDataScaling (DATA_LOG);

    RVector grad(raster->SLen()*2);
    GetGradient (mesh, raster, FWS, mua, mus, ref, freq, data, sd,
        qvec, mvec, grad);
    CopyVector (&plhs[0], grad);
}
