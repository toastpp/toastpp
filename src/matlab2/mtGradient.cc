// ========================================================================
// Implementation of class MatlabToast
// Gradient-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

#ifdef TOAST_THREAD_MATLAB_GRADIENT
#include "task.h"
#endif

using namespace std;

// =========================================================================
// Prototypes
// =========================================================================

void AddDataGradientReal (QMMesh *mesh, Raster *raster, const RFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, RVector *dphi,
    const RCompRowMatrix &mvec, RVector &grad);

void AddDataGradientCplx (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad);

void AddDataGradient2 (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad);

void GetGradientReal (QMMesh *mesh, Raster *raster, RFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref,
    const RVector &data, const RVector &sd,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    RVector &grad, RVector *phi=0, RVector *proj=0);

void GetGradientCplx (QMMesh *mesh, Raster *raster, CFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref, double freq,
    const RVector &data, const RVector &sd,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    RVector &grad, CVector *phi=0, RVector *proj=0);

// =========================================================================

void MatlabToast::Gradient (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // See if we need to do a real or complex computation
    double freq = mxGetScalar (prhs[7]);

    if (!freq) GradientReal (nlhs, plhs, nrhs, prhs);
    else       GradientCplx (nlhs, plhs, nrhs, prhs);
}

// =========================================================================

void MatlabToast::GradientReal (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n = mesh->nlen();

    Raster *raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    RVector grad(raster->SLen()*2);

    // source vectors
    RCompRowMatrix qvec;
    CopyTMatrix (qvec, prhs[2]);

    // measurement vectors
    RCompRowMatrix mvec;
    CopyTMatrix (mvec, prhs[3]);

    // nodal optical parameters
    RVector mua (n, mxGetPr (prhs[4]));
    RVector mus (n, mxGetPr (prhs[5]));
    RVector ref (n, mxGetPr (prhs[6]));

    // data & sd vectors (truncate phase parts if present)
    RVector data (mesh->nQM, mxGetPr (prhs[8]));
    RVector sd (mesh->nQM, mxGetPr (prhs[9]));
    
    char solver[128] = "direct";
    double tol = 1e-10;
    RVector *phi = 0;
    RVector *proj = 0;

    // read optional parameters as tupels
    for (i = 10; i < nrhs; i++) {
        char label[128];
        mxGetString (prhs[i], label, 128);

	if (!strcasecmp(label, "Method")) {     // linear solver method
	    mxGetString (prhs[++i], solver, 128);
	} else if (!strcasecmp(label, "Tolerance")) { // linear solver tolerance
  	    tol = mxGetScalar (prhs[++i]);
	} else if (!strcasecmp(label, "Fields")) { // fields for all sources
	    i++;
	    int q, nq = mesh->nQ, nlen = mesh->nlen();
	    mwSize m = mxGetM(prhs[i]);
	    mwSize n = mxGetN(prhs[i]);
	    double *pr = mxGetPr(prhs[i]);
	    xASSERT(m == nlen && n == nq, "Parameter phi wrong dimension");
	    phi = new RVector[nq];
	    for (q = 0; q < nq; q++) {
	        phi[q].New (nlen);
		for (j = 0; j < nlen; j++)
		    phi[q][j] = *pr++;
	    }
	} else if (!strcasecmp(label,"Projections")) {
	    proj = new RVector(mesh->nQM, mxGetPr (prhs[++i]));
	} else if (!strcasecmp(label,"Unwrap")) {
	    // N/A
	} else {
	    mexErrMsgTxt("Error parsing arguments");
	}
    }

    RFwdSolver FWS (mesh, solver, tol);
    FWS.SetDataScaling (DATA_LOG);

    GetGradientReal (mesh, raster, FWS, mua, mus, ref, data, sd,
        qvec, mvec, grad, phi, proj);
    CopyVector (&plhs[0], grad);

    if (phi) delete []phi;
    if (proj) delete proj;
}

// =========================================================================

void MatlabToast::GradientCplx (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    // mesh
    QMMesh *mesh = (QMMesh*)GetMesh(prhs[0]);
    ASSERTARG(mesh, 1, "Mesh not found");

    int n = mesh->nlen();

    Raster *raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    RVector grad(raster->SLen()*2);

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
    
    char solver[128] = "direct";
    double tol = 1e-10;
    bool unwrap = false;
    CVector *phi = 0;
    RVector *proj = 0;

    // read optional parameters as tupels
    for (i = 10; i < nrhs; i++) {
        char label[128];
        mxGetString (prhs[i], label, 128);

	if (!strcasecmp(label, "Method")) {     // linear solver method
	    mxGetString (prhs[++i], solver, 128);
	} else if (!strcasecmp(label, "Tolerance")) { // linear solver tolerance
  	    tol = mxGetScalar (prhs[++i]);
	} else if (!strcasecmp(label, "Fields")) { // fields for all sources
	    i++;
	    int q, nq = mesh->nQ, nlen = mesh->nlen();
	    mwSize m = mxGetM(prhs[i]);
	    mwSize n = mxGetN(prhs[i]);
	    double *pr = mxGetPr(prhs[i]);
	    double *pi = mxGetPi(prhs[i]);
	    xASSERT(m == nlen && n == nq, "Parameter phi wrong dimension");
	    phi = new CVector[nq];
	    for (q = 0; q < nq; q++) {
	        phi[q].New (nlen);
		for (j = 0; j < nlen; j++)
		    phi[q][j] = std::complex<double> (*pr++, *pi++);
	    }
	} else if (!strcasecmp(label,"Projections")) {
	    proj = new RVector(mesh->nQM*2, mxGetPr (prhs[++i]));
	} else if (!strcasecmp(label,"Unwrap")) {
	    unwrap = mxIsLogicalScalarTrue (prhs[++i]);	    
	} else {
	    mexErrMsgTxt("Error parsing arguments");
	}
    }

	int nth;
#ifdef TOAST_THREAD_MATLAB_GRADIENT
	nth = Task::GetThreadCount();
#else
	nth = 1;
#endif

    CFwdSolver FWS (mesh, solver, tol, nth);
    FWS.SetDataScaling (DATA_LOG);
    FWS.SetPhaseUnwrap (unwrap);

    GetGradientCplx (mesh, raster, FWS, mua, mus, ref, freq, data, sd,
		 qvec, mvec, grad, phi, proj);
    CopyVector (&plhs[0], grad);

    if (phi) delete []phi;
    if (proj) delete proj;
}

// ==========================================================================
// Implementation
// ==========================================================================

#ifndef TOAST_THREAD_MATLAB_GRADIENT

// ==========================================================================

void AddDataGradientReal (QMMesh *mesh, Raster *raster, const RFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, RVector *dphi,
    const RCompRowMatrix &mvec, RVector &grad)
{
    int i, j, q, m, n, idx, ofs_mod;
    double term;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    RVector wqa (mesh->nlen());
    //RVector wqb (mesh->nlen());
    RVector dgrad (slen);
    ofs_mod = 0;
    RVector grad_cmua (grad, 0, slen);
    RVector grad_ckappa (grad, slen, slen);

    for (q = 0; q < mesh->nQ; q++) {

	// expand field and gradient
        RVector cdfield (glen);
        RVector *cdfield_grad = new RVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh->nQMref[q];

	RVector y_mod (data, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector ypm_mod (proj, ofs_mod, n);
	RVector b_mod(n);
	b_mod = (y_mod-ypm_mod)/s_mod;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	RVector cproj(n);
	cproj = ProjectSingle (mesh, q, mvec, dphi[q]);
	wqa = 0.0;
	//wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const RVector qs = mvec.Row(m);
	    double rp = cproj[idx];
	    double dn = 1.0/(rp*rp);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * term*rp*dn;

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	RVector wphia (mesh->nlen());
	FWS.CalcField (wqa, wphia);

	RVector cafield(glen);
	RVector *cafield_grad = new RVector[dim];
	raster->Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster->Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua -= dgrad;

	// diffusion contribution
	// multiply complex field gradients
	RVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster->Map_GridToSol (gk, dgrad);
	grad_ckappa -= dgrad;

	delete []cdfield_grad;
	delete []cafield_grad;

	ofs_mod += n; // step to next source
    }
}

// ==========================================================================

void AddDataGradientCplx (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
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
    //RVector wqb (mesh->nlen());
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
	cproj = ProjectSingle (mesh, q, mvec, dphi[q]);
	wqa = std::complex<double>(0,0);
	//wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].real();
	    double ip = cproj[idx].imag();
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * std::complex<double> (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = -2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]);
	    wqa += qs * std::complex<double> (-term*ip*dn, -term*rp*dn);

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

	// DEBUG
	//if (q == 6)
	//  return;

	delete []cdfield_grad;
	delete []cafield_grad;

	ofs_mod += n; // step to next source
	ofs_arg += n;
    }
}

// ==========================================================================

#else // multithreaded version

// ==========================================================================

struct AddDataGradientReal_Threaddata {
    const QMMesh *mesh;
    const Raster *raster;
    const RFwdSolver *fws;
    const RCompRowMatrix *mvec;
    const RVector *proj;
    const RVector *data;
    const RVector *sd;
    RVector *dphi;
    RVector *grad;
};

void AddDataGradientReal_engine (task_data *td)
{
    int i, j, q, m, n, idx, ofs_mod;
    int itask = td->proc;
    int ntask = td->np;
    AddDataGradientReal_Threaddata *thdata =
        (AddDataGradientReal_Threaddata*)td->data;
    const QMMesh *mesh = thdata->mesh;
    const Raster *raster = thdata->raster;
    const RFwdSolver *fws = thdata->fws;
    const RCompRowMatrix *mvec = thdata->mvec;
    const RVector &proj = *thdata->proj;
    const RVector &data = *thdata->data;
    const RVector &sd = *thdata->sd;
    RVector *dphi = thdata->dphi;
    int nq = mesh->nQ;
    int q0 = (itask*nq)/ntask;
    int q1 = ((itask+1)*nq)/ntask;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    double term;
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    RVector *grad = thdata->grad;
    RVector grad_cmua_loc (slen);
    RVector grad_ckap_loc (slen);
    RVector wqa (mesh->nlen());
    RVector wqb (mesh->nlen());
    RVector dgrad (slen);
    ofs_mod = 0;          // data offset for Mod data
    for (i = 0; i < q0; i++) {
	ofs_mod += mesh->nQMref[i];
    }
    
    for (q = q0; q < q1; q++) {

        // expand field and gradient
        RVector cdfield (glen);
        RVector *cdfield_grad = new RVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh->nQMref[q];

	RVector b_mod(n);
	RVector s_mod(n);
	for (i = 0; i < n; i++) {
	    s_mod[i] = sd[ofs_mod+i];
 	    b_mod[i] = (data[ofs_mod+i]-proj[ofs_mod+i])/s_mod[i];
	}

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	RVector cproj(n);

	cproj = fws->ProjectSingle (q, *mvec, dphi[q], DATA_LIN);

	wqa = 0.0;
	wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const RVector qs = mvec->Row(m);
	    double rp = cproj[idx];
	    double dn = 1.0/(rp*rp);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * term*rp*dn;

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	RVector wphia (mesh->nlen());
	fws->CalcField (wqa, wphia);

	RVector cafield(glen);
	RVector *cafield_grad = new RVector[dim];
	raster->Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster->Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua_loc -= dgrad;

	// diffusion contribution
	// multiply complex field gradients
	RVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster->Map_GridToSol (gk, dgrad);
	grad_ckap_loc -= dgrad;

	ofs_mod += n; // step to next source

	delete []cdfield_grad;
	delete []cafield_grad;
    }

    // assemble into global gradient vector
    Task::UserMutex_lock();
    {
        RVector grad_cmua(*grad, 0, slen);    // mua part of grad
	RVector grad_ckap(*grad, slen, slen); // kappa part of grad
	grad_cmua += grad_cmua_loc;
	grad_ckap += grad_ckap_loc;
    }
    Task::UserMutex_unlock();
}

void AddDataGradientReal (QMMesh *mesh, Raster *raster, const RFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, RVector *dphi,
    const RCompRowMatrix &mvec, RVector &grad)
{
    AddDataGradientReal_Threaddata thdata = {
        mesh,
        raster,
        &FWS,
        &mvec,
        &proj,
        &data,
        &sd,
        dphi,
        &grad
    };

    Task::Multiprocess (AddDataGradientReal_engine, &thdata);
}

// ==========================================================================

struct AddDataGradientCplx_Threaddata {
    const QMMesh *mesh;
    const Raster *raster;
    const CFwdSolver *fws;
    const CCompRowMatrix *mvec;
    const RVector *proj;
    const RVector *data;
    const RVector *sd;
    CVector *dphi;
    RVector *grad;
};

void AddDataGradientCplx_engine (task_data *td)
{
    int i, j, q, m, n, idx, ofs_mod, ofs_arg;
    int itask = td->proc;
    int ntask = td->np;
    AddDataGradientCplx_Threaddata *thdata =
        (AddDataGradientCplx_Threaddata*)td->data;
    const QMMesh *mesh = thdata->mesh;
    const Raster *raster = thdata->raster;
    const CFwdSolver *fws = thdata->fws;
    const CCompRowMatrix *mvec = thdata->mvec;
    const RVector &proj = *thdata->proj;
    const RVector &data = *thdata->data;
    const RVector &sd = *thdata->sd;
    CVector *dphi = thdata->dphi;
    int nq = mesh->nQ;
    int q0 = (itask*nq)/ntask;
    int q1 = ((itask+1)*nq)/ntask;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    double term;
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    RVector *grad = thdata->grad;
    RVector grad_cmua_loc (slen);
    RVector grad_ckap_loc (slen);
    CVector wqa (mesh->nlen());
    RVector wqb (mesh->nlen());
    CVector dgrad (slen);
    ofs_mod = 0;          // data offset for Mod data
    ofs_arg = mesh->nQM;  // data offset for Arg data
    for (i = 0; i < q0; i++) {
	ofs_mod += mesh->nQMref[i];
	ofs_arg += mesh->nQMref[i];
    }
    
    for (q = q0; q < q1; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh->nQMref[q];

	RVector b_mod(n);
	RVector b_arg(n);
	RVector s_mod(n);
	RVector s_arg(n);
	for (i = 0; i < n; i++) {
	    s_mod[i] = sd[ofs_mod+i];
	    s_arg[i] = sd[ofs_arg+i];
 	    b_mod[i] = (data[ofs_mod+i]-proj[ofs_mod+i])/s_mod[i];
	    b_arg[i] = (data[ofs_arg+i]-proj[ofs_arg+i])/s_arg[i];
	}

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);

	cproj = fws->ProjectSingle (q, *mvec, dphi[q], DATA_LIN);

	wqa = std::complex<double>(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const CVector qs = mvec->Row(m);
	    double rp = cproj[idx].real();
	    double ip = cproj[idx].imag();
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * std::complex<double> (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = -2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]);
	    wqa += qs * std::complex<double> (-term*ip*dn, -term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh->nlen());
	fws->CalcField (wqa, wphia, 0, itask);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster->Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster->Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua_loc -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster->Map_GridToSol (gk, dgrad);
	grad_ckap_loc -= Re(dgrad);

	ofs_mod += n; // step to next source
	ofs_arg += n;

	delete []cdfield_grad;
	delete []cafield_grad;
    }

    // assemble into global gradient vector
    Task::UserMutex_lock();
    {
        RVector grad_cmua(*grad, 0, slen);    // mua part of grad
	RVector grad_ckap(*grad, slen, slen); // kappa part of grad
	grad_cmua += grad_cmua_loc;
	grad_ckap += grad_ckap_loc;
    }
    Task::UserMutex_unlock();
}

void AddDataGradientCplx (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad)
{
    AddDataGradientCplx_Threaddata thdata = {
        mesh,
        raster,
        &FWS,
        &mvec,
        &proj,
        &data,
        &sd,
        dphi,
        &grad
    };

    Task::Multiprocess (AddDataGradientCplx_engine, &thdata);
}

#endif

// =========================================================================
// Note: This does not work properly

void AddDataGradient2 (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad)
{
    int i, j, q, m, idx, ofs_mod, ofs_arg;
    double term;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    CCompRowMatrix wqa (mesh->nQ, mesh->nlen());
    CVector wqa_row(mesh->nlen());
    RVector wqb (mesh->nlen());
    CVector dgrad (slen);
    ofs_mod = 0;
    ofs_arg = mesh->nQM;
    RVector grad_cmua (grad, 0, slen);
    RVector grad_ckappa (grad, slen, slen);

    for (q = 0; q < mesh->nQ; q++) {
        int n = mesh->nQMref[q];

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
	wqa_row = std::complex<double>(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].real();
	    double ip = cproj[idx].imag();
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa_row += qs * std::complex<double> (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = -2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]);
	    wqa_row += qs * std::complex<double> (-term*ip*dn, -term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}
	wqa.SetRow (q, wqa_row);
    }

    CVector *wphia = new CVector[mesh->nQ];
    for (q = 0; q < mesh->nQ; q++) wphia[q].New (mesh->nlen());
    FWS.CalcFields (wqa, wphia);

    for (q = 0; q < mesh->nQ; q++) {
        int n = mesh->nQMref[q];

	// expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster->Map_MeshToGrid (wphia[q], cafield);
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

// =========================================================================

void GetGradientReal (QMMesh *mesh, Raster *raster, RFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref,
    const RVector &data, const RVector &sd,
    const RCompRowMatrix &qvec, const RCompRowMatrix &mvec,
    RVector &grad, RVector *phi, RVector *proj)
{
    const double c0 = 0.3;
    int i, n = mesh->nlen();
    int nQ = mesh->nQ;
    RVector *phi_local = 0;
    RVector proj_local(data.Dim());

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

    // Calculate fields
    FWS.Allocate ();
    FWS.Reset (msol);

    if (!phi || !proj) {
        if (!phi) {
	    phi_local = new RVector[nQ];
	    for (i = 0; i < nQ; i++)
	        phi_local[i].New (n);
	    FWS.CalcFields (qvec, phi_local);
	    phi = phi_local;
	}
	if (!proj) {
	    proj_local = FWS.ProjectAll_real (mvec, phi);
	    proj = &proj_local;
	}
    } else if (phi && FWS.LinSolver() == LSOLVER_DIRECT) {
	// HACK: SuperLU appears to fail occasionally unless it starts with
	// solving for a 'regular' source. Maybe a problem with condition
	// of linear system when RHS is an adjoint source?
	RVector tmp(n);
	FWS.CalcField(qvec.Row(0), tmp);
    }

    AddDataGradientReal (mesh, raster, FWS, *proj, data, sd, phi, mvec, grad);

    if (phi_local) delete []phi_local;
}

// =========================================================================

void GetGradientCplx (QMMesh *mesh, Raster *raster, CFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref, double freq,
    const RVector &data, const RVector &sd,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    RVector &grad, CVector *phi, RVector *proj)
{
    const double c0 = 0.3;
    int i, n = mesh->nlen();
    int nQ = mesh->nQ;
    CVector *phi_local = 0;
    RVector proj_local(data.Dim());

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

    // Calculate fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);

    if (!phi || !proj) {
        if (!phi) {
			phi_local = new CVector[nQ];
			for (i = 0; i < nQ; i++)
				phi_local[i].New (n);
			FWS.CalcFields (qvec, phi_local);
			phi = phi_local;
		}
		if (!proj) {
		    proj_local = FWS.ProjectAll_real (mvec, phi);
			proj = &proj_local;
		}
    } else if (phi && FWS.LinSolver() == LSOLVER_DIRECT) {
		// HACK: SuperLU appears to fail occasionally unless it starts with
		// solving for a 'regular' source. Maybe a problem with condition
		// of linear system when RHS is an adjoint source?
		CVector tmp(n);
		FWS.CalcField(qvec.Row(0), tmp);
    }

    AddDataGradientCplx (mesh, raster, FWS, *proj, data, sd, phi, mvec, grad);

    if (phi_local) delete []phi_local;
}

