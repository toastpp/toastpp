// ========================================================================
// Implementation of class MatlabToast
// Qvec/Mvec related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

#ifdef TOAST_THREAD_MATLAB_QMVEC
#include "task.h"
#endif

using namespace std;

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
    double phi0, phi1, a, f0, f1;
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

#ifdef TOAST_THREAD_MATLAB_QMVEC
struct Qvec_Threaddata {
    QMMesh *mesh;
    SourceMode qtype;
    SRC_PROFILE qprof;
    double qwidth;
    CCompRowMatrix *qvec;
};

void Qvec_engine (task_data *td)
{
    int i;
    int itask = td->proc;
    int ntask = td->np;
    Qvec_Threaddata *thdata = (Qvec_Threaddata*)td->data;
    QMMesh *mesh = thdata->mesh;
    int n = mesh->nlen();
    int nq = mesh->nQ;
    int q0 = (itask*nq)/ntask;
    int q1 = ((itask+1)*nq)/ntask;
    int dq = q1-q0;
    SourceMode qtype = thdata->qtype;
    SRC_PROFILE qprof = thdata->qprof;
    double qwidth = thdata->qwidth;
    CCompRowMatrix *qvec = thdata->qvec;

    CCompRowMatrix qvec_part(dq, n);

    for (i = q0; i < q1; i++) {
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
	qvec_part.SetRow (i-q0, q);
    }

    Task::UserMutex_lock();
    qvec->SetRows (q0, qvec_part);
    Task::UserMutex_unlock();
}
#endif

void MatlabToast::Qvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char typestr[256] = "";
    char profstr[256] = "";
    double w = 0.0;

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    int idx, n = mesh->nlen(), nQ = mesh->nQ;

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

#ifndef TOAST_THREAD_MATLAB_QMVEC
    for (int i = 0; i < nQ; i++) {
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
#else
    Qvec_Threaddata thdata = {
	mesh,
	qtype,
	qprof,
	qwidth,
	&qvec
    };
    Task::Multiprocess (Qvec_engine, &thdata);
#endif

    // return source vectors as matrix columns to matlab
    CopyTMatrix (&plhs[0], qvec);
}

// =========================================================================

#ifdef TOAST_THREAD_MATLAB_QMVEC
struct Mvec_Threaddata {
    QMMesh *mesh;
    SRC_PROFILE mprof;
    double mwidth;
    RVector *ref;
    CCompRowMatrix *mvec;
};

void Mvec_engine (task_data *td)
{
    const double c0 = 0.3;
    int i, j;
    int itask = td->proc;
    int ntask = td->np;
    Mvec_Threaddata *thdata = (Mvec_Threaddata*)td->data;
    QMMesh *mesh = thdata->mesh;
    RVector *ref = thdata->ref;
    int n = mesh->nlen();
    int nm = mesh->nM;
    int m0 = (itask*nm)/ntask;
    int m1 = ((itask+1)*nm)/ntask;
    int dm = m1-m0;
    SRC_PROFILE mprof = thdata->mprof;
    double mwidth = thdata->mwidth;
    CCompRowMatrix *mvec = thdata->mvec;

    CCompRowMatrix mvec_part(dm, n);

    for (i = m0; i < m1; i++) {
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
	if (ref) {
	    RVector &rref = *ref;
	    for (j = 0; j < n; j++) m[j] *= c0/(2.0*rref[j]*A_Keijzer(rref[j]));
	}
	mvec_part.SetRow (i-m0, m);
    }

    Task::UserMutex_lock();
    mvec->SetRows (m0, mvec_part);
    Task::UserMutex_unlock();
}
#endif

void MatlabToast::Mvec (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];

    QMMesh *mesh = (QMMesh*)GETMESH_SAFE(0);

    SRC_PROFILE mprof;
    mxGetString (prhs[1], cbuf, 256);
    if      (!strcasecmp (cbuf, "Gaussian"))  mprof = PROF_GAUSSIAN;
    else if (!strcasecmp (cbuf, "Cosine"))    mprof = PROF_COSINE;
    else if (!strcasecmp (cbuf, "TrigBasis")) mprof = PROF_COMPLETETRIG;
    else mexErrMsgTxt ("Invalid measurement profile");

    double mwidth = mxGetScalar (prhs[2]);

    int n, nM;

    n = mesh->nlen();
    nM = mesh->nM;
    CCompRowMatrix mvec;
    RVector ref(n);
    bool apply_c2a = true;

    if (nrhs >= 4) {
	int len = mxGetM(prhs[3])*mxGetN(prhs[3]);
	if ((len != 1 && len != n) || !mxIsDouble(prhs[3])) {
	    char cbuf[256];
	    sprintf (cbuf, "Mvec: parameter 3: expected double scalar or double vector of length %d", n);
	    mexErrMsgTxt (cbuf);
	}
	if (len == 1) {
	    double ref_homog = mxGetScalar(prhs[3]);
	    if (ref_homog) ref = mxGetScalar(prhs[3]);
	    else apply_c2a = false;
	} else
	    CopyVector (ref, prhs[3]);
    } else {
	mexErrMsgTxt("Mvec: no refractive index values supplied");
    }
    
    // build the measurement vectors
    mvec.New (nM, n);
#ifndef TOAST_THREAD_MATLAB_QMVEC
    const double c0 = 0.3;
    for (int i = 0; i < nM; i++) {
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
	if (apply_c2a)
	    for (int j = 0; j < n; j++) m[j] *= c0/(2.0*ref[j]*A_Keijzer(ref[j]));
	mvec.SetRow (i, m);
    }
#else
    Mvec_Threaddata thdata = {
	mesh,
	mprof,
	mwidth,
	(apply_c2a ? &ref : 0),
	&mvec
    };
    Task::Multiprocess (Mvec_engine, &thdata);
#endif

    // return source vectors as matrix columns to matlab
    CopyTMatrix (&plhs[0], mvec);
}
