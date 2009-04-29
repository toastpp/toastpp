// =========================================================================
// toastQvec
// Generate a matrix of source column vectors
//
// RH parameters:
//     1: mesh handle
//     2: source type (string: 'Neumann','Isotropic')
//     3: source profile (string: 'Gaussian','Cosine','TrigBasis')
//     4: source radius [mm] (double)
// LH parameters:
//     1: source matrix (sparse, sources in column vectors)
// =========================================================================

#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "source.h"
#include "util.h"

// =========================================================================
// local prototypes

CVector CompleteTrigSourceVector (const Mesh &mesh, int order);

// =========================================================================
// MAIN 

int CalcQvec (const QMMesh *mesh, SourceMode qtype, SRC_PROFILE qprof,
    double qwidth, const RVector &refind, mxArray **res) 
{
    int i, j, n, nQ;

    n = mesh->nlen();
    nQ = mesh->nQ;
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
	for (j = 0; j < n; j++) {
	    const double c0 = 0.3;
	    double c = c0/refind[j];
	    double A = A_Keijzer(refind[j]);
	    q[j] *= c/(2.0*A);
	}
	qvec.SetRow (i, q);
    }

    // return source vectors as matrix columns to matlab
    CopyTMatrix (res, qvec);

    return 0;
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

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int idx, n = mesh->nlen();

    char typestr[256] = "";
    char profstr[256] = "";
    double w = 0.0;
    RVector refind = mesh->plist.N(); // use mesh parameter by default
    
    if (nrhs >= 1 && mxIsStruct (prhs[1])) {

	// read source parameters from structure
	mxArray *field;
	field = mxGetField (prhs[1], 0, "type");
	if (field) mxGetString (field, typestr, 256);
	field = mxGetField (prhs[1], 0, "prof");
	if (field) mxGetString (field, profstr, 256);
	field = mxGetField (prhs[1], 0, "width");
	if (field) w = mxGetScalar (field);
	field = mxGetField (prhs[1], 0, "refind");
	if (field) CopyVector (refind, field);

    } else if (nrhs >= 2) {

	// read source parameters from function parameters
	if (mxIsChar(prhs[1])) mxGetString (prhs[1], typestr, 256);
	if (mxIsChar(prhs[2])) mxGetString (prhs[2], profstr, 256);
	if (nrhs >= 3) {
	    idx = 3;
	    if (mxIsNumeric(prhs[idx]) && mxGetNumberOfElements(prhs[idx])==1){
		w = mxGetScalar (prhs[idx]);
		idx++;
	    }
	    if (nrhs >= idx && mxIsNumeric(prhs[idx]) &&
		mxGetNumberOfElements(prhs[idx]) > 1) {
		CopyVector (refind, prhs[idx]);
		idx++;
	    }
	}
		
    } else {

	mexErrMsgTxt ("toastQvec: Invalid function parameters");

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

    CalcQvec (mesh, qtype, qprof, qwidth, refind, &plhs[0]);
}
