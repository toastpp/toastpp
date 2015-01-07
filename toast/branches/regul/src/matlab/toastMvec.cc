// =========================================================================
// toastMvec
// Generate a matrix of measurement column vectors
//
// RH parameters:
//     1: mesh handle
//     2: measurement profile (string: 'Gaussian','Cosine','TrigBasis')
//     3: measurement radius [mm] (double)
// LH parameters:
//     1: boundary measurement matrix (sparse, measurements in column vectors)
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

int CalcMvec (const QMMesh *mesh, SRC_PROFILE mprof, double mwidth, const RVector *ref,
    mxArray **res) 
{
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
	if (ref) {
	    const RVector &rref = *ref;
	    for (j = 0; j < n; j++) m[j] *= c0/(2.0*rref[j]*A_Keijzer(rref[j]));
	}
	mvec.SetRow (i, m);
    }

    // return source vectors as matrix columns to matlab
    CopyTMatrix (res, mvec);

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

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char cbuf[256];

    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int n = mesh->nlen();

    SRC_PROFILE mprof;
    mxGetString (prhs[1], cbuf, 256);
    if      (!strcasecmp (cbuf, "Gaussian"))  mprof = PROF_GAUSSIAN;
    else if (!strcasecmp (cbuf, "Cosine"))    mprof = PROF_COSINE;
    else if (!strcasecmp (cbuf, "TrigBasis")) mprof = PROF_COMPLETETRIG;
    else mexErrMsgTxt ("Invalid measurement profile");

    double mwidth = mxGetScalar (prhs[2]);
    RVector ref(n);
    bool apply_c2a = true;

    if (nrhs >= 4) {
	int len = mxGetM(prhs[3])*mxGetN(prhs[3]);
	if ((len != 1 && len != n) || !mxIsDouble(prhs[3])) {
	    char cbuf[256];
	    sprintf (cbuf, "Mvec: parameter 4: expected double scalar or double vector of length %d", n);
	    mexErrMsgTxt (cbuf);
	}
	if (len == 1) {
	    double ref_homog = mxGetScalar(prhs[3]);
	    if (ref_homog) ref = mxGetScalar(prhs[3]);
	    else apply_c2a = false;
	} else
	    CopyVector (ref, prhs[3]);
    } else {
	mexWarnMsgTxt("Mvec: Using nodal refractive index values stored with mesh");
	ref = mesh->plist.N();
    }
    
    CalcMvec (mesh, mprof, mwidth, (apply_c2a ? &ref:0), &plhs[0]);
}
