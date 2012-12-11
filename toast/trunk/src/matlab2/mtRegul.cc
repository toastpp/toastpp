// ========================================================================
// Implementation of class MatlabToast
// Regularisation-related methods
// ========================================================================

#include "matlabtoast.h"
#include "toastmex.h"

using namespace std;
using namespace toast;

// =========================================================================
// Prototypes
// =========================================================================

void ExtractKappa (const mxArray *prm, void **kappa, bool *istensor);

// =========================================================================
// Matlab interface
// =========================================================================

// =========================================================================

void MatlabToast::Regul (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Regularisation *reg = 0;    // regularisation instance
    const mxArray *regp = 0;    // parameter structure
    mxArray *field, *subfield;  // struct field pointers
    int prm = 0;
    char cbuf[256];

    // some global default settings
    char rtype[256] = "None";   // regularisation method
    double tau = 1e-4;          // regularisation hyperparameter
    void *kapref = 0;           // reference diffusivity
    bool istensor = false;      // reference diffusivity in tensor format?
    RVector kaprefimg;          // image for reference diffusivity
    bool brefimg = false;       // reference diffusivity from image?
    double sdr = 0.0;           // reference std
    double fT = 0.0;            // reference PM threshold

    // raster
    Raster *raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    // initial image
    RVector x0;
    CopyVector (x0, prhs[2]);

    plhs[0] = mxCreateNumericMatrix (1,1, mxUINT32_CLASS, mxREAL);
    unsigned int *ptr = (unsigned int*)mxGetData (plhs[0]);

    if (mxIsStruct (prhs[0])) { // scan parameters from struct

	regp = prhs[0]; // parameter structure

	// regularisation type
	field = mxGetField (regp,0,"method");
	if (field)
	    mxGetString (field, rtype, 256);
	else {
	    cerr << "Warning: toastRegul: missing required field: method\n";
	    cerr << "Disabling regularisation" << endl;
	}
	if (strcasecmp (rtype,"None")) {
	    field = mxGetField (regp,0,"tau");
	    if (field)
		tau = mxGetScalar (field);
	    else {
		cerr << "Warning: toastRegul: missing required field: tau\n";
		cerr << "Substituting default value: " << tau << endl;
	    }
	    field = mxGetField(regp,0,"prior");
	    if (field) {
		if (subfield = mxGetField (field, 0, "refimg")) {
		    CopyVector (kaprefimg, subfield);
		    brefimg = true;
		} else if (subfield = mxGetField (field,0,"refname")) {
		    mxGetString (subfield, cbuf, 256);
		    if (cbuf[0]) {
			ifstream ifs(cbuf);
			ifs >> kaprefimg;
			brefimg = true;
		    }
		}
		sdr = mxGetScalar(mxGetField(field,0,"smooth"));
		fT = mxGetScalar(mxGetField(field,0,"threshold"));
	    }
	    
	}   

    } else { // scan parameters from argument list

	// regularisation type
	prm = 3;
	mxGetString (prhs[0], rtype, 256);

	// hyperparameter tau
	tau = mxGetScalar (prhs[prm++]);

	// read generic parameters from argument list
	int prm0 = prm;
	while (prm0 < nrhs) {
	    mxGetString (prhs[prm0++], cbuf, 256);
	    if (!strcasecmp (cbuf, "KapRefImage")) {
		CopyVector (kaprefimg, prhs[prm0++]);
		brefimg = true;
	    } else if (!strcasecmp (cbuf, "KapRefScale")) {
		sdr = mxGetScalar (prhs[prm0++]);
	    } else if (!strcasecmp (cbuf, "KapRefPMThreshold")) {
		fT = mxGetScalar (prhs[prm0++]);
	    } else if (!strcasecmp (cbuf, "KapRefTensor")) {
		istensor = (mxGetScalar (prhs[prm0++]) != 0.0);
	    } else if (!strcasecmp (cbuf, "KapRef")) {
		ExtractKappa (prhs[prm0++], &kapref, &istensor);
		brefimg = false;
	    }
	}

    }
	
    if (!strcasecmp (rtype, "DIAG")) {
	    
	reg = new Tikhonov0 (tau, &x0, &x0);
	    
    } else if (!strcasecmp (rtype, "LAPLACIAN")) {
	    
	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("Regularisation: LAPLACIAN (obsolete)\n");
	    mexPrintf ("--> tau: %f\n", tau);
	}

	reg = new Tikhonov1 (tau, &x0, raster);
	    
    } else if (!strcasecmp (rtype, "TV")) {
	// total variation (TV)
	    
	// default settings
	double beta = 0.01;

	if (regp) { // read TV parameters from structure
	    field = mxGetField (regp,0,"tv");
	    if (field) field = mxGetField (field,0,"beta");
	    if (field) beta = mxGetScalar(field);
	} else { // read TV parameters from argument list
	    while (prm < nrhs) {
		mxGetString (prhs[prm++], cbuf, 256);
		if (!strcasecmp (cbuf, "Beta")) {
		    beta = mxGetScalar (prhs[prm++]);
		    break;
		}
	    }
	}
	    
	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("Regularisation: total variation (TV)\n");
	    mexPrintf ("--> tau: %f\n", tau);
	    mexPrintf ("--> beta: %f\n", beta);
	    mexPrintf ("--> diffusivity %s\n", (brefimg ?
		"from image" : kapref ? "from external array" : "none"));
	    if (brefimg || kapref) {
		mexPrintf ("--> diffusivity field format: %s\n",
			   (istensor ? "tensor" : "scalar"));
		mexPrintf ("--> diffusivity scale: %f\n", sdr);
		mexPrintf ("--> diffusivity PM threshold: %f\n", fT);
	    }
	}
	
	if (brefimg) {
	    reg = new TV (tau, beta, &x0, raster, kaprefimg, sdr, fT,
			  istensor);
	} else {
	    reg = new TV (tau, beta, &x0, raster, kapref, istensor);
	}
	
    } else if (!strcasecmp (rtype, "TK0")) {
	// 0-th order Tikhonov

	// default settings
	RVector xs (x0.Dim());
	xs = 1;
	if (regp) { // read TK0 parameters from structure
	    field = mxGetField(regp,0,"tk0");
	    if (field) field = mxGetField(field,0,"xs");
	    if (field) CopyVector (xs, field);
	} else { // read TK0 parameters from argument list
	    while (prm < nrhs) {
		mxGetString (prhs[prm++], cbuf, 256);
		if (!strcasecmp (cbuf, "Xs")) {
		    CopyVector (xs, prhs[prm++]);
		    break;
		}
	    }
	}

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("Regularisation: 0th order Tikhonov (TK0)\n");
	    mexPrintf ("--> tau: %f\n", tau);
	}

	reg = new Tikhonov0 (tau, &x0, &xs);

    } else if (!strcasecmp (rtype, "TK1")) {
	// 1-st order Tikhonov (Laplacian)

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("Regularisation: 1st order Tikhonov (TK1)\n");
	    mexPrintf ("--> tau: %f\n", tau);
	    mexPrintf ("--> diffusivity %s\n", (brefimg ?
		 "from image" : kapref ? "from external array" : "none"));
	    if (brefimg || kapref) {
		mexPrintf ("--> diffusivity field format: %s\n",
			   (istensor ? "tensor" : "scalar"));
		mexPrintf ("--> diffusivity scale: %f\n", sdr);
		mexPrintf ("--> diffusivity PM threshold: %f\n", fT);
	    }
	}

	if (brefimg) {
	    reg = new TK1 (tau, &x0, raster, kaprefimg, sdr, fT, istensor);
	} else {
	    reg = new TK1 (tau, &x0, raster, kapref, istensor);
	}
	
    } else if (!strcasecmp (rtype, "HUBER")) {
	// Huber function
	
	// default settings
	double eps = 0.01;
	
	if (regp) { // read Huber parameters from structure
	    field = mxGetField (regp,0,"huber");
	    if (field) field = mxGetField(field,0,"eps");
	    if (field) eps = mxGetScalar(field);
	} else { // read Huber parameters from argument list
	    while (prm < nrhs) {
		mxGetString (prhs[prm++], cbuf, 256);
		if (!strcasecmp (cbuf, "Eps")) {
		    eps = mxGetScalar (prhs[prm++]);
		    break;
		}
	    }
	}

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("Regularisation: Huber\n");
	    mexPrintf ("--> tau: %f\n", tau);
	    mexPrintf ("--> eps: %f\n", eps);
	    mexPrintf ("--> diffusivity %s\n", (brefimg ?
		 "from image" : kapref ? "from external array" : "none"));
	    if (brefimg || kapref) {
		mexPrintf ("--> diffusivity field format: %s\n",
			   (istensor ? "tensor" : "scalar"));
		mexPrintf ("--> diffusivity scale: %f\n", sdr);
		mexPrintf ("--> diffusivity PM threshold: %f\n", fT);
	    }
	}

	if (brefimg) {
	    reg = new Huber (tau, eps, &x0, raster, kaprefimg, sdr, fT,
			     istensor);
	} else {
	    reg = new Huber (tau, eps, &x0, raster, kapref, istensor);
	}
	
    } else if (!strcasecmp (rtype, "PM")) {
	// Perona-Malik (PM) regularisation
	
	// default settings
	double T = 1.0;
	
	// read specific parameters from argument list
	while (prm < nrhs) {
	    mxGetString (prhs[prm++], cbuf, 256);
	    if (!strcasecmp (cbuf, "T")) {
		T = mxGetScalar (prhs[prm++]);
	    }
	}

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    cerr << "Regularisation: Perona-Malik (PM)" << endl;
	    cerr << "--> tau=" << tau << endl;
	    cerr << "--> T=" << T << endl;
	    cerr << "--> diffusivity "
		 << (brefimg ? "from image" :
		 kapref ? "from external array" : "none") << endl;
	    if (brefimg || kapref) {
		cerr << "--> diffusivity field format: "
		     << (istensor ? "tensor" : "scalar") << endl;
		cerr << "--> diffusivity scale: " << sdr << endl;
		cerr << "--> diffusivity PM threshold: " << fT << endl;
	    }
	}

	if (brefimg) {
	    reg = new PM (tau, T, &x0, raster, kaprefimg, sdr, fT, istensor);
	} else {
	    reg = new PM (tau, T, &x0, raster, kapref, istensor);
	}
	
    } else if (!strcasecmp (rtype, "QPM")) {
	// Quadratic (?) Perona-Malik (QPM) regularisation
	
	// default settings
	double T = 1.0;
	
	// read specific parameters from argument list
	while (prm < nrhs) {
	    mxGetString (prhs[prm++], cbuf, 256);
	    if (!strcasecmp (cbuf, "T")) {
		T = mxGetScalar (prhs[prm++]);
	    }
	}

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    cerr << "Regularisation: Quadratic Perona-Malik (QPM)" << endl;
	    cerr << "--> tau=" << tau << endl;
	    cerr << "--> T=" << T << endl;
	    cerr << "--> diffusivity "
	     << (brefimg ? "from image" :
		 kapref ? "from external array" : "none") << endl;
	    if (brefimg || kapref) {
		cerr << "--> diffusivity field format: "
		     << (istensor ? "tensor" : "scalar") << endl;
		cerr << "--> diffusivity scale: " << sdr << endl;
		cerr << "--> diffusivity PM threshold: " << fT << endl;
	    }
	}

	if (brefimg) {
	    reg = new QPM (tau, T, &x0, raster, kaprefimg, sdr, fT, istensor);
	} else {
	    reg = new QPM (tau, T, &x0, raster, kapref, istensor);
	}
	
    } else if (!strcasecmp (rtype, "Tukey")) {
	// Tukey regularisation
	
	// default settings
	double T = 1.0;
	
	// read specific parameters from argument list
	while (prm < nrhs) {
	    mxGetString (prhs[prm++], cbuf, 256);
	    if (!strcasecmp (cbuf, "T")) {
		T = mxGetScalar (prhs[prm++]);
	    }
	}
	
	if (verbosity >= 1) {
	    // echo regularisation parameters
	    cerr << "Regularisation: Tukey" << endl;
	    cerr << "--> tau=" << tau << endl;
	    cerr << "--> T=" << T << endl;
	    cerr << "--> diffusivity "
		 << (brefimg ? "from image" :
		     kapref ? "from external array" : "none") << endl;
	    if (brefimg || kapref) {
		cerr << "--> diffusivity field format: "
		     << (istensor ? "tensor" : "scalar") << endl;
		cerr << "--> diffusivity scale: " << sdr << endl;
		cerr << "--> diffusivity PM threshold: " << fT << endl;
	    }
	}

	if (brefimg) {
	    reg = new Tukey (tau, T, &x0, raster, kaprefimg, sdr, fT,
			     istensor);
	} else {
	    reg = new Tukey (tau, T, &x0, raster, kapref, istensor);
	}
	
    } else {
	
	*ptr = 0;  // "no regularisation"
	return;
	
    }
    
    Regularisation **tmp = new Regularisation*[nreg+1];
    if (nreg) {
	memcpy (tmp, reglist, nreg*sizeof(Regularisation*));
	delete []reglist;
    }
    reglist = tmp;
    reglist[nreg++] = reg;

    *ptr = nreg;
}

// =========================================================================

void MatlabToast::ClearRegul (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    unsigned int regid;

    if (mxIsUint32(prhs[0])) {
	regid = *(unsigned int*)mxGetData(prhs[0]);
	if (regid == 0) // "no regularisation" flag
	    return;     // nothing to do
	regid--;
	if (regid >= nreg) {
	    if (nreg)
		mexPrintf ("ClearRegul: regularisation index out of range (expected 1..%d).\n", nreg);
	    else
		mexPrintf ("ClearRegul: regularisation index out of range (no meshes registered).\n");
	    return;
	}
    } else {
	mexPrintf ("ClearRegul: Invalid regularisation index format (expected uint32).\n");
	return;
    }

    if (reglist[regid]) {
	delete reglist[regid];
	reglist[regid] = 0;
    } else {
	mexPrintf ("ClearRegul: regularisation already cleared.\n");
    }
    
}

// =========================================================================

void MatlabToast::RegulValue (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    double prior = (reg ? reg->GetValue (x) : 0.0);
    plhs[0] = mxCreateDoubleScalar (prior);
}

// =========================================================================

void MatlabToast::RegulGradient (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    RVector grad(x.Dim());
    if (reg) grad = reg->GetGradient (x);

    CopyVector (&plhs[0], grad);    
}

// =========================================================================

void MatlabToast::RegulHDiag (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    RVector diag(x.Dim());
    if (reg) diag = reg->GetHessianDiag (x);
    CopyVector (&plhs[0], diag);
}

// =========================================================================

RCompRowMatrix BuildRHessian (Regularisation *reg, const RVector &x)
{
    int nprm = reg->GetNParam();
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

void MatlabToast::RegulHess (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);
    RCompRowMatrix H (BuildRHessian (reg, x));
    CopyMatrix (&plhs[0], H);    
}

// =========================================================================

void MatlabToast::RegulHess1f (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    // operand vector
    RVector f;
    CopyVector (f, prhs[2]);

    RVector y = reg->GetHess1f (x, f);
    CopyVector (&plhs[0], y);
}

// =========================================================================

void MatlabToast::RegulKappa (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // solution vector
    RVector x;
    CopyVector (x, prhs[1]);

    RVector kappa = reg->GetKappa (x);
    CopyVector (&plhs[0], kappa);
}

// =========================================================================
// Local functions
// =========================================================================

void ExtractKappa (const mxArray *prm, void **kappa, bool *istensor)
{
    mwSize ndim = mxGetNumberOfDimensions (prm);

    if (ndim <= 2) { // scalar field -> isotropic diffusivity
	cerr << "kref: scalar field" << endl;
	RVector *k = new RVector;
	CopyVector (*k, prm);
	*kappa = (void*)k;
	*istensor = false;
    } else { // tensor field -> anisotropic diffusivity
	cerr << "kref: tensor field" << endl;
	const mwSize *dim = mxGetDimensions (prm);
	RDenseMatrix *k = new RDenseMatrix[dim[0]];
	for (mwSize i = 0; i < dim[0]; i++) {
	    mxArray *ti = mxGetCell (prm, i);
	    CopyMatrix (k[i], ti);
	    *kappa = (void*)k;
	    *istensor = true;
	}
    }
}

