// ========================================================================
// Implementation of class MatlabToast
// Regularisation-related methods
// ========================================================================

#include "matlabtoast.h"
#include "mexutil.h"

using namespace std;

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
    RVector x0;
    Raster *raster = 0;

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
	    field = mxGetField(regp,0,"x0");
	    if (field)
	        CopyVector (x0, field);
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
        GETSTRING_SAFE(0, rtype, 256);

	// hyperparameter tau
	tau = mxGetScalar (prhs[3]);

	// generic parameters from key/value list
	if (nrhs > 3 && mxIsCell(prhs[4])) {
		int nprm = mxGetM(prhs[4])*mxGetN(prhs[4]);

	    for (prm = 0; prm < nprm; prm+=2) {
			field = mxGetCell (prhs[4], prm);
	        AssertArg_Char (field, __func__, 5);
		mxGetString (field, cbuf, 256);
		if (!strcasecmp (cbuf, "KapRefImage")) {
		    CopyVector (kaprefimg, mxGetCell (prhs[4], prm+1));
		    brefimg = true;
		} else if (!strcasecmp (cbuf, "KapRefScale")) {
		    sdr = mxGetScalar (mxGetCell (prhs[4], prm+1));
		} else if (!strcasecmp (cbuf, "KapRefPMThreshold")) {
		    fT = mxGetScalar (mxGetCell (prhs[4], prm+1));
		} else if (!strcasecmp (cbuf, "KapRefTensor")) {
		    istensor = (mxGetScalar (mxGetCell (prhs[4], prm+1))
				!= 0.0);
		} else if (!strcasecmp (cbuf, "KapRef")) {
		    ExtractKappa (mxGetCell (prhs[4], prm+1), &kapref,
				  &istensor);
		    brefimg = false;
		    // WARNING: race condition: relies on KapRefTensor
		    // being read before KapRef
		}
	    }
	}
    }

    // raster
    raster = GetBasis(prhs[1]);
    ASSERTARG(raster, 2, "Basis not found");

    // initial image
    CopyVector (x0, prhs[2]);

    mexPrintf ("Regularisation:\n");

    if (!strcasecmp (rtype, "DIAG")) {
	    
	reg = new Tikhonov0 (tau, &x0, &x0);
	    
    } else if (!strcasecmp (rtype, "LAPLACIAN")) {
	    
	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("--> Type............LAPLACIAN (obsolete)\n");
	    mexPrintf ("--> tau.............%f\n", tau);
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
	} else { // read TV parameters from key/value list
	    for (prm = 0; field = mxGetCell (prhs[4], prm); prm+=2) {
	        AssertArg_Char (field, __func__, 5);
		mxGetString (field, cbuf, 256);
		if (!strcasecmp (cbuf, "Beta")) {
		    beta = mxGetScalar (mxGetCell (prhs[4], prm+1));
		    break;
		}
	    }
	}
	    
	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("--> Type............TV (total variation)\n");
	    mexPrintf ("--> tau.............%f\n", tau);
	    mexPrintf ("--> beta............%f\n", beta);
	    mexPrintf ("--> diffusivity.....%s\n", (brefimg ?
		"from image" : kapref ? "from external array" : "none"));
	    if (brefimg || kapref) {
		mexPrintf ("--> diff. format....%s\n",
			   (istensor ? "tensor" : "scalar"));
		mexPrintf ("--> diff. scale.....%f\n", sdr);
		mexPrintf ("--> diff PM thresh..%f\n", fT);
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
	RVector tauvec;
	xs = 1;
	if (regp) { // read TK0 parameters from structure
	    field = mxGetField(regp,0,"tk0");
	    if (field) field = mxGetField(field,0,"xs");
	    if (field) CopyVector (xs, field);
	} else { // read TK0 parameters from key/value list
	    for (prm = 0; field = mxGetCell (prhs[4], prm); prm+=2) {
		AssertArg_Char (field, __func__, 5);
		mxGetString (field, cbuf, 256);
		if (!strcasecmp (cbuf, "Xs")) {
		    CopyVector (xs, mxGetCell (prhs[4], prm+1));
		} else if (!strcasecmp (cbuf, "Tauvec")) {
		    CopyVector (tauvec, mxGetCell (prhs[4], prm+1));
		    xASSERT(tauvec.Dim() == x0.Dim(), "TauVec: wrong length");
		}
	    }
	}

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("Regularisation: 0th order Tikhonov (TK0)\n");
	    mexPrintf ("--> tau: %f\n", tau);
	}

	reg = new Tikhonov0 (tau, &x0, &xs);
	if (tauvec.Dim() > 0) {
	    ((Tikhonov0*)reg)->SetTau(tauvec);
	    cerr << "Using Tauvec!" << endl;
	}

    } else if (!strcasecmp (rtype, "TK1")) {
	// 1-st order Tikhonov (Laplacian)

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("--> Type............1st order Tikhonov (TK1)\n");
	    mexPrintf ("--> tau.............%f\n", tau);
	    mexPrintf ("--> diffusivity.....%s\n", (brefimg ?
		 "from image" : kapref ? "from external array" : "none"));
	    if (brefimg || kapref) {
		mexPrintf ("--> diff. format....%s\n",
			   (istensor ? "tensor" : "scalar"));
		mexPrintf ("--> diff. scale.....%f\n", sdr);
		mexPrintf ("--> diff PM thresh..%f\n", fT);
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
	    for (prm = 0; field = mxGetCell (prhs[4], prm); prm+=2) {
	        AssertArg_Char (field, __func__, 5);
		mxGetString (field, cbuf, 256);
		if (!strcasecmp (cbuf, "Eps")) {
		    eps = mxGetScalar (mxGetCell (prhs[4], prm+1));
		    break;
		}
	    }
	}

	if (verbosity >= 1) {
	    // echo regularisation parameters
	    mexPrintf ("--> Type............Huber\n");
	    mexPrintf ("--> tau.............%f\n", tau);
	    mexPrintf ("--> eps.............%f\n", eps);
	    mexPrintf ("--> diffusivity.....%s\n", (brefimg ?
		"from image" : kapref ? "from external array" : "none"));
	    if (brefimg || kapref) {
		mexPrintf ("--> diff. format....%s\n",
			   (istensor ? "tensor" : "scalar"));
		mexPrintf ("--> diff. scale.....%f\n", sdr);
		mexPrintf ("--> diff PM thresh..%f\n", fT);
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
	
    } else if (!strcasecmp (rtype, "MRF")) {

	if (verbosity >= 1) {
	    cerr << "Regularisation: MRF" << endl;
	    cerr << "--> tau=" << tau << endl;
	}
	reg = new MRF (tau, &x0, raster);

    } else {
	
	reg = 0;  // "no regularisation"
	return;
	
    }
    
    plhs[0] = mxCreateNumericMatrix (1,1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle (reg);

    cerr << "Regul is " << reg->GetName() << endl;
}

// =========================================================================

void MatlabToast::ClearRegul (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Regularisation *reg = GetRegul (prhs[0]);
    delete reg;
    if (verbosity >= 1)
        mexPrintf ("<Regularisation object deleted>\n");
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

void MatlabToast::RegulSetLocalScaling (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    // regularisation
    Regularisation *reg = GetRegul(prhs[0]);

    // scaling image stack (concatenated for all parameters)
    RVector scale_all;
    CopyVector (scale_all, prhs[1]);

    reg->SetLocalScaling (scale_all);
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

