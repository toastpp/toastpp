#include "stoastlib.h"
#include <fstream>
#include <iomanip>
#include <string.h>
#include "util.h"
#include "solver_cw.h"
#include "source.h"
#include "timing.h"
#include <time.h>
#ifdef TOAST_MPI
#include "mpi.h"
#endif // TOAST_MPI

//#define OUTPUTFIELDS

using namespace std;
using namespace toast;

// ==========================================================================
// Global variables

enum PCG_PRECON {           // preconditioners for PCG solver
    PCG_PRECON_NONE,        //   no preconditioner
    PCG_PRECON_DIAGJTJ,     //   diagonal of JTJ
    PCG_PRECON_SPARSEJTJ,   //   sparse JTJ
    PCG_PRECON_FULLJTJ      //   complete JTJ
} g_pcg_precon = PCG_PRECON_NONE;

enum LM_PRECON {            // preconditioners for LM solver
    LM_PRECON_NONE,         //   no preconditioner
    LM_PRECON_HDIAG,        //   diagonal of Hessian
    LM_PRECON_CH,           //   Cholesky factorisation using explicit Hessian
    LM_PRECON_ICH,          //   Incomplete Cholesky with sparse Hessian
    LM_PRECON_GMRES         //   "matrix-less" Krylov subspace method (GMRES)
} g_lm_precon = LM_PRECON_NONE;

enum PARAM_SCALE {          // parameter scaling method
    PARAM_SCALE_NONE,       //   no scaling
    PARAM_SCALE_AVG,        //   scale with average of initial distribution
    PARAM_SCALE_LOG,        //   log scaling
    PARAM_SCALE_BOUNDLOG    //   x -> ln ((x-a)(b-x))
} g_pscale = PARAM_SCALE_NONE;

double g_lm_gmres_tol;      // convergence criterion for LM GMRES precon

double g_lsolver_tol;       // convergence criterion for linear solver

double g_param_cmuamin, g_param_cmuamax;
double g_param_ckappamin, g_param_ckappamax;
double g_refind;
char g_meshname[256];
int g_nimsize;
int g_imgfmt = IMGFMT_NIM;

// ==========================================================================
// PCG_RESET_INTERVAL: (PCG only) number of iterations between
// preconditioner resets

#define PCG_RESET_INTERVAL 10

// ==========================================================================
// RESCALE_HESSIAN: (PCG and LBFGS only): rescale initial Hessian with
// <y,s>/<y,y> at each iteration

//#define RESCALE_HESSIAN

// ==========================================================================
//#define RESCALE_JTJ
#define JTJ_SCALE 1.0

// Lower cutoff value for diag(JTJ). Only used in CG solver with
// diag(JTJ) preconditioner
#define DJTJ_LIMIT 1e-8

#define LINESEARCH_IN_LM

//#define NEWTON_CG

// ==========================================================================
// preconditioner method for LM solver

#define LM_CHOLESKY 0
#define PURE_CG 1
#define DIAG_JTJ 2
#define LM_GMRES 3

#define LM_PRECON LM_CHOLESKY


//#define COLUMN_NORMALISATION

//#define USE_KAPPA_IN_MUA_PMDF

//#define SPARSE_HESSIAN

// This parameter is used by LM solver and is usually essential!
#define LAPLACIAN_PRIOR

#define OUTPUT_INITIAL_PROJECTION

//#define LINSOLVE_TOL 1e-10;
#define ITMAX 1000
#define MAXREGION 100

// Multi-Wavelength Related Parameters
#define MAX_NOFWLENGTH 50

// Verbose timing output
//#define DO_PROFILE

//struct tms tme;
double clock0;

#ifdef DO_PROFILE
double solver_time = 0.0;
#endif

typedef enum {
    ISOTROPIC,
    NEUMANN
} SourceType;

// =========================================================================
// global parameters

SourceMode srctp = SRCMODE_NEUMANN;   // source type
int bWriteJ    = 0;
int bWriteGrad = 0;
//int bOutputNIM = 1;
//int bOutputRaw = 0;
int bOutputUpdate = 0;
int bOutputGradient = 0;
double avg_cmua = 1.0, avg_ckappa = 1.0;
ParamParser pp;

// =========================================================================
// local prototypes

void OutputProgramInfo ();
void SelectPCGOptions ();
void SelectLMOptions ();
void SelectMesh (char *meshname, QMMesh &mesh);
void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp);
void SelectMeasurementProfile (int &mtype, double &mwidth);
void SelectInitialParams (const Mesh &mesh, MWsolution &msol);
void SelectInitialReferenceParams (const Mesh &mesh, Solution &msol,
    int whichWavel);
void SelectData (int nqm, RVector &data);
void SelectBasis (IVector &gdim, IVector &bdim);
int  SelectImageFormat ();
PARAM_SCALE SelectParamScaling ();
int  SelectSDMode ();

void AssembleRealSystemMatrix (const Mesh &mesh, RGenericSparseMatrix &F,
    const RVector &cmua, const RVector &ckappa, const RVector &c2a);
void AssembleComplexSystemMatrix (const Mesh &mesh, CGenericSparseMatrix &CF,
    const RVector &cmua, const RVector &ckappa, const RVector &c2a,
    const double omega);
void CalcFieldE (const RCompRowMatrix &F, const RPreconditioner *precon,
    const RVector &qvec, RVector &x);
void CalcFieldE (const RCompRowMatrix &FL, const RVector &Fd,
    const RVector &qvec, RVector &x);
double OF_value (const RVector &data, const RVector &proj,
    const RVector &sd);
void OF_gradient_add_Mod (const Raster &raster, const CFwdSolver &FWS,
    const RVector &data, const RVector &proj, const RVector &sd, CVector *dphi,
    CCompRowMatrix &mvec, RVector &grad);
void OF_gradient_add_Arg (const Raster &raster, const CFwdSolver &FWS,
    const RVector &data, const RVector &proj, const RVector &sd, CVector *dphi,
    CCompRowMatrix &mvec, RVector &grad);

void OF_gradient_add_ModArg (const Raster &raster, const CFwdSolver &FWS,
    const RVector &data, const RVector &proj, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad, bool logparam = false,
    const Solution *sol=0);
// adds gradient for mod and arg data types to `grad', given pixel raster
// `raster', forward solver `FWS', data, projection and sd, fields `dphi',
// and measurement vector `mvec'. `sol' is only required if log parameter
// rescaling is required. `sol' must be given in user pixel basis

void OF_gradient_add (const RVector &data, const RVector &proj,
    const RVector &sd, const RDenseMatrix &J, RVector &grad);

void ReadDataFile (char *fname, RVector &data);
void ColumnScale (RDenseMatrix &A, const RVector &scale, int m0=-1, int m1=-1);
//void RowScale (RDenseMatrix &A, const RVector &rscale, int n0=-1, int n1=-1);
void BiScale (RSymMatrix &A, const RVector &rscale);
void LMsolve (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    ObjectiveFunction &OF, int itmax, const RVector &data, const RVector &sd,
    Solution &msol, Solution &bsol, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, double omega, bool logparam, double ftol);
void LMsolveUnder (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    ObjectiveFunction &OF, int itmax, const RVector &data, const RVector &sd,
    Solution &msol, Solution &bsol, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, double omega, bool logparam, double ftol);
void BFGSsolve (CFwdSolver &FWS, const Raster &raster, ObjectiveFunction &OF,
    int itmax, const RVector &data, const RVector &sd, Solution &msol,
    Solution &bsol, const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    double omega, bool logparam, double ftol);
bool LineSearch (CFwdSolver &FWS, const Raster &raster, ObjectiveFunction &OF,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &data, const RVector &sd, double omega, const Solution &s0,
    const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, RVector &proj, bool &proj_valid, bool logparam);
bool LineSearchWithPrior (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, bool logparam, const RVector &x0, const RVector& xs,
    const double &tau, const double tau1, const double tau2,
    const RCompRowMatrix &LTL);
void test_cholesky();
void test_biscale();

// =========================================================================
// MAIN 
int main (int argc, char *argv[])
{
    clock0 = tic();

    if (argc > 1 && pp.Open (argv[1]))
        cout << "Reading parameters from " << argv[1] << endl;
    if (argc > 2) {
        pp.LogOpen (argv[2]);
	cout << "Writing log to " << argv[2] << endl;
    } else {
	pp.LogOpen ("gridbasis.out");
	cout << "Writing log to gridbasis.out" << endl;
    }
    OutputProgramInfo ();

#ifdef TOAST_MPI
    // Initialise MPI
    int mpi_rank, mpi_size;
    int mpi_status = MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
    LOGOUT_1PRM("Initialising MPI session. MPI status = %d", mpi_status);
    LOGOUT_2PRM("Processor %d of %d", mpi_rank, mpi_size);

    cerr << "processor " << mpi_rank << " of " << mpi_size << endl;
#endif // TOAST_MPI

    const double c0 = 0.3;

    char meshname[256], cbuf[256];
    int    qprof, mprof;   // source/measurement profile (0=Gaussian, 1=Cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]
    int sdmode;
    int i, j, n, dim, nQ, nM, nQM, glen, blen, slen;
    IVector gdim, bdim;
    RVector gsize;
    QMMesh mesh;
    Raster *raster;
    Point bbmin, bbmax;
    bool logparam = false;

    SelectMesh (meshname, mesh);
    SelectSourceProfile (qprof, qwidth, srctp);
    SelectMeasurementProfile (mprof, mwidth);
    strcpy (g_meshname, meshname);
    n   = mesh.nlen();
    g_nimsize = n;
    dim = mesh.Dimension();
    nQ  = mesh.nQ;
    nM  = mesh.nM;
    nQM = mesh.nQM;
    mesh.BoundingBox (bbmin, bbmax);
    
    gsize.New(dim);
    gdim.New(dim);
    bdim.New(dim);
    SelectBasis (gdim, bdim);

    for (i = 0, glen = blen = 1; i < dim; i++) {
        gsize[i] = bbmax[i]-bbmin[i];
        glen *= gdim[i];
	blen *= bdim[i];
    }

    g_pscale  = SelectParamScaling ();
    Solver_CW *solver = Solver_CW::Create (&pp); // the nonlinear solver

    RCompRowMatrix qvec, mvec;
    RVector *dphi, *aphi;
    int cmd;
    RFwdSolverMW FWS (pp);
    FWS.WriteParams (pp);
    g_lsolver_tol = FWS.GetLinSolverTol();
    RVector data;
    bool useref = false;

    // --------------------------------------------------------------
    // multi-wavelength data information

    int nofMuachromo; // number of chromophores contributing to mua
    int nofwavel = 0; // number of measurement wavelengths
    double val;
    RVector wlength(MAX_NOFWLENGTH);

    // wavelengths
    if (!pp.GetString("WAVELENGTH", cbuf)) {
	cout << "Enter the wavelengths, ex: 690 750 830\n>> " << flush;
	cin.getline (cbuf, 256);
    }
    pp.PutString ("WAVELENGTH", cbuf);
    istringstream iss(cbuf);
    while (iss >> val) {
	wlength[nofwavel++] = val;
	xASSERT(nofwavel <= MAX_NOFWLENGTH, Max wavelength count exceeded);
    }

    // chromophores
    if (!pp.GetInt("NOFMUACHROMOPHORES", nofMuachromo)) {
	cout << "Number of mua chromophores?\n>> " << flush;
	cin >> nofMuachromo;
    }
    pp.PutInt ("NOFMUACHROMOPHORES", nofMuachromo);

    // extinction coefficients
    RDenseMatrix extcoef(nofwavel, nofMuachromo);
    double ecoeff;
    for (i = 0; i < nofwavel; i++) {
	for  (j = 0; j < nofMuachromo; j++) {
	    char exttype[30];
	    sprintf (exttype,"EXTINC_WAVEL_%d_CHROMO_%d", int(wlength[i]),j+1);
    
	    if (!pp.GetReal (exttype, extcoef(i,j))) {
		cout << "Enter the extinction value for Chromophore " << j+1
		     << " at wavelength " << wlength[i] << "\n>> " << flush;
		cin >> ecoeff;
		extcoef(i, j) = ecoeff;
	    }
	    pp.PutReal (exttype, extcoef(i,j) );
	}
    }
    
    // multi-wavelength solution in mesh basis
    int nprm = nofMuachromo; // chromophores
    nprm += 2;               // scattering parameters A and b
    nprm += 1;               // refractive index
    nprm += nofwavel;        // background mua at each wavelength
    MWsolution msol(nprm, n, nofwavel, extcoef, wlength);

    SelectInitialParams (mesh, msol);

    SelectData (nQM * nofwavel, data);
    sdmode = SelectSDMode();

    // Generate sd vector by assuming sd=data ('normalised')
    RVector sd(data);

    // build the source vectors
    qvec.New (nQ, n);
    for (i = 0; i < nQ; i++) {
	RVector q(n);
	switch (qprof) {
	case 0:
	    q = QVec_Point (mesh, mesh.Q[i], srctp);
	    break;
	case 1:
	    q = QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp);
	    break;
	case 2:
	    q = QVec_Cosine (mesh, mesh.Q[i], qwidth, srctp);
	    break;
	}
	qvec.SetRow (i, q);
    }

    // build the measurement vectors
    mvec.New (nM, n);
    for (i = 0; i < nM; i++) {
	RVector m(n);
	switch (mprof) {
	case 0:
	    m = QVec_Point (mesh, mesh.M[i], SRCMODE_NEUMANN);
	    break;
	case 1:
	    m = QVec_Gaussian (mesh, mesh.M[i], mwidth, SRCMODE_NEUMANN);
	    break;
	case 2:
	    m = QVec_Cosine (mesh, mesh.M[i], mwidth, SRCMODE_NEUMANN);
	    break;
	}
	for (j = 0; j < n; j++) m[j] *= mesh.plist[j].C2A();
	mvec.SetRow (i, m);
    }


    // build the field vectors
    dphi = new RVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new RVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // allocate system matrix
    cout << endl << "Allocating system matrix" << endl;
    FWS.Allocate (mesh);

    cout << "Generating intermediate pixel grid" << endl;
    
    if (blen == glen) // use single grid basis
	raster = new Raster_Pixel (gdim, gdim, &mesh);
    else              // use intermediate basis
	raster = new Raster_Pixel (gdim, bdim, &mesh);

    glen = raster->GLen();
    blen = raster->BLen();
    slen = raster->SLen();

    // solution in sparse user basis
    Solution bsol(nprm, raster->SLen());
    for (i = 0; i < nprm; i++)
	bsol.SetActive (i, msol.IsActive(i));
    raster->Map_MeshToSol (msol, bsol, true);

    // output image format
    g_imgfmt = SelectImageFormat();

    // reference data
    if (!pp.GetBool ("REFDATA", useref)) {
	cout << "\nUse reference data\n[1|0] >> ";
	cin >> cmd;
	useref = (cmd != 0);
    }
    pp.PutBool ("REFDATA", useref);

    // debugging flags
    if (!pp.GetInt ("DBG_WRITEJACOBIAN", bWriteJ)) {
	cout << "Output Jacobian during reconstruction (1/0)?\n>> " << flush;
	cin >> bWriteJ;
    }
    pp.PutInt ("DBG_WRITEJACOBIAN", bWriteJ);

    if (!pp.GetInt ("DBG_WRITEGRAD", bWriteGrad)) {
	cout << "Output gradient during reconstruction (1/0)?\n>> " << flush;
	cin >> bWriteGrad;
    }
    pp.PutInt ("DBG_WRITEGRAD", bWriteGrad);

    // set parameter scaling
    Scaler *pscaler;
    switch (g_pscale) {
    case PARAM_SCALE_NONE:
        pscaler = new NullScaler;
	break;
    case PARAM_SCALE_AVG: {
        RVector scale = bsol.GetActiveParams();
	for (i = 0; i < bsol.nActive(); i++) {
	    double sum = 0;
	    for (j = 0; j < slen; j++)
	        sum += scale[i*slen+j];
	    for (j = 0; j < slen; j++)
	        scale[i*slen+j] = slen/sum;
	}
	pscaler = new ConstScaler (scale);
	} break;
    case PARAM_SCALE_LOG:
        pscaler = new LogScaler;
	break;
    case PARAM_SCALE_BOUNDLOG:
	RVector xmin(slen), xmax(slen);
	for (i = 0; i < slen/2; i++) {
	    xmin[i] = g_param_cmuamin;
	    xmax[i] = g_param_cmuamax;
	    xmin[i+slen/2] = g_param_ckappamin;
	    xmax[i+slen/2] = g_param_ckappamax;
	}
	pscaler = new BoundLogScaler (xmin, xmax);
	break;
    }

    cout << "  Mapped solution" << endl;
    double v1, v2;
    v1 = vmin (bsol.GetParam(OT_CMUA));
    v2 = vmax (bsol.GetParam(OT_CMUA));
    cout << "    CMUA:   " << v1 << " to " << v2 << endl;
    v1 = vmin (bsol.GetParam(OT_CKAPPA));
    v2 = vmax (bsol.GetParam(OT_CKAPPA));
    cout << "    CKAPPA: " << v1 << " to " << v2 << endl;

    // map bsol back to msol to make sure msol contains a valid
    // basis representation of the solution
    raster->Map_SolToMesh (bsol, msol, true);

    cout << "  Re-mapped mesh:\n";
    v1 = vmin (msol.GetParam(OT_CMUA));
    v2 = vmax (msol.GetParam(OT_CMUA));
    cout << "    CMUA:   " << v1 << " to " << v2 << endl;
    v1 = vmin (msol.GetParam(OT_CKAPPA));
    v2 = vmax (msol.GetParam(OT_CKAPPA));
    cout << "    CKAPPA: " << v1 << " to " << v2 << endl;

    // convert chromophores, Scatter Params and N to cmua, ckappa, n
    msol.RegisterChange();

    switch (g_imgfmt) {
    case IMGFMT_NIM:
	for (i = 0; i < msol.nParam(); i++) {
	    char fname[256];
	    if (msol.IsActive(i)) {
		if (i < msol.nParam() - 3) 
		    sprintf (fname,"reconChromophore_%d.nim",i+1);
		if (i == msol.nParam() - 3)
		    sprintf (fname,"reconScatPrefactor_A.nim");
		if (i == msol.nParam() - 2) 
		    sprintf (fname,"reconScatPower_b.nim");
		WriteNimHeader (meshname, n, fname, "N/A");  
		msol.WriteImgGeneric (0, fname, i);
	    }
	}
	break;
    case IMGFMT_RAW:
	for (i = 0; i < bsol.nParam(); i++) {
	    if (bsol.IsActive(i)) {
		char fname[256];
		RVector img(blen);
		sprintf (fname, "reconParam_%d.raw", i+1);
		raster->Map_SolToBasis(bsol.GetParam(i), img);
		WriteRimHeader (raster->BDim(), fname);
		Solution::WriteImgGeneric (0, fname, img);
	    }
	}
	break;
    }

    if (bOutputUpdate) {
	WriteNimHeader (meshname, n, "update_mua.nim", "MUA");
	WriteNimHeader (meshname, n, "update_mus.nim", "MUS");
    }
    if (bOutputGradient) {
	WriteNimHeader (meshname, n, "gradient_mua.nim", "MUA");
	WriteNimHeader (meshname, n, "gradient_mus.nim", "MUS");
    }

    if (bWriteGrad) {
        WriteRimHeader (raster->BDim(), "grad_mua.raw");
	WriteRimHeader (raster->BDim(), "grad_mus.raw");
    }

    if (useref) {
        RVector refdata(nQM*nofwavel);
	RVector proj(nQM*nofwavel);
	bool refEqualinitial;

	if (!pp.GetString ("REFDATA_MOD", cbuf)) {
	    cout << "Reference file: ";
	    cin  >> cbuf;
	}
	ReadDataFile (cbuf, refdata);
	pp.PutString ("REFDATA_MOD", cbuf);

	if (!pp.GetBool ("REFEQUALINITIAL", refEqualinitial)) {
	    cout << "\nAre reference optical properties equal to initial "
		 << "guess\n[1|0] >> ";
	    cin >> cmd;
	    refEqualinitial = (cmd != 0);
	}
	pp.PutBool ("REFEQUALINITIAL", refEqualinitial);
	
	data -= refdata;
	
	for (i = 0; i < nofwavel; i++) {
	    cout << "Assembling and pre-processing system matrix at "
		 << "wavelength " << wlength[i] <<endl;
	  
	    if (refEqualinitial) {
		// use if ref state = initial guess 
		FWS.Reset (*msol.swsol[i], 0); 
	    } else {
		// reference optical properties are different to initial guess
		Solution refsol(OT_NPARAM, n);
		SelectInitialReferenceParams (mesh, refsol, int(wlength[i]));
		// now to be consistent with msol first map to a basis sol
		// then map it back!
		Solution dummybsol(OT_NPARAM, raster->SLen());
		raster->Map_MeshToSol (refsol, dummybsol, true);
		raster->Map_SolToMesh (dummybsol, refsol, true);
		FWS.Reset (refsol, 0);
	    }
          	  
	    cout << "Generating model baseline" << endl;
	    // this returns Phi
	    FWS.CalcFields (qvec, dphi/*, qcoup_lnmod*/);
	    // this returns Log(Phi) projected on Detectors
	    RVector proj_i(proj, i*nQM, nQM);
	    proj_i = FWS.ProjectAll (mvec, dphi/*, mcoup_lnmod*/);
	}
	data += proj;
    }

    // ==================================================================
    // Echo some essential parameters
    cout    << endl << "------------" << endl;
    logfile << endl << "------------" << endl;
    cout    << "Solution basis size: " << slen << endl;
    logfile << "Solution basis size: " << slen << endl; 
    cout    << "Basis grid:          " << bdim << endl;
    logfile << "Basis grid:          " << bdim << endl;
    if (blen != glen) {
	cout    << "Intermediate grid:      " << gdim << endl;
	logfile << "Intermediate grid:      " << gdim << endl;
    }

    // Select data scaling method
    switch (sdmode) {
    case 1:   // no scaling
	for (i = 0; i < sd.Dim(); i++)
	    sd[i] = 1.0;
	break;
    case 2:   // scale with individual data
	for (i = 0; i < sd.Dim(); i++)
	    sd[i] = data[i];
	break;
    case 3: { // scale with averages over data types
	cout << "Generating model baseline" << endl;
	RVector proj(nQM*2);
	FWS.CalcFields (qvec, dphi);
	proj = FWS.ProjectAll (mvec, dphi);
	sd = proj;
	double avd1 = 0.0, avd2 = 0.0;
	for (i = 0; i < nQM; i++) {
	    avd1 += sd[i]*sd[i];
	    avd2 += sd[i+nQM]*sd[i+nQM];
	}
	avd1 = sqrt(avd1);
	avd2 = sqrt(avd2);
	for (i = 0; i < nQM; i++) {
	    sd[i] = avd1;
	    sd[i+nQM] = avd2;
	}
	cout << "Averages: amplitude " << avd1 << ", phase " << avd2
	     << endl << endl;
	} break;		
    case 4: { // scale with difference averages over data types
	cout << "Generating model baseline" << endl;
	RVector proj(nQM*nofwavel);
	for (i = 0; i < nofwavel; i++) {
	    FWS.Reset (*msol.swsol[i],0);
	    FWS.CalcFields (qvec, dphi);
	    RVector proj_i(proj, i*nQM, nQM);
	    proj_i = FWS.ProjectAll (mvec, dphi);
	}
	sd = (data - proj);
	double avd1 = 0.0, avd2 = 0.0;
	for (i = 0; i < sd.Dim(); i++) {
	    avd1 += sd[i]*sd[i];
	}
	avd1 = sqrt(avd1);
	for (i = 0; i < sd.Dim(); i++) {
	    sd[i] = avd1;
	}
	cout << "Averages: amplitude " << avd1 << endl << endl;
	} break;
    }

    // set up objective function
    ObjectiveFunction OF (raster);
    OF.SetData (&data);
    OF.SetSD (&sd);

    cout << "Ref data  " << (useref ? "YES":"NO") << endl;
    cout << endl;

    // ==================================================================
    // Start the solver
    solver->Solve (FWS, *raster, pscaler, OF, data, sd, bsol, msol, qvec,mvec);
    delete solver;

#ifdef UNDEF
    switch (g_solver) {
	//case SOLVER_PCG:
        //cout << "Entering conjugate gradient solver ..." << endl;
	//CGsolve (FWS, raster, OF, ITMAX, data, sd, msol, bsol, qvec,
	//    mvec, omega, logparam, 1e-8);
	//break;
	//case SOLVER_LM:
        //cout << "Entering Levenberg-Marquardt solver ..." << endl;
	//LMsolve (FWS, raster, pscaler, OF, ITMAX, data, sd, msol, bsol, qvec,
	//    mvec, omega, logparam, 1e-6);
	//break;
    case SOLVER_LM_UNDER:
        cout << "Entering underdetermined LM solver ..." << endl;
	LMsolveUnder (FWS, raster, pscaler, OF, ITMAX, data, sd, msol, bsol,
            qvec, mvec, omega, logparam, 1e-6);
	break;
    case SOLVER_BFGS:
        cout << "Entering BFGS solver ..." << endl;
	BFGSsolve (FWS, raster, OF, ITMAX, data, sd, msol, bsol, qvec,
	    mvec, omega, logparam, 1e-6);
	break;
    case SOLVER_LINEAR:
	cout << "Entering linear solver ..." << endl;
	LinSolve (FWS, raster, pscaler, OF, data, sd, msol, bsol, qvec,
	    mvec, omega, logparam, 1e-6);
	break;
    default:
        cerr << "Nonlinear solver not defined. Aborting." << endl;
	exit (1);
    }
#endif

    // cleanup
    delete []dphi;
    delete []aphi;

    //times (&tme);
    //double total_time = (double)(tme.tms_utime-clock0)/(double)HZ;
    double total_time = toc(clock0);
    LOGOUT ("Final timings:");
    LOGOUT_1PRM ("Total: %f", total_time);
#ifdef DO_PROFILE
    LOGOUT_1PRM ("Solver: %f", solver_time);
#endif

#ifdef TOAST_MPI
    MPI_Finalize();
#endif

    return 0;
}                                                                              

void OutputProgramInfo ()
{
    char cbuf[256], cwd[250], *host;
    time_t tme = time(0);

    pp.Lineout ("+-------------------------------------------------+");
    pp.Lineout ("| Reconstruction of parameter distribution of the |");
    pp.Lineout ("| diffusion equation from frequency-domain data   |");
    pp.Lineout ("+-------------------------------------------------+");
    pp.Lineout (VERSION_STRING);
    sprintf (cbuf, "Executed %s", ctime(&tme));
    if (host = getenv("HOST"))
        sprintf (cbuf+strlen(cbuf), "on host %s ", host);
    sprintf (cbuf+strlen(cbuf), "(PID %d)", getpid());
    pp.Lineout (cbuf);
    if (getcwd (cwd, 250)) {
        sprintf (cbuf, "CWD: %s", cwd);
	pp.Lineout (cbuf);
    }
    pp.Lineout ("===================================================");
}

// ============================================================================

void SelectPCGOptions()
{
    char cbuf[256];
    bool def = false;

    // parse from definition file
    if (pp.GetString ("PCG_PRECON", cbuf)) {
        if (!strcasecmp (cbuf, "NONE"))
	    g_pcg_precon = PCG_PRECON_NONE, def = true;
	else if (!strcasecmp (cbuf, "DIAGJTJ"))
	    g_pcg_precon = PCG_PRECON_DIAGJTJ, def = true;
	else if (!strcasecmp (cbuf, "SPARSEJTJ"))
	    g_pcg_precon = PCG_PRECON_SPARSEJTJ, def = true;
	else if (!strcasecmp (cbuf, "FULLJTJ"))
	    g_pcg_precon = PCG_PRECON_FULLJTJ, def = true;
    }

    // ask user
    while (!def) {
        int cmd;
	cout << "\nSelect PCG preconditioner:\n";
	cout << "(1) None\n";
	cout << "(2) diagonal of JTJ\n";
	cout << "(3) sparse JTJ\n";
	cout << "(4) full JTJ\n";
	cout << ">> ";
	cin >> cmd;
	switch (cmd) {
	case 1: g_pcg_precon = PCG_PRECON_NONE, def = true; break;
	case 2: g_pcg_precon = PCG_PRECON_DIAGJTJ, def = true; break;
	case 3: g_pcg_precon = PCG_PRECON_SPARSEJTJ, def = true; break;
	case 4: g_pcg_precon = PCG_PRECON_FULLJTJ, def = true; break;
	}
    }

    // write back
    switch (g_pcg_precon) {
    case PCG_PRECON_NONE:      pp.PutString ("PCG_PRECON", "NONE"); break;
    case PCG_PRECON_DIAGJTJ:   pp.PutString ("PCG_PRECON", "DIAGJTJ"); break;
    case PCG_PRECON_SPARSEJTJ: pp.PutString ("PCG_PRECON", "SPARSEJTJ"); break;
    case PCG_PRECON_FULLJTJ:   pp.PutString ("PCG_PRECON", "FULLJTJ"); break;
    }
}

// ============================================================================

void SelectLMOptions()
{
    char cbuf[256];
    bool def = false;

    // 1. === PRECONDITIONER ===

    // parse from definition file
    if (pp.GetString ("LM_PRECON", cbuf)) {
        if (!strcasecmp (cbuf, "NONE"))
	    g_lm_precon = LM_PRECON_NONE,  def = true;
	else if (!strcasecmp (cbuf, "HDIAG"))
	    g_lm_precon = LM_PRECON_HDIAG, def = true;
	else if (!strcasecmp (cbuf, "CH"))
	    g_lm_precon = LM_PRECON_CH,    def = true;
	else if (!strcasecmp (cbuf, "ICH"))
	    g_lm_precon = LM_PRECON_ICH,   def = true;
	else if (!strcasecmp (cbuf, "GMRES"))
	    g_lm_precon = LM_PRECON_GMRES, def = true;
    }
    // ask user
    while (!def) {
        int cmd;
	cout << "\nSelect LM preconditioner:\n";
	cout << "(1) None\n";
	cout << "(2) Diagonal of Hessian\n";
	cout << "(3) Cholesky factorisation of full Hessian\n";
	cout << "(4) Incomplete CH factorisation of sparse Hessian\n";
	cout << "(5) GMRES solver with implicit Hessian\n";
	cout << ">> ";
	cin >> cmd;
	switch (cmd) {
	case 1: g_lm_precon = LM_PRECON_NONE,  def = true; break;
	case 2: g_lm_precon = LM_PRECON_HDIAG, def = true; break;
	case 3: g_lm_precon = LM_PRECON_CH,    def = true; break;
	case 4: g_lm_precon = LM_PRECON_ICH,   def = true; break;
	case 5: g_lm_precon = LM_PRECON_GMRES, def = true; break;
	}
    }
    // write back
    switch (g_lm_precon) {
    case LM_PRECON_NONE:  pp.PutString ("LM_PRECON", "NONE");  break;
    case LM_PRECON_HDIAG: pp.PutString ("LM_PRECON", "HDIAG"); break;
    case LM_PRECON_CH:    pp.PutString ("LM_PRECON", "CH");    break;
    case LM_PRECON_ICH:   pp.PutString ("LM_PRECON", "ICH");   break;
    case LM_PRECON_GMRES: pp.PutString ("LM_PRECON", "GMRES"); break;
    }

    // 2. === GMRES CONVERGENCE CRITERION ===

    if (g_lm_precon == LM_PRECON_GMRES) {
        if (!pp.GetReal ("LM_GMRES_TOL", g_lm_gmres_tol) ||
	  g_lm_gmres_tol <= 0.0) do {
	    cout << "\nSelect LM GMRES convergence criterion (>0):\n";
	    cout << ">> ";
	    cin >> g_lm_gmres_tol;
	} while (g_lm_gmres_tol <= 0.0);
	pp.PutReal ("LM_GMRES_TOL", g_lm_gmres_tol);
    }
}

// ============================================================================

void SelectMesh (char *meshname, QMMesh &mesh)
{
    char qmname[256];

    if (!pp.GetString ("MESHFILE", meshname)) {
        cout << "\nMesh file name:\n>> ";
	cin >> meshname;
    }
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();

    if (!pp.GetString ("QMFILE", qmname)) {
        cout << "\nQM file name:\n>> ";
	cin >> qmname;
    }
    ifstream qmf (qmname);
    mesh.LoadQM (qmf);

    // write back
    pp.PutString ("MESHFILE", meshname);
    pp.PutString ("QMFILE", qmname);
}

// ============================================================================

void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp)
{
    char cbuf[256];
    int cmd;

    bool typeok = false;
    if (pp.GetString ("SOURCETYPE", cbuf)) {
	if (!strcasecmp (cbuf, "NEUMANN")) {
	    srctp = SRCMODE_NEUMANN;
	    typeok = true;
	} else if (!strcasecmp (cbuf, "ISOTROPIC")) {
	    srctp = SRCMODE_ISOTROPIC;
	    typeok = true;
	}
    }
    while (!typeok) {
	cout << "\nSource type:\n";
	cout << "(1) Neumann boundary source\n";
	cout << "(2) Isotropic point source\n";
	cout << "[1|2] >> ";
	cin  >> cmd;
	switch (cmd) {
	    case 1: srctp = SRCMODE_NEUMANN;   typeok = true; break;
	    case 2: srctp = SRCMODE_ISOTROPIC; typeok = true; break;
	}
    }
    pp.PutString ("SOURCETYPE",
        srctp == SRCMODE_NEUMANN ? "NEUMANN" : "ISOTROPIC");

    qtype = -1;
    if (pp.GetString ("SOURCEPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    qtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    qtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    qtype = 2;
	}
    }
    while (qtype < 0) {
	cout << "\nSource profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> qtype;
	qtype -= 1;
    }
    if (qtype > 0 && !pp.GetReal ("SOURCEWIDTH", qwidth)) {
	switch (qtype) {
	case 1:
	    cout << "\nSource 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nSource support radius [mm]:\n>> ";
	    break;
	}
	cin >> qwidth;
    }
    switch (qtype) {
    case 0:
	pp.PutString ("SOURCEPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("SOURCEPROFILE", "GAUSSIAN");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    case 2:
	pp.PutString ("SOURCEPROFILE", "COSINE");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    }
}

// ============================================================================

void SelectMeasurementProfile (int &mtype, double &mwidth)
{
    char cbuf[256];
    mtype = -1;
    if (pp.GetString ("MEASUREMENTPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    mtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    mtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    mtype = 2;
	}
    }
    while (mtype < 0) {
	cout << "\nMeasurement profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> mtype;
	mtype -= 1;
    }
    if (mtype > 0 && !pp.GetReal ("MEASUREMENTWIDTH", mwidth)) {
	switch (mtype) {
	case 1:
	    cout << "\nMeasurement 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nMeasurement support radius [mm]:\n>> ";
	    break;
	}
	cin >> mwidth;
    }
    switch (mtype) {
    case 0:
	pp.PutString ("MEASUREMENTPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("MEASUREMENTPROFILE", "GAUSSIAN");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    case 2:
	pp.PutString ("MEASUREMENTPROFILE", "COSINE");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    }
}

// ============================================================================

int ScanRegions (const Mesh &mesh, int *nregnode)
{
    int i, reg, nreg;
    for (i = 0; i < MAXREGION; i++) nregnode[i] = 0;
    for (i = 0; i < mesh.nlen(); i++) {
	reg = mesh.nlist[i].Region();
	if (reg >= 0 && reg < MAXREGION) nregnode[reg]++;
    }
    for (nreg = i = 0; i < MAXREGION; i++)
	if (nregnode[i]) nreg++;
    return nreg;
}

void SelectInitialParams (const Mesh &mesh, MWsolution &msol)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    int nparam = msol.nParam();
    RVector param[nparam];
    int i, j, k, n, p, nreg, nregnode[MAXREGION];

    for (p = 0; p < nparam; p++) {
      char resetstr[256];
      char reconstr[256];

      if (p < msol.nmuaChromo) {
	  sprintf (resetstr,"CHROMOPHORE_%d",p+1);
	  sprintf (reconstr, "RECON_CHROMOPHORE_%d",p+1);
      } else if (p == msol.nmuaChromo) {
	  sprintf (resetstr,"SCATTERING_PREFACTOR_A");
	  sprintf (reconstr, "RECON_SCATTERING_PREFACTOR_A");
      } else if (p == msol.nmuaChromo+1) {
	  sprintf (resetstr,"SCATTERING_POWER_B");
	  sprintf (reconstr, "RECON_SCATTERING_POWER_B");
      } else if (p == msol.nmuaChromo+2) {
	  sprintf (resetstr,"RESET_N");
      } else {
	  sprintf (resetstr, "BACKGROUND_MUA_%d",
		   (int)msol.wlength[p-msol.nmuaChromo-3]);
      }

      param[p].New(mesh.nlen());
	if (pp.GetString (resetstr, cbuf)) {
	    pp.PutString (resetstr, cbuf);
	    //    if (!strcasecmp (cbuf, "MESH")) {
	    //	param[p] = mesh.plist.Param(prmtp[p]);
	    //} else
	    if (!strncasecmp (cbuf, "HOMOG", 5)) {
		sscanf (cbuf+5, "%lf", &prm);
		param[p] = prm;
	    } else if (!strncasecmp (cbuf, "REGION_HOMOG", 12)) {
		valstr = strtok (cbuf+12, " \t");
		for (n = 0; n < MAXREGION && valstr; n++) {
		    sscanf (valstr, "%lf", reg_prm+n);
		    valstr = strtok (NULL, " \t");
		}
		nreg = ScanRegions (mesh, nregnode);
		for (i = k = 0; k < n && i < MAXREGION; i++) {
		    if (nregnode[i]) {
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = reg_prm[k];
			k++;
		    }
		}	     
	    } else if (!strncasecmp (cbuf, "NIM", 3)) {
		ReadNim (cbuf+4, param[p]);
	    }
	} else {
	    cout << "\nSelect initial distribution for " << resetstr
		 << endl;
	    // cout << "(1) Use values stored in mesh\n";
	    cout << "(1) Global homogeneous\n";
	    cout << "(2) Homogeneous in regions\n";
	    cout << "(3) Nodal image file (NIM)\n";
	    cout << "[1|2|3] >> ";
	    cin >> resettp;
	    switch (resettp) {
	      //case 1:
	      //	param[p] = mesh.plist.Param(prmtp[p]);
	      //	strcpy (cbuf, "MESH");
	      //	break;
	    case 1:
		cout << "\nGlobal value:\n>> ";
		cin >> prm;
		param[p] = prm;
		sprintf (cbuf, "HOMOG %f", prm);
		break;
	    case 2:
		nreg = ScanRegions (mesh, nregnode);
		strcpy (cbuf, "REGION_HOMOG");
		cout << "\nFound " << nreg << " regions\n";
		for (i = 0; i < MAXREGION; i++) {
		    if (nregnode[i]) {
			cout << "Value for region " << i << " (" << nregnode[i]
			     << " nodes):\n>> ";
			cin >> prm;
			sprintf (cbuf+strlen(cbuf), " %f", prm);
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = prm;
		    }
		}
		break;
	    case 3:
		cout << "\nNIM file name:\n>> ";
		strcpy (cbuf, "NIM ");
		cin >> cbuf+4;
		ReadNim (cbuf+4, param[p]);
		break;
	    }
	    pp.PutString (resetstr, cbuf);
	}

	if (p < msol.nmuaChromo+2) {
	  bool active;
	    if (!pp.GetBool (reconstr, active)) {
		char c;
		do {
		    cout << "\nActivate " << resetstr
			 << " for reconstruction?\n[y|n] >> ";
		    cin >> c;
		} while (c != 'y' && c != 'n');
		active = (c == 'y');
	    }
	    pp.PutBool (reconstr, active);
	    msol.SetActive (p, active);
	}
	msol.SetParam (p, param[p]);
    }
    g_refind = mean (param[msol.nmuaChromo+2]); // assuming homogeneous n
    
    //msol.SetParam (OT_CMUA, param[0]*c0/param[2]);
    //msol.SetParam (OT_CKAPPA, c0/(3*param[2]*(param[0]+param[1])));
    //for (i = 0; i < param[2].Dim(); i++)
    //	param[2][i] = c0/(2*param[2][i]*A_Keijzer(param[2][i]));
    //msol.SetParam (OT_C2A, param[2]);   
}

// ===========================================================================

void SelectInitialReferenceParams (const Mesh &mesh, Solution &msol,
    int whichWavel)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    RVector param[3];
    int i, j, k, n, p, nreg, nregnode[MAXREGION];
    const char *rootstr[3] = {"RESET_REF_MUA", "RESET_REF_MUS", "RESET_REF_N"};
    const ParameterType prmtp[3] = {PRM_MUA, PRM_MUS, PRM_N};
    for (p = 0; p < 3; p++) {
        char resetstr[256];
	sprintf (resetstr,"%s_%d", rootstr[p], whichWavel);

	param[p].New(mesh.nlen());
	if (pp.GetString (resetstr, cbuf)) {
	    pp.PutString (resetstr, cbuf);
	    if (!strcasecmp (cbuf, "MESH")) {
		param[p] = mesh.plist.Param(prmtp[p]);
	    } else if (!strncasecmp (cbuf, "HOMOG", 5)) {
		sscanf (cbuf+5, "%lf", &prm);
		param[p] = prm;
	    } else if (!strncasecmp (cbuf, "REGION_HOMOG", 12)) {
		valstr = strtok (cbuf+12, " \t");
		for (n = 0; n < MAXREGION && valstr; n++) {
		    sscanf (valstr, "%lf", reg_prm+n);
		    valstr = strtok (NULL, " \t");
		}
		nreg = ScanRegions (mesh, nregnode);
		for (i = k = 0; k < n && i < MAXREGION; i++) {
		    if (nregnode[i]) {
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = reg_prm[k];
			k++;
		    }
		}	     
	    } else if (!strncasecmp (cbuf, "NIM", 3)) {
		ReadNim (cbuf+4, param[p]);
	    }
	} else {
	    cout << "\nSelect initial distribution for " << resetstr
		 << endl;
	    cout << "(1) Use values stored in mesh\n";
	    cout << "(2) Global homogeneous\n";
	    cout << "(3) Homogeneous in regions\n";
	    cout << "(4) Nodal image file (NIM)\n";
	    cout << "[1|2|3|4] >> ";
	    cin >> resettp;
	    switch (resettp) {
	    case 1:
		param[p] = mesh.plist.Param(prmtp[p]);
		strcpy (cbuf, "MESH");
		break;
	    case 2:
		cout << "\nGlobal value:\n>> ";
		cin >> prm;
		param[p] = prm;
		sprintf (cbuf, "HOMOG %f", prm);
		break;
	    case 3:
		nreg = ScanRegions (mesh, nregnode);
		strcpy (cbuf, "REGION_HOMOG");
		cout << "\nFound " << nreg << " regions\n";
		for (i = 0; i < MAXREGION; i++) {
		    if (nregnode[i]) {
			cout << "Value for region " << i << " (" << nregnode[i]
			     << " nodes):\n>> ";
			cin >> prm;
			sprintf (cbuf+strlen(cbuf), " %f", prm);
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = prm;
		    }
		}
		break;
	    case 4:
		cout << "\nNIM file name:\n>> ";
		strcpy (cbuf, "NIM ");
		cin >> cbuf+4;
		ReadNim (cbuf+4, param[p]);
		break;
	    }
	    pp.PutString (resetstr, cbuf);
	}
    }
    g_refind = mean (param[2]); // assuming homogeneous n
    //msol.SetRefind (g_refind);
    msol.SetParam (OT_CMUA, param[0]*c0/param[2]);
    msol.SetParam (OT_CKAPPA, c0/(3.0*param[2]*(param[0]+param[1])));
    for (i = 0; i < param[2].Dim(); i++)
	param[2][i] = c0/(2.0*param[2][i]*A_Keijzer(param[2][i]));
    msol.SetParam (OT_C2A, param[2]);
}

// ============================================================================

void SelectData (int nqm, RVector &data)
{
    char cbuf[256];
    bool def = false;
    RVector pdata;

    pp.PutString ("DATATYPE", "INTENSITY");
    data.New (nqm);
    if (!pp.GetString ("DATA_INTENSITY", cbuf)) {
	cout << "\nData file (log intensity)\n>> ";
	cin >> cbuf;
    }
    ReadDataFile (cbuf, data);
    pp.PutString ("DATA_INTENSITY", cbuf);
}

// ============================================================================

void SelectBasis (IVector &gdim, IVector &bdim)
{
    char cbuf[256];
    int i, n[3], dim = gdim.Dim();

    if (!pp.GetString ("BASIS", cbuf) ||
      sscanf (cbuf, "%d%d%d", n+0, n+1, n+2) < dim) {
        cout << "\nDimension of recon pixel basis [x y";
	cout << (dim > 2 ? " z]":"]") << endl << ">> ";
	cin >> n[0] >> n[1];
	if (dim > 2) cin >> n[2];
    }
    for (i = 0; i < dim; i++) bdim[i] = n[i];

    if (!pp.GetString ("GRID", cbuf) ||
      sscanf (cbuf, "%d%d%d", n+0, n+1, n+2) < dim) {
	char cmd;
	cout << "\nUse intermediate grid for mapping between mesh and basis?";
	cout << "\n[y|n] >> ";
	cin >> cmd;
	if (toupper (cmd == 'Y')) {
	    cout << "\nDimension of intermediate grid [x y";
	    cout << (dim > 2 ? " z]":"]") << endl << ">> ";
	    cin >> n[0] >> n[1];
	    if (dim > 2) cin >> n[2];
	} else {
	    for (i = 0; i < dim; i++) n[i] = bdim[i];
	}
    }
    for (i = 0; i < dim; i++) gdim[i] = n[i];

    // write back
    if (dim > 2) {
        sprintf (cbuf, "%d %d %d", gdim[0], gdim[1], gdim[2]);
	pp.PutString ("GRID", cbuf);
	sprintf (cbuf, "%d %d %d", bdim[0], bdim[1], bdim[2]);
	pp.PutString ("BASIS", cbuf);
    } else {
        sprintf (cbuf, "%d %d", gdim[0], gdim[1]);
	pp.PutString ("GRID", cbuf);
	sprintf (cbuf, "%d %d", bdim[0], bdim[1]);
	pp.PutString ("BASIS", cbuf);
    }
}

// ============================================================================

PARAM_SCALE SelectParamScaling ()
{
    char cbuf[256];
    bool def = false;
    PARAM_SCALE ps;

    if (pp.GetString ("PARAM_SCALE", cbuf)) {
        if (!strcasecmp (cbuf, "NONE"))
	    ps = PARAM_SCALE_NONE, def = true;
	else if (!strcasecmp (cbuf, "AVG"))
	    ps = PARAM_SCALE_AVG,  def = true;
	else if (!strcasecmp (cbuf, "LOG"))
	    ps = PARAM_SCALE_LOG,  def = true;
	else if (!strcasecmp (cbuf, "LOG_BOUNDED"))
	    ps = PARAM_SCALE_BOUNDLOG, def = true;
    }
    while (!def) {
        int cmd;
	cout << "\nSelect parameter scaling method:\n";
	cout << "(0) None\n";
	cout << "(1) Initial parameter averages\n";
	cout << "(2) Log parameters\n";
	cout << "(3) Bounded log parameters\n";
	cout << "[0|1|2|3] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: ps = PARAM_SCALE_NONE,     def = true; break;
	case 1: ps = PARAM_SCALE_AVG,      def = true; break;
	case 2: ps = PARAM_SCALE_LOG,      def = true; break;
	case 3: ps = PARAM_SCALE_BOUNDLOG, def = true; break;
	}
    }
    switch (ps) {
    case PARAM_SCALE_NONE:     pp.PutString ("PARAM_SCALE", "NONE"); break;
    case PARAM_SCALE_AVG:      pp.PutString ("PARAM_SCALE", "AVG");  break;
    case PARAM_SCALE_LOG:      pp.PutString ("PARAM_SCALE", "LOG");  break;
    case PARAM_SCALE_BOUNDLOG: pp.PutString ("PARAM_SCALE", "LOG_BOUNDED");
	break;
    }

    if (ps == PARAM_SCALE_BOUNDLOG) {
	if (!pp.GetReal ("PARAM_SCALE_CMUAMIN", g_param_cmuamin) ||
	    !pp.GetReal ("PARAM_SCALE_CMUAMAX", g_param_cmuamax) ||
	    g_param_cmuamin < 0 || g_param_cmuamax < 0 ||
	    g_param_cmuamin >= g_param_cmuamax) {
	    cout << "\nCMUA parameter range [1/ps]:\n";
	    cout << "<min> <max> >> ";
	    cin >> g_param_cmuamin >> g_param_cmuamax;
	}
	if (!pp.GetReal ("PARAM_SCALE_CKAPPAMIN", g_param_ckappamin) ||
	    !pp.GetReal ("PARAM_SCALE_CKAPPAMAX", g_param_ckappamax) ||
	    g_param_ckappamin < 0 || g_param_ckappamax < 0 ||
	    g_param_ckappamin >= g_param_ckappamax) {
	    cout << "\nCKAPPA parameter range [mm^2/ps]:\n";
	    cout << "<min> <max> >> ";
	    cin >> g_param_ckappamin >> g_param_ckappamax;
	}
	pp.PutReal ("PARAM_SCALE_CMUAMIN", g_param_cmuamin);
	pp.PutReal ("PARAM_SCALE_CMUAMAX", g_param_cmuamax);
	pp.PutReal ("PARAM_SCALE_CKAPPAMIN", g_param_ckappamin);
	pp.PutReal ("PARAM_SCALE_CKAPPAMAX", g_param_ckappamax);
    }

    return ps;
}

// ============================================================================

int SelectImageFormat ()
{
    static char *fmtstr[4] = {"RAW", "NIM", "PGM", "PPM"};
    char cbuf[256];
    int fmt = 0;

    if (pp.GetString ("IMAGEFORMAT", cbuf)) {
	if      (!strcasecmp (cbuf, "RAW")) fmt = IMGFMT_RAW;
	else if (!strcasecmp (cbuf, "NIM")) fmt = IMGFMT_NIM;
	else if (!strcasecmp (cbuf, "PGM")) fmt = IMGFMT_PGM;
	else if (!strcasecmp (cbuf, "PPM")) fmt = IMGFMT_PPM;
    }
    while (fmt < 1 || fmt > 4) {
	cout << "\nOutput image format:\n";
	cout << "(1) Raw data\n";
	cout << "(2) NIM (nodal image\n";
	cout << "(3) PGM (portable grayscale map) - 2D only\n";
	cout << "(4) PPM (portable pixmap) - 2D only\n";
	cout << "[1|2|3|4] >> ";
	cin  >> fmt;
    }
    pp.PutString ("IMAGEFORMAT", fmtstr[fmt-1]);
    return fmt;
}

// ============================================================================

int SelectSDMode ()
{
    static char *modestr[4] = {"NONE", "DATA", "AVG_DATA", "AVG_DIFFDATA"};
    char cbuf[256];
    int mode = 0;

    if (pp.GetString ("DATA_SCALING", cbuf)) {
	if      (!strcasecmp (cbuf, "NONE"))         mode = 1;
	else if (!strcasecmp (cbuf, "DATA"))         mode = 2;
	else if (!strcasecmp (cbuf, "AVG_DATA"))     mode = 3;
	else if (!strcasecmp (cbuf, "AVG_DIFFDATA")) mode = 4;
    }
    while (mode < 1 || mode > 4) {
	cout << "\nData scaling method:\n";
	cout << "(1) none\n";
	cout << "(2) data\n";
	cout << "(3) data type average\n";
	cout << "(4) data type difference average\n";
	cout << "[1|2|3|4] >> ";
	cin >> mode;
    }
    pp.PutString ("DATA_SCALING", modestr[mode-1]);
    return mode;
}

// ============================================================================

void AssembleRealSystemMatrix (const Mesh &mesh, RGenericSparseMatrix &F,
    const RVector &cmua, const RVector &ckappa, const RVector &c2a)
{
    AddToSysMatrix (mesh, F,  &cmua,   ASSEMBLE_PFF);
    AddToSysMatrix (mesh, F,  &ckappa, ASSEMBLE_PDD);
    AddToSysMatrix (mesh, F,  &c2a,    ASSEMBLE_BNDPFF);
}

void AssembleComplexSystemMatrix (const Mesh &mesh, CGenericSparseMatrix &CF,
    const RVector &cmua, const RVector &ckappa, const RVector &c2a,
    const double omega)
{
    AddToSysMatrix (mesh, CF, &cmua,   ASSEMBLE_PFF);
    AddToSysMatrix (mesh, CF, &ckappa, ASSEMBLE_PDD);
    AddToSysMatrix (mesh, CF, &c2a,    ASSEMBLE_BNDPFF);
    AddToSysMatrix (mesh, CF, omega,   ASSEMBLE_iCFF);
}

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
		bool bfront = (z > 0 && (!mask || (mask[idx-stridez] >= 0)));
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
// instantiation
template void ImageGradient (const IVector &gdim, const RVector &gsize,
    const RVector &im, RVector *grad, const int *mask);
template void ImageGradient (const IVector &gdim, const RVector &gsize,
    const CVector &im, CVector *grad, const int *mask);

// calculate intensity field for source vector qvec
// using iterative solver
void CalcFieldE (const RCompRowMatrix &F, const RPreconditioner *precon,
    const RVector &qvec, RVector &x)
{
    double tol = g_lsolver_tol;
    IterativeSolve (F, qvec, x, tol, precon);
}

// calculate intensity field for source vector qvec
// using direct solver
void CalcFieldE (const RCompRowMatrix &FL, const RVector &Fd,
    const RVector &qvec, RVector &x)
{
    CholeskySolve (FL, Fd, qvec, x);
}

double OF_value (const RVector &data, const RVector &proj,
    const RVector &sd)
{
    int i, dim = data.Dim();
    double temp, sum = 0.0;
    for (i = 0; i < dim; i++) {
        temp = (data[i] - proj[i]) / sd[i];
	sum += temp*temp;
    }
    return sum;
}

void OF_gradient_add_Mod (const Raster &raster, const CFwdSolver &FWS,
    const RVector &data, const RVector &proj, const RVector &sd, CVector *dphi,
    CCompRowMatrix &mvec, RVector &grad)
{
    // Calculate gradient of objective function (Mod) using adjoint method
    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs, nQ = mesh.nQ;
    int glen = raster.GLen();
    int blen = raster.BLen();
    int dim  = raster.Dim();
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    CVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    CVector dgrad (blen);
    ofs = 0;  // data offset for Mod data

    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster.Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh.nQMref[q];
	RVector y (data, ofs, n);
	RVector s (sd, ofs, n);
	RVector ypm (proj, ofs, n);
	RVector b(n);
	b = (y-ypm)/s;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	Project_cplx (mesh, q, dphi[q], cproj);
	wqa = complex(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    CVector qs = mvec.Row(m);
	    double rp = cproj[idx].re;
	    double ip = cproj[idx].im;
	    double dn = 1.0/(rp*rp + ip*ip);
	    double term = -2.0 * b[idx] / (ype[idx]*s[idx]);
	    wqa += qs * complex (term*rp*dn, -term*ip*dn);
	    wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);
	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster.Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster.Map_GridToBasis (cdfield * cafield, dgrad);
	RVector grad_mua(grad, 0, blen);
	grad_mua -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster.Map_GridToBasis (gk, dgrad);
	RVector grad_kappa (grad, blen, blen);
	grad_kappa -= Re(dgrad);
	delete []cdfield_grad;
	delete []cafield_grad;
	ofs += n; // step to next source
    }
}

void OF_gradient_add_Arg (const Raster &raster, const CFwdSolver &FWS,
    const RVector &data, const RVector &proj, const RVector &sd, CVector *dphi,
    CCompRowMatrix &mvec, RVector &grad)
{
    // Calculate gradient of objective function (phase component)
    // using adjoint method

    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs, nQ = mesh.nQ;
    int glen = raster.GLen();
    int blen = raster.BLen();
    int dim  = raster.Dim();
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    CVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    CVector dgrad (blen);
    ofs = mesh.nQM;  // data offset for Arg data

    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster.Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh.nQMref[q];
	RVector y (data, ofs, n);
	RVector s (sd, ofs, n);
	RVector ypm (proj, ofs, n);
	RVector b(n);
	b = (y-ypm)/s;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	Project_cplx (mesh, q, dphi[q], cproj);
	wqa = complex(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    CVector qs = mvec.Row(m);
	    double rp = cproj[idx].re;
	    double ip = cproj[idx].im;
	    double dn = 1.0/(rp*rp + ip*ip);
	    double term = -2.0 * b[idx] / (ype[idx]*s[idx]);
	    wqa += qs * complex (term*rp*dn, -term*ip*dn);
	    wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);
	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster.Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster.Map_GridToBasis (cdfield * cafield, dgrad);
	RVector grad_mua(grad, 0, blen);
	grad_mua -= Im(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster.Map_GridToBasis (gk, dgrad);
	RVector grad_kappa (grad, blen, blen);
	grad_kappa -= Im(dgrad);
	delete []cdfield_grad;
	delete []cafield_grad;
	ofs += n; // step to next source
    }
}

void OF_gradient_add_ModArg (const Raster &raster, const CFwdSolver &FWS,
    const RVector &data, const RVector &proj, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad, bool logparam,
    const Solution *sol)
{
#ifdef COMPUTE_FLOPS
    ResetFlops();
#endif
    double tk = tic();
    // Calculate gradient of objective function (Mod) using adjoint method
    const QMMesh &mesh = *FWS.meshptr;
    int i, j, q, m, n, idx, ofs_mod, ofs_arg, nQ = mesh.nQ;
    int glen = raster.GLen();
    int slen = raster.SLen();
    int dim  = raster.Dim();
    double term;
    const IVector &gdim = raster.GDim();
    const RVector &gsize = raster.GSize();
    const int *elref = raster.Elref();
    CVector wqa (mesh.nlen());
    RVector wqb (mesh.nlen());
    CVector dgrad (slen);
    ofs_mod = 0;         // data offset for Mod data
    ofs_arg = mesh.nQM;  // data offset for Arg data
    RVector grad_cmua(grad, 0, slen);       // mua part of grad
    RVector grad_ckappa (grad, slen, slen); // kappa part of grad
    
    double tm_mesh2grid = 0.0;
    double tm_grid2sol = 0.0;
    double tm_gradient = 0.0;
    double tm_innerloop = 0.0;

    for (q = 0; q < mesh.nQ; q++) {

        // expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	tic();
	raster.Map_MeshToGrid (dphi[q], cdfield);
	tm_mesh2grid += toc();
	tic();
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);
	tm_gradient += toc();

        n = mesh.nQMref[q];

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
	Project_cplx (mesh, q, dphi[q], cproj);
	wqa = complex(0,0);
	wqb = 0.0;

	tic();
	for (m = idx = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
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
	tm_innerloop += toc();

	// adjoint field and gradient
	CVector wphia (mesh.nlen());
	FWS.CalcField (wqa, wphia);
	//WriteImage (Re(wphia), q, "wphia_re.nim");
	//WriteImage (Im(wphia), q, "wphia_im.nim");

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	tic();
	raster.Map_MeshToGrid (wphia, cafield);
	tm_mesh2grid += toc();
	tic();
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);
	tm_gradient += toc();

	// absorption contribution
	tic();
	raster.Map_GridToSol (cdfield * cafield, dgrad);
	tm_grid2sol += toc();
	grad_cmua -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	tic();
	raster.Map_GridToSol (gk, dgrad);
	tm_grid2sol += toc();
	grad_ckappa -= Re(dgrad);
	delete []cdfield_grad;
	delete []cafield_grad;
	ofs_mod += n; // step to next source
	ofs_arg += n;
    }

    if (logparam && sol) {  // log parameter scaling
	grad_cmua   *= exp (sol->GetParam (OT_CMUA));
	grad_ckappa *= exp (sol->GetParam (OT_CKAPPA));
    }
#ifdef DO_PROFILE
    cerr << "  Gradient timings" << endl;
    cerr << "    mesh to grid: " << tm_mesh2grid << endl;
    cerr << "    grid to sol:  " << tm_grid2sol << endl;
    cerr << "    img gradient: " << tm_gradient << endl;
    cerr << "    meas loop:    " << tm_innerloop << endl;
    cerr << "    total:        " << toc(tk) << endl;
#endif
#ifdef COMPUTE_FLOPS
    cerr << "  Flops" << endl;
    cerr << "    Additions:      " << FlopsAdd() << endl;
    cerr << "    Multiplications:" << FlopsMul() << endl;
#endif
}

void OF_gradient_add (const RVector &data, const RVector &proj,
		      const RVector &sd, const RDenseMatrix &J,
		      RVector &grad)
{
    // Calculate gradient of objective function, given full Jacobian matrix
    int i, n = J.nRows();

    xASSERT(data.Dim() == n, Invalid length of data vector);
    xASSERT(sd.Dim() == n, Invalid length of sd vector);
    xASSERT(proj.Dim() == n, Invalid length of projection vector);
    xASSERT(grad.Dim() == J.nCols(), Invalid length of gradient vector);

    for (i = 0; i < n; i++) {
        RVector pmdf = J.Row(i);
	const double y  = data[i];
	const double s  = sd[i];
	const double yp = proj[i];
	double backproject = (y-yp)/s;
	grad += pmdf * (-2.0*backproject);
    }
}

// norm of 2 vectors
double error_norm (const RVector &a, const RVector &b)
{
    double temp, sum = 0.0;
    for (int i = 0; i < a.Dim(); i++) {
        temp = a[i]-b[i];
	sum += temp*temp;
    }
    return sum;
}

// norm of 2 vectors scaled by sd
double error_norm (const RVector &a, const RVector &b, const RVector &sd)
{
    double temp, sum = 0.0;
    for (int i = 0; i < a.Dim(); i++) {
        temp = (a[i]-b[i])/sd[i];
	sum += temp*temp;
    }
    return sum;
}

#ifdef UNDEF
// returns the gradient of the objective function
void GetGradient (const QMMesh &mesh, CFwdSolver &FWS, const CVector *dphi,
    const CVector *aphi, const RVector &data, const RVector &proj,
    const RVector &sd, RVector &grad)
{
    int q, m, idx1, idx2;
    int n = mesh.nlen(), ofs = mesh.nQM;
    double idenom, y, yp, s;
    RVector res(n);

    for (q = idx1 = 0, idx2 = ofs; q < mesh.nQ; q++) {
        for (m = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q,m)) continue;
	    CVector pmdf = PMDF_mua (dphi[q], aphi[m]);
	    switch (FWS.datatype) {
	    case MEAS_FRE_FIM:
		y   = data[idx1];
		yp  = proj[idx1];
		s   = sd[idx1];
		grad += Re(pmdf) * (-2.0*(y-yp)/(s*s));
		y   = data[idx2];
		yp  = proj[idx2];
		s   = sd[idx2];
		grad += Im(pmdf) * (-2.0*(y-yp)/(s*s));
		break;
	    case MEAS_FMOD_FARG:
	        break;
	    }
	    idx1++, idx2++;
	}
    }
}
#endif

void ReadDataFile (char *fname, RVector &data)
{
    int i, n = data.Dim();
    char c;
    ifstream ifs (fname);
    do {
        ifs >> c;
    } while (ifs.good() && c != '[');
    if (!ifs.good())
	xERROR(Data file not found or invalid);
    for (i = 0; i < n; i++)
        ifs >> data[i];
}

void WriteJacobian (const RMatrix *J, const Raster &raster,
    const QMMesh &mesh)
{
    int q, m, i, j, idx, xofs, yofs, dim = raster.Dim();
    double vmin, vmax, scal;

    if (dim != 2) { // for now, just dump the whole thing
	for (i = 0; i < J->nRows(); i++) {
	    char cbuf[256];
	    sprintf (cbuf, "pmdf_%03d.dat", i);
	    ofstream ofs (cbuf);
	    ofs << J->Row(i) << endl;
	}
	return; // this function only handles 2D problems
    }

    IVector bdim = raster.BDim();
    int blen = raster.BLen();
    int slen = raster.SLen();
    int imgw = mesh.nM*bdim[0];
    int imgh = mesh.nQ*bdim[1];
    IVector idim(2);
    idim[0] = imgw, idim[1] = imgh;
    RVector img(imgw*imgh);
    RVector Ji(J->nCols());

    // 1: real, mua
    for (q = idx = 0; q < mesh.nQ; q++) {
        yofs = q*bdim[1];
        for (m = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    xofs = m*bdim[0];
	    Ji = J->Row(idx);
	    RVector sol(Ji, 0, slen);
	    RVector p(blen);
	    raster.Map_SolToBasis (sol, p);
	    ImageScale (p, vmin, vmax);
	    scal = 1.0/(vmax-vmin);
	    for (j = 0; j < bdim[1]; j++)
	        for (i = 0; i < bdim[0]; i++)
		    img[(yofs+j)*imgw + (xofs+i)] = (p[j*bdim[0]+i]-vmin)*scal;
	    idx++;
	}
    }
    vmin = 0.0, vmax = 1.0;
    WritePPM (img, idim, &vmin, &vmax, "j_re_mua.ppm");
    cout << "Jacobian (real, mua) written to j_re_mua.ppm" << endl;

    // 2: real, kappa
    if (J->nCols() >= slen*2) { // J contains kappa
        for (q = idx = 0; q < mesh.nQ; q++) {
	    yofs = q*bdim[1];
	    for (m = 0; m < mesh.nM; m++) {
	        if (!mesh.Connected (q, m)) continue;
		xofs = m*bdim[0];
		Ji = J->Row(idx);
		RVector sol(Ji, slen, slen);
		RVector p(blen);
		raster.Map_SolToBasis (sol, p);
		ImageScale (p, vmin, vmax);
		scal = 1.0/(vmax-vmin);
		for (j = 0; j < bdim[1]; j++)
		    for (i = 0; i < bdim[0]; i++)
		        img[(yofs+j)*imgw + (xofs+i)] =
			    (p[j*bdim[0]+i]-vmin)*scal;
		idx++;
	    }
	}
	vmin = 0.0, vmax = 1.0;
	WritePPM (img, idim, &vmin, &vmax, "j_re_kappa.ppm");
	cout << "Jacobian (real, kappa) written to j_re_kappa.ppm" << endl;
    }

    if (J->nRows() >= mesh.nQM*2) {
        // 3: imag, mua
        for (q = 0, idx = mesh.nQM; q < mesh.nQ; q++) {
	    yofs = q*bdim[1];
	    for (m = 0; m < mesh.nM; m++) {
	        if (!mesh.Connected (q, m)) continue;
		xofs = m*bdim[0];
		Ji = J->Row(idx);
		RVector sol(Ji, 0, slen);
		RVector p(blen);
		raster.Map_SolToBasis (sol, p);
		ImageScale (p, vmin, vmax);
		scal = 1.0/(vmax-vmin);
		for (j = 0; j < bdim[1]; j++)
		    for (i = 0; i < bdim[0]; i++)
		        img[(yofs+j)*imgw + (xofs+i)] =
			    (p[j*bdim[0]+i]-vmin)*scal;
		idx++;
	    }
	}
	vmin = 0.0, vmax = 1.0;
	WritePPM (img, idim, &vmin, &vmax, "j_im_mua.ppm");
	cout << "Jacobian (imag, mua) written to j_im_mua.ppm" << endl;

	// 4: real, kappa
	if (J->nCols() >= slen*2) { // J contains kappa
	    for (q = 0, idx = mesh.nQM; q < mesh.nQ; q++) {
	        yofs = q*bdim[1];
		for (m = 0; m < mesh.nM; m++) {
		    if (!mesh.Connected (q, m)) continue;
		    xofs = m*bdim[0];
		    Ji = J->Row(idx);
		    RVector sol(Ji, slen, slen);
		    RVector p(blen);
		    raster.Map_SolToBasis (sol, p);
		    ImageScale (p, vmin, vmax);
		    scal = 1.0/(vmax-vmin);
		    for (j = 0; j < bdim[1]; j++)
		        for (i = 0; i < bdim[0]; i++)
			    img[(yofs+j)*imgw + (xofs+i)] =
			        (p[j*bdim[0]+i]-vmin)*scal;
		    idx++;
		}
	    }
	    vmin = 0.0, vmax = 1.0;
	    WritePPM (img, idim, &vmin, &vmax, "j_im_kappa.ppm");
	    cout << "Jacobian (imag, kappa) written to j_im_kappa.ppm" << endl;
	}
    }
}

void ColumnScale (RDenseMatrix &A, const RVector &cscale, int m0, int m1)
{
    // Scales the columns of matrix A given a vector of scaling factors.
    // A sub-range of columns can be specified by the limits m0 and m1
    // Dimension of cscale must be >= m1-m0

    int i, j;
    int n = A.nRows();
    int m = A.nCols();

    if (m0 < 0) m0 = 0;
    if (m1 < 0 || m1 > m) m1 = m;

    for (j = m0; j < m1; j++)
        for (i = 0; i < n; i++)
	    A(i,j) *= cscale[j-m0];
}

#ifdef UNDEF
void RowScale (RDenseMatrix &A, const RVector &rscale, int n0, int n1)
{
    // Scales the rows of matrix A given a vector of scaling factors.
    // A sub-range of rows can be specified by the limits n0 and n1
    // Dimension of rscale must be >= n1-n0


    int i, j;
    int n = A.nRows();
    int m = A.nCols();

    if (n0 < 0) n0 = 0;
    if (n1 < 0 || n1 > n) n1 = n;

    for (i = n0; i < n1; i++)
        for (j = 0; j < m; j++)
	    A(i,j) *= rscale[i-n0];
}
#endif

void BiScale (RSymMatrix &A, const RVector &cscale)
{
    // Scales the diagonal matrix A given a vector of scaling factors.

    int i, j;
    int n = A.nRows();
    int m = A.nCols();
    /* m should be equal to n ! */
    for (j = 0; j < m; j++)
      for (i = 0; i <=j ; i++) { /* its symmetric ! */
	//	cout << "A(" << i << "," << j << ") was " << A(i,j);
	    A(i,j) *= cscale[j]*cscale[i];
	//	    cout << " -> " << A(i,j) << endl;
      }
}

#ifdef UNDEF

// This is an implementation of preconditioned nonlinear CG
// from the Shewchuk paper, B5 (pg.53) but without Secant method

void CGsolve (CFwdSolver &FWS, const Raster &raster, ObjectiveFunction &OF,
    int itmax, const RVector &data, const RVector &sd, Solution &msol,
    Solution &bsol, const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    double omega, bool logparam, double ftol)
{
    const int reset_count = PCG_RESET_INTERVAL;

    // initialisations
    int i, j;
    const QMMesh &mesh = raster.mesh();
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh.nlen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int slen = raster.SLen();
    int n    = slen*2;
    bool pvalid;
    double delta_new, delta_old, delta_mid, delta_0, delta_d, beta;
    double of, of_value, of_prior, fmin, alpha = -1.0;
    double gamma = 1.0;
    RVector r(n), r0(n), s(n), d(n), M(n);
    RVector proj(ndat);
    RDenseMatrix J;
    RSymMatrix JTJ;
    RCompRowMatrix JTJ_L;
    RVector JTJ_d;
    RVector x0 (bsol.GetActiveParams()); // initial state
    RVector mu0(n/2), kap0(n/2);
    //    x0.Relink(mu0,0,n/2); 
    //    x0.Relink(kap0,n/2,n/2); 
    for(i = 0; i < n/2; i++) kap0[i] = x0[i+n/2];
    //cout << "Initial kappa\n" << kap0 << endl;

    CVector *dphi = new CVector[mesh.nQ];
    for (i = 0; i < mesh.nQ; i++) dphi[i].New (nlen);
    CVector *aphi = new CVector[mesh.nM];
    for (i = 0; i < mesh.nM; i++) aphi[i].New (nlen);

    // Start of Shewchuk implementation

    int i_count = 0; // iteration counter
    int k_count = 0; // reset counter

    cout << "  Generating fields and gradient" << endl;
    FWS.Reset (msol, omega);
    FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
    FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
    ProjectAll (mesh, FWS, mvec, dphi, proj);

#ifdef OUTPUT_INITIAL_PROJECTION
    ofstream ofs1 ("init_fmod.fem");
    RVector proj1(proj, 0, ndat/2);
    ofs1 << proj1 << endl;
    ofstream ofs2 ("init_farg.fem");
    RVector proj2(proj, ndat/2, ndat/2);
    ofs2 << proj2 << endl;
#endif

    // r = -f'(x)
    OF.get_gradient (raster, FWS, proj, dphi, mvec, logparam, &bsol, r);
    r0 = r;

    if (bWriteGrad) {
        RVector r_mua (r, 0, slen);
	RVector r_mus (r, slen, slen);
	RVector rim(blen);
	raster.Map_SolToBasis (r_mua, rim);
	WriteImage (rim, 0, "grad_mua.raw");
	raster.Map_SolToBasis (r_mus, rim);
	WriteImage (rim, 0, "grad_mus.raw");
    }

    r = -r;
    of_value = OF.get_posterior (&proj);
    of_prior = OF.get_prior (&bsol);
    of = of_value + of_prior;

    if (g_pcg_precon != PCG_PRECON_NONE) {
        // calculate preconditioner M
        cout << "  Calculating Jacobian ..." << endl;
	J.New (ndat, n);
	GenerateJacobian (raster, dphi, aphi, FWS.datatype, J, sd, kap0);
	if (logparam) // Rescale Jacobian with parameters
	    ColumnScale (J, exp(bsol.GetActiveParams()));
	if (bWriteJ) WriteJacobian (&J, raster, mesh);
	switch (g_pcg_precon) {
	case PCG_PRECON_FULLJTJ:
	    // Full Hessian preconditioner setup
	    cout << "  Generating dense JTJ" << endl;
	    ATA_dense (raster, J, JTJ);
	    cout << "  Using preconditioner FULLJTJ" << endl;
	    break;
	case PCG_PRECON_SPARSEJTJ:
	    // Incomplete CH preconditioner setup
	    cout << "  Generating sparse JTJ" << endl;
	    ATA_sparse (raster, J, JTJ_L, JTJ_d);
	    cout << "  Using preconditioner SPARSEJTJ" << endl;
	    break;
	case PCG_PRECON_DIAGJTJ:
	    // Diagonal Hessian preconditioner setup
	    cout << "  Calculating diagonal of JTJ" << endl;
	    ATA_diag (J, M);
	    cout << "    Range " << vmin (M) << " to " << vmax (M) << endl;
#ifdef DJTJ_LIMIT
	    M.Clip (DJTJ_LIMIT, 1e50);
	    cout << "    Cutoff at " << DJTJ_LIMIT << endl;
#endif // DJTJ_LIMIT
	    cout << "  Using preconditioner DIAGJTJ" << endl;
	    break;
	}
	cerr << "  Precon reset interval: " << reset_count << endl;
    } else {
	cout << "  Using preconditioner NONE" << endl;
    }

    times (&tme);
    cout << "Iteration 0  CPU "
	 << (double)(tme.tms_utime-clock0)/(double)HZ 
	 << "  OF " << of
	 << "  (prior " << of_prior << " )" << endl;

    // apply preconditioner
    switch (g_pcg_precon) {
    case PCG_PRECON_NONE:
        s = r;
	break;
    case PCG_PRECON_DIAGJTJ:
        s = r/M; // element-wise division
	break;
    case PCG_PRECON_SPARSEJTJ:
        CholeskySolve (JTJ_L, JTJ_d, r, s);
	break;
    case PCG_PRECON_FULLJTJ:
        s = CHsubst (JTJ, r);
	break;
    }

    d = s;
    delta_new = r & d;                 // r^t M^-1 r
    delta_0 = delta_new;

    while (i_count < itmax && delta_new > ftol*ftol*delta_0) {
        delta_d = d & d;

	if (alpha < 0.0) { // initialise step length
	    alpha = of / l2norm (d);
	    cout << "  Initial step length reset to " << alpha << endl;
	}
	// line search. this replaces the Secant method of the Shewchuk code
	if (LineSearch (FWS, raster, OF, qvec, mvec, data, sd, omega, bsol,
			d, of, fmin, alpha, msol, proj, pvalid, logparam)) {

	    // update: x = x + alpha d
	    bsol += d*alpha;
	    raster.Map_SolToMesh (logparam ? exp(bsol) : bsol, msol, true);
	    msol.UnscaleParams ();
	    if (bOutputNIM) {
	        msol.WriteImg_mua (i_count+1, "recon_mua.nim");
		msol.WriteImg_mus (i_count+1, "recon_mus.nim");
	    }
	    if (bOutputRaw) {
	        Solution rsol(blen);
		raster.Map_SolToBasis (logparam ? exp(bsol):bsol, rsol, true);
		rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
		rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
	    }

	} else {
  	    cout << "  ** Line search failed. Resetting." << endl;
	    k_count = reset_count-1; // force reset
	    i_count--; // don't increment iteration counter
	}

	// r = -f'(x)
	cout << "  Generating fields and gradient" << endl;
	if (!pvalid) { // need to re-generate fields and projections
	    FWS.Reset (msol, omega);
	    for (i = 0; i < mesh.nQ; i++) dphi[i].Clear();
	    FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
	    ProjectAll (mesh, FWS, mvec, dphi, proj);
	}
	FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
	OF.get_gradient (raster, FWS, proj, dphi, mvec, logparam, &bsol, r);
	if (bWriteGrad) {
	    RVector r_mua (r, 0, slen);
	    RVector r_mus (r, slen, slen);
	    RVector rim(blen);
	    raster.Map_SolToBasis (r_mua, rim);
	    WriteImage (rim, i_count+1, "grad_mua.raw");
	    raster.Map_SolToBasis (r_mus, rim);
	    WriteImage (rim, i_count+1, "grad_mus.raw");
	}

#ifdef RESCALE_HESSIAN
	RVector x1(bsol.GetActiveParams());
	RVector S(x1-x0);
	RVector Y(r-r0);
	gamma = (Y&S) / (Y&Y);
	cout << "Hessian scale " << gamma << endl;
	x0 = x1;
	r0 = r;
#endif

	r = -r;
	of_value = OF.get_posterior (&proj);
	of_prior = OF.get_prior (&bsol);
	of = of_value + of_prior;
	delta_old = delta_new;
	delta_mid = r & s;

	k_count++;

	if (g_pcg_precon != PCG_PRECON_NONE && k_count == reset_count) {
	    // re-calculate preconditioner and reset CG
	    cout << "  Calculating Jacobian ..." << endl;
	    GenerateJacobian (raster, dphi, aphi, FWS.datatype, J, sd, kap0);
	    if (logparam) // Rescale Jacobian with parameters
	        ColumnScale (J, exp(bsol.GetActiveParams()));
	    switch (g_pcg_precon) {
	    case PCG_PRECON_FULLJTJ:
	        cout << "  Generating dense JTJ ..." << endl;
		ATA_dense (raster, J, JTJ);
		break;
	    case PCG_PRECON_SPARSEJTJ:
  	        cout << "  Generating sparse JTJ ..." << endl;
		ATA_sparse (raster, J, JTJ_L, JTJ_d);
		break;
	    case PCG_PRECON_DIAGJTJ:
	        cout << "  Calculating diagonal of JTJ ..." << endl;
		ATA_diag (J, M);
		cout << "    Range " << vmin (M) << " to " << vmax (M) << endl;
#ifdef DJTJ_LIMIT
		M.Clip (DJTJ_LIMIT, 1e50);
		cout << "    Cutoff at " << DJTJ_LIMIT << endl;
#endif // DJTJ_LIMIT
	    }
	}

	// apply preconditioner
	switch (g_pcg_precon) {
	case PCG_PRECON_NONE:
	    s = r;
	    break;
	case PCG_PRECON_DIAGJTJ:
	    s = (r/M)*gamma; // element-wise division
	    break;
	case PCG_PRECON_SPARSEJTJ:
	    CholeskySolve (JTJ_L, JTJ_d, r, s);
	    break;
	case PCG_PRECON_FULLJTJ:
	    s = CHsubst (JTJ, r);
	    break;
	}

	delta_new = r & s;
	beta = (delta_new - delta_mid) / delta_old;
	if (k_count == reset_count || beta <= 0.0) {
	    d = s;
	    k_count = 0;
	} else {
	    d = s + d * beta;
	}
	i_count++;
	times (&tme);
	cout << "Iteration " << i_count << "  CPU "
	     << (double)(tme.tms_utime-clock0)/(double)HZ
	     << "  OF " << of
	     << "  (prior " << of_prior << " )" << endl;

#ifdef DO_PROFILE
	cout << "  Solver time: " << solver_time << endl;
#endif
    }
}
#endif // UNDEF

#ifdef UNDEF // old BFGS scheme
void BFGSsolve (CFwdSolver &FWS, const Raster &raster, int itmax,
    const RVector &data, const RVector &sd, Solution &msol, Solution &bsol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, bool logparam,
    double ftol)
{
    int i, j;
    const QMMesh &mesh = raster.mesh();
    int ndat = data.Dim();
    int nlen = mesh.nlen();
    int slen = raster.SLen();
    int n    = slen*2;
    double of, fmin, alpha;
    RVector proj(ndat);
    RVector g(n);
    RVector p(n);
    RVector s(n), y(n);

    CVector *dphi = new CVector[mesh.nQ];
    for (i = 0; i < mesh.nQ; i++) dphi[i].New (nlen);

    int i_count = 0;

    // find initial guess for minimizer
    RVector x0 (bsol.GetActiveParams());

    // find initial objective function
    FWS.Reset (msol, omega);
    FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
    ProjectAll (mesh, FWS, mvec, dphi, proj);
    of = OF_value (data, proj, sd);

    times (&tme);
    cout << "Iteration 0  CPU "
	 << (double)(tme.tms_utime-clock0)/(double)HZ 
	 << "  OF " << of << endl;

#if BFGS_HESSIAN == BFGS_HESSIAN_IDENTITY
    // initial guess for Hessian: identity
    cout << "  Allocating Hessian:  " << n << "x" << n << " ("
	 << (4*n*(n+1))/1048576 << " Mbytes)" << endl;
    RSymMatrix H(n,n);
    for (i = 0; i < n; i++) H(i,i) = 1.0;
    cout << "  Initial Hessian: Identity" << endl;
#elif BFGS_HESSIAN == BFGS_HESSIAN_DIAG
    cout << "  Allocating Hessian:  " << n << "x" << n << " ("
	 << (4*n*(n+1))/1048576 << " Mbytes)" << endl;
    RSymMatrix H(n,n);
    CVector *aphi = new CVector[mesh.nM];
    for (i = 0; i < mesh.nM; i++) aphi[i].New (nlen);
    FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
    cout << "  Allocating Jacobian: " << ndat << "x" << n << " ("
	 << (8*ndat*n)/1048576 << " Mbytes)" << endl;
    RDenseMatrix J (ndat, n);
    RVector M(n);
    RVector kap0(n/2);
    for(i = 0; i < n/2; i++) kap0[i] = x0[i+n/2];
    GenerateJacobian (raster, dphi, aphi, FWS.datatype, J, sd, kap0);
    if (logparam) // Rescale Jacobian with parameters
        ColumnScale (J, exp(bsol.GetActiveParams()));
    ATA_diag (J, M);
    J.New(0,0);
    for (i = 0; i < n; i++) H(i,i) = M[i];
    cout << "  Initial Hessian: diag(JTJ)" << endl;
#elif BFGS_HESSIAN == BFGS_MEMORYLESS
    // nothing to do here
#else
    cerr << "  ** Invalid Hessian initialisation. Aborting." << endl;
    exit (1);
#endif

    // initial gradient
    OF_gradient_add_ModArg (raster, FWS, data, proj, sd, dphi, mvec, g,
        logparam, &bsol);

    // initial guess of step length
    alpha = of / l2norm (g);

    // begin Quasi-Newton iterations
    while (i_count < itmax) {

#if BFGS_HESSIAN == BFGS_MEMORYLESS
        if (!i_count) {
	    p = -g;
	} else {
	    double a1 = 1.0 / (s & y);
	    double a2 = (y & y) * a1 + 1.0;
	    double a3 = a1 * a2;
	    RVector Hi(n);
	    for (i = 0; i < n; i++) {
	        for (j = 0; j < n; j++)
		    Hi[j] = s[i]*s[j]*a3 - (y[i]*s[j] + s[i]*y[j])*a1;
		Hi[i] += 1.0;
		p[i] = -(Hi & g);
	    }
	}
#else
        RSymMatrix HF(H); // terrible waste of memory
        CHdecomp (HF, true);
	p = -CHsubst (HF, g);
#endif

	// line search
	LineSearch (FWS, raster, OF, qvec, mvec, data, sd, omega, bsol, p, of,
		    fmin, alpha, msol, proj, pvalid, logparam);

	// update approximate solution
	bsol += p*alpha;
	raster.Map_SolToMesh (logparam ? exp(bsol) : bsol, msol, true);

	if (bOutputNIM) {
	    msol.WriteImg_mua (i_count+1, "recon_mua.nim");
	    msol.WriteImg_mus (i_count+1, "recon_mus.nim");
	}
	if (bOutputRaw) {
	    Solution rsol(blen);
	    raster.Map_SolToBasis (logparam ? exp(bsol):bsol, rsol, true);
	    rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
	    rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
	}

	// new gradient
	FWS.Reset (msol, omega);
	FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
	ProjectAll (mesh, FWS, mvec, dphi, proj);
	RVector g1(n);
	OF_gradient_add_ModArg (raster, FWS, data, proj, sd, dphi, mvec, g1,
				logparam, &bsol);
	of = OF_value (data, proj, sd);

	RVector x1(bsol.GetActiveParams());
	s = x1-x0;
	y = g1-g;

#if BFGS_HESSIAN != BFGS_MEMORYLESS
	// update Hessian
	RVector h1(n), h2(n);
	H.Ax(s, h1);
	double den1 = 1.0/(h1 & s);
	double den2 = 1.0/(y & s);
	for (i = 0; i < n; i++) {
	    double h = 0.0;
	    for (j = 0; j < n; j++) h += s[j] * H(j,i);
	}
	for (i = 0; i < n; i++)
  	    for (j = 0; j <= i; j++) // symmetric, so only need lower triangle
	        H(i,j) += -h1[i]*h2[j]*den1 + y[i]*y[j]*den2;
#endif
	i_count++;
	x0 = x1;
	g  = g1;
	times (&tme);
	cout << "Iteration " << i_count << "  CPU "
	     << (double)(tme.tms_utime-clock0)/(double)HZ
	     << "  OF " << of << endl;
    }
    
}
#endif

void BFGSsolve (CFwdSolver &FWS, const Raster &raster, ObjectiveFunction &OF,
    int itmax, const RVector &data, const RVector &sd, Solution &msol,
    Solution &bsol, const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    double omega, bool logparam, double ftol)
{
    int i, j, k;
    const QMMesh *mesh = FWS.MeshPtr();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int n    = slen*2;
    bool pvalid;
    double of, fmin, alpha, sum, vt;
    RVector proj(ndat);
    RVector g(n);
    RVector p(n);
    RVector s(n), y(n);

    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);

    int i_count = 0;

    // find initial guess for minimizer
    RVector x0 (bsol.GetActiveParams());

    // find initial objective function
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    proj = FWS.ProjectAll_real (mvec, dphi);
    of = OF_value (data, proj, sd);

    cout << "Iteration 0  CPU "
	 << toc(clock0) << "  OF " << of << endl;

    cout << "  Allocating inverse Hessian:  " << n << "x" << n << " ("
	 << (4*n*(n+1))/1048576 << " Mbytes)" << endl;
    RSymMatrix HI(n,n);

#if BFGS_HESSIAN == BFGS_HESSIAN_IDENTITY
    // initial guess for Hessian: identity
    for (i = 0; i < n; i++) HI(i,i) = 1.0;
    cout << "  Initial Hessian: Identity" << endl;
#elif BFGS_HESSIAN == BFGS_HESSIAN_DIAG
    CVector *aphi = new CVector[mesh->nM];
    for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);
    FWS.CalcFields (*mesh, mesh->nM, mvec, aphi);
    cout << "  Allocating Jacobian: " << ndat << "x" << n << " ("
	 << (8*ndat*n)/1048576 << " Mbytes)" << endl;
    RDenseMatrix J (ndat, n);
    RVector M(n);
    RVector kap0(n/2);
    for(i = 0; i < n/2; i++) kap0[i] = x0[i+n/2];
    GenerateJacobian (raster, dphi, aphi, FWS.datatype, J, sd, kap0);
    if (logparam) // Rescale Jacobian with parameters
        ColumnScale (J, exp(bsol.GetActiveParams()));
    ATA_diag (J, M);
    J.New(0,0);
    for (i = 0; i < n; i++) HI(i,i) = 1.0/M[i];
    cout << "  Initial Hessian: diag(JTJ)" << endl;
#else
    cerr << "  ** Invalid Hessian initialisation. Aborting." << endl;
    exit (1);
#endif

    // initial gradient
    OF_gradient_add_ModArg (raster, FWS, data, proj, sd, dphi, mvec, g,
        logparam, &bsol);

    // initial guess of step length
    alpha = of / l2norm (g);

    // begin Quasi-Newton iterations
    while (i_count < itmax) {

        HI.Ax (-g, p);

	// line search
	LineSearch (FWS, raster, OF, qvec, mvec, data, sd, omega, bsol, p, of,
		    fmin, alpha, msol, proj, pvalid, logparam);

	// update approximate solution
	bsol += p*alpha;
	raster.Map_SolToMesh (logparam ? exp(bsol) : bsol, msol, true);

	switch (g_imgfmt) {
	case IMGFMT_NIM:
	    msol.WriteImg_mua (i_count+1, "recon_mua.nim");
	    msol.WriteImg_mus (i_count+1, "recon_mus.nim");
	    break;
	case IMGFMT_RAW:
	    Solution rsol(OT_NPARAM, blen);
	    raster.Map_SolToBasis (logparam ? exp(bsol):bsol, rsol, true);
	    rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
	    rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
	    break;
	}

	// new gradient
	FWS.Reset (msol, omega);
	FWS.CalcFields (qvec, dphi);
	proj = FWS.ProjectAll_real (mvec, dphi);
	RVector g1(n);
	OF_gradient_add_ModArg (raster, FWS, data, proj, sd, dphi, mvec, g1,
				logparam, &bsol);
	of = OF_value (data, proj, sd);

	RVector x1(bsol.GetActiveParams());
	s = x1-x0;
	y = g1-g;

	// update inverse Hessian
	// Byrd et. al. 1996 (p.2)

	cout << "  Updating inverse Hessian" << endl;
	ofstream ofs("HI.dat"); // temporary storage of new inverse Hessian
	RVector Vj(n), HVj(n), VTHVj(n);
	double rho = 1.0 / (y & s);
	for (j = 0; j < n; j++) {
  	    cout << j << endl;
  	    // build j-th column of V = y s^T
	    for (i = 0; i < n; i++)
	        Vj[i] = -rho*y[i]*s[j];
	    Vj[j] += 1.0;
	    // build j-th column of HV
	    for (i = 0; i < n; i++) {
	        for (sum = 0.0, k = 0; k < n; k++)
		    sum += HI(i,k)*Vj[k];
		HVj[i] = sum;
	    }
	    // build j-th column of V^T HV
	    for (i = 0; i < n; i++) {
	        for (sum = 0.0, k = 0; k < n; k++) {
		    vt = -y[k]*s[i];
		    if (i == k) vt += 1.0;
		    sum += vt*HVj[k];
		}
		VTHVj[i] = sum;
	    }
	    // add rho * s s^T
	    for (i = 0; i < n; i++)
	        VTHVj[i] += rho * s[i]*s[j];
	    // write column j to disk
	    ofs << VTHVj << endl;
	}
	ofs.close();
	// read new inverse Hessian back into H
	ifstream ifs("HI.dat");
	for (j = 0; j < n; j++) {
	    ifs >> VTHVj;
	    for (i = j; i < n; i++) // only need lower triangle
	        HI(i,j) = VTHVj[i];
	}
	ifs.close();
	//RVector h1(n), h2(n);
	//H.Ax(s, h1);
	//double den1 = 1.0/(h1 & s);
	//double den2 = 1.0/(y & s);
	//for (i = 0; i < n; i++) {
	//    double h = 0.0;
	//    for (j = 0; j < n; j++) h += s[j] * H(j,i);
	//}
	//for (i = 0; i < n; i++)
  	//    for (j = 0; j <= i; j++) // symmetric, so only need lower triangle
	//        H(i,j) += -h1[i]*h2[j]*den1 + y[i]*y[j]*den2;

	i_count++;
	x0 = x1;
	g  = g1;
	cout << "Iteration " << i_count << "  CPU "
	     << toc(clock0) << "  OF " << of << endl;
    }
    
}

// ==========================================================================
// This data structure defines the Hessian implicitly

struct HESS_DATA {
    const RDenseMatrix *J;           // Jacobian
    const RCompRowMatrix *Lap;       // Laplacian for prior
    const RVector *dg;               // diagonal additive component
    const double *tau, *tau1, *tau2; // scaling factors
    const double *lambda;            // diagonal scaling factor
};

// ==========================================================================

RSymMatrix Hess_full (const HESS_DATA &data)
{
    // Computes the explicit Hessian from its components, defined as
    // (J^T J + P + lambda D)
    // where J and P are full, D is the diagonal of J^T J + P, and
    // lambda is a scalar

    // unpack the data
    const RDenseMatrix &J     = *data.J;
    const RCompRowMatrix &Lap = *data.Lap;
    const RVector &dg         = *data.dg;
    const double tau          = *data.tau;
    const double tau1         = *data.tau1;
    const double tau2         = *data.tau2;
    const double lambda       = (data.lambda ? *data.lambda : 0.0);
    int m = J.nRows(), n = J.nCols(), n2 = Lap.nRows();

    RSymMatrix hess (ATA (J));

    // form prior P
#ifdef LAPLACIAN_PRIOR
    int i, j, k, nz, *colidx = new int[n2];
    double *val = new double[n2];
    for (i = 0; i < n2; i++) {
        nz = Lap.SparseRow (i, colidx, val);
	for (k = 0; k < nz && (j=colidx[k]) <= i; k++) {
	    hess(i,j)       += val[k] * tau1;
	    hess(i+n2,j+n2) += val[k] * tau2;
	}
    }
    delete []colidx;
    delete []val;
#else
    hess.AddDiag (tau/sqr(dg));
#endif
    hess.AddDiag (hess.Diag() * lambda);
    return hess;
}

// ==========================================================================

RVector Hess_diag (const HESS_DATA &data)
{
    // Computes the diagonal of the Hessian given implicitly as
    // (J^T J + P + lambda D)
    // where J and P are full, D is the diagonal of J^T J + P, and
    // lambda is a scalar

    // unpack the data
    const RDenseMatrix &J     = *data.J;
    const RCompRowMatrix &Lap = *data.Lap;
    const RVector &dg         = *data.dg;
    const double tau          = *data.tau;
    const double tau1         = *data.tau1;
    const double tau2         = *data.tau2;
    const double lambda       = (data.lambda ? *data.lambda : 0.0);
    int m = J.nRows(), n = J.nCols(), n2 = Lap.nRows();

    RVector diag(n);
    diag = ATA_diag (J);

    // form prior P
#ifdef LAPLACIAN_PRIOR
    // in this version, D is the Laplacian
    int i;
    RVector dd = Lap.Diag();
    for (i = 0; i < n2; i++) {
        diag[i]    += dd[i] * tau1;
	diag[i+n2] += dd[i] * tau2;
    }
#else
    // In this version, D is a diagonal prior
    diag += tau/sqr(dg) * y;
#endif
    diag *= 1+lambda;
    return diag;
}

// ==========================================================================

RVector JTJx_clbk (const RVector &x, void *context)
{
    // Computes (J^T J + P + lambda D) x
    // where J and P are full, D is the diagonal of J^T J + P, and
    // lambda is a scalar

    // unpack the data
    HESS_DATA *data           = (HESS_DATA*)context;
    const RDenseMatrix &J     = *data->J;
    const RCompRowMatrix &Lap = *data->Lap;
    const RVector &dg         = *data->dg;
    const double tau          = *data->tau;
    const double tau1         = *data->tau1;
    const double tau2         = *data->tau2;
    const double lambda       = (data->lambda ? *data->lambda : 0.0);
    int m = J.nRows(), n = J.nCols(), n2 = Lap.nRows();

    RVector Px(n);

    // form Px
#ifdef LAPLACIAN_PRIOR
    int i, j, nz, *colidx = new int[n2];
    double *val = new double[n2];
    for (i = 0; i < n2; i++) {
        nz = Lap.SparseRow (i, colidx, val);
	for (j = 0; j < nz; j++) {
	    Px[i]    += val[j] * tau1 * x[colidx[j]];
	    Px[i+n2] += val[j] * tau2 * x[colidx[j]+n2];
	}
    }
    delete []colidx;
    delete []val;
#else
    Px = tau/sqr(dg) * x;
#endif
    return ATx (J, Ax (J, x)) + Px + Hess_diag(*data) * x;
}

// ==========================================================================

void LMsolve (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    ObjectiveFunction &OF, int itmax, const RVector &data, const RVector &sd,
    Solution &msol, Solution &bsol, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, double omega, bool logparam, double ftol)
{
    bool Use_precon = true;

    // initialisations
    int i, inr, nrmax = 40;
    const QMMesh *mesh = FWS.MeshPtr();
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int n    = slen*2;
    double of, fmin, errstart, err0, err00, err1;
    double lambda = 10, alpha = -1.0, tau = 1e-5; //5e-3;
    double tau1, tau2;

    RVector r(n), s(n), d(n);
    RVector h(n), dold(n);

    RVector proj(ndat);
    RDenseMatrix J(ndat,n);
    RCompRowMatrix Graph(slen,slen);

    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    // initial solution (scaled)

    RVector x(x0), x1(x0);
    // current solution, trial solution (scaled);

    HESS_DATA hdata = {&J, &Graph, &x, &tau, &tau1, &tau2, &lambda};
    // set the components for the implicit definition of the Hessian

    tau1 = tau/* /x[0] */; tau2 = tau/* /x[n/2] */;
    cout << "tau 1 " << tau1 << " tau2 " << tau2 << endl;

    // allocate field vectors
    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);
    CVector *aphi = new CVector[mesh->nM];
    for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);

    int i_count = 0; // iteration counter
    cout << "Resetting forward solver\n";
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);
    proj = FWS.ProjectAll_real (mvec, dphi);
    err0 = errstart = of = OF_value (data, proj, sd) + tau; // estimate 
                                         // initial effect of prior as 1 ! 
    err00 = 0; // this is previous error

    //    test_biscale();
    cout << "Starting error " << errstart << endl;

    cout << "Building neighbour graph" << endl;
    int *rowptr, *colidx, nzero, row, col, i0, i1;
    raster.NeighbourGraph (rowptr, colidx, nzero);

    double *valptr = new double[nzero];
    for (row = 0; row < slen; row++) {
        i0 = rowptr[row]; i1 = rowptr[row+1];
	for (i = i0; i < i1; i++) {
	    col = colidx[i];
	    valptr[i] = (col == row ? i1-i0-1 : -1); // Discrete Laplacian
	}
    }
    Graph.Initialise (rowptr, colidx, valptr);
    Graph.PrintFillinGraph ("fillin.pgm", 600, true, false);
    delete []rowptr;
    delete []colidx;
    delete []valptr;

    // LM iteration loop
    for (inr = 0; inr < nrmax &&
	   err0 > ftol*errstart && fabs(err00-err0) > ftol; inr++) {

	RVector kap = bsol.GetParam (OT_CKAPPA);

	cout << "Calculating Jacobian ..." << endl;
	GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
            FWS.GetDataScaling(), J);

	// NOTE: MUST PREFORM DATA SCALING HERE!

	pscaler->ScaleJacobian (bsol.GetActiveParams(), J);
	if (bWriteJ) WriteJacobian (&J, raster, *mesh);

	// calculate Gradient
	RVector b = (data-proj)/sd;             // negative residual
	J.ATx (b,r);                            // gradient
	cout << "gradient norm " << l2norm(r) << endl;
	
	// penalty_gradient_tikh0_rescale_add(bsol.GetActiveParams(),x0,r,tau);
#ifdef LAPLACIAN_PRIOR
	penalty_gradient_tikh1_add (x, x0, r, tau1, tau2, Graph);
#else
	penalty_gradient_tikh0_add (x, x0, x0, r, tau);
#endif
	cout << "gradient norm with penalty " << l2norm(r) << endl;

	bool KeepGoing = true;
	i_count = 0;
	while (i_count < itmax && KeepGoing) {
	    cout << "LM iteration " << i_count << endl;

	    // solve Hh = r

	    switch (g_lm_precon) {
	    case LM_PRECON_NONE:    // assume H = I
	        h = r;
		break;
	    case LM_PRECON_HDIAG:   // assume H = diag(H)
	        h = r / Hess_diag (hdata);
		break;
	    case LM_PRECON_CH: {    // solve with Cholesky factorisation
	        RSymMatrix hess = Hess_full (hdata);  // need explicit Hessian
		CHdecomp (hess);
		h = CHsubst(hess, r);
	        } break;
	    case LM_PRECON_ICH:     // Incomplete Cholesky factorisation
	        cerr << "LM ICH preconditioner not yet implemented" << endl;
		exit(1);
	    case LM_PRECON_GMRES: {  // GMRES solver with implicit Hessian
	        double tol = g_lm_gmres_tol;
		GMRES (JTJx_clbk, &hdata, r, h, tol);
		//GMRES (Hess_full (hdata), r, h, tol); // sanity check
	        } break;
	    }
	    d = h;
	    double stepsize = l2norm(h)/l2norm(r); // Newton step/SD step
	    stepsize = 0.1;
	    cout << "\nstep size " << stepsize  << endl;

#ifdef NEWTON_CG
	    double beta;
	    if(inr >  0) {
	        beta = (h&r)/hrold;
		d += dold*beta;
	    }
	    dold = d;
	    hrold = h&r;
#endif

#ifdef LINESEARCH_IN_LM
	    if (alpha < 0.0) alpha = stepsize; // initialise step length
	    if (!LineSearchWithPrior (FWS, raster, pscaler, qvec, mvec, data,
		  sd, omega, d, err0, fmin, alpha, msol, logparam, x0, x, tau,
                  tau1, tau2, Graph)) {
	        lambda *= 4;
		cout << "No decrease in line search \n";
		continue;
	    }
	    err1 = fmin;
	    x1 = x + d*alpha;
#else
	    alpha = 1.0;
	    x1 = x + d*alpha;  // trial solution	
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(x1), msol);
	    FWS.Reset (msol, omega);
	    FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
	    FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
	    ProjectAll (mesh, FWS, mvec, dphi, proj);
	    err1 = OF_value (data, proj, sd);
#endif
	    if(err1 < err0) {  // accept update
	        err00 = err0;
		err0 = err1;
		x = x1;
		bsol.SetActiveParams (pscaler->Unscale(x));
		raster.Map_SolToMesh (bsol, msol);
		FWS.Reset (msol, omega);
		FWS.CalcFields (qvec, dphi);
		FWS.CalcFields (mvec, aphi);
		proj = FWS.ProjectAll_real (mvec, dphi);

		switch (g_imgfmt) {
		case IMGFMT_NIM:
		    msol.WriteImg_mua (i_count+1, "recon_mua.nim");
		    msol.WriteImg_mus (i_count+1, "recon_mus.nim");
		    break;
		case IMGFMT_RAW:
		    Solution rsol(OT_NPARAM, blen);
		    raster.Map_SolToBasis (bsol, rsol, true);
		    rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
		    rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
		    break;
		}
		lambda /= 4;
		KeepGoing = false;
	    } else {          // reject update
	        lambda *= 4;
	    }
	    cout << "error 0 " << err0 << " error 1 " << err1 << " lambda "
		 << lambda << " alpha " << alpha << " tau " << tau << endl;
	    i_count++;
	}
    } // end of NR loop;
}

// ==========================================================================

void LMsolveUnder (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    ObjectiveFunction &OF, int itmax, const RVector &data, const RVector &sd,
    Solution &msol, Solution &bsol, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, double omega, bool logparam, double ftol)
{
    bool Use_precon = true;

    // initialisations
    int i,j, inr, nrmax = 40;
    const QMMesh *mesh = FWS.MeshPtr();
    int dim  = raster.Dim();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int n    = slen*2;
    double of, fmin, errstart, err0, err00, err1;
    double lambda =10, alpha = -1.0,tau = 5e-3;
    double tau1, tau2;

    RVector r(n), s(n), d(n), M(n);
    RVector h(n), dold(n);
    RVector S(n), SS(n), IS(n);          // solution rescaling matrices
    RVector R(ndat), RR(ndat), IR(ndat); // data rescaling matrices
    RVector C(ndat); // holds diagonal JJT

    RVector proj(ndat);
    RDenseMatrix J;
    RSymMatrix Hess0, Hess;
    RVector x0 (bsol.GetActiveParams()); // initial state
    //    cout << "Initial state\n" << x0 << endl;
    RVector mu0(n/2), kap0(n/2);
    //    x0.Relink(mu0,0,n/2); 
    //    x0.Relink(kap0,n/2,n/2); 
    for(i = 0; i < n/2; i++) kap0[i] = x0[i+n/2];
    cout << "Initial kappa\n" << kap0 << endl;
    RVector x (bsol.GetActiveParams()); // current state
    RVector mu(n/2), kap(n/2);

    Solution bsol1 (OT_NPARAM, slen); // updated state
    RVector xs(bsol.GetActiveParams());       // assume variance of solution 
                                              // is initial state
    //    RVector mus(n/2), kaps(n/2);
    //    xs.Relink(mus,0,n/2); 
    //    xs.Relink(kaps,n/2,n/2); 

    tau1 = tau/xs[0]; tau2 = tau/xs[n/2];
    cout << "tau 1 " << tau1 << " tau2 " << tau2 << endl;

    // allocate field vectors
    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);
    CVector *aphi = new CVector[mesh->nM];
    for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);

    int i_count = 0; // iteration counter
    int k_count = 0; // reset counter
    cout << "Reseting forward solver\n";
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);
    proj = FWS.ProjectAll_real (mvec, dphi);
    err0 = errstart = of = OF_value (data, proj, sd) + tau; // estimate 
                                         // initial effect of prior as 1 ! 
    err00 = 0; // this is previous error

    //    test_biscale();
    cout << "Starting error " << errstart << endl;

    cout << "Building neighbour graph" << endl;
    int *rowptr, *colidx, nzero, row, col, i0, i1;
    raster.NeighbourGraph (rowptr, colidx, nzero);

    //raster.MapGraph2Sol (browptr, bcolidx, bnzero,rowptr, colidx, nzero);
    // MapGraphToSol: no longer required (done in NeighbourGraph)

    double *valptr = new double[nzero];
    for (row = 0; row < slen; row++) {
        i0 = rowptr[row]; i1 = rowptr[row+1];
	for (i = i0; i < i1; i++) {
	    col = colidx[i];
	    if(col > slen)
	      cout << "neighbour " << col << " is outside support\n";
//	    else
//	      cout << "(" << row << "," << col << ")";
	    valptr[i] = (col == row ? i1 -i0 - 1:-1); // Discrete Laplacian
	}
//	cout << endl;
    }
    RCompRowMatrix Graph(slen, slen, rowptr, colidx, valptr);
    Graph.PrintFillinGraph ("fillin.pgm", 600, true, false);
    delete []rowptr;
    delete []colidx;
    delete []valptr;

    for(inr = 0; inr < nrmax  && err0 > ftol*errstart && fabs(err00-err0) > ftol; inr++) {
      x = bsol.GetActiveParams();
         // calculate Hessian
      cout << "Calculating Jacobian ..." << endl;
      J.New (ndat, n);
      for(i = 0; i < n/2; i++) kap[i] = x[i+n/2];
      GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
          FWS.GetDataScaling(), J);

      // NOTE MUST PERFORM DATA SCALING HERE!

      //      cout << "Jacobian\n" << J << endl;
      //      ATA_diag (J, M);
      //      cout << "column sum squares of J\n" << M << endl;
      for(j = 0; j < n; j++) {
	S[j] = xs[j];                 // form a scaling matrix
	IS[j] = 1/xs[j];              // and its inverse  
      }

      if (bWriteJ) WriteJacobian (&J, raster, *mesh);
#ifdef ROW_NORMALISATION
      cout << "Calculating diagonal of JJT ..." << endl;
      AAT_diag (J, C);
      for(j = 0; j < ndat; j++) R[j] = (C[j] = (C[j] ? 1/sqrt(C[j]) : 1.0));
      cout << "rescaling Jacobian \n";
      J.RowScale (M);
#endif
      cout << "Calculating JJT ..." << flush;
      Hess0 = AAT(J);
      cout << " Finished!" << endl;
      //      cout << "Diagonal of Hessian\n";
      //      for(i = 0; i< n ; i++)cout << Hess0(i,i) << " ";


      for(i = 0; i < ndat ; i++)    {  // add prior Hessian to Hess
	    //       	  cout << Hess(i,i) << " -> ";
	Hess0(i,i) += tau; //*** THINK : should it have sd in it ??? ***//
       	  cout << Hess0(i,i) << " ";
	  R[i] = (C[i] = 1/sqrt(Hess0(i,i)));  // form rescaler
      }
      cout << "\nDone diagonal of JJT ..." << endl;
      cout << R << endl;
      /* Rescale the Hessian */
      BiScale(Hess0,C);
      cout << "\nDone rescaling " << endl;
      // calculate Gradient
      r.Clear();
      RVector b = (data-proj)/sd;             // negative residual
      bool KeepGoing = true;
      // test_cholesky();

      while (i_count < itmax && KeepGoing) {
        cout << "LM iteration " << i_count << endl;
        Hess = Hess0;
	for(i = 0; i < ndat ; i++)
          {
	    //       	  cout << Hess(i,i) << " -> ";
       	  Hess(i,i) += lambda ;
       	  // cout << Hess(i,i) << " ";
        }


        CHdecomp (Hess);
	r = CHsubst(Hess, b);
	J.ATx(r,h);
	double stepsize= l2norm(h)/l2norm(r); // Newton step/SD step
        cout << "\nstep size " << stepsize  << endl;
        cout << "\nunscaling\n";
        for(j = 0; j < n; j++) {	
	  d[j] = h[j]*S[j]; /* unscale */
	  //	  if(Hess(j,j)==1.0) d[j] = 1.0;   /* prevent zero pixels...*/
	}


	//	/* leave out line search for now */
#ifdef LINESEARCH_IN_LM
	if (alpha < 0.0) { // initialise step length
//	      alpha = (logparam ? 1.0 : 0.1 * l2norm (bsol) / l2norm (r));
	  //alpha = stepsize;
	  alpha = 1.0;
	}
	if (!LineSearchWithPrior (FWS, raster, pscaler, qvec, mvec, data, sd,
          omega,d,err0, fmin, alpha, msol, logparam,x0,xs,tau,tau1,tau2,Graph))
	  {
	    lambda *= 4;
	    cout << "No decrease in line search \n";
	    continue;
	  }
	
	err1 = fmin;
	bsol1 = bsol +d*alpha;
#else
      	alpha = 1.0;
	//	alpha = (logparam ? 1.0 : 0.1 * l2norm (bsol) / l2norm (r));

       	// update: x = x + alpha d
       
	//	cout << d << endl;
	//if(logparam) {
	//  cout << "unscaling log\n";
	//  bsol1 = bsol;
	//  bsol1 *= d;
	//  bsol1 += bsol;
        //}
	//else
	/* leave out test. Just do line search! */

	bsol1 = bsol + d*alpha; // trial solution
	
	raster.Map_SolToMesh (logparam ? exp(bsol1): bsol1, msol, true);

	FWS.Reset (msol, omega);
	FWS.CalcFields (mesh, mesh.nQ, qvec, dphi);
	FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
	ProjectAll (mesh, FWS, mvec, dphi, proj);
	err1 = OF_value (data, proj, sd);
#endif
	if(err1 < err0) {
	  err00 = err0;
	  err0 = err1;
	  bsol = bsol1; // trial solution	
	  lambda /= 4;
   	  raster.Map_SolToMesh (logparam ? exp(bsol): bsol, msol, true);

	  switch (g_imgfmt) {
	  case IMGFMT_NIM:
	      msol.WriteImg_mua (i_count+1, "recon_mua.nim");
	      msol.WriteImg_mus (i_count+1, "recon_mus.nim");
	      break;
	  case IMGFMT_RAW:
	      Solution rsol(OT_NPARAM, blen);
	      raster.Map_SolToBasis (logparam ? exp(bsol):bsol, rsol, true);
	      rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
	      rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
	      break;
	  }

	  FWS.Reset (msol, omega);
	  FWS.CalcFields (qvec, dphi);
	  FWS.CalcFields (mvec, aphi);
	  proj = FWS.ProjectAll_real (mvec, dphi);
	  KeepGoing = false;
	}
	else {
	  lambda *= 4;
	}
	cout << "error 0 " << err0 << " error 1 " << err1 << " lambda " << lambda << " alpha " << alpha << " tau " << tau << endl;
	i_count++;
      }
    } // end of NR loop;
}


bool LineSearch (CFwdSolver &FWS, const Raster &raster, ObjectiveFunction &OF,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &data, const RVector &sd, double omega, const Solution &s0,
    const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, RVector &proj, bool &proj_valid, bool logparam)
{
    cout << "  Linesearch Start" << endl;
    const int MAXIT = 16;
    double x0 = 0.0, xm, fm, fb, xb = lambda;
    proj_valid = false;
    Solution stepsol(s0);

    stepsol = s0+grad*xb;
    raster.Map_SolToMesh (stepsol, meshsol);
    while (!meshsol.Valid()) {
        cout << "    ** Invalid nodal parameters. Reducing step size" << endl;
	xb *= 0.5;
	stepsol = s0+grad*xb;
	raster.Map_SolToMesh (stepsol, meshsol);
    }
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fb = OF.get_value (&proj, &stepsol);

    cout << "    Step  " << xb << "  OF " << fb << endl;

    if (fb < f0) { // increase interval
        xm = xb; fm = fb;
	xb *= 2.0;
	stepsol = s0+grad*xb;
	raster.Map_SolToMesh (stepsol, meshsol);
	if (!meshsol.Valid()) {
   	    cout << "    ** Invalid nodal parameters. Truncating step size."
		 << endl;
	    fb = fm*10.0; // invalidate this step
	} else {
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fb = OF.get_value (&proj, &stepsol);
	    cout << "    Step  " << xb << "  OF " << fb << endl;
	}
	while (fb < fm) {
	    x0 = xm; f0 = fm;
	    xm = xb; fm = fb;
	    xb *= 2.0;
	    stepsol = s0+grad*xb;
	    raster.Map_SolToMesh (stepsol, meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fb = OF.get_value (&proj, &stepsol);
	    cout << "    Step  " << xb << "  OF " << fb << endl;
	}
    } else { // decrease interval
        xm = 0.5*xb;
	stepsol = s0+grad*xm;
	raster.Map_SolToMesh (stepsol, meshsol);
	proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	fm = OF.get_value (&proj, &stepsol);
	cout << "    Step  " << xm << "  OF " << fm << endl;
	int itcount = 0;
	while (fm > f0) {
  	    if (++itcount > MAXIT) return false;
	    xb = xm; fb = fm;
	    xm = 0.5*xb;
	    stepsol = s0+grad*xm;
	    raster.Map_SolToMesh (stepsol, meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fm = OF.get_value (&proj, &stepsol);
	    cout << "    Step  " << xm << "  OF " << fm << endl;
	}
    }
    // quadratic interpolation
    double a = ((f0-fb)/(x0-xb) - (f0-fm)/(x0-xm)) / (xb-xm);
    double b = (f0-fb)/(x0-xb) - a*(x0+xb);
    lambda = -b/(2.0*a);
    stepsol = s0+grad*lambda;
    raster.Map_SolToMesh (stepsol, meshsol);
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fmin = OF.get_value (&proj, &stepsol);
    if (fmin > fm) {  // interpolation didn't give improvement
        lambda = xm, fmin = fm;
    } else proj_valid = true;
    cout << "    Final " << lambda << "  OF " << fmin << endl;
    cout << "  Linesearch End" << endl;
    return true;
}

bool LineSearchWithPrior (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &grad, double f0, double &fmin, double &lambda,
    Solution &meshsol, bool logparam, const RVector &p0, const RVector &p,
    const double &tau, const double tau1, const double tau2,
    const RCompRowMatrix &LTL)
{
    const int MAXIT = 6;
    double x0 = 0.0, xm, fm, fmd, fmp, fb, fbd, fbp, fmind, fminp,xb = lambda;
    RVector proj(data.Dim());
    
    RVector p1 = p+grad*xb;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fbd = OF_value (data, proj, sd);

#ifdef LAPLACIAN_PRIOR
    fbp = penalty_tikh1 (p1, p0, tau1, tau2, LTL);
    fb = fbd + fbp;
#else
    fbp = penalty_tikh0 (p1, p0, p);
    fb = fbd + tau*fbp;
#endif

    cout << "  Lsearch: step " << xb << ", of=" << fbd << " prior " << fbp
	 << " total " << fb << endl;

    if (fb < f0) { // increase interval
        xm = xb; fm = fb;
	xb *= 2.0;
	p1 = p+grad*xb;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	fbd = OF_value (data, proj, sd);

#ifdef LAPLACIAN_PRIOR
	fbp = penalty_tikh1 (p1, p0, tau1, tau2, LTL);
	fb = fbd + fbp;
#else
	fbp = penalty_tikh0 (p1, p0, p);
	fb = fbd + tau*fbp;
#endif

	cout << "  Lsearch: step " << xb << ", of=" << fbd << " prior " << fbp
	     << " total " << fb << endl;
	while (fb < fm) {
	    x0 = xm; f0 = fm;
	    xm = xb; fm = fb;
	    xb *= 2.0;
	    p1 = p+grad*xb;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    	    fbd = OF_value (data, proj, sd);

#ifdef LAPLACIAN_PRIOR
	    fbp = penalty_tikh1 (p1, p0, tau1, tau2, LTL);
	    fb = fbd + fbp;
#else
	    fbp = penalty_tikh0 (p1, p0, p);
	    fb = fbd + tau*fbp;
#endif

	    cout << "  Lsearch: step " << xb << ", of=" << fbd << " prior "
		 << fbp << " total " << fb << endl;
	}
    } else { // decrease interval
        xm = 0.5*xb;
	p1 = p+grad*xm;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    	fmd = OF_value (data, proj, sd);

#ifdef LAPLACIAN_PRIOR
	fmp = penalty_tikh1 (p1, p0, tau1, tau2, LTL);
	fm = fmd + fmp;
#else
	fmp = penalty_tikh0 (p1, p0, p);
	fm = fmd + tau*fmp;
#endif

	cout << "  Lsearch: step " << xm << ", of=" << fmd << " prior " << fmp
	     << " total " << fm << endl;
	int itcount = 0;
	while (fm > f0) {
  	    if (++itcount > MAXIT) return false;
	    xb = xm; fb = fm;
	    xm = 0.5*xb;
	    p1 = p+grad*xm;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fmd = OF_value (data, proj, sd);

#ifdef LAPLACIAN_PRIOR
	    fmp = penalty_tikh1 (p1, p0, tau1, tau2, LTL);
	    fm = fmd + fmp;
#else
	    fmp = penalty_tikh0 (p1, p0, p);
	    fm = fmd + tau*fmp;
#endif

 	    cout << "  Lsearch: step " << xm << ", of=" << fmd << " prior "
		 << fmp << " total " << fm << endl;
	}
    }
    // quadratic interpolation
    double a = ((f0-fb)/(x0-xb) - (f0-fm)/(x0-xm)) / (xb-xm);
    double b = (f0-fb)/(x0-xb) - a*(x0+xb);
    lambda = -b/(2.0*a);
    p1 = p+grad*lambda;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fmind = OF_value (data, proj, sd);

#ifdef LAPLACIAN_PRIOR
    fminp = penalty_tikh1 (p1, p0, tau1, tau2, LTL);
    fmin = fmind + fminp;
#else
    fminp = penalty_tikh0 (p1, p0, p);
    fmin = fmind + tau*fminp;
#endif

    if (fmin > fm) {  // interpolation didn't give improvement
        lambda = xm, fmin = fm;
	cout << "  Lsearch (no interpolation) : final step " << lambda
	     << ", total of + prior =" << fmin << " tau " << tau << endl;
    } else {
        cout << "  Lsearch (with interpolation) : final step " << lambda
	     << ", of=" << fmind << " prior " << fminp << " total " << fmin
	     << " tau " << tau << endl;
    }
    // restimate tau 
    //    tau = fmind/fminp;
 
    return true;
}

void test_cholesky()
{
  /* Test... */
      RDenseMatrix A(2,3);
      for(int i = 0; i < 2; i++)
	for(int j = 0; j < 3; j++)
	  A(i,j) = 1 + i*j;
      RSymMatrix H = ATA(A);
      cout << "A spd matrix... \n" << H << endl;
      CHdecomp (H);
      cout << "A cholesky decomposition of H ..." << H << endl;
}
void test_biscale()
{
  /* Test... */
      int i,j;
      RDenseMatrix A(2,3);
      RVector S(3);
      A(0,0) = 3; A(1,0) = 4;
      A(0,1) = 5; A(1,1) = 12;
      A(0,2) = 1; A(1,2) = 1;
      for(i = 0; i < 2; i++) {
	for(j = 0; j < 3; j++) {
	  S[j] += A(i,j)*A(i,j);
	}
      }
      cout << "original matrix ...\n" << A << endl;
      cout << "column sum squares " << S << endl;
      for(j = 0; j < 3; j++) S[j] = 1/sqrt(S[j]);
      cout << "scale vector " << S << endl;
      RSymMatrix H = ATA(A);
      cout <<  "ATA before biscaling... \n" << H << endl;
      BiScale(H,S);
      cout << "result of biscaling... \n" << H << endl;
      ColumnScale(A,S);
      cout << "rescaled matrix ...\n" << A << endl;
      RSymMatrix H2 = ATA(A);
      cout << "ATA of rescaled matrix... \n" << H2 << endl;
}
