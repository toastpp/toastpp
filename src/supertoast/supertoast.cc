#include "stoastlib.h"
#include <iomanip>
#include <string.h>
#include "util.h"
#include "solver.h"
#include "source.h"
#include "timing.h"
#include "supertoast.h"
#include "supertoast_util.h"
#include <time.h>
#if defined(WIN32) || defined(WIN64)
#else
#include <unistd.h>
#endif

//latest revision
using namespace std;

// Max number of mesh regions
#define MAXREGION 100

// Verbose timing output
//#define DO_PROFILE

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
int g_imgfmt = IMGFMT_NIM;

double clock0, wt0;

#ifdef DO_PROFILE
double solver_time = 0.0;
#endif

// =========================================================================
// global parameters

SourceMode srctp = SRCMODE_NEUMANN;   // source type
int bWriteJ    = 0;
int bWriteGrad = 0;
int bOutputUpdate = 1;
int bOutputGradient = 1;
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
void SelectInitialParams (const Mesh &mesh, Solution &msol);
void SelectData (DataScale dscale, int nqm, double &freq, RVector &data);
void SelectBasis (IVector &gdim, IVector &bdim);
int  SelectImageFormat ();
PARAM_SCALE SelectParamScaling ();
int  SelectSDMode ();

// =========================================================================
// MAIN

int main (int argc, char *argv[])
{
    clock0 = tic();
    wt0 = walltic();

    if (argc > 1 && pp.Open (argv[1]))
        cout << "Reading parameters from " << argv[1] << endl;
    if (argc > 2) {
        pp.LogOpen (argv[2]);
	cout << "Writing log to " << argv[2] << endl;
    } else {
	pp.LogOpen ("supertoast.out");
	cout << "Writing log to supertoast.out" << endl;
    }
    OutputProgramInfo ();

#ifdef TOAST_THREAD
    int nth = 0;
    pp.GetInt("NTHREAD", nth);
    Task_Init (nth);
#endif

    const double c0 = 0.3;

    char meshname[256], cbuf[256];
    double omega;
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
    Solver *solver = Solver::Create (&pp); // the nonlinear solver

    CCompRowMatrix qvec, mvec;
    CVector *dphi, *aphi;
    int cmd;
    CFwdSolver FWS (&mesh, pp);
    FWS.WriteParams (pp);
    g_lsolver_tol = FWS.GetLinSolverTol();
    RVector data, data1, data2;
    bool useref = false;

    // solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    SelectInitialParams (mesh, msol);

    SelectData (FWS.GetDataScaling(), nQM, omega, data);
    sdmode = SelectSDMode();
    omega *= 2.0*Pi*1e-6; // convert from MHz to rad/ps

    data1.Relink (data, 0, nQM);
    data2.Relink (data, nQM, nQM);

    // Generate sd vector by assuming sd=data ('normalised')
    RVector sd(data);
#ifdef USE_DATA_AVERAGE_AS_SCALING
    // now try just the average of the data
    double avd1 =0.0, avd2 = 0.0;
    for (i = 0; i < nQM; i++){
      avd1 += data1[i];
      avd2 += data2[i];
    }
    avd1 /= nQM;
    avd2 /= nQM;
    cout << "Average amplitude : " << avd1 << " phase " << avd2 << endl;
    for (i = 0; i < nQM; i++){
      sd[i] = avd1;
      sd[i+nQM] = avd2;
    }
#endif

    // build the source vectors
    qvec.New (nQ, n);
    for (i = 0; i < nQ; i++) {
	CVector q(n);
	switch (qprof) {
	case 0:
	    SetReal (q, QVec_Point (mesh, mesh.Q[i], srctp));
	    break;
	case 1:
	    SetReal (q, QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp));
	    break;
	case 2:
	    SetReal (q, QVec_Cosine (mesh, mesh.Q[i], qwidth, srctp));
	    break;
	}
	qvec.SetRow (i, q);
    }

    // build the measurement vectors
    mvec.New (nM, n);
    RVector c2a = msol.GetParam(OT_C2A);
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
	CVector cm(n);
	SetReal (cm, m*c2a);
	//for (j = 0; j < n; j++) m[j] *= mesh.plist[j].C2A();
	mvec.SetRow (i, cm);
    }


    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new CVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // allocate system matrix
    cout << endl << "Allocating system matrix" << endl;
    FWS.Allocate ();

    cout << "Generating intermediate pixel grid" << endl;
    
    if (blen == glen) // use single grid basis
	raster = new Raster_Pixel (gdim, gdim, &mesh);
    else              // use intermediate basis
	raster = new Raster_Pixel (gdim, bdim, &mesh);

    glen = raster->GLen();
    blen = raster->BLen();
    slen = raster->SLen();

    // solution in sparse user basis
    Solution bsol(OT_NPARAM, raster->SLen());
    bsol.SetActive (OT_CMUA, true);
    bsol.SetActive (OT_CKAPPA, true);
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

    cout << "  Original solution:\n";
    cout << "    CMUA:   " << vmin (msol.GetParam(OT_CMUA)) << " to "
	 << vmax (msol.GetParam(OT_CMUA)) << endl;
    cout << "    CKAPPA: " << vmin (msol.GetParam(OT_CKAPPA)) << " to "
	 << vmax (msol.GetParam(OT_CKAPPA)) << endl;

    cout << "  Mapped solution" << endl;
    cout << "    CMUA:   " << vmin (bsol.GetParam(OT_CMUA)) << " to "
	 << vmax (bsol.GetParam(OT_CMUA)) << endl;
    cout << "    CKAPPA: " << vmin (bsol.GetParam(OT_CKAPPA)) << " to "
	 << vmax (bsol.GetParam(OT_CKAPPA)) << endl;

    // map bsol back to msol to make sure msol contains a valid
    // basis representation of the solution
    //  raster->Map_SolToMesh (bsol, msol, true);

    //cout << "  Re-mapped solution:\n";
    //cout << "    CMUA:   " << vmin (msol.GetParam(OT_CMUA)) << " to "
    // << vmax (msol.GetParam(OT_CMUA)) << endl;
    //cout << "    CKAPPA: " << vmin (msol.GetParam(OT_CKAPPA)) << " to "
    // << vmax (msol.GetParam(OT_CKAPPA)) << endl;

    cout << "Assembling and pre-processing system matrix" << endl;
    FWS.Reset (msol, omega);

    if (g_imgfmt != IMGFMT_RAW) {
        WriteNimHeader (meshname, n, "recon_mua.nim", "MUA");
	WriteNimHeader (meshname, n, "recon_mus.nim", "MUS");
	// write out the initial images
	msol.WriteImg_mua (0, "recon_mua.nim");
	msol.WriteImg_mus (0, "recon_mus.nim");
    }
    if (g_imgfmt != IMGFMT_NIM) {
        Solution rsol (OT_NPARAM, raster->BLen());
	raster->Map_SolToBasis (bsol, rsol, true);
	WriteRimHeader (raster->BDim(), "recon_mua.raw");
	WriteRimHeader (raster->BDim(), "recon_mus.raw");
	rsol.WriteImg_mua (0, "recon_mua.raw");
	rsol.WriteImg_mus (0, "recon_mus.raw");
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
        RVector refdata(nQM);
	RVector proj(nQM*2);
	char file_re[32], file_im[32];
	char name_re[32], name_im[32];
        switch (FWS.GetDataScaling()) {
	case DATA_LIN:
	    strcpy (file_re, "REFDATA_REAL");
	    strcpy (file_im, "REFDATA_IMAG");
	    strcpy (name_re, "real component");
	    strcpy (name_im, "imaginary component");
	    break;
	case DATA_LOG:
	    strcpy (file_re, "REFDATA_MOD");
	    strcpy (file_im, "REFDATA_ARG");
	    strcpy (name_re, "log amplitude");
	    strcpy (name_im, "phase");
	    break;
	}

	if (!pp.GetString (file_re, cbuf)) {
	    cout << "Reference file: " << name_re << ": ";
	    cin  >> cbuf;
	}
	ReadDataFile (cbuf, refdata);
	pp.PutString (file_re, cbuf);
	data1 -= refdata;
	
	if (!pp.GetString (file_im, cbuf)) {
	    cout << "Reference file: " << name_im << ": ";
	    cin  >> cbuf;
	}
	ReadDataFile (cbuf, refdata);
	pp.PutString (file_im, cbuf);
	data2 -= refdata;
	cout << "Generating model baseline" << endl;
	FWS.CalcFields (qvec, dphi);
	proj = FWS.ProjectAll_real (mvec, dphi);
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
	proj = FWS.ProjectAll_real (mvec, dphi);
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
	RVector proj(nQM*2);

	FWS.CalcFields (qvec, dphi);
	proj = FWS.ProjectAll_real (mvec, dphi);

	sd = (data - proj);
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
    }

    // set up objective function
    ObjectiveFunction OF (raster);
    OF.SetData (&data);
    OF.SetSD (&sd);

    cout << "Ref data  " << (useref ? "YES":"NO") << endl;
    cout << endl;

    // ==================================================================
    // Start the solver
    LOGOUT("**** SOLVER started: CPU %f", toc(clock0));
    solver->Solve (FWS, *raster, pscaler, OF, data, sd, bsol, msol, qvec, mvec,
        omega);
    LOGOUT("**** SOLVER finished: CPU %f", toc(clock0));
    delete solver;

    // cleanup
    delete []dphi;
    delete []aphi;

    //times (&tme);
    //double total_time = (double)(tme.tms_utime-clock0)/(double)HZ;
    double total_time = toc(clock0);
    double total_wallc = walltoc(wt0);

#ifdef DO_PROFILE
    LOGOUT("Solver: %f", solver_time);
#endif
    LOGOUT("Total timings: %f real, %f CPU", total_wallc, total_time);

    return 0;
}                                                                              

// ============================================================================

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
    if ((host = getenv("HOST")))
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

void SelectInitialParams (const Mesh &mesh, Solution &msol)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    RVector param[3];
    int i, j, k, n, p, nreg, nregnode[MAXREGION];
    const char *resetstr[3] = {"RESET_MUA", "RESET_MUS", "RESET_N"};
    const ParameterType prmtp[3] = {PRM_MUA, PRM_MUS, PRM_N};
    for (p = 0; p < 3; p++) {

	param[p].New(mesh.nlen());
	if (pp.GetString (resetstr[p], cbuf)) {
	    pp.PutString (resetstr[p], cbuf);
	    if (!strcasecmp (cbuf, "MESH")) {
		xERROR("This option is no longer supported");
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
	    cout << "\nSelect initial distribution for " << resetstr[p]
		 << endl;
	    cout << "(1) Global homogeneous\n";
	    cout << "(2) Homogeneous in regions\n";
	    cout << "(3) Nodal image file (NIM)\n";
	    cout << "[1|2|3] >> ";
	    cin >> resettp;
	    switch (resettp) {
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
		cin >> (cbuf+4);
		ReadNim (cbuf+4, param[p]);
		break;
	    }
	    pp.PutString (resetstr[p], cbuf);
	}
    }
    g_refind = mean (param[2]); // assuming homogeneous n
    msol.SetParam (OT_CMUA, param[0]*c0/param[2]);
    msol.SetParam (OT_CKAPPA, c0/(3.0*param[2]*(param[0]+param[1])));
    msol.SetParam (OT_N, param[2]);
    for (i = 0; i < param[2].Dim(); i++)
	param[2][i] = c0/(2.0*param[2][i]*A_Keijzer(param[2][i]));
    msol.SetParam (OT_C2A, param[2]);
}

// ============================================================================

void SelectData (DataScale dscale, int nqm, double &freq, RVector &data)
{
    char cbuf[256];
    bool def = false;
    RVector pdata;

    switch (dscale) {
    case DATA_LIN:
	data.New (2*nqm);
	if (!pp.GetString ("DATA_REAL", cbuf)) {
	    cout << "\nData file 1 (real):\n>> ";
	    cin >> cbuf;
	}
	pdata.Relink (data,0,nqm);
	ReadDataFile (cbuf, pdata);
	pp.PutString ("DATA_REAL", cbuf);
	if (!pp.GetString ("DATA_IMAG", cbuf)) {
	    cout << "\nData file 2 (imag)\n>> ";
	    cin >> cbuf;
	}
	pdata.Relink (data,nqm,nqm);
	ReadDataFile (cbuf, pdata);
	pp.PutString ("DATA_IMAG", cbuf);
	break;
    case DATA_LOG:
	data.New (2*nqm);
	if (!pp.GetString ("DATA_MOD", cbuf)) {
	    cout << "\nData file 1 (mod):\n>> ";
	    cin >> cbuf;
	}
	pdata.Relink (data,0,nqm);
	ReadDataFile (cbuf, pdata);
	pp.PutString ("DATA_MOD", cbuf);
	if (!pp.GetString ("DATA_ARG", cbuf)) {
	    cout << "\nData file 2 (arg)\n>> ";
	    cin >> cbuf;
	}
	pdata.Relink (data,nqm,nqm);
	ReadDataFile (cbuf, pdata);
	pp.PutString ("DATA_ARG", cbuf);
	break;
    }

    if (!pp.GetReal ("FREQ", freq) || freq < 0.0) do {
	cout << "\nSource modulation frequency [MHz]:\n>> ";
	cin >> freq;
    } while (freq < 0.0);
    pp.PutReal ("FREQ", freq);
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
    static const char *fmtstr[4] = {"RAW", "NIM", "RAW+NIM"};
    char cbuf[256];
    int fmt = 0;

    if (pp.GetString ("IMAGEFORMAT", cbuf)) {
	if      (!strcasecmp (cbuf, "RAW"))     fmt = IMGFMT_RAW;
	else if (!strcasecmp (cbuf, "NIM"))     fmt = IMGFMT_NIM;
	else if (!strcasecmp (cbuf, "RAW+NIM")) fmt = IMGFMT_RAW_NIM;
	//else if (!strcasecmp (cbuf, "PGM")) fmt = IMGFMT_PGM;
	//else if (!strcasecmp (cbuf, "PPM")) fmt = IMGFMT_PPM;
    }
    while (fmt < 1 || fmt > 3) {
	cout << "\nOutput image format:\n";
	cout << "(1) Raw data\n";
	cout << "(2) NIM (nodal image\n";
	cout << "(3) Both (RAW+NIM)\n";
	//cout << "(3) PGM (portable grayscale map) - 2D only\n";
	//cout << "(4) PPM (portable pixmap) - 2D only\n";
	cout << "[1|2|3] >> ";
	cin  >> fmt;
    }
    pp.PutString ("IMAGEFORMAT", fmtstr[fmt-1]);
    return fmt;
}

// ============================================================================

int SelectSDMode ()
{
    static const char *modestr[4] =
        {"NONE", "DATA", "AVG_DATA", "AVG_DIFFDATA"};
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

// ==========================================================================

static bool CheckRange (const Solution &sol)
{
    bool inrange = true;

#ifdef LIMIT_RANGE
    const double MIN_CMUA = 0;
    const double MAX_CMUA = 0.1;
    const double MIN_CKAPPA = 0;
    const double MAX_CKAPPA = 10;

    double vmin, vmax;
    sol.Extents (OT_CMUA, vmin, vmax);
    if (vmin < MIN_CMUA || vmax > MAX_CMUA) {
	cerr << "WARNING: " << vmin << " < CMUA < " << vmax
	     << " in trial solution" << endl;
	inrange = false;
    }
    sol.Extents (OT_CKAPPA, vmin, vmax);
    if (vmin < MIN_CKAPPA || vmax > MAX_CKAPPA) {
	cerr << "WARNING: " << vmin << " < CKAPPA < " << vmax
	     << " in trial solution" << endl;
	inrange = false;
    }
#endif

    return inrange;
}

// ==========================================================================
// The callback function for obtaining the objective function during
// line search

double of_clbk (const RVector &x, double *of_sub, void *context)
{
    OF_CLBK_DATA *ofdata = (OF_CLBK_DATA*)context;
    CFwdSolver *fws = ofdata->fws;
    const Raster *raster = ofdata->raster;
    const Scaler *pscaler = ofdata->pscaler;
    Solution *meshsol = ofdata->meshsol;
    const Regularisation *reg = ofdata->reg;
    const CCompRowMatrix *qvec = ofdata->qvec;
    const CCompRowMatrix *mvec = ofdata->mvec;
    double omega = ofdata->omega;
    const RVector *data = ofdata->data;
    const RVector *sd = ofdata->sd;

    raster->Map_ActiveSolToMesh (pscaler->Unscale(x), *meshsol);
    if (!CheckRange (*meshsol)) return -1.0; // error

    RVector proj = fws->ProjectAll_real (*qvec, *mvec, *meshsol, omega);
    if (visnan (proj)) return -1.0; // error

    double fd = ObjectiveFunction::get_value (*data, proj, *sd);
    double fp = reg->GetValue (x);
    if (of_sub) {
	of_sub[0] = fd;
	of_sub[1] = fp;
    }
    return fd + fp;
}
