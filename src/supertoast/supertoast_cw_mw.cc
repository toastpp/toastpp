#include "stoastlib.h"
#include <fstream>
#include <iomanip>
#include <string.h>
#include "util.h"
#include "solver_cw.h"
#include "source.h"
#include "supertoast_cw_mw.h"
#include "supertoast_util.h"
#include "timing.h"
#include <time.h>
#include <unistd.h>

#ifdef TOAST_MPI
#include "mpi.h"
#endif // TOAST_MPI

using namespace std;

// Max number of mesh regions
#define MAXREGION 100

// Multi-Wavelength Related Parameters
#define MAX_NOFWLENGTH 50

// Verbose timing output
//#define DO_PROFILE

#define LIMIT_RANGE

// ==========================================================================
// Global variables

enum PARAM_SCALE {          // parameter scaling method
    PARAM_SCALE_NONE,       //   no scaling
    PARAM_SCALE_AVG,        //   scale with average of initial distribution
    PARAM_SCALE_LOG,        //   log scaling
    PARAM_SCALE_LINLOG,     //   lin+log scaling
    PARAM_SCALE_BOUNDLOG    //   x -> ln ((x-a)(b-x))
} g_pscale = PARAM_SCALE_NONE;

double g_lsolver_tol;       // convergence criterion for linear solver

double g_param_cmuamin, g_param_cmuamax;
double g_param_ckappamin, g_param_ckappamax;
double g_refind;
char g_meshname[256];
char g_prefix[256] = "\0";
int g_imgfmt = IMGFMT_NIM;

SourceMode srctp = SRCMODE_NEUMANN;   // source type
int bOutputUpdate = 0;
int bOutputGradient = 0;
double avg_cmua = 1.0, avg_ckappa = 1.0;
ParamParser pp;

double clock0, wclock0;

#ifdef DO_PROFILE
double solver_time = 0.0;
#endif

// =========================================================================
// local prototypes

void OutputProgramInfo ();
void SelectMesh (char *meshname, QMMesh &mesh);
void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp);
void SelectMeasurementProfile (int &mtype, double &mwidth);
void SelectInitialParams (const Mesh &mesh, MWsolution &msol,
    const RVector &wlength);
void SelectInitialReferenceParams (const Mesh &mesh, Solution &msol,
    int whichWavel);
void SelectData (DataScale dscale, int nqm, int nlambda, const RVector &wlength,
    RVector &data);
void SelectRefdata (DataScale dscale, int nqm, int nlambda,
    const RVector &wlength, RVector &data, bool &refEqualinitial);
void SelectBasis (IVector &gdim, IVector &bdim);
int  SelectImageFormat ();
PARAM_SCALE SelectParamScaling ();
int  SelectSDMode ();

// =========================================================================
// MAIN

int main (int argc, char *argv[])
{
    clock0 = tic();
    wclock0 = walltic();

    if (argc > 1 && pp.Open (argv[1]))
        cout << "Reading parameters from " << argv[1] << endl;
    if (argc > 2) {
        pp.LogOpen (argv[2]);
	cout << "Writing log to " << argv[2] << endl;
    } else {
	pp.LogOpen ("supertoast_cw_mw.out");
	cout << "Writing log to supertoast_cw_mw.out" << endl;
    }
    OutputProgramInfo ();
    pp.GetString("PREFIX", g_prefix);
    pp.PutString("PREFIX", g_prefix);
    
#ifdef TOAST_MPI
    // Initialise MPI
    int mpi_rank, mpi_size;
    int mpi_status = MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
    LOGOUT("Initialising MPI session. MPI status = %d", mpi_status);
    LOGOUT("Processor %d of %d", mpi_rank, mpi_size);

    cerr << "processor " << mpi_rank << " of " << mpi_size << endl;
#endif // TOAST_MPI

#ifdef TOAST_THREAD
    int nth = 0;
    pp.GetInt("NTHREAD", nth);
    Task_Init (nth);
    pp.PutInt("NTHREAD", nth);
#endif

    const double c0 = 0.3;

    char meshname[256], cbuf[256];
    double wscale;
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
    Solver_CW *solver = Solver_CW::Create (&pp); // the nonlinear solver

    RCompRowMatrix qvec, mvec;
    RVector *dphi, *aphi;
    int cmd;
    RFwdSolverMW FWS (&mesh, pp);
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
	xASSERT(nofwavel <= MAX_NOFWLENGTH, "Max wavelength count exceeded");
    }

    if (!pp.GetReal ("WAVELENGTH_SCALE", wscale))
	wscale = 1.0;
    pp.PutReal ("WAVELENGTH_SCALE", wscale);

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
    MWsolution msol(nprm, n, nofwavel, extcoef, wlength*wscale);

    SelectInitialParams (mesh, msol, wlength);
    SelectData (FWS.GetDataScaling(), nQM, nofwavel, wlength, data);
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
    RVector c2a = msol.GetC2A();
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
	mvec.SetRow (i, m*c2a);
    }


    // build the field vectors
    dphi = new RVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new RVector[nM];
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

#ifdef UNDEF
    RVector tmp1(n); tmp1 = 1;
    RVector tmp1m(slen);
    raster->Map_MeshToSol(tmp1,tmp1m);
    ofstream ofs1("mesh2sol.dat");
    ofs1 << tmp1m << endl;

    RVector tmp2(slen); tmp2=1;
    RVector tmp2m(n);
    raster->Map_SolToMesh (tmp2,tmp2m);
    ofstream ofs2 ("sol2mesh.dat");
    ofs2 << tmp2m << endl;
    exit(0);
#endif

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
    case PARAM_SCALE_LINLOG: {
	RVector scale = bsol.GetActiveParams();
	for (i = 0; i < bsol.nActive(); i++) {
	    double sum = 0;
	    for (j = 0; j < slen; j++)
	        sum += scale[i*slen+j];
	    for (j = 0; j < slen; j++)
	        scale[i*slen+j] = slen/sum;
	}
	pscaler = new LinLogScaler (scale);
	} break;
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

    if (g_imgfmt != IMGFMT_RAW) {
	for (i = 0; i < msol.nParam(); i++) {
	    char fname[256];
	    if (msol.IsActive(i)) {
		if (i < msol.nmuaChromo) 
		    sprintf (fname,"%sreconChromophore_%d.nim",g_prefix, i+1);
		if (i == msol.nmuaChromo)
		    sprintf (fname,"%sreconScatPrefactor_A.nim",g_prefix);
		if (i == msol.nmuaChromo + 1) 
		    sprintf (fname,"%sreconScatPower_b.nim",g_prefix);
		WriteNimHeader (meshname, n, fname, "N/A");  
		msol.WriteImgGeneric (0, fname, i);
	    }
	}
    }
    if (g_imgfmt != IMGFMT_NIM) {
	Solution gsol(nprm, raster->BLen());
	raster->Map_SolToBasis (bsol, gsol, true);
	for (i = 0; i < msol.nParam(); i++) {
	    if (msol.IsActive(i)) {
		char fname[256];
		if (i < msol.nmuaChromo) 
		    sprintf (fname,"%sreconChromophore_%d.raw",g_prefix,i+1);
		if (i == msol.nmuaChromo)
		    sprintf (fname,"%sreconScatPrefactor_A.raw",g_prefix);
		if (i == msol.nmuaChromo + 1) 
		    sprintf (fname,"%sreconScatPower_b.raw",g_prefix);
		WriteRimHeader (raster->BDim(), fname);
		gsol.WriteImgGeneric (0, fname, i);
	    }
	}
    }

    if (bOutputUpdate) {
	WriteNimHeader (meshname, n, "update_mua.nim", "MUA");
	WriteNimHeader (meshname, n, "update_mus.nim", "MUS");
    }
    if (bOutputGradient) {
	WriteNimHeader (meshname, n, "gradient_mua.nim", "MUA");
	WriteNimHeader (meshname, n, "gradient_mus.nim", "MUS");
    }

    if (useref) {
        RVector refdata(nQM*nofwavel);
	RVector proj(nQM*nofwavel);
	bool refEqualinitial;

	SelectRefdata (FWS.GetDataScaling(), nQM, nofwavel, wlength, refdata,
	    refEqualinitial);
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
	RVector proj = FWS.ProjectAll_wavel (qvec, mvec, msol, 0);
	sd = proj;
	for (i = 0; i < nofwavel; i++) {
	    double avd = 0.0;
	    for (j = 0; j < nQM; j++) avd += sd[i*nQM+j]*sd[i*nQM+j];
	    avd = sqrt(avd);
	    for (j = 0; j < nQM; j++) sd[i*nQM+j] = avd;
	}
	} break;		
    case 4: { // scale with difference averages over data types
	cout << "Generating model baseline" << endl;
	RVector proj = FWS.ProjectAll_wavel (qvec, mvec, msol, 0);
	sd = (data - proj);
	for (i = 0; i < nofwavel; i++) {
	    double avd = 0.0;
	    for (j = 0; j < nQM; j++) avd += sd[i*nQM+j]*sd[i*nQM+j];
	    avd = sqrt(avd);
	    for (j = 0; j < nQM; j++) sd[i*nQM+j] = avd;
	}
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

    // cleanup
    delete []dphi;
    delete []aphi;

    double total_time = toc(clock0);
    double total_wtime = walltoc(wclock0);

    LOGOUT ("Final timings:");
    LOGOUT("Total: CPU: %f, Wall: %f", total_time, total_wtime);
#ifdef TOAST_THREAD
    LOGOUT("Multiprocessing: CPU: %f, Wall %f", Task::GetThreadCPUTiming(),
	   Task::GetThreadWallTiming());
#endif

#ifdef DO_PROFILE
    LOGOUT("Solver: %f", solver_time);
#endif

#ifdef TOAST_MPI
    MPI_Finalize();
#endif

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

// ============================================================================

void SelectInitialParams (const Mesh &mesh, MWsolution &msol,
    const RVector &wlength)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    int nparam = msol.nParam();
    RVector *param = new RVector[nparam];
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
		   (int)wlength[p-msol.nmuaChromo-3]);
      }

      param[p].New(mesh.nlen());
	if (pp.GetString (resetstr, cbuf)) {
	    pp.PutString (resetstr, cbuf);
	    if (!strncasecmp (cbuf, "HOMOG", 5)) {
	        if (sscanf (cbuf+5, "%lf", &prm) != 1)
		    xERROR("Parse error on initial parameters!");
		param[p] = prm;
	    } else if (!strncasecmp (cbuf, "REGION_HOMOG", 12)) {
		valstr = strtok (cbuf+12, " \t");
		for (n = 0; n < MAXREGION && valstr; n++) {
		    if (sscanf (valstr, "%lf", reg_prm+n) != 1)
		        xERROR("Parse error on initial parameters!");
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
    
    delete []param;
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
		xERROR("This option is no longer supported");
	    } else if (!strncasecmp (cbuf, "HOMOG", 5)) {
	        if (sscanf (cbuf+5, "%lf", &prm) != 1)
		    xERROR("Parse error on reference parameters!");
		param[p] = prm;
	    } else if (!strncasecmp (cbuf, "REGION_HOMOG", 12)) {
		valstr = strtok (cbuf+12, " \t");
		for (n = 0; n < MAXREGION && valstr; n++) {
		    if (sscanf (valstr, "%lf", reg_prm+n) != 1)
		        xERROR("Parse error on reference parameters!");
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

void SelectData (DataScale dscale, int nqm, int nlambda, const RVector &wlength,
    RVector &data)
{
    char cbuf[256], tag[256];
    RVector adata;
    int i, len = nqm*nlambda;

    data.New (len);

    for (i = 0; i < nlambda; i++) {
        adata.Relink (data, i*nqm, nqm);

	switch (dscale) {
	case DATA_LIN:
  	    sprintf (tag, "DATA_REAL_WAVEL_%0.0f", wlength[i]);
	    if (!pp.GetString (tag, cbuf)) {
	        cout << "\nData file for " << tag << ":\n>> ";
		cin >> cbuf;
	    }
	    pp.PutString (tag, cbuf);
	    ReadDataFile (cbuf, adata);
	    break;
	case DATA_LOG:
  	    sprintf (tag, "DATA_MOD_WAVEL_%0.0f", wlength[i]);
	    if (!pp.GetString (tag, cbuf)) {
	        cout << "\nData file for " << tag << ":\n>> ";
		cin >> cbuf;
	    }
	    pp.PutString (tag, cbuf);
	    ReadDataFile (cbuf, adata);
	    break;
	}
    }
}

// ============================================================================

void SelectRefdata (DataScale dscale, int nqm, int nlambda,
    const RVector &wlength, RVector &data, bool &refEqualinitial)
{
    char cbuf[256], tag[256];
    RVector adata;
    int i, cmd, len = nqm*nlambda;

    data.New (len);

    for (i = 0; i < nlambda; i++) {
        adata.Relink (data, i*nqm, nqm);

	switch (dscale) {
	case DATA_LIN:
  	    sprintf (tag, "REFDATA_REAL_WAVEL_%0.0f", wlength[i]);
	    if (!pp.GetString (tag, cbuf)) {
	        cout << "\nReference data file for " << tag << ":\n>> ";
		cin >> cbuf;
	    }
	    pp.PutString (tag, cbuf);
	    ReadDataFile (cbuf, adata);
	    break;
	case DATA_LOG:
  	    sprintf (tag, "REFDATA_MOD_WAVEL_%0.0f", wlength[i]);
	    if (!pp.GetString (tag, cbuf)) {
	        cout << "\nReference data file for " << tag << ":\n>> ";
		cin >> cbuf;
	    }
	    pp.PutString (tag, cbuf);
	    ReadDataFile (cbuf, adata);
	    break;
	}
    }

    if (!pp.GetBool ("REFEQUALINITIAL", refEqualinitial)) {
        cout << "\nAre reference optical properties equal to initial "
	     << "guess\n[1|0] >> ";
	cin >> cmd;
	refEqualinitial = (cmd != 0);
    }
    pp.PutBool ("REFEQUALINITIAL", refEqualinitial);
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
	else if (!strcasecmp (cbuf, "LINLOG"))
	    ps = PARAM_SCALE_LINLOG, def = true;
	else if (!strcasecmp (cbuf, "LOG_BOUNDED"))
	    ps = PARAM_SCALE_BOUNDLOG, def = true;
    }
    while (!def) {
        int cmd;
	cout << "\nSelect parameter scaling method:\n";
	cout << "(0) None\n";
	cout << "(1) Initial parameter averages\n";
	cout << "(2) Log parameters\n";
	cout << "(3) Lin+log scaling\n";
	cout << "(4) Bounded log parameters\n";
	cout << "[0|1|2|3|4] >> ";
	cin >> cmd;
	switch (cmd) {
	case 0: ps = PARAM_SCALE_NONE,     def = true; break;
	case 1: ps = PARAM_SCALE_AVG,      def = true; break;
	case 2: ps = PARAM_SCALE_LOG,      def = true; break;
	case 3: ps = PARAM_SCALE_LINLOG,   def = true; break;
	case 4: ps = PARAM_SCALE_BOUNDLOG, def = true; break;
	}
    }
    switch (ps) {
    case PARAM_SCALE_NONE:     pp.PutString ("PARAM_SCALE", "NONE"); break;
    case PARAM_SCALE_AVG:      pp.PutString ("PARAM_SCALE", "AVG");  break;
    case PARAM_SCALE_LOG:      pp.PutString ("PARAM_SCALE", "LOG");  break;
    case PARAM_SCALE_LINLOG:   pp.PutString ("PARAM_SCALE", "LINLOG");  break;
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
	//else if (!strcasecmp (cbuf, "PGM"))     fmt = IMGFMT_PGM;
	//else if (!strcasecmp (cbuf, "PPM"))     fmt = IMGFMT_PPM;
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

static bool CheckRange (const MWsolution &sol)
{
    bool inrange = true;
    
#ifdef LIMIT_RANGE
    const double MIN_CMUA = 0;
    const double MAX_CMUA = 1.0;
    const double MIN_CKAPPA = 0;
    const double MAX_CKAPPA = 10;

    double vmin, vmax;
    int i;
    int nofwavel = sol.nofwavel;

    for (i = 0; i < nofwavel; i++) {

	sol.swsol[i]->Extents (OT_CMUA, vmin, vmax);
	if (vmin < MIN_CMUA || vmax > MAX_CMUA) {
	    cerr << "WARNING: " << vmin << " < CMUA < " << vmax
		 << " in trial solution" << endl;
	    inrange = false;
	}
	sol.swsol[i]->Extents (OT_CKAPPA, vmin, vmax);
	if (vmin < MIN_CKAPPA || vmax > MAX_CKAPPA) {
	    cerr << "WARNING: " << vmin << " < CKAPPA < " << vmax
		 << " in trial solution" << endl;
	    inrange = false;
	}
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
    RFwdSolverMW *fws = ofdata->fws;
    const Raster *raster = ofdata->raster;
    const Scaler *pscaler = ofdata->pscaler;
    MWsolution *meshsol = ofdata->meshsol;
    const Regularisation *reg = ofdata->reg;
    const RCompRowMatrix *qvec = ofdata->qvec;
    const RCompRowMatrix *mvec = ofdata->mvec;
    const RVector *data = ofdata->data;
    const RVector *sd = ofdata->sd;

    raster->Map_ActiveSolToMesh (pscaler->Unscale(x), *meshsol);
    meshsol->RegisterChange();
    if (!CheckRange (*meshsol)) return -1.0; // error

    RVector proj = fws->ProjectAll_wavel (*qvec, *mvec, *meshsol, 0);
    if (visnan (proj)) return -1.0; // error
    
    double fd = ObjectiveFunction::get_value (*data, proj, *sd);
    double fp = reg->GetValue (x);
    if (of_sub) {
	of_sub[0] = fd;
	of_sub[1] = fp;
    }
    return fd + fp;
}
