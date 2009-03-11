/*****************************************************************************
 * toast.cc             Martin Schweiger, Simon Arridge              2.09.96 *
 *                                                                           *
 * Time resolved optical absorption and scattering tomography                *
 * Reconstruction program                                                    *
 *                                                                           *
 * Version 2 - S.Arridge 12.09.96                                            *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 * Notes :                                                                   *
 * 1. Added stack to record history of error norm reductions                 *
 * 2. rearrange structure to create and destruct kernels at each NR step     *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#ifdef __BORLANDC__
#include <strstrea.h>
#include <conio.h>
#include <process.h>

typedef unsigned pid_t;
#else
#include <sstream>
#include <unistd.h>

#endif
#include <time.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include <toast.h>

using namespace toast;

#define TINY_CHANGE 1e-4
//#define DO_OWN_PARAM_SCALE
// ***************************************************************************
// local prototypes and definitions

char *ReconTypeString (ReconType recontp);
void SwapSolution (State **st1, Projection **pr1, State **st2,
    Projection **pr2);
void GetImageData (char *fname, RVector &params, ParameterType &which);
void DisplayInfo (void);
char Choice_MainMenu ();
void AskReconType (ReconType &recontp);
Basis *AskBasis (Mesh &mesh, bool *fixnode);
void AskSolverType (KernelType &kerneltp);
void AskGradientMethod (GradientEvalType &gradtype);
void AskIterations (int &niter);
void AskConvCriterion (double &convcrit);
void AskMTypeStatus (Environment &env, bool *active);

RVector *CreateParameterScaling(const Environment &env, const QMMesh &qmmesh, 
    const Basis *basis, State &state, ReconType recontp);
RDenseMatrix ***CreateDataFilter(const Environment &env, const QMMesh &qmmesh, 
    const Basis *basis, State &state, const SD& sd, ReconType recontp, 
    RVector *S, RVector** &W);

// ***************************************************************************
// global variables

char clear[]  = "\033[2J\033[0;0f";
char normal[] = "\033[0m";
char bold[]   = "\033[1m";

// ***************************************************************************
// main program

int main (int argc, char *argv[])
{
    time_t t1 = time(NULL);

    int i;
    Environment env;		// environment containing global options
    QMMesh qmmesh;		// FEM mesh with qm description
    Mesh mesh2;                 // FEM mesh for reconstruction *** SRA ***
    Basis *basis;               // solution basis
    State *base_state;          // pointer to `confirmed' state
    Projection *base_proj;      // associated projections
    SD *sd;                     // standard deviation
    Data *refdata;              // reference data
    bool *fixnode;              // list of fixed nodes
    int nfixed;                 // number of fixed nodes

// check expiry date

    CHECK_EXPIRED ();
    SetVersion ("15");

// command line interpretation

    env.ParseCmdLine (argc, argv);
    env.mesh = &qmmesh;
    env.mesh2 = &mesh2;         // *** SRA ***
    //env.verify ();
    if (env.INTERACTIVE)
	cout << clear << bold << "iTOAST - interactive TOAST" << normal
	     << endl << "Initialising ..." << endl;

#ifdef TOAST_PARALLEL
    cerr << "Starting thread pool\n";
    g_tpool = new ThreadPool (Task::GetThreadCount());
    cerr << "Finished thread pool\n";
    // global, defined in task.cc
#endif

// open log file and echo parameters

    LogfileOpen (env.logfile);
    logfile << env;
    logfile << endl << "====== TOAST log ======" << endl;
    logfile << VERSION_STRING;
    OutputProgramInfo (logfile);
    logfile << "=======================" << endl;

    LOGOUT_OPEN ("Initialisations");

// read mesh

    {
	LOGOUT_OPEN ("Read mesh");
	ifstream ifs (env.meshfile);
	ifs >> qmmesh;
	xASSERT (ifs.good(), Problem reading mesh file.);
	LOGOUT_1PRM ("%d elements", qmmesh.elen());
	LOGOUT_1PRM ("%d nodes", qmmesh.nlen());
	Point min, max;
	qmmesh.BoundingBox (min, max);
	char cbuf[256];
	std::ostringstream oss(cbuf);
	oss << "Bounding box: " << min << ' ' << max << '\0';
	LOGOUT(oss.str().c_str());
	LOGOUT_EXIT; // read mesh
    }

// rescale mesh

    if (env.MSCALE) {
	logfile << "Rescaling mesh by factor " << env.MeshScale << endl;
	qmmesh.ScaleMesh (env.MeshScale);
    }

// set mesh boundary type

    switch (env.BndCondTp) {
    case BC_DIRICHLET:
	LOGOUT("Using boundary condition DIRICHLET");
	qmmesh.BndType = MESH_DIRICHLET;
	for (i = 0; i < qmmesh.nlist.Len(); i++)
	    if (qmmesh.nlist[i].isBnd())
		qmmesh.nlist[i].SetBndTp (BND_DIRICHLET);
	break;
    case BC_ROBIN:
	LOGOUT("Using boundary condition ROBIN");
	qmmesh.BndType = MESH_ROBIN;
	for (i = 0; i < qmmesh.nlist.Len(); i++)
	    if (qmmesh.nlist[i].isBnd())
		qmmesh.nlist[i].SetBndTp (BND_ROBIN);
	break;
    }

// mesh setup: generate element matrices etc.

    LOGOUT_OPEN("Initialise mesh");
    qmmesh.Setup ();
    LOGOUT_EXIT; // intialise mesh

// read QM description

    LOGOUT_OPEN("Setup QM data");
    if (env.qmfile[0]) {
	ifstream ifs (env.qmfile);
	qmmesh.LoadQM (ifs);
    } else
	qmmesh.SetupRegularQM (env.nQnM);
    LOGOUT_2PRM("%d sources, %d detectors", qmmesh.nQ, qmmesh.nM);
    LOGOUT_1PRM("%d combinations", qmmesh.nQM);
    LOGOUT_EXIT; // setup qm

// set up data map
    //    cout << "Initialising Data Map\n";
    env.InitDataMap (qmmesh);
    //    cout << "Initialised Data Map\n";

// reset optical parameters of the mesh

    LOGOUT_OPEN ("Parameter reset");
    for (i = 0; i < 3; i++) {
	ParameterType prmtp;
	char pstr[10];
	switch (i) {
	    case 0: prmtp =  PRM_MUA, strcpy (pstr, "MUA"); break;
	    case 1: prmtp =  PRM_KAPPA, strcpy (pstr, "P2"); break;
	    case 2: prmtp =  PRM_N, strcpy (pstr, "N"); break;
	    default: prmtp = PRM_ZERO; break; // never reached
	}
	RVector params(qmmesh.plist.Len());
	switch (env.ResetTp[i]) {
	case RESET_MESH:
	    qmmesh.plist.Param (prmtp, params);
	    LOGOUT_1PRM ("Global: %s from mesh", pstr);
	    break;
	case RESET_HOMOG:
	case RESET_PROBE:
	    params = env.Reset[i];
	    if (i == P2 && !env.P2ResetIsKappa) prmtp = PRM_MUS;
	    LOGOUT_2PRM ("Global: %s %f", pstr, env.Reset[i]);
	    break;
	case RESET_NIM: {
	    ParameterType which;
	    GetImageData (env.ResetNIMfile[i], params, which);
	    if (prmtp == PRM_MUA)
	      xASSERT (which == PRM_MUA,
		  Image file does not match parameter type);
	    else if (prmtp == PRM_KAPPA)
	      xASSERT (which == PRM_KAPPA || which == PRM_CKAPPA ||
		  which == PRM_MUS, Image file does not match parameter type);
	    LOGOUT_2PRM ("Global: %s %s", pstr, env.ResetNIMfile[i]);
	    }
	    break;
	}
	for (int region = 0; region < MAXREG; region++) {
	    int nd;
	    if (!env.DoResetReg[region]) continue;
	    switch (env.ResetTpReg[region][i]) {
		case RESET_MESH:
		    for (nd = 0; nd < qmmesh.nlist.Len(); nd++)
			if (qmmesh.nlist[nd].Region() == region)
			    params[nd] = qmmesh.plist[nd].Param (prmtp);
		    LOGOUT_2PRM ("Reg. %d: %s from mesh", region, pstr);
		    break;
		case RESET_HOMOG:
		    for (nd = 0; nd < qmmesh.nlist.Len(); nd++)
			if (qmmesh.nlist[nd].Region() == region)
			    params[nd] = env.ResetReg[region][i];
		    //xASSERT (env.P2ResetIsKappaReg[region] ==
		    //qmmesh.elist.isKappa,
		    //P2 type mismatch during parameter reset.);
		    LOGOUT_3PRM ("Reg. %d: %s %f", region, pstr,
				 env.ResetReg[region][i]);
		    break;
	        case RESET_PROBE:
		    xERROR(RESET_PROBE not implemented for region reset);
		    break;
		case RESET_NIM:
		    xERROR(RESET_NIM not implemented for region reset);
		    break;
	    }
	}
	qmmesh.plist.SetParam (prmtp, params);
    }
    if (env.PARAM_SCALE) {
	env.InitParamScale(qmmesh);  // init scales from parameter averages
	LOGOUT_OPEN("Scale parameters");
	LOGOUT_1PRM("P1 scale: %f", env.ParamScale[0]);
	LOGOUT_1PRM("P2 scale: %f", env.ParamScale[1]);
	LOGOUT_1PRM("P3 scale: %f", env.ParamScale[2]);
	LOGOUT_EXIT;
    }
    LOGOUT_EXIT;  // Parameter reset

    //#define DUMP_PARAMS
#ifdef DUMP_PARAMS
    {
	ofstream ofs("param.dbg");
	ofs << qmmesh.plist << endl;
    }
#endif

// write out modified mesh if requested

    if (env.dbg_writemesh) {
        char cbuf[256];
	strcpy (cbuf, env.rootname);
	strcat (cbuf, "_reset.opt");
	ofstream ofs (cbuf);
	ofs << qmmesh << endl;
    }

// scan for fixed nodes

    LOGOUT_OPEN("Scan for fixed nodes");
    fixnode = new bool[qmmesh.nlist.Len()];
    for (i = 0; i < qmmesh.nlist.Len(); i++) {
        fixnode[i] = false;

	// fixed boundary?
	if (env.HOLD_BOUNDARY) {
	    if (qmmesh.nlist[i].isBnd()) {
	        fixnode[i] = true;
	    } else if (env.HOLD_BOUNDARY_WIDTH > 0) {
	        if (env.HOLD_BOUNDARY_WIDTH >= qmmesh.BoundaryDistance (i))
		    fixnode[i] = true;
	    }
	}
	// fixed regions?
	if (env.FixedRegion (qmmesh.nlist[i].Region()))
	    fixnode[i] = true;
    }
    for (i = nfixed = 0; i < qmmesh.nlist.Len(); i++)
        if (fixnode[i]) nfixed++;
    LOGOUT_1PRM ("Fixed nodes: %d", nfixed);
    LOGOUT_EXIT; // scan for fixed nodes

// set up inverse solution basis
    //    NodeBasis nb(qmmesh, fixnode);
    LOGOUT_OPEN("Setup solution basis");
    switch (env.BasisTp) {
    case BASIS_FWDMESH:
      //        cout << "Creating Node Basis\n";
	basis = new NodeBasis (qmmesh, fixnode);
	break;
    case BASIS_REGION:
	basis = new RegionBasis (qmmesh, fixnode);
	break;
    case BASIS_FUZZY:
        basis = new FuzzyBasis (qmmesh, fixnode, env.FuzzyBasisFile);
	break;
    case BASIS_PIXEL:
	basis = new PixelBasis (qmmesh, env.SolXdim, env.SolYdim, 1,
	    fixnode);
	break;
    case BASIS_PIXEL_3D:
	basis = new PixelBasis3D (qmmesh, env.SolXdim, env.SolYdim,
	    env.SolZdim, fixnode);
	break;
    case BASIS_CUBICPIXEL_2D:
        basis = new CubicPixelBasis2D (qmmesh, env.SolXdim, env.SolYdim,
	    fixnode);
	break;
    case BASIS_CUBICPIXEL_3D:
        basis = new CubicPixelBasis3D (qmmesh, env.SolXdim, env.SolYdim,
	    env.SolZdim, fixnode);
	break;
    case BASIS_FOURIER:
	basis = new FourierBasis (qmmesh, env.SolXdim, env.SolYdim, 1,
	    fixnode);
	break;
	//    case BASIS_WAVELET:
	//	basis = new WaveletBasis (qmmesh, env.SolXdim, env.SolYdim, 1,
	//	    fixnode);
	//	break;
    case BASIS_GAUSSBLOB:
        basis = new GaussBlobBasis (qmmesh, env.BlobSigma, env.BlobSup,
	    env.SolXdim, env.SolYdim, fixnode);
	break;
    case BASIS_BESSELBLOB:
        basis = new BesselBlobBasis (qmmesh, env.BlobSigma, env.BlobSup,
	    env.SolXdim, env.SolYdim, 1, fixnode);
	break;
    case BASIS_HANNINGBLOB:
        basis = new HanningBlobBasis (qmmesh, env.BlobSup,
	    env.SolXdim, env.SolYdim, fixnode);
	break;
    case BASIS_RAMPBLOB:
        basis = new RampBlobBasis (qmmesh, env.BlobSup,
	    env.SolXdim, env.SolYdim, fixnode);
	break;
    case BASIS_SPLINEBLOB:
        basis = new SplineBlobBasis (qmmesh, env.BlobSup,
	    env.SolXdim, env.SolYdim, fixnode);
	break;
    case BASIS_BESSELBLOB_3D:
        basis = new BesselBlobBasis3D (qmmesh, env.BlobSigma, env.BlobSup,
	    env.SolXdim, env.SolYdim, env.SolZdim, fixnode);
	break;
    case BASIS_SECONDMESH:
        // need to read in the second mesh, then set up basis
        // *** to be implemented ***
       // read mesh

        {
	LOGOUT_OPEN ("Read 2nd mesh");
	ifstream ifs (env.mesh2file);
	ifs >> mesh2;
	xASSERT (ifs.good(), Problem reading 2nd mesh file.);
	LOGOUT_1PRM ("%d elements", mesh2.elen());
	LOGOUT_1PRM ("%d nodes", mesh2.nlen());
	Point min, max;
	mesh2.BoundingBox (min, max);
	char cbuf[256];
	std::ostringstream oss(cbuf);
	oss << "Bounding box: " << min << ' ' << max;
        mesh2.Setup();         // possibly this is overkill 
	LOGOUT(oss.str().c_str());
	LOGOUT_EXIT;
        }
        basis = new SecondMeshBasis (qmmesh, mesh2, fixnode);
        break;
    default:
	xERROR(Unknown basis type);
	basis = 0;
	break;
    }
    LOGOUT_1PRM("Basis type: %s", basis->IdString());
    LOGOUT_1PRM("Total Size: %d", basis->Order());
    if (env.PRECALC_BASISMATRIX) {
        if (env.basisfile[0] && basis->LoadBasisMatrices (env.basisfile)) {
	    LOGOUT_1PRM("Read basis matrices from %s", env.basisfile);
	} else {
	    LOGOUT("Pre-calculating basis matrices ...");
	    double sup_avg = basis->PrecomputeBasisMatrices (true, true);
	    int ncd, ncr;
	    for (i = ncd = ncr = 0; i < basis->Order(); i++) {
	        const RGenericSparseMatrix *bsm = (env.ReconTp == RECON_KAPPA ?
		     basis->BasisFDD(i) : basis->BasisFFF(i));
		if (bsm) {
		    if (bsm->StorageType() == MATRIX_COORD) ncd++;
		    else                                    ncr++;
		}
	    }
	    LOGOUT_1PRM("Supported nodes per BF (avg): %f", sup_avg);
	    LOGOUT_2PRM("BF storage: %d coord, %d compressed row", ncd, ncr);
	    if (env.basisfile[0]) {
	        basis->SaveBasisMatrices (env.basisfile);
		LOGOUT_1PRM("Basis matrices written to %s", env.basisfile);
	    }
	}
    }
    LOGOUT_EXIT; // solution basis

// set up regularisation

    LOGOUT("Setting up regularisation ...");
    Regularisation *reg = Regularisation::Create (env, basis);

// generate neighbour list if required
#ifdef UNDEF
    if ((env.BasisTp == BASIS_FWDMESH && env.FilterTp == FILTER_MEDIAN) ||
	env.MARKOV == true) {
	LOGOUT("Generating element neighbour list ...");
	qmmesh.SetupNeighbourList ();
    }
#endif
// convert parameters to log if requested
#ifdef UNDEF
    if (env.LOG_PARAM) {
	logfile << "Transforming data to log (Experimental!) ..." << flush;
	for (i = 0; i < qmmesh.plist.Len(); i++) {
	    qmmesh.plist[i].SetMua (log (qmmesh.plist[i].Mua()));
	    qmmesh.plist[i].SetKappa (log (qmmesh.plist[i].Kappa()));
	}
	logfile << ", finished." << endl;
    }
#endif
// read data from disk

    LOGOUT("Reading data");
    Data data (env, qmmesh);
    data.Read ();
    if (env.USING_REFERENCE_DATA) {
	LOGOUT("Reading reference data");
	refdata = new Data(env, qmmesh, true);
	refdata->Read ();
    } else {
        refdata = 0; // not used
    }

// create a state and initialise it from the mesh as a default

    LOGOUT_OPEN("Setup states and updates");
    base_state = new State (env, qmmesh, basis);
    LOGOUT_EXIT; // states & updates

// initialise projections and standard deviation

    LOGOUT_OPEN("Setup projections and sd's");
    base_proj = new Projection (env, qmmesh);
    base_proj->Update (base_state->FWS);
    sd = new SD (env, qmmesh);
    if (env.READ_SD) {
        LOGOUT ("Reading SD");
	sd->Read();
    } else
        sd->SetDefaults (data);
    if (env.USING_REFERENCE_DATA) {
	data -= *refdata;
	data += *base_proj;
    }     
    // must do renormalise _after_ reference data, if required!
    if (env.ERROR_NORM_SCALE) sd->Renormalise (data, *base_proj);
    LOGOUT_EXIT; // projections & sds

// write out regularisation image for initial distribution
    if (env.dbg_output_regular && reg->Active()) {
        LOGOUT("Writing regularisation map to reg.nim");
	reg->WriteToImage (*base_state, "reg");
    }

// global parameter optimisation

    bool globalOptimise = false;
    RVector param(3);
    for (i = 0; i < 3; i++) {
	if (env.ResetTp[i] == RESET_PROBE) 
	    param[i] = env.Reset[i], globalOptimise = true;
	else param[i] = -1.0;
    }
    if (globalOptimise) {
        GlobalOptimise (env, qmmesh, data, *sd, param);
	base_state->InitFromMesh (qmmesh);
 	base_proj->Update (base_state->FWS);
	sd->Update (base_state->FWS);
	if (env.ERROR_NORM_SCALE) sd->Renormalise (data, *base_proj);
    }

    Update update (env, basis);

#ifdef CHECK_NORMS_HACK
    RVector tmp(data);
    cerr << "constructed tmp OK\n";
    tmp -= *base_proj;
        cerr << tmp << endl;
        cerr << *base_proj;
	cerr << *sd;
    cerr << "Initial error norm " << l2normsq(tmp) << endl;
    tmp /= *sd;
    cerr << "\tChi-squared " << l2normsq( tmp) << endl;
    cerr << "finished error\n";
#endif
    //sd->Update (base_state->FWS);
    //if (env.ERROR_NORM_SCALE) sd->Renormalise (data, *base_proj);

    // adjust Markov hyper-parameter
    if (reg->Active() && env.MARKOV == true) {
	if (env.ReconTp != RECON_KAPPA && env.markov_hyper_mua < 0.0) {
	    double mhp = reg->AdjustMarkovHyper (*base_state, *sd, MUA);
	    LOGOUT_1PRM("Markov hyper mua adjusted to %g", mhp);
	}
	if (env.ReconTp != RECON_MUA && env.markov_hyper_kappa < 0.0) {
	    double mhp = reg->AdjustMarkovHyper (*base_state, *sd, KAPPA);
	    LOGOUT_1PRM("Markov hyper kappa adjusted to %g", mhp);
	}
    }

    base_state->OpenImage ();
    if (env.dbg_output_prior) {
	logfile << "Writing initial guess to image file(s) ..." << flush;
	base_state->WriteImage ();
	logfile << ", finished." << endl;
    }

    base_state->RegisterChange();
    // this forces a map of parameters from basis to nodal, so we know
    // that the nodal distribution corresponds with basis distribution
    // even if initial image cannot be correctly represented in the basis

    // interactive solver loop

    bool terminate = false;
    int maxiter = env.MaxNR;
    char cmd;
    double convcrit = env.StopCriterion;
    GradientEvalType grad_method = env.GradientMethod;
    ReconType recontp = env.ReconTp;
    Solver *solver;
    RVector *IHess = 0;
    RDenseMatrix ***Cfilt =0;
    RVector **Wfilt = 0;
    if(env.FiltBackProp || env.ScaledGradient) {
      IHess = CreateParameterScaling(env, qmmesh, basis, *base_state, recontp);

      cerr << "Outputing IHess" << endl;
      Mesh *tmpmeshptr = env.mesh;
      int sdim = IHess->Dim();
      NodalImage nim (basis, tmpmeshptr, env.meshfile, false);
      nim.Open ("IHess");
      if  (recontp == RECON_BOTH) {
	  RVector g1(*(IHess),0,sdim/2);
	  RVector g2(*(IHess),sdim/2,sdim/2);
	  nim.Write(g1);
	  nim.Write(g2);
      }
      else
	  nim.Write (*(IHess));
      if(!env.ScaledGradient) {
         IHess = 0;
      }
      else {
        cerr << "Using non zero IHess\n";      
#ifndef DO_OWN_PARAM_SCALE
        cerr << "Deleting IHess\n";
        base_state->SetPscale(*IHess);
        delete IHess;
        IHess = 0;       // from now on, GSolver uses IHess inside state
#endif
      }
      Cfilt = CreateDataFilter(env, qmmesh, basis, *base_state, *sd, recontp, 
	      IHess, Wfilt);
      Wfilt = 0; // disable it !
    }

    if(!env.FiltBackProp)
      Cfilt = 0;
    else
      cerr << "Using data space filters\n";

    if (env.KernelTp == KERNEL_CG_FR ||
	env.KernelTp == KERNEL_CG_PR ||
	env.KernelTp == KERNEL_TRUNCATED_NEWTON ||
	env.KernelTp == KERNEL_BLOCK_KACZMARZ ||
	env.KernelTp == KERNEL_QUASI_NEWTON) {
      	solver = CreateGradientSolver (env, env.KernelTp, Cfilt, Wfilt);      
    } else {
	solver = new JacobianSolver (env, qmmesh, data, *sd, *reg);
    }

    LOGOUT_EXIT;  // Initialisations

    do {
	if (env.INTERACTIVE) cmd = Choice_MainMenu();
	else                 cmd = 'G', terminate = true;
	switch (cmd) {
	case '1':
	    AskReconType (recontp);
	    if (recontp != env.ReconTp) {
		LOGOUT_OPEN("Reset recon type");
		LOGOUT_1PRM("From: %s", ReconTypeString (env.ReconTp));
		LOGOUT_1PRM("To:   %s", ReconTypeString (recontp));
		env.ReconTp = recontp;
		base_state->ResetReconType (recontp, qmmesh);
		update.ResetReconType (recontp);
		base_state->OpenImage (env.rootname, true);
		LOGOUT_EXIT;
	    }
	    break;
	case '2': {
	    Basis *nbasis = AskBasis (qmmesh, fixnode);
	    LOGOUT_OPEN("Map basis");
	    LOGOUT_1PRM("Mapping from: %s", basis->IdString());
	    LOGOUT_1PRM("Mapping to:   %s", nbasis->IdString());
	    LOGOUT1("Mapping base state");
	    base_state->ResetBasis (nbasis, true);
	    LOGOUT1("Mapping update");
	    update.ResetBasis (nbasis, true);
	    delete basis;
	    basis = nbasis;
	    LOGOUT_EXIT;
	} break;
	case '3':
	    AskSolverType (env.KernelTp);
	    delete solver;
	    if (env.KernelTp == KERNEL_CG_FR ||
		env.KernelTp == KERNEL_CG_PR ||
		env.KernelTp == KERNEL_TRUNCATED_NEWTON ||
		env.KernelTp == KERNEL_QUASI_NEWTON) {
		solver = CreateGradientSolver (env, env.KernelTp);
	    } else {
		solver = new JacobianSolver (env, qmmesh, data, *sd, *reg);
	    }
	    break;
	case '4':
	    AskGradientMethod (grad_method);
	    break;
	case '5':
	    AskIterations (maxiter);
	    break;
	case '6':
	    AskConvCriterion (convcrit);
	    break;
	case '7':
	    {
		int i;
		bool *active = new bool[env.nmeasure()];
		for (i = 0; i < env.nmeasure(); i++)
		    active[i] = data.MeasureActive(i);
		AskMTypeStatus (env, active);
		for (i = 0; i < env.nmeasure(); i++) {
		    data.ActivateMeasure (i, active[i]);
		    base_proj->ActivateMeasure (i, active[i]);
		}
		delete []active;
	    }
	    break;
	case 'G': {
	    if (env.INTERACTIVE)
		cout << "Reconstruction running ..." << endl;
	    
	    ObjectiveFunction OF (env, qmmesh, basis, *base_state, data,
				  *sd, *reg);
	    OF.SetGradientMethod (grad_method);
	    solver->Solve (OF, *base_state, update, maxiter, convcrit, IHess);
	} break;
	case 'Q':
	    terminate = true;
	    break;
	}
    } while (!terminate);

    // clean-up

    delete solver;
    delete base_state;
    delete base_proj;
    delete basis;
    delete sd;
    delete reg;
    delete []fixnode;

    logfile << "Total CPU time: " << env.runtime() << endl;
#define DEBUG_TIMING
#ifdef DEBUG_TIMING
    extern double cgtime;
    logfile << "CPU time spent in CG solver: " << cgtime << endl;
#endif
    time_t t2 = time(NULL);
    LOGOUT_1PRM ("Wallclock timer: %0.2f", difftime (t2, t1));
    LOGOUT_1PRM ("Total linear solver iterations: %d", IterCount);

#ifdef __BORLANDC__
    while (!kbhit());
#endif
    return 0;
}

char *ReconTypeString (ReconType recontp)
{
    static char *str[3] = {"MUA", "KAPPA", "BOTH"};
    return str[recontp];
}

// ***************************************************************************

void SwapSolution (State **st1, Projection **pr1, State **st2,
    Projection **pr2)
{
    State *st_tmp = *st1; *st1 = *st2; *st2 = st_tmp;
    Projection *pr_tmp = *pr1; *pr1 = *pr2; *pr2 = pr_tmp;
}

// ***************************************************************************

void GetImageData (char *fname, RVector &params, ParameterType &which)
{
    char cbuf[256], ctype[10];
    int len, ok, i;
    double val;
    ifstream ifs (fname);
    do {
	ifs.getline (cbuf, 256);
	if (!strncmp (cbuf, "ImageSize", 9)) {
	    sscanf (cbuf+11, "%d", &len);
	    xASSERT (len >= params.Dim(),
		Image size of eim file too small for mesh.);
	} else if (!strncmp (cbuf, "SolutionType", 12)) {
	    sscanf (cbuf+14, "%s", ctype);
	    if (!strncmp (ctype, "MUA", 3))
		which = PRM_MUA;
	    else if (!strncmp (ctype, "MUS", 3))
	        which = PRM_MUS;
	    else if (!strncmp (ctype, "KAPPA", 5))
		which = PRM_KAPPA;
	    else if (!strncmp (ctype, "CKAPPA", 6))
	        which = PRM_CKAPPA;
	    else {
		cerr << "Unknown recon type in image file." << endl;
		exit (1);
	    }
	}
    } while (strncmp (cbuf, "EndHeader", 9));
    do {
	ifs.getline (cbuf, 256);
	if ((ok = (!strncmp (cbuf, "Image", 5)))) {
	    for (i = 0; i < len; i++) {
		ifs >> val;
		if (i < params.Dim()) params[i] = val;
	    }
	    ifs.getline (cbuf, 256);	// consume eol
	}
    } while (ok);
}

char Choice_MainMenu ()
{
    char choice;
    cout << clear << bold << "iTOAST - Main Menu" << normal << endl << endl;
    cout << "Modify parameters" << endl;
    cout << "(1) Recon parameters" << endl;
    cout << "(2) Recon basis" << endl;
    cout << "(3) Inverse solver" << endl;
    cout << "(4) Gradient method" << endl;
    cout << "(5) Number of iterations" << endl;
    cout << "(6) Convergence criterion" << endl;
    cout << "(7) Measurement types" << endl;
    cout << endl;
    cout << "(G) Go" << endl;
    cout << "(Q) Quit" << endl << endl;
    cout << "[1-7|G|Q] ? ";
    do {
	cin >> choice;
	choice = toupper(choice);
    } while (strchr ("1234567GQ", choice) == 0);
    return choice;
}

void AskReconType (ReconType &recontp)
{
    int choice;

    cout << clear << bold << "iTOAST - Select recon parameters" << normal
	 << endl << endl;
    cout << "1 = MUA" << endl;
    cout << "2 = KAPPA" << endl;
    cout << "3 = BOTH" << endl << endl;
    cout << "[1-3] ? ";
    do { cin >> choice; } while (choice < 1 || choice > 3);
    switch (choice) {
    case 1:
	recontp = RECON_MUA;
	break;
    case 2:
	recontp = RECON_KAPPA;
	break;
    case 3:
	recontp = RECON_BOTH;
	break;
    }
}

Basis *AskBasis (Mesh &mesh, bool *fixnode)
{
    int choice;
    int xdim, ydim, zdim;
    Basis *newbasis = 0;

    cout << clear << bold << "iTOAST - Select recon basis" << normal
	 << endl << endl;
    cout << "1 = Node basis" << endl;
    cout << "2 = Pixel basis" << endl;
    cout << "3 = Fourier basis" << endl << endl;
    cout << "[1-3] ? ";
    do { cin >> choice; } while (choice < 1 || choice > 3);
    switch (choice) {
    case 1:
	newbasis = new NodeBasis (mesh, fixnode);
	break;
    case 2:
	cout << "Grid dimension:" << endl;
	cout << "xdim = ";
	cin >> xdim;
	cout << "ydim = ";
	cin >> ydim;
	if (mesh.nlist[0].Dim() == 2)
	    newbasis = new PixelBasis (mesh, xdim, ydim, 1, fixnode);
	else {
	    cout << "zdim = ";
	    cin >> zdim;
	    newbasis = new PixelBasis3D (mesh, xdim, ydim, zdim, fixnode);
	}
	break;
    case 3:
	cout << "Number of terms:" << endl;
	cout << "xdim = ";
	cin >> xdim;
	cout << "ydim = ";
	cin >> ydim;
	newbasis = new FourierBasis (mesh, xdim, ydim, 1, fixnode);
	cout << "Apply Gaussian filter to Fourier coefficients (1|0) ? ";
	cin >> choice;
	if (choice == 1) {
	    double sigmax, sigmay;
	    cout << "Gauss-width x: ";
	    cin >> sigmax;
	    cout << "Gauss-width y: ";
	    cin >> sigmay;
	    ((FourierBasis*)newbasis)->SetAttenuation (sigmax, sigmay);
	}
	break;
    }
    return newbasis;
}

void AskSolverType (KernelType &kerneltp)
{
    char cmd[256];
    int choice;
    do {
	cout << clear << bold << "iTOAST - Select inverse solver" << normal
	     << endl << endl;
	cout << "Jacobian-based solvers" << endl;
	cout << "0 = SVD" << endl;
	cout << "1 = Gauss-Jordan" << endl;
	cout << "2 = ART" << endl;
	cout << "3 = Block-ART" << endl;
	cout << "4 = Augmented ART" << endl;
	cout << "5 = Augmented Block-ART" << endl << endl;
	cout << "Gradient-based solvers" << endl;
	cout << "6 = Conjugate gradient (Fletcher-Reeves)" << endl;
	cout << "7 = Conjugate gradient (Polak-Ribiere)" << endl;
	cout << "8 = Truncated Newton" << endl;
	cout << "9 = Quasi-Newton" << endl << endl;
	cout << "Current value   = " << (int)kerneltp << endl;
	cout << "New value [0-9] = ";
	cin.getline (cmd, 256);
	choice = cmd[0]-'0';
    } while (choice < 0 || choice > 9);
    kerneltp = (KernelType)choice;
}

void AskGradientMethod (GradientEvalType &gradtype)
{
    char cmd[256];
    int choice;
    do {
	cout << clear << bold << "iTOAST - Select gradient method" << normal
	     << endl << endl;
	cout << "1 = explicit perturbation" << endl;
	cout << "2 = basis system matrix" << endl;
	cout << "3 = adjoint PMDF" << endl;
	cout << "4 = adjoint direct" << endl << endl;
	cout << "Current value   = " << (gradtype+1) << endl;
	cout << "New value [1-4] = ";
	cin.getline (cmd, 256);
	choice = cmd[0]-'0';
    } while (choice < 1 || choice > 4);
    gradtype = (GradientEvalType)(choice-1);
}

void AskIterations (int &niter)
{
    cout << clear << bold << "iTOAST - Select number of iterations" << normal
	 << endl << endl;
    cout << "Current value = " << niter << endl;
    cout << "New value     = ";
    cin >> niter;
}

void AskConvCriterion (double &convcrit)
{
    cout << clear << bold << "iTOAST - Select relative convergence criterion"
	 << normal << endl << endl;
    cout << "Current value = " << convcrit << endl;
    cout << "New value     = ";
    cin  >> convcrit;
}

void AskMTypeStatus (Environment &env, bool *active)
{
    int choice;

    do {
	cout << clear << bold << "iTOAST - Activate/deactivate data sets"
	     << normal << endl << endl;
	for (int i = 0; i < env.nmeasure(); i++) {
	    cout << '(' << (i+1) << ")\t" << env.mtp_fname[i] << endl;
	    cout << '\t' << "[Status: " << bold
		 << (active[i] ? "ACTIVE" : "INACTIVE") << normal
		 << ']' << endl << endl;
	}
	cout << "Switch activation status for [1-" << (env.nmeasure())
	     << ", 0 to quit]: ";
	cin >> choice;
	if (choice > 0 && choice <= env.nmeasure())
	    active[choice-1] = !active[choice-1];
    } while (choice > 0);
}

RDenseMatrix ***CreateDataFilter(const Environment &env, const QMMesh &qmmesh, 
       const Basis *basis, State &state, const SD& sd, ReconType recontp, 
       RVector *S, RVector** &Wfilt)
{
  //  #define CHECK_INVERSE
    int q, i,t,k,nt = env.nmeasure(),nq=qmmesh.nQ,nm, sdim, ofs;
    double mu = 2;   // try multiplying diagonal for now !
    RDenseMatrix ***Cfilt = new RDenseMatrix ** [nt];
    Wfilt = new RVector * [nt];
    for (t = 0; t < nt ; t++) {
        Cfilt[t] = new RDenseMatrix * [nq];
	sdim = basis->Order() * (recontp == RECON_BOTH ? 2 : 1);
	Wfilt[t] = new RVector(sdim);
	*(Wfilt[t]) = 0;
    }
    for (t = 0; t < nt; t++) {
        cerr << "Filter for datatype " << *(env.mtp[t]) << endl;
        for (q = 0; q < qmmesh.nQ; q++) {
	  //	    cerr << "    source " << q << endl;
	    //	    cerr << "calling Mcpmdf\n";
	    Mpmdf J(env, qmmesh, basis, state, recontp, env.mtp[t], q);
	    //	    cerr << "making AAt\n";
	    ofs = env.DMap().TQindex (t, q);
	    if (env.XSCALE) {       // include s.d. in PMDF
	      for (k = 0; k < J.Dim(ROW); k++) { 
		//   cerr << "Dividing PMDF by " << sd[ofs+k] << endl;
		  J.SetRow (k, J.Row(k) / sd[ofs+k]);
              }
            }
#ifdef DO_OWN_PARAM_SCALE
	    if(S){       // parameter rescale
	      //	      cerr << "Parameter Rescaling in Filter\n";
	      for (k = 0; k < J.Dim(ROW); k++) 
		for(i = 0 ; i < J.Dim(COL); i++)
		     J(k,i) *= (*S)[i];
              
            }
#endif
	    RSymMatrix AAt = AAT (J);		// AAt is for sure symmetric
	    nm = AAt.Dim(ROW);
	    //	    cerr << "Add identity regulariser\n";
	    //	    double trace = 0.0;
       	    for (i = 0; i < nm; i++) { 
		AAt(i,i) *= mu;		// add diagonal regulariser
		//		trace += AAt[i][i];
	    }
	    //	    cout << "matrix trace = " << trace ;
	    /*     create a dummy matrix to test inversion procedures */
	    /*	    
	    for(k=0; k< nm; k++) 
	      for(i=0; i<=k; i++)  AAt[k][i] = (k==i?mu:0.0);
	    */
#ifdef CHECK_INVERSE
		    /* check ! */
	     RDenseMatrix tmp(nm,nm); 
	    for(k=0; k< nm; k++) 
	      for(i=0; i<=k; i++) tmp(k,i) = tmp(i,k)= AAt(k,i);
#endif  
	    /* consider storing Cholesky decomposition */
	    //	    cerr << "Cholesky decompose\n";
	    CHdecomp (AAt);			// Cholesky decompose
	    
	    //	    SymMatrix AI = inverse(AAt); // inverse not available..
	    Cfilt[t][q] = new RDenseMatrix(nm,nm);
	    // cerr << "inverting row : ";
	    for(k = 0; k < nm; k++) {      // inverse by substitution!
	      // cerr << k << " ";
	      RVector rhs(nm); rhs = 0.0; rhs[k] = 1.0;
	      RVector x =  CHsubst (AAt, rhs); // Cholesky substitute
	      Cfilt[t][q]->SetRow (k, x);
	    }
	    // cerr << endl;
	    //	    for (trace = 0.0, i = 0; i < nm; i++) { 
	    //		trace += (*(Cfilt[t][q]))[i][i];
	    //	    }
	    //cout << " inverse trace = " << trace << endl;

#ifdef CHECK_INVERSE
	    cerr << *(Cfilt[t][q]) * tmp << endl;
#endif
	    for( k = 0; k < nm ; k++)     // sum pmdfs for normalisation
#ifdef L2COLUMNNORM
	      for(i = 0; i < sdim; i++)
		(*(Wfilt[t]))[i] += sqr(J[k][i]);
#else
	      *(Wfilt[t]) += J.Row(k);
#endif
        }                                 // end loop on sources
	for(i = 0; i < sdim; i++)
#ifdef L2COLUMNNORM
	 (*(Wfilt[t]))[i] = 1.0/sqrt((*(Wfilt[t]))[i]);
#else
	(*(Wfilt[t]))[i] = 1.0/((*(Wfilt[t]))[i]);
#endif
	cerr << "Outputing Wfilt" << endl;
	Mesh *tmpmeshptr = env.mesh;
	NodalImage nim (basis, tmpmeshptr, env.meshfile, false);
	nim.Open ("wfilt");
	if  (recontp == RECON_BOTH) {
	  RVector g1(*(Wfilt[t]),0,sdim/2);
	  RVector g2(*(Wfilt[t]),sdim/2,sdim/2);
	  nim.Write(g1);
	  nim.Write(g2);
        }
	else
	  nim.Write (*(Wfilt[t]));
    }
    return Cfilt;
}

RVector *CreateParameterScaling(const Environment &env, const QMMesh &qmmesh, 
       const Basis *basis, State &state,  ReconType recontp)
{
  int ndim1  = basis->Order();
  int ndim = (recontp == RECON_BOTH ? 2*ndim1 : ndim1);
  int hh;
  RVector *IHess = new RVector (ndim);

  (*IHess) = 1.0;          // set to identity;

#define NORMAL_LOG
//#define RADIAL_HACK
//#define PMDF_SHAPE
#ifdef NORMAL_LOG
       for(hh = 0; hh < ndim; hh++)  
       (*IHess)[hh] *= state[hh];
       // MS 17.3.00
       //RVector prm(ndim);
       //state.GetParam (prm);
       //*IHess *= prm;
#endif
#ifdef PMDF_SHAPE
        int q,m;

	RVector wq(state.FWS.soldim());
	for (q = 0; q < qmmesh.nQ; q++) {
	    wq += state.FWS.qvec[q];
	} // end loop q
	RVector dphi = state.FWS.Solve_E (wq, 0);

	RVector wm(state.FWS.soldim());
	for (m = 0; m < qmmesh.nM; m++) {
	    wm += state.FWS.mvec[m];
	} // end loop m
	RVector wphi = state.FWS.Solve_E (wm, 0);

	Update Result (env, basis);
	if (env.ReconTp != RECON_KAPPA) {
	  Result.param[0] = basis->IntFG (dphi, wphi);
	}
	if (env.ReconTp != RECON_MUA) {
	  Result.param[1] = basis->IntGradFGradG (dphi, wphi);
	}
	cerr << "Outputing PMDF SCALE" << endl;
	Mesh *tmpmeshptr = env.mesh;
	NodalImage nim (basis, tmpmeshptr, env.meshfile, false);
	nim.Open ("pmdf_scale");
	if  (recontp == RECON_BOTH) {
	  RVector g1(Result,0,ndim1);
	  RVector g2(Result,ndim1,ndim1);
	  nim.Write(g1);
	  nim.Write(g2);
        }
	else
	  nim.Write (Result);
        for(hh = 0; hh < ndim; hh++)  
	 (*IHess)[hh] /= Result[hh];
#endif
#ifdef RADIAL_HACK

       Point bbmin, bbmax;
       cout << "mesh is " << qmmesh.Dimension() << "-dimensional\n";
       /*
       if(qmmesh.Dimension() == 3) {
	 bbmin.New(3); bbmax.New(3);
       }
       else {
	 bbmin.New(2); bbmax.New(2);
       }
       */
       qmmesh.BoundingBox (bbmin, bbmax);
       double size = sqrt(sqr(bbmax[0] - bbmin[0]) + sqr(bbmax[1] - bbmin[1]));
       cerr << "Hessian dimension " << ndim << endl;
       cerr << "Mesh dimension " << qmmesh.nlist.Len() << endl;
       cerr << "Mesh size " << size << endl;
       if(recontp != RECON_BOTH)
         for(int hh = 0; hh < ndim; hh++) {       // ad hoc radial scaling
	   //   cout << "calling basis->cdist2(hh) \n";
	   double d = basis->cdist2(hh);
	   // cout << hh << " " << d << endl;
	   (*IHess)[hh] *= exp(-sqr((3*d-size)/size));
         }
       else
         for(int hh = 0; hh < ndim/2; hh++) {       // ad hoc radial scaling
	   //cout << "calling basis->cdist2(hh) \n";
	   double d = basis->cdist2(hh);
	   // cout << hh << " " << d << endl;
	   (*IHess)[hh] *= exp(-sqr((3*d-size)/size));
	   (*IHess)[hh+ndim/2] *= exp(-sqr((3*d-size)/size));
         }
#endif
       return IHess;
}










