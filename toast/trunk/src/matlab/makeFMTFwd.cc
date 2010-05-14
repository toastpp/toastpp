// Interface:
// RH-1: mesh handle
// RH-2: basis handle
// RH-3: regularisation handle
// RH-4: parameter file name (string)
// LH-1: FMT solver handle

#include "mex.h"
#include "felib.h"
#include "stoastlib.h"
#include "util.h"
#include "projector.h"
#include "muaSolver.h"
#include "matrixFreeSolver.h"
#include "MLEMSolver.h"

#define MAXREGION 100

using namespace std;

// =========================================================================
// local prototypes

void positionSourcesAndDetectors(int nQ, int qprof, double qwidth,
    SourceMode srctp, QMMesh *mesh, RCompRowMatrix * qvec, RVector * norms);

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFMTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    
    const double c0 = 0.3;
    double tol = 1e-8;
    double lin_tol;
    SourceMode srctp = SRCMODE_NEUMANN;   // source type
    char meshname[256];
    char cbuf[256];
    double freq, omega;    // modulation frequency (MHz and cycles/ps)
    int qprof, mprof;      // source/measurement profile
                           // (0=point, 1=Gaussian, 2=cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]

    // get mesh pointer from handle
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar(prhs[0]));

    // open parameter file
    ParamParser pp;
    mxGetString (prhs[3], cbuf, 256);
    if (pp.Open (cbuf))
        mexPrintf ("Reading parameters from %s", cbuf);
    else
        mexErrMsgTxt ("Parameter file not found");

    int nQ, nM, nQM;
    RCompRowMatrix * qvec = new RCompRowMatrix;
    int i, j;
    double t0, t1;
    Point bbmin, bbmax;

    RFwdSolver *FWS = new RFwdSolver (pp); // memory leak
    lin_tol = FWS->GetLinSolverTol();

    //SelectMesh (pp, meshname, *mesh);
    mesh->MarkBoundary ();
    SelectSourceProfile (pp, qprof, qwidth, srctp);
    SelectMeasurementProfile (pp, mprof, mwidth);
    int n = mesh->nlen();
    int dim = mesh->Dimension();
    nQ = mesh->nQ;
    nM = mesh->nM;
    nQM = mesh->nQM;
    mesh->BoundingBox (bbmin, bbmax);
    cout << "BBox min = "<<bbmin<<" BBox max = "<<bbmax<<" diff = "<<(bbmax-bbmin)<<endl;

    //SelectData (pp, freq);
    //omega = freq * 2.0*Pi*1e-6;

    // build the source vectors 
    RVector * srcNorms = new RVector[nQ];
    RVector * detNorms = new RVector[nM];
    RVector * srcNormsP = 0, * detNormsP = 0;
    if (SelectNormals(pp, srcNorms, detNorms, nQ, nM))
    {
	    srcNormsP = srcNorms;
	    detNormsP = detNorms;
    }
    positionSourcesAndDetectors(nQ, qprof, qwidth, srctp, mesh,qvec,srcNormsP);
    delete[] srcNorms;

    // reset initial parameter estimates
    Solution sol(OT_NPARAM, n);
    SelectInitialParams (pp, *mesh, sol);

    // allocate system matrix
    cout << endl << "Allocating system matrix" << endl;
    FWS->Allocate (*mesh);

    cout << endl << endl << "----------" << endl;

    // solve for excitation fields
    cout << "Generating test data, phi_e and phi_f" << endl;

    Camera ** camPList = new Camera*[nQ];
    Projector ** projPList = new Projector*[nQ];
    SelectProjectors (pp, *mesh, projPList, camPList, detNormsP);
    delete[] detNorms;
    int imageW = camPList[0]->getImageWidth(), imageH = camPList[0]->getImageHeight();
    int nImagePts = imageW*imageH;
    
    // Get grid size
    Raster *rast = (Raster*)Handle2Ptr (mxGetScalar(prhs[1]));
    IVector gDim = rast->GDim();
    int gSize = rast->GLen();

    // Load in ct_data for prior, and labels to generate fluorescence
    IVector priorDim(3), labelDim(3);// priorDim[0]=328; priorDim[1]=888; priorDim[2]=197;
    char priorImgFile[200], labelImgFile[200];
    RVector * priorDataRaw=0, * priorData=0, * labelRaw=0, * labelData=0;
    if (SelectPriorImage(pp, priorImgFile, priorDim))
    {
	int priorSize = priorDim[0]*priorDim[1]*priorDim[2];
	priorDataRaw = new RVector(priorSize);
	priorData = new RVector(gSize);
	ReadData(priorImgFile/*"/cs/research/vision/home0/green/trudge/digimouse/ct_data/ct_clip_328x888x197.img"*/, *priorDataRaw);
	// Sub sample to grid resolution
	double rx = (double)(priorDim[0])/gDim[0], ry = (double)(priorDim[1])/gDim[1], rz = (double)(priorDim[2])/gDim[2];
	cout << "Down-sampling ct data by "<<rx<<", "<<ry<<", "<<rz<<endl;
	double ii=0.0, jj=0.0, kk=0.0;
	for (int k=0; kk<priorDim[2] && k<gDim[2]; ++k, kk+=rz)
	{
	    jj=0.0;
	    for (int j=0; jj<priorDim[1] && j<gDim[1]; ++j, jj+=ry)
	    {
		ii=0.0;
		for (int i=0; ii<priorDim[0] && i<gDim[0]; ++i, ii+=rx)
		{
		    int gridIdx = i + j*gDim[0] + k*gDim[0]*gDim[1];
		    int priorIdx = ii + (int)(jj+0.5)*priorDim[0] + (int)(kk+0.5)*priorDim[0]*priorDim[1];
		    (*priorData)[gridIdx] = (*priorDataRaw)[priorIdx];
		}
	    }
	}
	cout << "CTDataRaw max "<<vmax(*priorDataRaw) << " min " <<vmin(*priorDataRaw)<<endl;
	cout << "CTData max "<<vmax(*priorData) << " min " <<vmin(*priorData)<<endl;

	delete priorDataRaw;
    }

    if (SelectLabelImage(pp, labelImgFile, labelDim))
    {
	int labelSize = labelDim[0]*labelDim[1]*labelDim[2];
	labelRaw = new RVector(labelSize);
	labelData = new RVector(gSize);
	ReadData(labelImgFile /*"/home/trudge/digimouse/atlas/atlas_clip_328x888x197.img"*/, *labelRaw);
	// Sub sample to grid resolution
	double rx = (double)(labelDim[0])/gDim[0], ry = (double)(labelDim[1])/gDim[1], rz = (double)(labelDim[2])/gDim[2];
	cout << "Down-sampling label data by "<<rx<<", "<<ry<<", "<<rz<<endl;
	double ii=0.0, jj=0.0, kk=0.0;
	for (int k=0; kk<labelDim[2] && k<gDim[2]; ++k, kk+=rz)
	{
	    jj=0.0;
	    for (int j=0; jj<labelDim[1] && j<gDim[1]; ++j, jj+=ry)
	    {
		ii=0.0;
		for (int i=0; ii<labelDim[0] && i<gDim[0]; ++i, ii+=rx)
		{
		    int gridIdx = i + j*gDim[0] + k*gDim[0]*gDim[1];
		    int labelIdx = ii + (int)(jj+0.5)*labelDim[0] + (int)(kk+0.5)*labelDim[0]*labelDim[1];
		    (*labelData)[gridIdx] = (*labelRaw)[labelIdx];
		}
	    }
	}
	cout << "labelDataRaw max "<<vmax(*labelRaw) << " min " <<vmin(*labelRaw)<<endl;
	cout << "labelData max "<<vmax(*labelData) << " min " <<vmin(*labelData)<<endl;

	delete labelRaw;
    }

    // Generate fluorescence and absorption distribution
    // Create a sphere (r=5) of fluorophore at (0,10,0)
    RVector cnt1(3), cnt2(3), fluoData, fluo(gSize, 0.0), mMua(n, 0.01), mua(gSize), testFieldG(gSize);
    rast->Map_MeshToGrid(mMua, mua);
    cnt1[0]=0.0; cnt1[1]=10.0; cnt1[2]=0.0;
    cnt2[0]=0.0; cnt2[1]=-10.0; cnt2[2]=0.0;
    for (int k=0; k<gDim[2]; ++k)
    {
	for (int j=0; j<gDim[1]; ++j)
	{
	    for (int i=0; i<gDim[0]; ++i)
	    {
		RVector pos(3); pos[0]=i-gDim[0]/2; pos[1]=j-gDim[1]/2; pos[2]=k-gDim[2]/2;	
		RVector d1 = pos-cnt1, d2 = pos-cnt2;
		d1[2] = 0.0; d2[2] = 0.0;
		//if (l2norm(d)<5.0) fluo[i + j*gDim[0] + k*gDim[0]*gDim[1]]=1.0;
		int idx = i + j*gDim[0] + k*gDim[0]*gDim[1];
		testFieldG[idx] = l2norm(d1);
		if (labelData)
		{
		    if ((*labelData)[idx]==1)
		    {
			fluo[idx] = 1.0;
		    }
		    if ((*labelData)[idx]==18)
		    {
			mua[idx] = 0.01;
		    }
		}
		else
		{
		    // sphere radius 5 at (5,0,0)
		    if ( l2norm(d1)<5.0 || l2norm(d2)<5.0 /*|| (fabs(pos[1])<=20.0 && fabs(pos[2])<=10.0)*/) 
		    {
			if (fabs(pos[2])<10.0)
			{
			    fluo[idx]=1.0; //0.5639;
			    mua[idx]=0.01;
			}
		    }
		}
	    }
	}
    }
    if (labelData) delete labelData;
    RVector testFieldM;
    rast->Map_GridToMesh(testFieldG, testFieldM); 	


    cout << "Sum of fluo vector is "<<sum(fluo)<<endl;
    RVector sfluo;
    rast->Map_GridToSol(fluo, sfluo); 	
    rast->Map_GridToMesh(mua, mMua);
    cout << "Assembling and pre-processing system matrix" << endl;
    FWS->Reset (sol, 0.0 /*omega*/);

    // Setup regularisation
    Regularisation * regul = (Regularisation*)Handle2Ptr (mxGetScalar(prhs[2]));
    double tau = regul->GetTau();
    bool linReg = false; // TODO
    //SelectRegularisation (pp, *mesh, rast, &regul, linReg, priorData, tau);

    // Create solvers
    bool solveMua, solveFluo, muaLogSolve;
    int muaIters;
    double fluoTol;
    IVector dataWin(4);
    int dataWinW, dataWinH, dataWinSize;
    SelectSolvers (pp, solveMua, solveFluo, muaIters, muaLogSolve, fluoTol, dataWin);
    if (solveMua)
    {
	sol.SetParam(OT_CMUA, mMua);
    }
    MuaSolver muaSolver(*FWS, *mesh, regul, rast, nQ, *qvec, projPList, sol);
    FluoSolverType fst = SelectFluoSolver(pp);
    FluoSolver * fSolver;
    if (fst == FSOLVER_MF)
    {
	if (dataWin[1] == 0) dataWin[1] = imageW;
	if (dataWin[3] == 0) dataWin[3] = imageH;
	dataWinW = dataWin[1] - dataWin[0];
	dataWinH = dataWin[3] - dataWin[2];
	dataWinSize = dataWinW * dataWinH;
	fSolver = new MatrixFreeSolver(*FWS, *mesh, regul, tau, rast, nQ, *qvec,
				       projPList, dataWin);
	cout << "Created matrix-free solver" << endl;
    }else if (fst == FSOLVER_MLEM)
    {
	fSolver = new MLEMSolver(*FWS, *mesh, regul, rast, nQ, *qvec, projPList);
	cout << "Created MLEM solver" << endl;
    }



    plhs[0] = mxCreateScalarDouble (Ptr2Handle (fSolver));

    //for (int i = 0; i < nQ; i++) {
    //    delete camPList[i];
    //	  delete projPList[i];
    //}
    delete []camPList;
    delete []projPList;
    // note that these lists can be deleted here, because they have been
    // copied into the solver instances.

    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
}

// =========================================================================

void positionSourcesAndDetectors(int nQ, int qprof, double qwidth,
    SourceMode srctp, QMMesh *mesh, RCompRowMatrix * qvec, RVector * norms)
{ 
    int * bndellist, * bndsdlist;
    int nBndElems = mesh->BoundaryList (&bndellist, &bndsdlist);

    int nn = mesh->nlen();
    qvec->New (nQ, nn);
    //PinholeCamera * camList = new PinholeCamera[nQ];
    for (int i = 0; i < nQ; i++) {
	// Setup source vector
	RVector q(nn);
	cout << "Setting source at *mesh source "<<mesh->Q[i]<<endl;
	double dist;
	RVector n;
	if (norms)
	{
	    n=norms[i];
	}else{
	    n = -mesh->Q[i];
	}
	n /= l2norm(n);

	//mesh->Q[i] *= 2.0;
	RVector sourcePos = mesh->BndIntersect (mesh->Q[i] + l2norm(mesh->Q[i])*n, mesh->Q[i]); 
	//ProjectToBoundaryNode(*mesh, mesh->Q[i], n, nBndElems, bndellist, bndsdlist);
	cout << "Source intersects mesh at "<<sourcePos<<endl;

	switch (qprof) {
	    case 0:
		q = QVec_Point (*mesh, sourcePos, srctp);
		break;
	    case 1:
		q = QVec_Gaussian (*mesh, sourcePos, qwidth, srctp);
		break;
	    case 2:
		q = QVec_Cosine (*mesh, sourcePos, qwidth, srctp);
		break;
	}
	cout << "min(qvec["<<i<<"] = "<<vmin(q)<<endl;
	cout << "qvec row "<<i<<" has norm = "<<l2norm(q)<<endl;
	qvec->SetRow (i, q);
    }

    delete[] bndellist;
    delete[] bndsdlist;
}

// ============================================================================

#ifdef UNDEF



// =========================================================================
// MAIN 


int main (int argc, char *argv[]) {

    RVector reconMua(rast->SLen(), 0.001), gReconMua(rast->GLen(), 0.0);
    int iters = 0;
    if (solveMua)
    {
	cout <<"Solving for Mua..."<<endl;
	t0 = tic();
	muaSolver.Solve(excitData, reconMua, 0, 100, muaIters, 1e-8, muaLogSolve);
        rast->Map_SolToGrid(reconMua, gReconMua);
        rast->Map_GridToMesh(reconMua, mMua);
	sol.SetParam(OT_CMUA, mMua);
	cout <<"...done"<<endl;
	t1 = tic();
	cout << "Time: " << t1-t0 << " seconds" << endl << endl;
    }
//    FWS.Reset(sol, 0.0);    // reset, dont use recon mua for fluor recon

    // Solve for the test data into x 
    RVector sx(rast->SLen()), x(gSize);
    if (solveFluo)
    {
	cout << "Solving for fluorescence" << endl;
	t0 = tic();
	if (linReg)
	{
	    //sx = sfluo; // initial guess with correct distribution
	    iters = fSolver->Solve(data, sx, 0, 30, fluoTol);
	}else{
	    RVector diag(rast->SLen(), 1.0);
	    fSolver->fwdOperator(diag);
	    //diag += vmin(diag) + 1e4;
	    RPrecon_Diag * precon = (RPrecon_Diag *) (RPreconditioner::Create (PRECON_DIAG));
	    precon->ResetFromDiagonal(diag);
	    cout << "norm(diag) = " << l2norm(diag) << endl;
	    cout << "min(diag) = " << vmin(diag)<< ", max(diag) = " << vmax(diag) << endl;
	    cout << "min(abs(diag)) = "<< sqrtf(vmin(diag*diag)) << endl;
	    iters = fSolver->SolveNonLin(data, sx, 0 /*precon*/, 10, fluoTol);
	}
	rast->Map_SolToGrid(sx, x);	// map to full grid
	cout << "Number of iterations = "<<iters<<endl;
	t1 = tic();
	cout << "Time: " << t1-t0 << " seconds" << endl << endl;
	cout << "Max x = "<<vmax(x)<<", Min x = "<<vmin(x)<<endl;
    }

    // Generate comparison data from soln
    RVector soldata,solExcitData;
    fSolver->fwdOperator(sx, soldata);	// test data
    fSolver->getExcitationImages(solExcitData);

    double objFunc = l2normsq(soldata-data);
    if (regul)
    {
	// Check size of regularisation
	double regterm = regul->GetValue(sx);
	cout << "Regularisation term / obj func = "<<regterm/objFunc<<endl;
    }
    cout << "Obj func = "<<objFunc<<endl;
    cout << "NMSE = " << (objFunc/l2normsq(data)) <<endl;

    // output fields as NIM files
    int cmd=0;
//    cout << "Output nodal fields (1/0)? " << flush;
//    cin >> cmd;
    RVector projTestField;
    fSolver->testProj(projTestField);
    if (cmd) {
	RVector im(nImagePts), fld(gDim[0]*gDim[1]);
	IVector imgDim(2), flddim(2);
	double scale = 255/vmax(x);

//	if (priorData) WriteBinaryData(*priorData, "prior.dat");
	WriteBinaryData(x, "soln2.dat");
//	WriteBinaryData(fluo, "fluo.dat");
//	WriteBinaryData(gReconMua, "reconmua.dat");
//	WriteBinaryData(mua, "targetmua.dat");
	RVector logData = data;
	//logData.Clip(1e-25, vmax(data));
	//logData = log(logData);
	WriteBinaryData(soldata, "soldata.dat");
	WriteBinaryData(logData, "data.dat");
	WriteBinaryData(excitData, "excitData.dat");
	WriteBinaryData(simExcitData, "simExcitData.dat");
	WriteBinaryData(projTestField, "projTestField.nim");
	if (regul)
	{
	    RVector kappa, skappa = regul->GetKappa(sx);
	    rast->Map_SolToGrid(skappa, kappa);
	    WriteBinaryData(kappa, "kappa.dat");
	}
	imgDim[0] = dataWinW; imgDim[1] = dataWinH;
	flddim[0] = gDim[0]; flddim[1] = gDim[1];
	for (int i=0; i<gDim[2]; ++i)
	{
	    char fname[100];
	    sprintf(fname, "soln%d.pgm", i);
	    fld.Copy(x, 0, i*gDim[1]*gDim[0], gDim[1]*gDim[0]);
	    WritePGM(fld*scale, flddim, fname, false); 
	    //sprintf(fname, "fluo%d.pgm", i);
	    //fld.Copy(fluo, 0, i*gDim[1]*gDim[0], gDim[1]*gDim[0]);
	    //WritePGM(fld, flddim, fname, true); 
	}
	RVector * excitFields = fSolver->getExcitationFields();
	WriteBinaryData(excitFields[0], "excitField0.nim");
	for (int i=0; i<nQ; ++i)
	{
	    char fname[100];

	    sprintf(fname, "simexcit_img%d.pgm", i);
	    im.Copy(simExcitData, 0, i*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);

	    sprintf(fname, "logdata_img%d.pgm", i);
	    im.Copy(logData, 0, i*nImagePts, nImagePts);
	    //im.Clip(1e-25, 100.0);
	    cout << "max log data "<<vmax(im)<<" min log data "<<vmin(im)<<endl;
	    WritePGM(im, imgDim, fname);

	   /* sprintf(fname, "logsolndata_img%d.pgm", i);
	    im.Copy(soldata, 0, i*nImagePts, nImagePts);
	    im = im + vmin(im);
	    im.Clip(1e-25, 100.0);
	    cout << "max log soln "<<vmax(im)<<" min log soln "<<vmin(im)<<endl;
	    WritePGM(log(im), imgDim, fname);
*/
	    projPList[i]->projectFieldToImage(testFieldM, im);
	    sprintf(fname, "testfield_img%d.pgm", i);
	    WritePGM(im, imgDim, fname);    

	    sprintf(fname, "excit_data_img%d.pgm", i);
	    im.Copy(excitData, 0, i*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);
	    
	    sprintf(fname, "solnexcit_data_img%d.pgm", i);
	    im.Copy(solExcitData, 0, i*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);
	    
	    sprintf(fname, "data_img%d.pgm", i);
	    im.Copy(data, 0, i*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);

	    sprintf(fname, "solndata_img%d.pgm", i);
	    im.Copy(soldata, 0, i*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);

	    // Write excitation fields
	    sprintf(fname, "excitField%d.pgm", i);
	    WriteBinaryData(excitFields[i], fname);
	    
	}
    }


#ifdef TOAST_MPI
    MPI_Finalize();
#endif

    for (int i=0; i<nM; ++i)
    {
	delete projPList[i], camPList[i];
    }
    delete[] projPList;
    delete[] camPList;
    if (regul) delete regul;
    delete mesh, qvec;
    if (priorData)
	delete priorData;
    delete fSolver;
    return 0;
}                                                                              

// ============================================================================

#endif
