/*
---
Some utility functions for reading/writing data and parameters
from file or console input
---
*/

#include "util.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "camera.h"
#include "projector.h"
#ifdef MESA_SUPPORT
#include "GLProjector.h"
#endif

using namespace std;

const char *colormap = "../../scales/fire2.pal";

bool SelectLabelImage(ParamParser &pp, char * labelFname, IVector & gDim)
{
    if (!pp.GetString ("LABELFILE", labelFname)) 
    {
	return false;
    } 
    if (!pp.GetInt ("LABELSIZE_X", gDim[0])) 
    {
	cout << "\nLabel x size:\n>> ";
	cin >> gDim[0];
    }
    if (!pp.GetInt ("LABELSIZE_Y", gDim[1])) 
    {
	cout << "\nLabel y size:\n>> ";
	cin >> gDim[1];
    }
    if (!pp.GetInt ("LABELSIZE_Z", gDim[2])) 
    {
	cout << "\nLabel z size:\n>> ";
	cin >> gDim[2];
    }
    return true;
}

FluoSolverType SelectFluoSolver(ParamParser &pp)
{
    char fsolvername[100];
    if (!pp.GetString ("FSOLVER", fsolvername)) 
    {
	return FSOLVER_NONE;
    }
    if (!strcasecmp (fsolvername, "MATRIXFREE")) {
	return FSOLVER_MF;
    }
    else if (!strcasecmp (fsolvername, "MLEM")) {
	return FSOLVER_MLEM;
    }
    else {
	return FSOLVER_NONE;
    }
}

bool SelectPriorImage(ParamParser &pp, char* priorFname, IVector & gDim)
{
    if (!pp.GetString ("PRIORFILE", priorFname)) 
    {
	return false;
    } 
    if (!pp.GetInt ("PRIORSIZE_X", gDim[0])) 
    {
	cout << "\nPrior x size:\n>> ";
	cin >> gDim[0];
    }
    if (!pp.GetInt ("PRIORSIZE_Y", gDim[1])) 
    {
	cout << "\nPrior y size:\n>> ";
	cin >> gDim[1];
    }
    if (!pp.GetInt ("PRIORSIZE_Z", gDim[2])) 
    {
	cout << "\nPrior z size:\n>> ";
	cin >> gDim[2];
    }
    return true;
}

typedef enum { CAMTYPE_PINHOLE, CAMTYPE_ORTHO, CAMTYPE_NONE } CameraType;
typedef enum { PROJTYPE_MESAGL, PROJTYPE_NONE } ProjectorType;
void SelectProjectors (ParamParser &pp, QMMesh & mesh, Projector ** projPtrs, Camera ** camPtrs, RVector * norms)
{
    // Get camera type
    CameraType ctp;
    char str[100];
    if (!pp.GetString ("CAMERATYPE", str)) {
	cout << "\nCamera type:\n>> ";
	cin >> str;
    }
    double f;
    if (!strcasecmp (str, "PINHOLE")) {
	ctp = CAMTYPE_PINHOLE;
	// Get focal length
	if (!pp.GetReal ("FLEN", f)) {
	    cout << "\nCamera focal length:\n>> ";
	    cin >> f;
	}
    }else{
	if (!strcasecmp (str, "ORTHO")) {
	    ctp = CAMTYPE_ORTHO;
	}
	else{
	    cout<<"Camera type "<<str<<" not supported";
	    return;
	}
    }
    // Get projector type
    ProjectorType projtp;
    if (!pp.GetString ("PROJECTORTYPE", str)) {
	cout << "\nProjector type:\n>> ";
	cin >> str;
    }
    if (!strcasecmp (str, "MESAGL")) {
	projtp = PROJTYPE_MESAGL;
    }else{
	cout<<"Projector type "<<str<<" not supported";
	return;
    }

    // Get image dimensions
    int w,h;
    if (!pp.GetInt ("IMAGEW", w)) {
	cout << "\nCamera image width:\n>> ";
	cin >> w;
    }
    if (!pp.GetInt ("IMAGEH", h)) {
	cout << "\nCamera image height:\n>> ";
	cin >> h;
    }
    double pixSize;
    if (!pp.GetReal ("PIXELSIZE", pixSize)) {
	cout << "\nPixel size:\n>> ";
	cin >> h;
    }

    // Get gantry axis
    RVector axis(3);
    if (    (!pp.GetReal ("GANTRYAXIS_X", axis[0])) ||
	    (!pp.GetReal ("GANTRYAXIS_Y", axis[1])) ||
	    (!pp.GetReal ("GANTRYAXIS_Z", axis[2])) )
    {
	cout << "\nGantry rotation axis (i j k):\n>> ";
	cin >> axis;
    }
    axis /= l2norm(axis);

    // Get gantry radius
    double gantryR;
    if (!pp.GetReal ("GANTRYRADIUS", gantryR))
    {
	cout << "\nGantry radius:\n>> ";
	cin >> gantryR;
    }

    // Create cameras - assume:
    //	    1. Detector positions are set in mesh
    //	    2. Each detector looks at the axis
    //	    3. camera pointer array is allocated
    int * bndellist, * bndsdlist;
    int nBndElems=0;
 //   if (projtp==PROJTYPE_MESAGL)
 //   {
    nBndElems = mesh.BoundaryList (&bndellist, &bndsdlist);
  //  }

    // Set up shape function integral matrix
/*    RCompRowMatrix B(mesh.nlist.Len(), mesh.nlist.Len());
    for (int be=0; be<nBndElems; ++be)
    {
	int elemIdx = bndellist[be];
	int sideIdx = bndsdlist[be];
	Element * el = mesh.elist[elemIdx];
	int nSideNodes = el->nSideNode(sideIdx);
	RVector ei = el->BndIntF();
	for (int sn=0; sn<nSideNodes; ++sn)
	{
	    int elemNodeIdx = el->SideNode(sideIdx, sn);
	    int meshNodeIdx = el->Node[elemNodeIdx];
	    // Sum integrals of shape functions for all elements node is in
	    B(meshNodeIdx, meshNodeIdx) += ei[sideIdx];
	}
    }
*/
    // Create projectors
    for (int m=0; m<mesh.nM; ++m)
    {
	// Setup camera
	Point cpos = mesh.M[m]; 
	RVector z;
	if (norms)
	{
	    z = norms[m];
	}else{
	    z = dot(cpos, axis)*axis - cpos; // Look at axis
	}
	cpos = cpos + z * (1.0 - gantryR / l2norm(z)); // force camera onto gantry radius
	RVector x(3), y(3);
	x = axis;   // image x-axis is along rotation axis of gantry
	y = cross3(x, z);  // image y-axis is perp. to x and z
	y /= l2norm(y);
	z /= l2norm(z);
	cpos = cpos - y * 2.1428; // manual adjustment for optical axis offset in image
	cout <<"Cam pos = "<<cpos<<endl;
	cout << "x = "<<x<<" y= "<<y<<" z = "<<z<<endl;

	// Create camera
	if (ctp==CAMTYPE_PINHOLE)
	{
	    camPtrs[m] = new PinholeCamera(w, h, f, cpos, x, y, z);
	}else{
	    if (ctp==CAMTYPE_ORTHO)
	    {
		camPtrs[m] = new OrthoCamera(w, h, pixSize, cpos, x, y, z);
	    }
	}
	if (projtp==PROJTYPE_MESAGL) {
#ifdef MESA_SUPPORT
	    // note: currently only GLProjector is implemented, so MESA_SUPPORT
	    // is required
	    projPtrs[m] = new GLProjector(camPtrs[m], &mesh, nBndElems,
					  bndellist, bndsdlist);
#else
	    xERROR("Mesa-3D support required but not provided.");
#endif
	} else {
	    xERROR("Unsupported projector type requested.");
	}
    }
}

RVector ProjectToBoundaryNode(const QMMesh & mesh, const RVector & pos, const RVector & ray, int nBndElems, int * bndellist, int * bndsdlist)
{
    //double dmin = 1000000.0;
    //int idx=-1;
    RVector sourcePos;
 /*   for (int nidx=0; nidx<mesh.nlen(); ++nidx)
    {
        Node & n = mesh.nlist[nidx];
	if (n.isBnd())
	{
	    RVector n_minus_pos = n - pos;
	    double dperp = l2norm(cross3(n_minus_pos, ray));
	    double dparal = dot(n_minus_pos, ray);
	    double cos = dot(n,pos);
	    if (dperp<dmin && cos>0.0)
	    {
		dmin = dperp;
		idx = nidx;
	    }
	}
    }*/


    // Find intersection with mesh
    RDenseMatrix jacMat; // not used, just for argument
    for (int e=0; e<nBndElems; ++e)
    {
	int eidx = bndellist[e];
	Element * elem = mesh.elist[eidx];
	int side = bndsdlist[e];
	RVector normal = elem->DirectionCosine(side, jacMat);
	RVector elemPos = mesh.nlist[elem->Node[elem->SideNode(side, 0)]];
	Point intersect = pos + ray * dot(elemPos-pos, normal) / dot(ray, normal);
	RVector shapeFun = elem->GlobalShapeF(mesh.nlist, intersect);

	if (dot(normal, ray)<0.0)
	{
	    if ( (shapeFun[0]>=0.0) && (shapeFun[1]>=0.0) && (shapeFun[2]>=0.0) && (shapeFun[3]>=0.0))
	    {
		sourcePos = intersect;
	    }
	}
    }

    return sourcePos;
}

/* Set ring of sources about axis, distance d along axis from origin 
void SetRingSources(QMMesh & mesh, int nqm, RVector axis, double d)
{
    Point *newQ = new Point[nqm];
    Point *newM = new Point[nqm];
    Point cnt(dim);
    double radius = 50.0;// Size (&cnt);
    for (int i = 0; i < nqm; i++) {
	newQ[i].New (dim);
	newM[i].New (dim);
	double angle = (2.0 * Pi * i) / nqm;
	newQ[i][0] = cnt[0] + radius * cos (angle);
	newQ[i][1] = cnt[1] + radius * sin (angle);
	newQ[i][2] = cnt[2];
	newM[i][0] = cnt[0] + radius * cos (angle + Pi);
	newM[i][1] = cnt[1] + radius * sin (angle + Pi);
	newM[i][2] = cnt[2];
	ProjectToBoundaryNode (mesh, newQ[i], newQ[i]);
    }
    mesh.SetupQM (newQ, nqm, newM, nqm);
}*/

double SelectNoise(ParamParser &pp)
{
    double noise;
    if (!pp.GetReal ("NOISE", noise)) {
	cout << "\nSimulation noise :\n>> ";
	cin >> noise;
    }
    return noise;
}

void SelectSolvers (ParamParser &pp, bool & solveMua, bool & solveFluo, int & iters, bool & logSolve, double & fluoTol, IVector & dataWin)
{
    if (pp.GetBool ("MUASOLVE", solveMua) && solveMua)
    {
	if (!pp.GetInt ("MUAITERS", iters)) {
	    cout << "\nMua solver number of iterations:\n>> ";
	    cin >> iters;
	}	
	if (!pp.GetBool ("MUALOGSOLVE", logSolve)) {
	    cout << "\nSolve for log(mua)? (1=yes, 0=no):\n>> ";
	    cin >> logSolve;
	}
    }

    pp.GetBool("FLUOSOLVE", solveFluo);
    pp.GetReal("FSOLVERTOL", fluoTol);

    pp.GetInt("DATAWIN_XMIN", dataWin[0]);
    pp.GetInt("DATAWIN_XMAX", dataWin[1]);
    pp.GetInt("DATAWIN_YMIN", dataWin[2]);
    pp.GetInt("DATAWIN_YMAX", dataWin[3]);

}

void SelectGrid (ParamParser &pp, QMMesh &mesh, 
			    Raster ** rast)
{
    IVector gDim(3);
    RDenseMatrix bb(2,3), * bbp;
    double val;

    bbp = &bb;

    if (!pp.GetInt ("GRIDSIZE_X", gDim[0])) {
	bbp = 0;
    }
    if (!pp.GetInt ("GRIDSIZE_Y", gDim[1])) {
	bbp = 0;
    }
    if (!pp.GetInt ("GRIDSIZE_Z", gDim[2])) {
	bbp = 0;
    }
    if (!pp.GetReal ("GRIDSIZE_BBMINX", val)) {
	bbp = 0;
    }
    bb.Set(0,0,val);
    if (!pp.GetReal ("GRIDSIZE_BBMINY", val)) {
	bbp = 0;
    }
    bb.Set(0,1,val);
    if (!pp.GetReal ("GRIDSIZE_BBMINZ", val)){
	bbp = 0;
    }
    bb.Set(0,2,val);
    if (!pp.GetReal ("GRIDSIZE_BBMAXX", val)) {
	bbp = 0;
    }
    bb.Set(1,0,val);
    if (!pp.GetReal ("GRIDSIZE_BBMAXY", val)){
	bbp = 0;
    }
    bb.Set(1,1,val);
    if (!pp.GetReal ("GRIDSIZE_BBMAXZ", val)) {
	bbp = 0;
    }
    bb.Set(1,2,val);
    if (!bbp)
    {
	cout << "No grid bounding box specified, using full mesh bounding box" << endl;
    }
    *rast = new Raster_Pixel(gDim, gDim, (Mesh*) &mesh, bbp);

}
			    
typedef enum {REGTYPE_TK1SIGMA, REGTYPE_TK1, REGTYPE_TV, REGTYPE_HUBER, REGTYPE_NONE} RegType;
void SelectRegularisation (ParamParser &pp, QMMesh &mesh, 
			    Raster * rast, Regularisation ** reg, bool & linReg,
			    const RVector * kref, double & tau)
{
    char cbuf[256];

    RegType rtp;
    if (!pp.GetString ("REGTYPE", cbuf)) {
	cout << "\nRegularisation type:\n>> ";
	cin >> cbuf;
    }
    if (!strcasecmp (cbuf, "TK1SIGMA")) {
	rtp = REGTYPE_TK1SIGMA;
    }else{
	if (!strcasecmp (cbuf, "TK1")) {
	    rtp = REGTYPE_TK1;
	}
    }
    if (!strcasecmp (cbuf, "HUBER")) {
	rtp = REGTYPE_HUBER;
    }
    if (!strcasecmp (cbuf, "TV")) {
	rtp = REGTYPE_TV;
    }
    if (!strcasecmp (cbuf, "NONE")) {
	rtp = REGTYPE_NONE;
    }

    bool nonLinReg = true;
    if (!pp.GetBool ("NONLINREG", nonLinReg)) {
	cout << "\nUsing linear regularisation";
	nonLinReg = false;
    }
    linReg = !nonLinReg;

    double fT=0.0, beta=0.0, eps=0.0, sd=3.0, sdr=3.0;
    if (!pp.GetReal ("REGTAU", tau)) {
	//cout << "\nRegularisation tau:\n>> ";
	//cin >> tau;
    }
    if (rtp != REGTYPE_NONE)
    {
	if (!pp.GetReal ("REGFT", fT)) {
	    //cout << "\nRegularisation fT:\n>> ";
	    //cin >> fT;
	}
	if (!pp.GetReal ("REGBETA", beta)) {
	    //cout << "\nRegularisation beta:\n>> ";
	    //cin >> beta;
	}
	if (!pp.GetReal ("REGEPSILON", eps)) {
	    //cout << "\nRegularisation epsilon:\n>> ";
	    //cin >> eps;
	}
 	if (!pp.GetReal ("REGSD", sd)) {
	    //cout << "\nRegularisation smoothing diameter:\n>> ";
	    //cin >> eps;
	    cout << "Using default smoothing diameter of 3" << endl;
	}
	if (!pp.GetReal ("REGREFSD", sdr)) {
	    //cout << "\nRegularisation ref image smoothing diameter:\n>> ";
	    //cin >> eps;
	    cout << "Using default ref image smoothing diameter of 3" << endl;
	}
   }

    RVector x0(rast->SLen());
    switch (rtp)
    {
	case REGTYPE_TK1SIGMA:
	    if (kref)
	    {
		*reg = new TK1Sigma (tau, sd, &x0,
				    rast, *kref, sdr, fT, true, false);
	    }else{
		cerr << "No prior data for TK1Sigma regularisation"<<endl;
		*reg = 0;
	    }
	    break;
	case REGTYPE_TK1:
	    *reg = new TK1(tau, &x0, rast);
	    break;
	case REGTYPE_TV:
	    *reg = new TV(tau, beta, &x0, rast);
	    break;
	case REGTYPE_HUBER:
	    *reg = new Huber(tau, eps, &x0, rast);
	    break;
	case REGTYPE_NONE:
	    *reg = 0;
	    break;
	default:
	    cerr << "Unknown regularisation type "<<cbuf<<endl;
	    *reg = 0;
	    break;
    }
}

bool SelectNormals(ParamParser &pp, RVector * srcNorms, RVector * detNorms, int nQ, int nM)
{
    char snname[256];

    if (!pp.GetString ("NORMALSFILE", snname)) {
        cout << "\nNo source normal vector file name: pointing sources at axis\n>> ";
	return false;
    }
    cout << "Loading normals file: "<<snname<<endl;
    ifstream ifs (snname);

    for (int i=0; i<nQ; ++i)
    {
	ifs >> srcNorms[i]; 
    }
    for (int i=0; i<nM; ++i)
    {
	ifs >> detNorms[i]; 
    }

    return true;
}

void SelectMesh (ParamParser &pp, char *meshname, QMMesh &mesh)
{
    char qmname[256];

    if (!pp.GetString ("MESHFILE", meshname)) {
        cout << "\nRegistered mesh file name:\n>> ";
	cin >> meshname;
    }
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();

    if (!pp.GetString ("QMFILE", qmname)) {
        cout << "\nQM file name:\n>> ";
	cin >> qmname;
    }
    cout << "Loading QM file: "<<qmname<<endl;
    ifstream qmf (qmname);
    mesh.LoadQM (qmf);
//    mesh.SetupRegularQM_XZ(12);

    // write back
    pp.PutString ("MESHFILE", meshname);
    pp.PutString ("QMFILE", qmname);
}

// ============================================================================

bool SelectDataImage(ParamParser &pp, RVector & fluoData, RVector & excitData, int imageX, int imageY, int nImages)
{
    char cbuf1[256], cbuf2[256];
    if (pp.GetString ("FLUODATAIMAGEFILE", cbuf1)) 
    {
	// load fluorescence data
	int size = imageX*imageY*nImages;
	fluoData.New(size);	
	ReadDoubleData(cbuf1, fluoData);

	if (pp.GetString ("EXCITDATAIMAGEFILE", cbuf2)) 
	{   
	    // load excitation data
	    excitData.New(size);
	    ReadDoubleData(cbuf2, excitData);
	    return true;
	}
    }
    return false;
}

void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp)
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

void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth)
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

void SelectData (ParamParser &pp, double &freq)
{
    if (!pp.GetReal ("FREQ", freq) || freq < 0.0) do {
	cout << "\nSource modulation frequency [MHz]:\n>> ";
	cin >> freq;
    } while (freq < 0.0);
    pp.PutReal ("FREQ", freq);
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

void SelectInitialParams (ParamParser &pp, const Mesh &mesh, Solution &msol)
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
	    cout << "(1) Use values stored in mesh\n";
	    cout << "(2) Global homogeneous\n";
	    cout << "(3) Homogeneous in regions\n";
	    cout << "(4) Nodal image file (NIM)\n";
	    cout << "[1|2|3|4] >> ";
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
    msol.SetParam (OT_CMUA,   param[0]*c0/param[2]);
    msol.SetParam (OT_CKAPPA, c0/(3.0*param[2]*(param[0]+param[1])));
    msol.SetParam (OT_N, param[2]);
    for (i = 0; i < param[OT_C2A].Dim(); i++)
	param[OT_C2A][i] = c0/(2*param[2][i]*A_Keijzer(param[OT_C2A][i]));
    msol.SetParam (OT_C2A, param[OT_C2A]);
}

void OpenNIM (const char *nimname, const char *meshname, int size)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = N/A" << endl;
    ofs << "ImageSize = " << size << endl;
    ofs << "EndHeader" << endl;
}

void WriteNIM (const char *nimname, const RVector &img, int size, int no)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << no << endl;
    for (int i = 0; i < size; i++)
        ofs << img[i] << ' ';
    ofs << endl;
}

bool ReadNim (char *nimname, RVector &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (nimname);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    do {
        ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
        ifs >> img[i];
    return true;
}

void WritePGM (const RVector &img, const IVector &gdim, char *fname, bool doScale)
{
    int i, ii, dim = gdim[0]*gdim[1];
    double imgmin = 1e100, imgmax = -1e100;
    unsigned char *pixmap = new unsigned char[dim];

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    for (i = 0; i < dim; i++) {
        if (img[i+ii] < imgmin) imgmin = img[i+ii];
	if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    double scale = 1.0;
    if (doScale) scale = 256.0/(imgmax-imgmin);
    for (i = 0; i < dim; i++) {
        int v = (int)((img[i+ii]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P5" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++) ofs << pixmap[i];
    ofs << endl;
    delete []pixmap;
}

void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname)
{
    typedef struct {
        unsigned char r,g,b;
    } RGB;

    int i, ii, dim = gdim[0]*gdim[1];
    unsigned char *pixmap = new unsigned char[dim];
    double imgmin, imgmax, scale;

    static RGB colmap[256];
    static bool have_colmap = false;

    if (!have_colmap) {
        int r, g, b;
        ifstream ifs (colormap);
	for (i = 0; i < 256; i++) {
	    ifs >> r >> g >> b;
	    colmap[i].r = (unsigned char)r;
	    colmap[i].g = (unsigned char)g;
	    colmap[i].b = (unsigned char)b;
	}
	have_colmap = true;
    }

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    // rescale image
    if (scalemin) {
        imgmin = *scalemin;
    } else {
	for (i = 0, imgmin = 1e100; i < dim; i++)
	    if (img[i+ii] < imgmin) imgmin = img[i+ii];
    }
    if (scalemax) {
        imgmax = *scalemax;
    } else {
      for (i = 0, imgmax = -1e100; i < dim; i++)
	  if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    scale = 256.0/(imgmax-imgmin);

    for (i = 0; i < dim; i++) {
        int v = (int)((img[i+ii]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P6" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++)
        ofs << colmap[pixmap[i]].r
	    << colmap[pixmap[i]].g
	    << colmap[pixmap[i]].b;
    ofs << endl;

    delete []pixmap;
}

void ReadDoubleData (char *fname, RVector &data)
{
    int i, n = data.Dim();
    ifstream ifs (fname);
    if (!ifs.is_open())
    {
	cout<<"Cannot open file "<<fname<<endl;
	return;
    }
    int nn = n*sizeof(double);
    char * cbuf = new char[nn];
    ifs.read(cbuf, nn);
    for (i = 0; i < n; i++)
    {
	double * valp = (double*)(cbuf + i*sizeof(double));
	data[i] = *valp;
    }
    ifs.close();
    delete[] cbuf;
}

void ReadData (char *fname, RVector &data)
{
    int i, n = data.Dim();
    ifstream ifs (fname);
    char * cbuf = new char[n];
    ifs.read(cbuf, n);
    for (i = 0; i < n; i++)
    {
	data[i] = (double)(cbuf[i]);
    }
    delete cbuf;
    ifs.close();
}

void WriteBinaryData(const RVector &data, const char *fname)
{
    ofstream ofs (fname, ofstream::binary);
    /*int i;
    for (i=0; i<data.Dim(); ++i)
    {
	double val = data[i];
	ofs<<val;
    }*/
    ofs.write((const char *)data.data_buffer(), data.Dim()*sizeof(double));
    ofs.close();
}

void WriteData (const RVector &data, char *fname)
{
    ofstream ofs (fname);
    ofs << setprecision(14);
    ofs << data << endl;
}

void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname)
{
    int q, m, i;
    ofstream ofs (fname);
    for (m = i = 0; m < mesh.nM; m++) {
	for (q = 0; q < mesh.nQ; q++) {
	    if (mesh.Connected (q,m)) ofs << data[i++];
	    else                      ofs << '-';
	    ofs << (q == mesh.nQ-1 ? '\n' : '\t');
	}
    }   
}


