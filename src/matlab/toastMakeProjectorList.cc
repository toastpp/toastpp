// Interface:
// RH-1: mesh handle
// RH-2: Image dimensions (2x1) 
// RH-3: Axis (3x1)
// RH-4: Camera type = {'ORTHO', 'PINHOLE'}
// RH-5: Camera parameters:
//				'flen' = focal length of pinhole 
//				'pixelsize' = scale of ortho projection in mm

#include "mex.h"
#include "util.h"
#include "toastmex.h"
#include "fdotlib.h"
//#include "felib.h"
//#include "stoastlib.h"
//#include "camera.h"
#include "GLProjector.h"


using namespace std;

typedef enum { CAMTYPE_PINHOLE, CAMTYPE_ORTHO, CAMTYPE_NONE } CameraType;
typedef enum { PROJTYPE_MESAGL, PROJTYPE_NONE } ProjectorType;


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
      std::ofstream   fout("makeFMTFwd.log");
      std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
     */
    size_t nr, nc;
    double * dpr;

    // Get mesh pointer
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar(prhs[0]));
    mesh->MarkBoundary ();
    int nQ, nM;
    nQ = mesh->nQ;
    nM = mesh->nM;

    // Image dimensions
    nr = mxGetM(prhs[1]); nc = mxGetN(prhs[1]);
    if (nr*nc < 2)
    {
	cerr << "Warning: Less than two image dimensions specified";
    }
    dpr = (double*)mxGetPr (prhs[1]);
    int w = (int)dpr[0];
    int h = (int)dpr[1];

    // Gantry axis - orients camera up
    RVector axis(3, 0.0);
    nr = mxGetM(prhs[2]);
    nc = mxGetN(prhs[2]);
    if (nr*nc < 3)
    {
	cerr << "Warning: axis had less than 3 elements - missing elements will be = 0" << endl;
    }
    dpr = mxGetPr(prhs[2]);
    for (int i=0; i<(nr*nc); ++i) axis[i] = dpr[i];
    double gantryR = l2norm(axis);
    axis /= gantryR;

    // Get camera type
    int pidx = 3;
    CameraType ctp;
    char str[100];
    double pixSize, f;
    mxGetString(prhs[pidx++], str, 100);
    if (!strcasecmp (str, "PINHOLE")) {
	ctp = CAMTYPE_PINHOLE;
	// Get focal length
	mxGetString (prhs[pidx++],str,100);
	if (!strcmp(str, "flen"))
	    f = mxGetScalar (prhs[pidx++]);
	else {
	    mexErrMsgTxt
		("toastMakeProjectorList: missing required field: flen");
	}
    } else if (!strcasecmp (str, "ORTHO")) {
	ctp = CAMTYPE_ORTHO;
    } else {
	char cbuf[256];
	sprintf (cbuf, "Camera type %s not supported", str);
	mexErrMsgTxt (cbuf);
    }

    // pixel size
    mxGetString (prhs[pidx++],str,100);
    if (!strcmp(str, "pixelsize"))
	pixSize = mxGetScalar (prhs[pidx++]);
    else {
	mexWarnMsgTxt
	    ("Warning: toastMakeProjectorList: missing required field: pixelsize set to 1");
	pixSize = 1.0;
    }

    RVector shift(2, 0.0);
    if (nrhs > pidx-1)
    {
	mxGetString (prhs[pidx++],str,100);
	if (!strcasecmp(str, "shift"))
	{
	    nr = mxGetM(prhs[pidx]);
	    nc = mxGetN(prhs[pidx]);
	    if (nr*nc < 2)
	    {
		mexErrMsgTxt
		    ("Must specify both x and y camera position shift");
	    } else{
		dpr = mxGetPr(prhs[pidx]);
		shift[0] = dpr[0];
		shift[1] = dpr[1];
	    }
	}
    }

    // Construct cameras and projectors
    Projector ** projPList = new Projector*[nM];
    Camera ** camPList = new Camera*[nM];

    // Get projector type - assume MESAGL for now
    ProjectorType projtp = PROJTYPE_MESAGL;

    // Create cameras - assume:
    //	    1. Detector positions are set in mesh
    //	    3. camera pointer array is allocated
    int * bndellist, * bndsdlist;
    int nBndElems=0;
    nBndElems = mesh->BoundaryList (&bndellist, &bndsdlist);

    // Create vector of handles for returning pointer list
    plhs[0] = mxCreateDoubleMatrix (nM, 1, mxREAL);
    dpr = mxGetPr (plhs[0]);

    // Create projectors
    for (int m=0; m<mesh->nM; ++m)
    {
	// Setup camera
	Point cpos = mesh->M[m]; 
	RVector z = mesh->MN[m];
	if (l2norm(z)==0)
	{
	    z = dot(cpos, axis)*axis - cpos; // Look at axis
	}
	cpos = cpos + z * (1.0 - gantryR / l2norm(z)); // force camera onto gantry radius
	RVector x(3), y(3);
	x = axis;   // image x-axis is along rotation axis of gantry
	y = cross3(x, z);  // image y-axis is perp. to x and z
	y /= l2norm(y);
	z /= l2norm(z);
	cpos = cpos + x*shift[0] + y*shift[1]; // - y * 2.1428; // manual adjustment for optical axis offset in image
//	cout <<"Cam pos = "<<cpos<<endl;
//	cout << "x = "<<x<<" y= "<<y<<" z = "<<z<<endl;

	// Create camera
	if (ctp==CAMTYPE_PINHOLE)
	{
	    camPList[m] = new PinholeCamera(w, h, f, cpos, x, y, z, pixSize);

	}else{
	    if (ctp==CAMTYPE_ORTHO)
	    {
		camPList[m] = new OrthoCamera(w, h, pixSize, cpos, x, y, z);
	    }
	}

	if (projtp==PROJTYPE_MESAGL)
	{
	    projPList[m] = new GLProjector(camPList[m], mesh, nBndElems, bndellist, bndsdlist);
	}

	dpr[m] = Ptr2Handle (projPList[m]);

    } 

    delete []camPList;
    delete []projPList;
    // note that these lists can be deleted here, because they have been
    // copied into the solver instances.

    // ...
    //    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
}

// =========================================================================


