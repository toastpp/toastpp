#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"

using namespace std;

bool ReadNim (char *name, int idx, RVector &img)
{
    char cbuf[256];
    int i, j = 0, imgsize = 0;

    ifstream ifs (name);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM") && strcmp (cbuf, "RIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    for (;;) {
	do {
	    ifs.getline (cbuf, 256);
	} while (ifs.good() && strncasecmp (cbuf, "Image", 5));
	if (!ifs.good()) break;
	for (i = 0; i < imgsize; i++)
	    ifs >> img[i];
	if (++j == idx) break;
    }
    return true;
}

bool ReadNimAll (char *name, RDenseMatrix &img)
{
    char cbuf[256];
    int i, j = 0, imgsize = 0;

    ifstream ifs (name);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM") && strcmp (cbuf, "RIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;

    img.New(imgsize,0);
    RDenseMatrix img_single(imgsize,1);
    for (;;) {
	do {
	    ifs.getline (cbuf, 256);
	} while (ifs.good() && strncasecmp (cbuf, "Image", 5));
	if (!ifs.good()) break;
	for (i = 0; i < imgsize; i++)
	    ifs >> img_single(i,0);
	img = cath(img,img_single);
    }
    return true;
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Reads a NIM (nodal image) from a file.
    //
    // RH parameters:
    //     1: NIM file name (string)
    //     2: image index (integer >= 1; 0 = last image)
    // LH parameters:
    //     1: image data (double array)

    char nimname[256];
    mxGetString (prhs[0], nimname, 256);
    bool read_all = false;

    int idx;
    if (nrhs < 2) {
	idx = 1;
    } else if (mxIsChar(prhs[1])) {
	char cbuf[256];
	mxGetString (prhs[1], cbuf, 256);
	if (!strcasecmp (cbuf,"all"))
	    read_all = true;
	else
	    mexErrMsgTxt ("Argument 2 must be index or 'all'");
    } else {
	idx = (int)mxGetScalar (prhs[1]);
    }

    if (read_all) {
	RDenseMatrix img;
	ReadNimAll (nimname, img);
	CopyMatrix (&plhs[0], img);
	mexPrintf ("NIM: size = %d x %d\n", img.nRows(), img.nCols());
    } else {
	RVector img;
	ReadNim (nimname, idx, img);

	plhs[0] = mxCreateDoubleMatrix (img.Dim(), 1, mxREAL);
	memcpy (mxGetPr (plhs[0]), img.data_buffer(), img.Dim()*sizeof(double));
	
	mexPrintf ("NIM: size = %d\n", img.Dim());
    }
}
