// =========================================================================
// toastReadMesh
// Reads a TOAST mesh from file to dynamic heap, initialises it
// and returns a pointer to MATLAB.
//
// RH parameters:
//     1: mesh file name (string)
// LH parameters:
//     1: mesh handle (pointer)
// =========================================================================


#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"
#include "toastmex.h"

using namespace std;

bool fileExists(const std::string& fileName)
{
  std::fstream fin;
  fin.open(fileName.c_str(),std::ios::in);
  if( fin.is_open() )
  {
    fin.close();
    return true;
  }
  fin.close();
  return false;
}


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char meshname[256];
    mxGetString (prhs[0], meshname, 256);
    
    if (fileExists(meshname) == 1) {
	QMMesh *mesh = new QMMesh;
	ifstream ifs (meshname);
	ifs >> *mesh;
	if (!ifs.good())
	    mexErrMsgTxt ("Mesh file not found or invalid format.\n");

	mesh->Setup();
	plhs[0] = mxCreateDoubleScalar (Ptr2Handle (mesh));
    } else {
	mexErrMsgTxt ("Mesh file not found or invalid format.\n");
    }
}
