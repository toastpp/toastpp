// =========================================================================
// toastReadQM
// Reads a QM (source-detector) specification for a TOAST mesh
// from a file
//
// RH parameters:
//     1: mesh handle (double)
//     2: QM file name (string)
// LH parameters:
//     <none>
// =========================================================================


#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "util.h"

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
    char qmname[256];
    mxGetString (prhs[1], qmname, 256);

    if (fileExists(qmname) == 1) {
        QMMesh *mesh = (QMMesh*)Handle2Ptr(mxGetScalar(prhs[0]));
	ifstream ifs(qmname);
	mesh->LoadQM (ifs);
    } else {
        mexErrMsgTxt ("QM file not found or invalid format.\n");
    }
}
