// Interface:
// RH-1: FDOTFwd handle
// LH-1: FEM system matrix

#include "mex.h"
#include "felib.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "util.h"
#include "FDOTFwd.h"

// =========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("makeFDOTFwd.log");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    */
    // get FDOTFwd pointer from handle
    FDOTFwd *fsolver = (FDOTFwd*)Handle2Ptr (mxGetScalar(prhs[0]));
    CopyMatrix( &plhs[0], fsolver->getSysmat() );    

    /*
    // ...
    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer 
    */
}
