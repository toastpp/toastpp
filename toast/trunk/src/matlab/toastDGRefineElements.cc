// =========================================================================
// toastDGRefineElements.cc
// Refines the given elements of a TOAST mesh 
//
// RH parameters:
//     1: mesh handle (pointer) 
//     2: element list (integer, m, where m is number of
//        elements). The list is
//        1-based.
// LH parameters:
//     1: new mesh handle (pointer)
// =========================================================================

#include "mex.h"
#include "toastmex.h"
#include "nonconformingMesh.h"
#include "util.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"
void Assert (bool cond, const char *msg)
{
    if (!cond) {
	char cbuf[256] = "toastDGRefineElements: ";
	strcat (cbuf, msg);
	mexErrMsgTxt (cbuf);
    }
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, k;
    RVector elems;
    CopyVector(elems, prhs[1]);

    int nel = elems.Dim();
    NonconformingMesh *mesh  = (NonconformingMesh *)Handle2Ptr (mxGetScalar (prhs[0]));;
    mexPrintf("Number of elements to refine : %d\n", nel);
    for(int i=0; i<nel; i++){
	Assert(mesh->elist[i]->Type() == ELID_TET4, "Unsupported element: Only 4-noded tetrahedra currently supported");
        mesh->RefineTetElem((int)(elems[i]-0.5));
        if((i+1)%1000 == 0){ fflush(stdout); mexPrintf("Refined %d elements\n", i+1);}
               
    }
	

    plhs[0] = mxCreateScalarDouble (Ptr2Handle (mesh));

   }
