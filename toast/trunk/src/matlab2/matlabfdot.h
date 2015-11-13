// -*-C++-*-
// ========================================================================
// Declaration of class MatlabFDOT
// MEX interface for TOAST-FDOT functions
// Fluorescence optical tomography solver
// ========================================================================

#ifndef __MATLABFDOT_H
#define __MATLABFDOT_H

#include "FDOTFwd.h"
#include "matlabtoast.h"

class MatlabFDOT {
public:
    MatlabFDOT();
    ~MatlabFDOT();

    void MakeFwd (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void DelFwd (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void FwdOp (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void FwdOpRatio (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void AdjOp (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void AdjOpRatio (int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]);
    void Excit (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void Sysmat (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

    void MakeProjectorList (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]);
    void ProjectToImage (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]);
    void ProjectToField (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]);
    
private:
    FDOTFwd *GetFDOTFwd (const mxArray *idx);

    FDOTFwd **fdotfwdlist;     // list of FDOT forward solvers
    unsigned int nfdotfwd;     // number of FDOT forward solvers
};

#endif // !__MATLABFDOT_H
