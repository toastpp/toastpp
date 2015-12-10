#ifndef MATRIX_FREE_SOLVER_H
#define MATRIX_FREE_SOLVER_H

#include "stoastlib.h"
#include "supermatrix.h"
#include "fwdsolver.h"
#include "qmmesh.h"
#include "projector.h"
#include <vector>
#include "fluoSolver.h"

using namespace std;

class FDOTLIB MatrixFreeSolver : public FluoSolver
{
    public:
	MatrixFreeSolver(   RFwdSolver * _FEMSolver, QMMesh & mesh, 
			    Regularisation * reg, double & tau_, Raster * rast,
			    int numSources, const RCompRowMatrix & qVecs,
			    Projector ** projList, IVector & dataWin);
	~MatrixFreeSolver(); 
	int Solve(const RVector & data, 
		    RVector & result, RPreconditioner * precon, 
		    int maxIters, double tol=1e-7);
	int SolveNonLin(const RVector & data, 
		    RVector & result, RPreconditioner * precon, 
		    int maxIters, double tol);
	void fwdOperator(RVector & x);
	void adjOperator(RVector & x);
	void fwdOperator(const RVector & x, RVector & result);
	void adjOperator(const RVector & x, RVector & result);
	void adjFwdOperator(RVector & x);
	//void projectFieldToImage(const RVector & phi, RVector & image);
	//void projectImageToField(const RVector & image, RVector & phi);
	void calcExcitationData();
	void getExcitationImages(RVector & img);
	RVector * getExcitationFields();
	void testProj(RVector & field);

    protected:
	RFwdSolver * FEMSolver;
	QMMesh & FEMMesh;
	Regularisation * regul;
	Raster * raster;
	int nQ, nNodes, nBndNodes, * bndNodes, nImagePts, dataWinSize;
	RVector * phi_e;
	RVector excitImg;
	Projector ** projectors;
	RCompRowMatrix regHess, qVecs;
	bool linReg;	
	double tau; // copy of regularisation parameter
};

// A helper function to pass into BiCGSTAB, since it cannot take member funcs
RVector mf_adjFwdCaller(const RVector& x, void * context);

#endif
