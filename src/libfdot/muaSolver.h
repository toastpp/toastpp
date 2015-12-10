#ifndef MUA_FREE_SOLVER_H
#define MUA_FREE_SOLVER_H

#include "stoastlib.h"
#include "supermatrix.h"
#include "fwdsolver.h"
#include "qmmesh.h"
#include "projector.h"
#include <vector>

using namespace std;

class FDOTLIB MuaSolver
{
    public:
	MuaSolver(   RFwdSolver & _FEMSolver, QMMesh & mesh, 
			    Regularisation * reg, Raster * rast,
			    int numSources, const RCompRowMatrix & qVecs_,
			    Projector ** projList, Solution & initSol);
	~MuaSolver(); 
	int Solve(const RVector & data, 
		    RVector & result, RPreconditioner * precon, 
		    int maxLinIters, int nonLinIters, double tol=1e-7, bool logSolve=false);

	void fwdOperator(RVector & x);
	void adjOperator(RVector & x);
	void fwdOperator(const RVector & x, RVector & result);
	void adjOperator(const RVector & x, RVector & result);
	void adjFwdOperator(RVector & x);

	void logFwdOperator(RVector & x, const RVector & linx);
	void logAdjOperator(RVector & x, const RVector & linx);
	void logFwdOperator(const RVector & x, const RVector & mLinx, RVector & result);
	void logAdjOperator(const RVector & x, const RVector & mLinx, RVector & result);
	void logAdjFwdOperator(RVector & x, const RVector & mLinx);

	const RVector & getLinearisationPoint() { return linx; }

	void projectFieldToImage(const RVector & phi, RVector & image);
	void projectImageToField(const RVector & image, RVector & phi);
	void getExcitationImages(RVector & img);

    protected:
	RFwdSolver & FEMSolver;
	QMMesh & FEMMesh;
	Regularisation * regul;
	Raster * raster;
	int nQ, nNodes, nBndNodes, * bndNodes, nImagePts;
	RVector * phi_e, linx;
	Projector ** projectors;
	RCompRowMatrix regHess, qVecs;
	Solution sol;
	IVector gDim;
};

// A helper function to pass into BiCGSTAB, since it cannot take member funcs
RVector mua_adjFwdCaller(const RVector& x, void * context);
TVector<double> mua_logAdjFwdCaller(const TVector<double>& logx, void * context);

#endif
