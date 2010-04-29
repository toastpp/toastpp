#ifndef FLUO_SOLVER_H
#define FLUO_SOLVER_H

#include "stoastlib.h"
#include "supermatrix.h"
#include "fwdsolver.h"
#include "qmmesh.h"
#include "projector.h"
#include <vector>

using namespace std;

class FDOTLIB FluoSolver
{
    public:
	FluoSolver()
	{}
	~FluoSolver()
	{}
	virtual int Solve(const RVector & data, 
		    RVector & result, const RPreconditioner * precon, 
		    int maxIters, double tol=1e-7) = 0;
	virtual int SolveNonLin(const RVector & data, 
		    RVector & result, const RPreconditioner * precon, 
		    int maxIters, double tol) = 0;

	virtual void fwdOperator(RVector & x) = 0;
	virtual void adjOperator(RVector & x) = 0;
	virtual void fwdOperator(const RVector & x, RVector & result) = 0;
	virtual void adjOperator(const RVector & x, RVector & result) = 0;
	virtual void adjFwdOperator(RVector & x) = 0;
	//virtual void projectFieldToImage(const RVector & phi, RVector & image) = 0;
	//virtual void projectImageToField(const RVector & image, RVector & phi) = 0;
	virtual void getExcitationImages(RVector & img){};
	virtual RVector * getExcitationFields() = 0;
	virtual void testProj(RVector & field)=0;
	RCompRowMatrix & getMask() {return mask;}
	RCompRowMatrix & getTranspMask() {return maskT;}
    protected:
	RCompRowMatrix mask, maskT;
};

#endif
