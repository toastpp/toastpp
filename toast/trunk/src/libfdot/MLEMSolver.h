#ifndef MLEM_SOLVER_H
#define MLEM_SOLVER_H

#include "stoastlib.h"
#include "supermatrix.h"
#include "fwdsolver.h"
#include "qmmesh.h"
#include "projector.h"
#include <vector>
#include "fluoSolver.h"

using namespace std;

class FDOTLIB MLEMSolver : public FluoSolver
{
    public:
	MLEMSolver(   RFwdSolver & _FEMSolver, QMMesh & mesh, 
			    Regularisation * reg, Raster * rast,
			    int numSources, const RCompRowMatrix & qVecs,
			    Projector ** projList);
	~MLEMSolver(); 
	int Solve(const RVector & data, 
		    RVector & result, const RPreconditioner * precon, 
		    int maxIters, double tol=1e-7);
	int SolveNonLin(const RVector & data, 
		    RVector & result, const RPreconditioner * precon, 
		    int maxIters, double tol)
	{
	    Solve(data, result, precon, maxIters, tol);
		return 0; // MS100428: added
	}
	void fwdOperator(RVector & x);
	void adjOperator(RVector & x);
	void fwdOperator(const RVector & x, RVector & result);
	void adjOperator(const RVector & x, RVector & result);
	void adjFwdOperator(RVector & x);
	//void projectFieldToImage(const RVector & phi, RVector & image);
	//void projectImageToField(const RVector & image, RVector & phi);
	void getExcitationImages(RVector & img);
	void testProj(RVector & field){}
	RVector * getExcitationFields();

    protected:
	RFwdSolver & FEMSolver;
	QMMesh & FEMMesh;
	Regularisation * regul;
	Raster * raster;
	int nQ, nNodes, nBndNodes, * bndNodes, nImagePts;
	RVector * phi_e;
	Projector ** projectors;
	RCompRowMatrix regHess, qVecs;
};

// A helper function to pass into BiCGSTAB, since it cannot take member funcs
RVector mf_adjFwdCaller(const RVector& x, void * context);

#endif
