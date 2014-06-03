#ifndef FDOTFWD_H
#define FDOTFWD_H

#include "fdotlib.h"
//#include "stoastlib.h"
//#include "supermatrix.h"
//#include "fwdsolver.h"
//#include "qmmesh.h"
#include "projector.h"
#include <vector>

using namespace std;

class FDOTLIB FDOTFwd 
{
    public:
	FDOTFwd(   RFwdSolver * _FEMSolver, QMMesh & mesh, 
			    Raster * rast,
			    int numSources, const RCompRowMatrix & qVecs,
			    Projector ** projList);
	~FDOTFwd(); 
	void fwdOperator(RVector & x, bool ratio=false, double epsilon=0.0);
	void adjOperator(RVector & x, bool ratio=false, double epsilon=0.0);
	void fwdOperator(const RVector & x, RVector & result, bool ratio=false, double epsilon=0.0);
	void adjOperator(const RVector & x, RVector & result, bool ratio=false, double epsilon=0.0);
	void adjOperator(RVector &b, int q, bool ratio=false, double epsilon=0.0);	
	void adjFwdOperator(RVector & x, bool ratio=false, double epsilon=0.0);
	void adjOperator(const RVector &x, RVector &result, int q, bool ratio, double epsilon=0.0);
	void calcExcitationData();
	void getExcitationImages(RVector & img);
	RVector * getExcitationFields();
	RCompRowMatrix & getSysmat();

    protected:
	RFwdSolver * FEMSolver;
	QMMesh & FEMMesh;
	Raster * raster;
	double voxelSize;
	int nQ, nNodes, nBndNodes, * bndNodes, nImagePts;
	RVector * phi_e;
	RVector excitImg;
	Projector ** projectors;
	RCompRowMatrix qVecs;
	RCompRowMatrix meshToGridMap, meshToGridMapT;
};

// A helper function to pass into BiCGSTAB, since it cannot take member funcs
RVector fdot_adjFwdCaller(const RVector& x, void * context);

#endif // !FDOTFWD_H
