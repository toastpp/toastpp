#include "FDOTFwd.h"
#include "matrix.h"
#include "dgmatrix.h"
#include "util.h"

FDOTFwd::FDOTFwd( RFwdSolver * _FEMSolver, QMMesh & mesh, 
				    Raster * rast,
				    int numSources, const RCompRowMatrix & qVecs_, 
				    Projector ** projList)
    :   FEMSolver(_FEMSolver), FEMMesh(mesh), 
	nQ(numSources), nNodes(mesh.nlen()), 
	nBndNodes(mesh.nbnd()),
	projectors(projList), raster(rast),
	qVecs(qVecs_)
{
    nImagePts = projectors[0]->getImageSize(); 

    phi_e = new RVector[nQ];
    for (int i=0; i<nQ; ++i)
    {
	phi_e[i].New(nNodes);
    }
   	
    calcExcitationData();
}

FDOTFwd::~FDOTFwd() 
{
    delete[] phi_e; 
    delete FEMSolver;
    delete projectors;
}

RCompRowMatrix & FDOTFwd::getSysmat()
{
    // Return system matrix to MATLAB
    return *(FEMSolver->F);
}

void FDOTFwd::calcExcitationData()
{
    for (int i=0; i<nQ; ++i) phi_e[i] = 0.0; // reset for initial guess of iterative solver

    FEMSolver->CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);
    RVector tmpImg(nImagePts);
    excitImg.New(0);
    for (int i=0; i<nQ; ++i)
    {
	projectors[i]->projectFieldToImage(phi_e[i], tmpImg);
	append(excitImg, tmpImg);
    }
}

void FDOTFwd::getExcitationImages(RVector & img)
{
    img = excitImg;
}

RVector * FDOTFwd::getExcitationFields()
{
    return phi_e;
}

RVector fdot_adjFwdCaller(const RVector& x, void * context)
{
    FDOTFwd * solver = (FDOTFwd*) context;
    RVector result = x;
    solver->adjFwdOperator(result);
    return result;
}

void FDOTFwd::adjFwdOperator(RVector & x, bool ratio, double epsilon)
{
    RVector ox = x;
    double xnorm = l2norm(x);
    if (xnorm==0.0)
    {
	cout<<"x = 0, returning zero vector"<<endl;
	return;
    }

    // Compute ATA(x) 
    fwdOperator(x, ratio, epsilon);
    adjOperator(x, ratio, epsilon);
}

void FDOTFwd::fwdOperator(RVector & x, bool ratio, double epsilon)
{
    RVector phi_f;
    RVector fld(nNodes), tmpFld(nNodes), tmpImg(nImagePts);

    // Map field from grid to node space
    //fld = *nodeMap * x;
    raster->Map_SolToMesh(x, fld);
    
    // Run fwd solver for each source
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate fwd field for source "<<i<<endl;
	RVector pf = fld * phi_e[i];
	if (l2norm(pf)==0.0)
	{
	    cout<<"pf = 0, returning zero vector"<<endl;
	    tmpImg.Clear();
	}
	else
	{
	    FEMSolver->CalcField (/*FEMMesh,*/ pf, tmpFld);
	    projectors[i]->projectFieldToImage(tmpFld, tmpImg);
	}
	append(phi_f, tmpImg);
	cout<<"Done"<<endl;
    }
    
    // in-place operation on x
    if (!ratio)
    {
	x = phi_f;
    }else{
	x = phi_f / (excitImg + epsilon) ;
    }
    cout << "fwdOperator: Range of images = (" << vmin(x) << ", "<< vmax(x) << ")"<<endl;
}

void FDOTFwd::adjOperator(RVector & b, bool ratio, double epsilon)
{
    RVector tmpImg(nImagePts),
	    tmpFld(nNodes),
	    adjPhi_f(nNodes), 
	    result(nNodes);
    if (ratio) b /= (excitImg + epsilon);
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate adj field for source "<<i<<endl;
	tmpImg.Copy(b, 0, i*nImagePts, nImagePts);	
	projectors[i]->projectImageToField(tmpImg, tmpFld);
	FEMSolver->CalcField (/*FEMMesh,*/ tmpFld, adjPhi_f);
	result += adjPhi_f*phi_e[i];
	cout<<"Done"<<endl;
    }

    raster->Map_MeshToSol(result, b);
}

// single-source version
void FDOTFwd::adjOperator(RVector &b, int q, bool ratio, double epsilon)
{
    RVector tmpImg(nImagePts),
	    tmpFld(nNodes),
	    adjPhi_f(nNodes), 
	    result(nNodes, 0.0);
    if (ratio) b /= (excitImg + epsilon);


    cout<<"Calculate adj field for source "<<q<<endl;
    tmpImg.Copy(b, 0, q*nImagePts, nImagePts);	
    projectors[q]->projectImageToField(tmpImg, tmpFld);
    FEMSolver->CalcField (/*FEMMesh,*/ tmpFld, adjPhi_f);
    result += adjPhi_f*phi_e[q];
    cout<<"Done"<<endl;
    

    raster->Map_MeshToSol(result, b);
}

void FDOTFwd::adjOperator(const RVector &x, RVector &result, int q, bool ratio, double epsilon)
{
    result = x;
    adjOperator (result, q, ratio, epsilon);
}

void FDOTFwd::adjOperator(const RVector & x, RVector & result, bool ratio, double epsilon)
{
    result = x;
    adjOperator(result, ratio, epsilon);
}

void FDOTFwd::fwdOperator(const RVector & x, RVector & result, bool ratio, double epsilon)
{
    result = x;
    fwdOperator(result, ratio, epsilon);
}


