#define FDOTLIB_IMPLEMENTATION
#include "fdotlib.h"

//#include "MLEMSolver.h"
//#include "matrix.h"
//#include "util.h"


MLEMSolver::MLEMSolver( RFwdSolver & _FEMSolver, QMMesh & mesh, 
				    Regularisation * reg, Raster * rast,
				    int numSources, const RCompRowMatrix & qVecs_, 
				    Projector ** projList)
    :   FEMSolver(_FEMSolver), FEMMesh(mesh), 
	nQ(numSources), nNodes(mesh.nlen()), 
	nBndNodes(mesh.nbnd()),
	projectors(projList), regul(reg), raster(rast),
	qVecs(qVecs_)
{
    nImagePts = projectors[0]->getImageSize(); 
    phi_e = new RVector[nQ];
    for (int i=0; i<nQ; ++i)
    {
	phi_e[i].New(nNodes);
    }
   	
    if (regul)
    {
	RVector x(rast->SLen());
	regul->SetHess1(regHess, x, 0);
    }	
}

MLEMSolver::~MLEMSolver() 
{
    delete[] phi_e; 
}

void MLEMSolver::getExcitationImages(RVector & img)
{
    FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);
    RVector tmpImg(nImagePts);
    img.New(0);
    for (int i=0; i<nQ; ++i)
    {
	projectors[i]->projectFieldToImage(phi_e[i], tmpImg);
	append(img, tmpImg);
    }
}

RVector * MLEMSolver::getExcitationFields()
{
    return phi_e;
}

RVector mlem_adjFwdCaller(const RVector& x, void * context)
{
    MLEMSolver * solver = (MLEMSolver*) context;
    RVector result = x;
    solver->adjFwdOperator(result);
    return result;
}


int MLEMSolver::Solve(const RVector & data, 
		    RVector & result, const RPreconditioner * precon, 
		    int maxIters, double tol)
{
    FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);
    RVector alpha(raster->SLen()), gAlpha(raster->GLen()), ones(nImagePts*nQ, 1.0);
    cout<<"Calculate adjoint field for data = 1"<<endl;
    adjOperator(ones, alpha);
    raster->Map_SolToGrid(alpha, gAlpha);
    WriteBinaryData(gAlpha, "alpha.dat");

    // Set initial guess to 1's
    result = 1.0;

    cout<<"Run Newton iterations"<<endl;
    RVector y(nImagePts*nQ), r(nImagePts*nQ), z(raster->SLen()), zz(raster->SLen());
    for (int i=0; i<1; ++i)
    {
	fwdOperator(result, y);
	char fname[100];
	sprintf(fname, "fwdOp_on_result%d.dat", i);
	WriteBinaryData(y, fname);
	for (int j=0; j<y.Dim(); ++j)
	{
	    if (y[j]==0.0)
	    {
		y[j] = 1.0;	// catch divide by zero
		if (data[i]!=0.0)
		{   
		    cout << "Setting y to 1.0, but data is non-zero!"<<endl;
		}
	    }
	}
	r = data/y;
	adjOperator( r, z );
	zz = z / alpha;
	cout << "Norm(zz) = "<<l2norm(zz)<<endl;
	result *= zz;
	cout << "Iteration "<<(i-1)<<", error = "<<l2norm(y-data)/l2norm(data)<<endl;
    }
    RVector ATr(raster->GLen()), ATr_over_alpha(raster->GLen());
    WriteBinaryData(r, "r.dat");
    raster->Map_SolToGrid(z, ATr);
    WriteBinaryData(ATr, "ATr.dat");
    raster->Map_SolToGrid(zz, ATr_over_alpha);
    WriteBinaryData(ATr_over_alpha, "ATr_over_alpha.dat");

	return 0; // MS100428: added
}

void MLEMSolver::adjFwdOperator(RVector & x)
{
    if (l2norm(x)==0.0)
    {
	cout<<"x = 0, returning zero vector"<<endl;
	return;
    }

    // Construct regularisation term 
    RVector reg(raster->SLen());
    if (regul)
    {
	//RCompRowMatrix hess;
	//regul->SetHess1(hess, x, 0);
	reg = regHess*x;
	cout << "Max reg = "<<vmax(reg)<<" Min reg = "<<vmin(reg)<<"  l2norm(reg) = "<<l2norm(reg)<<endl;
    }else{
	reg = x * 1e-8;
    }

    // Compute ATA(x) 
    cout<<"Max x = "<<vmax(x)<<", min x = "<<vmin(x)<<endl;    
    fwdOperator(x);
    cout<<"Max fwd = "<<vmax(x)<<", min fwd = "<<vmin(x)<<endl;
    adjOperator(x);
    cout<<"Max adj = "<<vmax(x)<<", min adj = "<<vmin(x)<<endl;

    // Add regularisation term
    x += reg;

}

void MLEMSolver::fwdOperator(RVector & x)
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
	cout<<"Max phi = "<<vmax(phi_e[i])<<", min phi = "<<vmin(phi_e[i])<<endl;
	RVector pf = fld * phi_e[i];
	if (l2norm(pf)==0.0)
	{
	    cout<<"pf = 0, returning zero vector"<<endl;
	    tmpImg.Clear();
	}
	else
	{
	    cout<<"Max pf = "<<vmax(pf)<<", min pf = "<<vmin(pf)<<endl;
	    FEMSolver.CalcField (/*FEMMesh,*/ pf, tmpFld);
	    cout<<"Max tmpFld = "<<vmax(tmpFld)<<", min tmpFld = "<<vmin(tmpFld)<<endl;	
	    projectors[i]->projectFieldToImage(tmpFld, tmpImg);
	}
	append(phi_f, tmpImg);
	cout<<"Done"<<endl;
    }
    
    // in-place operation on x
    x = phi_f;
}

void MLEMSolver::adjOperator(RVector & b)
{
    RVector tmpImg(nImagePts),
	    tmpFld(nNodes),
	    adjPhi_f(nNodes), 
	    result(nNodes),
	    gAdjPhi_f(raster->SLen()),
	    gPhi_e(raster->SLen());
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate adj field for source "<<i<<endl;
	tmpImg.Copy(b, 0, i*nImagePts, nImagePts);	
	projectors[i]->projectImageToField(tmpImg, tmpFld);
	char fname[100];
	sprintf(fname, "tmpFld%d.dat", i);
	WriteBinaryData(tmpFld, fname);
	FEMSolver.CalcField (/*FEMMesh,*/ tmpFld, adjPhi_f);
	raster->Map_MeshToGrid(adjPhi_f, gAdjPhi_f);
	raster->Map_MeshToGrid(phi_e[i], gPhi_e);
/*	if (i<3)
	{
	//    char fname[100];
	    sprintf(fname, "backproj%d.dat", i);
	    WriteBinaryData(gAdjPhi_f, fname);
	    sprintf(fname, "fwdfield%d.dat", i);
	    WriteBinaryData(gPhi_e, fname);
	}*/
	result += adjPhi_f*phi_e[i];
	cout<<"Done"<<endl;
    }

    // in-place operation on b, map into regular grid basis
    //b = *gridMap * result;
    raster->Map_MeshToSol(result, b);
}

void MLEMSolver::adjOperator(const RVector & x, RVector & result)
{
    result = x;
    adjOperator(result);
}

void MLEMSolver::fwdOperator(const RVector & x, RVector & result)
{
    result = x;
    fwdOperator(result);
}


