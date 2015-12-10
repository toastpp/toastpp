#define FDOTLIB_IMPLEMENTATION
#include "fdotlib.h"

//#include "muaSolver.h"
//#include "matrix.h"
//#include "util.h"


MuaSolver::MuaSolver( RFwdSolver & _FEMSolver, QMMesh & mesh, 
				    Regularisation * reg, Raster * rast,
				    int numSources, const RCompRowMatrix & qVecs_, 
				    Projector ** projList, Solution & initSol)
    :   FEMSolver(_FEMSolver), FEMMesh(mesh), 
	nQ(numSources), nNodes(mesh.nlen()), 
	nBndNodes(mesh.nbnd()),
	projectors(projList), regul(reg), raster(rast), sol(initSol),
	qVecs(qVecs_), gDim(rast->GDim())
{
    nImagePts = projectors[0]->getImageSize(); 
    phi_e = new RVector[nQ];
    for (int i=0; i<nQ; ++i)
    {
	phi_e[i].New(nNodes);
    }
    FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);
   	
    if (regul)
    {
	RVector x(rast->SLen());
	regul->SetHess1(regHess, x, 0);
    }	
}

MuaSolver::~MuaSolver() 
{
    delete[] phi_e; 
}

void MuaSolver::getExcitationImages(RVector & img)
{
    RVector tmpImg(nImagePts);
    img.New(0);
    for (int i=0; i<nQ; ++i)
    {
	projectors[i]->projectFieldToImage(phi_e[i], tmpImg);
	append(img, tmpImg);
    }
}

RVector mua_adjFwdCaller(const RVector& x, void * context)
{
    MuaSolver * solver = (MuaSolver*) context;
    RVector result = x;
    solver->adjFwdOperator(result);
    return result;
}

TVector<double> mua_logAdjFwdCaller(const TVector<double>& logx, void * context)
{
    MuaSolver * solver = (MuaSolver*) context;
    RVector result = logx;
    const RVector & linx = solver->getLinearisationPoint();
    solver->logAdjFwdOperator(result, linx);
    return result;
}

int MuaSolver::Solve(const RVector & data, 
		    RVector & result, RPreconditioner * precon, 
		    int maxLinIters, int nonLinIters, double tol, bool logSolve)
{
	TVector<double> y(nNodes), Fx(data.Dim()), delta(result.Dim());

    // initial guess is supplied in result vector 
    sol.SetParam(OT_CMUA, result);
    FEMSolver.Reset(sol, 0.0);

    // Compute excitation field for estimated mua
    FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);    
    getExcitationImages(Fx);

    cout<<"Initial error = " << l2norm(data-Fx) << endl;

    for (int i=0; i<nonLinIters; ++i)
    {
	linx = result; // set linearisation point

	// Get delta_mua
	if (logSolve)
	{
	    cout<<"Calculate adjoint field"<<endl;
	    logAdjOperator(data-Fx, linx, y); //difference data
	    cout<<"Run BiCGSTAB"<<endl;
	    BiCGSTAB<double> (mua_logAdjFwdCaller, (void*)this, y, delta, tol, precon, maxLinIters);
	    result *= exp(delta);    // log(result)+=delta ==> result*=exp(delta)
	}
	else
	{
	    cout<<"Calculate adjoint field"<<endl;
	    adjOperator(data-Fx, y);
	    cout<<"Run BiCGSTAB"<<endl;
	    BiCGSTAB (mua_adjFwdCaller, this, y, delta, tol, precon, maxLinIters);
	    result += delta;
	}


	RVector mResult;
	raster->Map_SolToMesh(result, mResult);
	sol.SetParam(OT_CMUA, mResult);
	FEMSolver.Reset(sol, 0.0);

	// Recompute excitation field for estimated mua
	FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);    
	getExcitationImages(Fx);
	cout<<"Iteration "<< i+1 << " error = " << l2norm(data-Fx) << endl;

	// Write data out
	// Excitation data
	char fname[100];
	RVector im(nImagePts);
	IVector imgDim = projectors[0]->getImageDim();
	for (int j=0; j<nQ; ++j)
	{
	    sprintf(fname, "solnexcit_data_iter%d_img%d.pgm", i, j);
	    im.Copy(Fx, 0, j*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);

	    sprintf(fname, "excit_err_data_iter%d_img%d.pgm", i, j);
	    im.Copy(data-Fx, 0, j*nImagePts, nImagePts);
	    WritePGM(im, imgDim, fname);
	}

	// Mua sample slices
	RVector soln, gResult;
	raster->Map_SolToGrid(result, gResult);
	im.New(gDim[0]*gDim[1]);
	for (int j=0; j<4; ++j)
	{
	    int slice = j*(gDim[2]-1)/3;
	    im.Copy(gResult, 0, slice*gDim[0]*gDim[1], gDim[0]*gDim[1]);
	    append(soln, im);
	}
	sprintf(fname, "mua_soln_iter%d_slices.dat", i);
	WriteBinaryData(soln, fname);
    }

	return 0; // MS100428: added
}


void MuaSolver::logAdjFwdOperator(RVector & logx, const RVector & linx)
{
    // Solve in log(x), linx is the linearisation point

    if (l2norm(logx)==0.0)
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
	reg = linx * (regHess * exp(logx));
	cout << "Max reg = "<<vmax(reg)<<" Min reg = "<<vmin(reg)<<"  l2norm(reg) = "<<l2norm(reg)<<endl;
    }else{
	reg = logx * 1e-8;
    }

    // Map linearisation point into mesh
    RVector mLinx(nNodes);
    raster->Map_SolToMesh(linx, mLinx);

    // Compute ATA(x) 
    cout<<"Max logx = "<<vmax(logx)<<", min logx = "<<vmin(logx)<<endl;    
    logFwdOperator(logx, mLinx);
    cout<<"Max fwd = "<<vmax(logx)<<", min fwd = "<<vmin(logx)<<endl;
    logAdjOperator(logx, mLinx);
    cout<<"Max adj = "<<vmax(logx)<<", min adj = "<<vmin(logx)<<endl;

    // Add regularisation term
    logx += reg;

}

void MuaSolver::adjFwdOperator(RVector & x)
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
	reg = x * 1e-3;
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

void MuaSolver::logFwdOperator(RVector & x, const RVector & mLinx)
{
    RVector phi_f;
    RVector fld(nNodes), tmpFld(nNodes), tmpImg(nImagePts);

    // Calculate excitation fields
    FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);

    // Map field from grid to node space
    //fld = *nodeMap * x;
    raster->Map_SolToMesh(x, fld);
    
    // Run fwd solver for each source
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate fwd field for source "<<i<<endl;
	cout<<"Max phi = "<<vmax(phi_e[i])<<", min phi = "<<vmin(phi_e[i])<<endl;
	RVector pf = -fld * phi_e[i];
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
	    tmpFld *= mLinx; // gradient of log
	    projectors[i]->projectFieldToImage(tmpFld, tmpImg);
	}
	append(phi_f, tmpImg);
	cout<<"Done"<<endl;
    }
    
    // in-place operation on x
    x = phi_f;
}

void MuaSolver::fwdOperator(RVector & x)
{
    RVector phi_f;
    RVector fld(nNodes), tmpFld(nNodes), tmpImg(nImagePts);

    // Calculate excitation fields
    FEMSolver.CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);

    // Map field from grid to node space
    //fld = *nodeMap * x;
    raster->Map_SolToMesh(x, fld);
    
    // Run fwd solver for each source
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate fwd field for source "<<i<<endl;
	cout<<"Max phi = "<<vmax(phi_e[i])<<", min phi = "<<vmin(phi_e[i])<<endl;
	RVector pf = -fld * phi_e[i];
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

void MuaSolver::logAdjOperator(RVector & b, const RVector & mLinx)
{
    RVector tmpImg(nImagePts),
	    tmpFld(nNodes),
	    adjPhi_f(nNodes), 
	    result(nNodes);
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate adj field for source "<<i<<endl;
	tmpImg.Copy(b, 0, i*nImagePts, nImagePts);	
	projectors[i]->projectImageToField(tmpImg, tmpFld);
	FEMSolver.CalcField (/*FEMMesh,*/ tmpFld, adjPhi_f);
	result -= adjPhi_f*phi_e[i];
	cout<<"Done"<<endl;
    }

    // in-place operation on b, map into regular grid basis
    result *= mLinx;	// grad of log
    raster->Map_MeshToSol(result, b);
}

void MuaSolver::adjOperator(RVector & b)
{
    RVector tmpImg(nImagePts),
	    tmpFld(nNodes),
	    adjPhi_f(nNodes), 
	    result(nNodes);
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate adj field for source "<<i<<endl;
	tmpImg.Copy(b, 0, i*nImagePts, nImagePts);	
	projectors[i]->projectImageToField(tmpImg, tmpFld);
	FEMSolver.CalcField (/*FEMMesh,*/ tmpFld, adjPhi_f);
	result -= adjPhi_f*phi_e[i];
	cout<<"Done"<<endl;
    }

    // in-place operation on b, map into regular grid basis
    raster->Map_MeshToSol(result, b);
}

void MuaSolver::adjOperator(const RVector & x, RVector & result)
{
    result = x;
    adjOperator(result);
}

void MuaSolver::fwdOperator(const RVector & x, RVector & result)
{
    result = x;
    fwdOperator(result);
}

void MuaSolver::logAdjOperator(const RVector & x, const RVector & mLinx, RVector & result)
{
    result = x;
    logAdjOperator(result, mLinx);
}

void MuaSolver::logFwdOperator(const RVector & x, const RVector & mLinx, RVector & result)
{
    result = x;
    logFwdOperator(result, mLinx);
}


