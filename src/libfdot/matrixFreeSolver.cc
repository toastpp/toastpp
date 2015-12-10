#define FDOTLIB_IMPLEMENTATION
#include "fdotlib.h"

//#include "matrixFreeSolver.h"
//#include "matrix.h"
//#include "dgmatrix.h"
//#include "util.h"

MatrixFreeSolver::MatrixFreeSolver( RFwdSolver * _FEMSolver, QMMesh & mesh, 
				    Regularisation * reg, double & tau_, Raster * rast,
				    int numSources, const RCompRowMatrix & qVecs_, 
				    Projector ** projList, IVector & dataWin)
    :   FEMSolver(_FEMSolver), FEMMesh(mesh), 
	nQ(numSources), nNodes(mesh.nlen()), 
	nBndNodes(mesh.nbnd()),
	projectors(projList), regul(reg), raster(rast),
	qVecs(qVecs_), tau(tau_)
{
    nImagePts = projectors[0]->getImageSize(); 

    /* Construct mask matrix */
    int dataWinW = dataWin[1] - dataWin[0];
    int dataWinH = dataWin[3] - dataWin[2];
    dataWinSize = dataWinW * dataWinH;
/*    int w = projectors[0]->getImageWidth();
    int h = projectors[0]->getImageHeight();
    mask.New(nImagePts, nImagePts, 1.0);
    
    RVector r(nImagePts,0.0);
    for (int i=dataWin[0], row=0; i<dataWin[1]; ++i)
    {
	for (int j=dataWin[2]; j<dataWin[3]; ++j)
	{
	    int k = i + j*w;
	    r[k] = 1.0;
	    mask.SetRow(row, r);
	    r[k] = 0.0;
	    ++row;
	}
    }
    
    maskT = transp(mask);*/

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
    calcExcitationData();
}

MatrixFreeSolver::~MatrixFreeSolver() 
{
    delete[] phi_e; 
    delete FEMSolver;
}

void MatrixFreeSolver::calcExcitationData()
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

void MatrixFreeSolver::testProj(RVector & field)
{
    RVector tmpImg(nImagePts); //, tmpFld;
    tmpImg.Copy(excitImg, 0, 0, nImagePts);	
    projectors[0]->projectImageToField(tmpImg, field);
    //raster->Map_MeshToGrid(tmpFld, field);
}

void MatrixFreeSolver::getExcitationImages(RVector & img)
{
    img = excitImg;
}

RVector * MatrixFreeSolver::getExcitationFields()
{
    return phi_e;
}

RVector mf_adjFwdCaller(const RVector& x, void * context)
{
    MatrixFreeSolver * solver = (MatrixFreeSolver*) context;
    RVector result = x;
    solver->adjFwdOperator(result);
    return result;
}


int MatrixFreeSolver::Solve(const RVector & data, 
		    RVector & result, RPreconditioner * precon, 
		    int maxIters, double tol)
{
    calcExcitationData();

    RVector y(nNodes);
    cout<<"Calculate adjoint field"<<endl;
    adjOperator(data, y);

    double ATynorm = l2norm(y);
    cout << "MatrixFreeSolver::Solve: ATy norm = "<<ATynorm<<endl;
    RVector g;
    fwdOperator(result, g);
    double dataErr = l2norm(g-data) / l2norm(data);
    double sysErr = l2norm(y - mf_adjFwdCaller(result, this)) / l2norm(y);
    RVector ATAx;
    adjOperator(g, ATAx);
    double sysErr2 = l2norm(y - ATAx) / l2norm(y);
    cout	<< "MatrixFreeSolver::Solve: initial" 
	<< " data error = " << dataErr 
	<< ", system error = " << sysErr 
	<< ", system error2 = " << sysErr2
	<< endl;

    linReg = true;


    cout<<"Run PCG"<<endl;
    int iters;
    iters = BiCGSTAB<double> (&mf_adjFwdCaller, this, y, result, tol, precon, maxIters);
    //GMRES (&mf_adjFwdCaller, this, y, result, tol, precon, 10, &iters, &tol);
    cout << "MatrixFreeSolver::Solve: PCG error = " << tol << " in " << iters << " iterations" << endl;

    fwdOperator(result, g);
    dataErr = l2norm(g-data) / l2norm(data);
    sysErr = l2norm(y - mf_adjFwdCaller(result, this)) / l2norm(y);
    double dataTerm = l2normsq(g-data), regTerm = tau * l2normsq(result);
    cout	<< "MatrixFreeSolver::Solve: final " 
	<< " data error = " << dataErr 
	<< " data term = " << dataTerm
	<< ", prior term = " << regTerm
	<< ", system error = " << sysErr 
	<< endl;

    cout << "GREP THIS: " << dataTerm << ", " << regTerm << endl;

    /*GMRES (TVector<MT> (*Av_clbk)(const TVector<MT> &v, void *context),
    void *context, const TVector<MT> &b, TVector<MT> &x, double tol,
    const TPreconditioner<MT> *precon, int restart, int *iter, double *res)*/
    //    return GMRES (&mf_adjFwdCaller, this, y, result, tol, precon, 10, &iters, &tol);

	return 0; // MS100428: added
}

int MatrixFreeSolver::SolveNonLin(const RVector & data, 
		    RVector & result, RPreconditioner * precon, 
		    int maxIters, double tol)
{
    //FEMSolver->CalcFields(/*FEMMesh, nQ,*/ qVecs, phi_e);
    calcExcitationData();
    RVector y(nNodes, 0.0), oy(nNodes);
    cout<<"Calculate adjoint field"<<endl;
    adjOperator(data, oy);

    linReg = false;
    RVector delta(result.Dim());
    int iters=0; 
    RVector g;
    double dataErr, sysErr;
    for (int i=0; i<3; ++i)
    {
	if (regul)
	{
	    regul->SetHess1(regHess, result, 0);
	    RVector regGrad = regul->GetGradient(result);
	    cout << "MatrixFreeSolver::Solve: norm(reg->Grad) = " << l2norm(regGrad) << endl;
	    y = oy - mf_adjFwdCaller(result, this) - regGrad;
	}

	double ATynorm = l2norm(y);
	cout << "MatrixFreeSolver::Solve: ATy norm = "<<ATynorm<<endl;
	fwdOperator(result, g);
	dataErr = l2norm(g-data) / l2norm(data);
	sysErr = ATynorm/l2norm(oy);
	cout	<< "MatrixFreeSolver::Solve: iteration " << i 
		<< " data error = " << dataErr 
		<< ", system error = " << sysErr 
		<< endl;
	if (sysErr < tol)
	{
	    cout << "MatrixFreeSolver::Solve: Stopping iterations early, system error < " << tol << endl;
	    break;
	}

	cout<<"MatrixFreeSolver::Solve: Run BiCGSTAB"<<endl;
	double ttol=1e-2;
	int iiters;
	delta = 0.0;
	iiters = BiCGSTAB<double> (&mf_adjFwdCaller, this, y, delta, ttol, precon, maxIters);
	//GMRES (&mf_adjFwdCaller, this, y, delta, ttol, precon, 2, &iiters, &ttol);
/*MATHLIB int GMRES (TVector<MT> (*Av_clbk)(const TVector<MT> &v, void *context),
    void *context, const TVector<MT> &b, TVector<MT> &x, double tol,
    const TPreconditioner<MT> *precon, int restart, int *iter, double *res)*/
	cout << "MatrixFreeSolver::Solve: BiCGSTAB delta error = "<<ttol<<", iterations = "<<iiters<<endl;
	iters += iiters;
	result += delta;
    }
    fwdOperator(result, g);
    dataErr = l2norm(g-data) / l2norm(data);
    sysErr = l2norm(y - mf_adjFwdCaller(result, this)) / l2norm(y);
    double dataTerm = l2normsq(g-data), regTerm = tau * l2normsq(result);
    cout	<< "MatrixFreeSolver::Solve: final " 
	<< " data error = " << dataErr 
	<< " data term = " << dataTerm
	<< ", prior term = " << regTerm
	<< ", system error = " << sysErr 
	<< endl;
    cout << "GREP THIS: " << dataTerm << ", " << regTerm << endl;
    return iters;
}

void MatrixFreeSolver::adjFwdOperator(RVector & x)
{
    RVector ox = x;
    double xnorm = l2norm(x);
    if (xnorm==0.0)
    {
	cout<<"x = 0, returning zero vector"<<endl;
	return;
    }

    // Construct regularisation term 
    RVector reg(raster->SLen());
/*    if (linReg)
    {*/
	if (regul)
	{
	    //RCompRowMatrix hess;
	    //regul->SetHess1(hess, x, 0);
	    reg = regHess * x /*1e-12*/;
	    cout << "MatrixFreeSolver::adjFwdOperator: Norm(x) = " << l2norm(x) << endl;
	    cout << "MatrixFreeSolver::adjFwdOperator: Max reg = "<<vmax(reg)<<" Min reg = "<<vmin(reg)<<"  l2norm(reg) = "<<l2norm(reg)<<endl;
	    /*cout << "Adding 1*x to reg term --> + lambda*I to system"<<endl;
	    reg = reg + x*1e-12;*/
	}else{
	    reg = x * tau;
	}
  /*  }*/

    // Compute ATA(x) 
    cout<<"Max x = "<<vmax(x)<<", min x = "<<vmin(x)<<endl;    
    fwdOperator(x);
    cout<<"Max fwd = "<<vmax(x)<<", min fwd = "<<vmin(x)<<endl;
    adjOperator(x);
    cout<<"Max adj = "<<vmax(x)<<", min adj = "<<vmin(x)<<endl;

    xnorm = l2norm(x);
    cout << "MatrixFreeSolver::adjFwdOperator: Norm ATA(x) without regularisation = " << xnorm << endl;

    // Add regularisation term
    x += reg;
    xnorm = l2norm(x);
    cout << "MatrixFreeSolver::adjFwdOperator: Norm ATA(x) + regularisation = " << xnorm << endl;
    double xTATAx = dot(ox, x);
    cout << "MatrixFreeSolver::adjFwdOperator: xTATAx = " << xTATAx << endl;
    if (xTATAx < 0.0) cout << "MatrixFreeSolver::adjFwdOperator: Negative xTATAx"<< endl;
}

void MatrixFreeSolver::fwdOperator(RVector & x)
{
    RVector phi_f;
    RVector fld(nNodes), tmpFld(nNodes), tmpImg(nImagePts);

    // Map field from grid to node space
    //fld = *nodeMap * x;
    raster->Map_SolToMesh(x, fld);
    cout << "fwdOperator: Range of field = (" << vmin(fld) << ", "<< vmax(fld) << ")"<<endl;
    
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
	    FEMSolver->CalcField (/*FEMMesh,*/ pf, tmpFld);
	    cout<<"Max tmpFld = "<<vmax(tmpFld)<<", min tmpFld = "<<vmin(tmpFld)<<endl;	
	    projectors[i]->projectFieldToImage(tmpFld, tmpImg);
	}
	append(phi_f, tmpImg);
	cout<<"Done"<<endl;
    }
    
    // in-place operation on x
    cout << "Min excitImg = "<<vmin(excitImg)<<endl;
    x = phi_f / (excitImg + 1e-9);
    cout << "fwdOperator: Range of images = (" << vmin(x) << ", "<< vmax(x) << ")"<<endl;
}

void MatrixFreeSolver::adjOperator(RVector & b)
{
    RVector tmpImg(nImagePts),
	    tmpFld(nNodes),
	    adjPhi_f(nNodes), 
	    result(nNodes);
    cout << "Min excitImg = "<<vmin(excitImg)<<endl;
    b /= (excitImg + 1e-9);
    for (int i=0; i<nQ; ++i)
    {
	cout<<"Calculate adj field for source "<<i<<endl;
	tmpImg.Copy(b, 0, i*nImagePts, nImagePts);	
	cout << "adjOperator: Range of image = (" << vmin(tmpImg) << ", "<< vmax(tmpImg) << ")"<<endl;
	projectors[i]->projectImageToField(tmpImg, tmpFld);
	cout << "adjOperator: Range of projected data = (" << vmin(tmpFld) << ", "<< vmax(tmpFld) << ")"<<endl;
	FEMSolver->CalcField (/*FEMMesh,*/ tmpFld, adjPhi_f);
	cout << "adjOperator: Range of resulting field = (" << vmin(adjPhi_f) << ", "<< vmax(adjPhi_f) << ")"<<endl;
	result += adjPhi_f*phi_e[i];
	cout<<"Done"<<endl;
    }

    // in-place operation on b, map into regular grid basis
    //b = *gridMap * result;
    cout << "adjOperator: Range of mesh basis vector = (" << vmin(result) << ", " << vmax(result) << ")"<<endl;
    raster->Map_MeshToSol(result, b);
    cout << "adjOperator: Range of solution basis vector = (" << vmin(b) << ", " << vmax(b) << ")"<<endl;
    //b = raster->SmoothImage(b, 10);
//    WriteBinaryData(b, "AT.dat");
}

void MatrixFreeSolver::adjOperator(const RVector & x, RVector & result)
{
    result = x;
    adjOperator(result);
}

void MatrixFreeSolver::fwdOperator(const RVector & x, RVector & result)
{
    result = x;
    fwdOperator(result);
}


