#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "regul.h"

using namespace std;

// ==========================================================================

Regularisation::Regularisation (const Raster *_raster, double _tau,
    const RVector *_x0)
    : raster(_raster), tau(_tau)
{
    //LOGOUT("Entering Regularisation constructor");
    if(!raster) {
	cerr << "No raster defined in Regularisation constructor\n";
	x0 = 0;
	rowptr = 0;
	colidx = 0;
    } else {
	slen = raster->SLen();
	glen = raster->GLen();
	dim  = raster->Dim();
	raster->NeighbourGraph (rowptr, colidx, nzero);
	CreateHessStruct1param(Hs1p);
	pdim = Hs1p.nCols();
	nset = _x0->Dim() / pdim;

	xASSERT (nset*pdim == _x0->Dim(), "Invalid parameter dimensions");
	x0 = new RVector (*_x0);
    }
    //LOGOUT("Leaving Regularisation constructor");
}

Regularisation::~Regularisation ()
{
    //LOGOUT("Entering Regularisation destructor");
    if (x0) delete x0;
    if (rowptr) delete []rowptr;
    if (colidx) delete []colidx;
    //LOGOUT("Leaving Regularisation destructor");
}

Regularisation *Regularisation::Create (ParamParser *pp, const RVector *_x0,
    const Raster *_raster, const RVector *_xs)
{
    char cbuf[256];
    int prior = -1;
    Regularisation *reg = NULL;

    if (pp->GetString ("PRIOR", cbuf) || pp->GetString ("LM_PRIOR", cbuf)) {
	if (!strcasecmp (cbuf, "NONE"))
	    prior = 0;
	else if (!strcasecmp (cbuf, "TV"))
	    prior = 1;
	else if (!strcasecmp (cbuf, "TK0"))
	    prior = 2;
	else if (!strcasecmp (cbuf, "TK1"))
	    prior = 3;
	else if (!strcasecmp (cbuf, "HUBER"))
	    prior = 4;
	else if (!strcasecmp (cbuf, "MRF"))
	    prior = 5;
    }
    while (prior < 0) {
	cout << "\nSelect a regularisation method:\n";
	cout << "(0) None\n";
	cout << "(1) Total variation (TV)\n";
	cout << "(2) Tikhonov 0th order (TK0)\n";
	cout << "(3) Tikhonov 1st order (TK1)\n";
	cout << "(4) Huber\n";
	cout << "(5) MRF\n";
	cout << "[0|1|2|3|4|5] >> ";
	cin >> prior;
    }

    switch (prior) {
    case 0: reg = new NullRegularisation(); break;
    case 1: reg = new TV (1, 1, _x0, _raster); break;
    case 2: reg = new Tikhonov0 (1, _x0, _xs); break;
    case 3: reg = new TK1 (1, _x0, _raster); break;
    case 4: reg = new Huber (1, 1, _x0, _raster); break;
    case 5: reg = new MRF(1, _x0, _raster); break;
    }
    if (reg) {
	reg->ReadParams (pp);
	reg->WriteParams (pp);
    }
    return reg;
}

void Regularisation::ReadParams (ParamParser *pp)
{
    // === Hyperparameter tau ===
    if (!pp->GetReal ("PRIOR_TAU", tau) && !pp->GetReal ("LM_TAU", tau)) {
	cout << "\nSelect regularisation hyperparameter tau (>= 0)\n";
	cout << ">> ";
	cin >> tau;
    }
}

void Regularisation::WriteParams (ParamParser *pp)
{
    pp->PutString ("PRIOR", GetName());

    // === Hyperparameter tau ===
    pp->PutReal ("PRIOR_TAU", tau);
}

void Regularisation::CreateHessStruct1param (RCompRowMatrix &Hess)
{
    // compute structure of sparse regularisation matrix as
    // Laplacian

    int row, col, i, i0, i1, cnt_idx, j;
    double valsum, *valptr;
    valptr = new double[nzero];

    for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];
#ifdef OLD_BIT
        int cntj;
	for (i = i0; i < i1; i++) { // find central pixel
	    col = colidx[i];
	    if (col == row) {
		cnt_idx = i;
		cntj = col;
		//		cntj = raster->Sol2Grid (col);
		break;
	    }
	}
#endif
	valsum = 0.0;
	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) {
	        j  = row;
	      //		j = raster->Sol2Grid (col);
		valptr[i] = -1;
		valsum += valptr[i];
	    }
	    else
		cnt_idx = i;      // find central pixel
	}
	valptr[cnt_idx] = -valsum;
	//	cout << valptr[cnt_idx] << " ";
    }
    Hess.New(slen,slen);
    Hess.Initialise (rowptr, colidx, valptr);
    //    cout << Hess << endl;

    delete []valptr;
}
#ifdef ROW_BY_ROW_HESSIAN
void Regularisation::CreateNonLinearHessian (const RVector &x, RCompRowMatrix &Hess) 
{
  //    cout << "Setting up Hessian" << endl;
    int i,j,k;


    int *ci = new int[slen];
    double *v = new double[slen];
    double *hv = new double[nzero];
    for(i = j = 0; i < slen; i++) {
      int nz;
      nz = GetHessianRow (x, i, ci, v);
      if (nz != rowptr[i+1]-rowptr[i])
	  cerr << "Inconsistent matrix structure in row " << i << endl;
      for (k = 0; k < nz; k++)
	  hv[j++] = v[k];
    }
    cout << "Allocate Hessian" << endl;
    Hess.New(slen,slen);
    Hess.Initialise (rowptr, colidx, hv);
    delete []ci;
    delete []v;
    delete []hv;
}
#endif
#ifdef REG_1ST_ORDER
// =========================================================================
// =--------------- Base class for Regularisation with Gradient -----------=
// =========================================================================
Regularisation_1stOrder::Regularisation_1stOrder (const Raster *_raster,const RVector *_x0) : Regularisation(_raster, _x0) {
    //LOGOUT("Entering Regularisation_1stOrder constructor");
    //LOGOUT("Leaving Regularisation_1stOrder constructor");
};


// =========================================================================
// =------------------------ Generic with PreSmoothing -----------------   =
// =========================================================================

#else
GenericSigma::GenericSigma (const Raster *_raster, const RVector *_x0,
    const int _npreg, double * _preg, double _tau, double _sdx, double _sdf,
    const void *_kref, bool _KRefTensFlag)
    : Regularisation (_raster, _tau, _x0), npreg(_npreg), sdx(_sdx), sdf(_sdf),
  kref( ( RVector *) (_KRefTensFlag ? 0 : _kref)),
  kreftens( ( RDenseMatrix *) (_KRefTensFlag ? _kref : 0)),
  KRefTensFlag(_KRefTensFlag)
#endif
{
    //LOGOUT("Entering GenericSigma constructor");
    //cout << "Entering GenericSigma constructor\n" ;

    if(npreg > 0) {
	preg = new double[npreg];
	for (int ip = 0; ip < npreg; ip++)
	    preg[ip] = _preg[ip];
    }
    
    if(sdx > 0.0) {
	LOGOUT("Smoothing Kappa ",sdx);
	cout << "Smoothing Kappa "<< sdx << endl;
    }
    if(sdf > 0.0) {
	LOGOUT("Smoothing image argument ",sdf);
	cout << "Smoothing image argument " << sdf << endl;
    }
    sdr = fT = 0.0;
    if( KRefTensFlag) {
	cout << "Building Neighbour Graph for Reference Tensor\n";
	raster->NeighbourGraph (NG);
      Df = new RCompRowMatrix [dim];
      Db = new RCompRowMatrix [dim];
      CreateDerivativeOperators ();    
    }
    KRefIsLocal = false;
    krefimg_name = NULL;
    //LOGOUT("Leaving GenericSigma constructor");
    //cout << "Leaving GenericSigma constructor\n";
}

GenericSigma::GenericSigma (const Raster *_raster, const RVector *_x0,
    int _npreg, double *_preg, double _tau, const RVector &_refimg,
    double _sdr, double _sdx, double _sdf, double _fT, bool _KRefTensFlag)
    : Regularisation (_raster, _tau, _x0), npreg(_npreg), sdx(_sdx),
      sdf(_sdf), sdr(_sdr), fT(_fT), KRefTensFlag (_KRefTensFlag)
{
    //cout << "Entering GenericSigma constructor\n" ;

    if(npreg > 0) {
	preg = new double[npreg];
	for (int ip = 0; ip < npreg; ip++)
	    preg[ip] = _preg[ip];
    }
    
    if(sdx > 0.0) {
	LOGOUT("Smoothing Kappa ",sdx);
	cout << "Smoothing Kappa "<< sdx << endl;
    }
    if(sdf > 0.0) {
	LOGOUT("Smoothing image argument ",sdf);
	cout << "Smoothing image argument " << sdf << endl;
    }
    if( KRefTensFlag) {
      cout << "Building Neighbour Graph for Reference Tensor\n";
      raster->NeighbourGraph (NG);
      Df = new RCompRowMatrix [dim];
      Db = new RCompRowMatrix [dim];
      CreateDerivativeOperators ();    

    }
 
    kref     = (KRefTensFlag ? NULL : MakeKref (_refimg, _sdr, _fT));
    kreftens = (KRefTensFlag ? MakeKrefTens (_refimg, _sdr, _fT) : NULL);
    KRefIsLocal = true;
    krefimg_name = NULL;

    //LOGOUT("Leaving GenericSigma constructor");
    //cout << "Leaving GenericSigma constructor\n";

}

GenericSigma::~GenericSigma ()
{
    if(npreg > 0)
	delete [] preg;

    if( KRefTensFlag){
	delete [] Df;
	delete [] Db;
    }

    if (KRefIsLocal) {
	if (kref) delete kref;
	if (kreftens) delete []kreftens;
    }

    if (krefimg_name)
	delete []krefimg_name;
}

void GenericSigma::ReadParams (ParamParser *pp)
{
    char cbuf[256];
    int cmd;

    Regularisation::ReadParams (pp);

    // === Image for reference kappa ===
    if (krefimg_name) {
	delete []krefimg_name;
	krefimg_name = NULL;
    }
    if (!pp->GetString ("PRIOR_KAPREFIMG", cbuf)) {
	for (cmd = -1; cmd < 0 || cmd > 1; ) {
	    cout << "\nImage for reference diffusivity:\n";
	    cout << "(0) None\n";
	    cout << "(1) From file\n";
	    cout << "[0|1] >> ";
	    cin >> cmd;
	}
	switch (cmd) {
	case 0:
	    cbuf[0] = '\0';
	    break;
	case 1:
	    cout << "Image data file name:\n>> ";
	    cin >> cbuf;
	    break;
	}
    }
    if (cbuf[0] && strcasecmp (cbuf, "NONE")) {
	krefimg_name = new char[strlen(cbuf)+1];
	strcpy (krefimg_name, cbuf);
    }

    // Build diffusivity fields from image
    if (krefimg_name) {
	ifstream ifs(krefimg_name);
	xASSERT (ifs.good(), "Diffusivity image not found");
	RVector krefimg;
	ifs >> krefimg;

	if (!pp->GetReal ("PRIOR_DIFFSCALE", sdr)) {
	    cout << "\nScaling factor for diffusivity image:\n>> ";
	    cin >> sdr;
	}
	if (!pp->GetReal ("PRIOR_PMTHRESHOLD", fT)) {
	    cout << "\nPM threshold value:\n>> ";
	    cin >> fT;
	}
	kref = (KRefTensFlag ? NULL : MakeKref (krefimg, sdr, fT));
	kreftens = (KRefTensFlag ? MakeKrefTens (krefimg, sdr, fT) : NULL);
    }
}

void GenericSigma::WriteParams (ParamParser *pp)
{
    Regularisation::WriteParams (pp);

    // === Image for reference kappa ===
    if (krefimg_name) {
	pp->PutString ("PRIOR_KAPREFIMG", krefimg_name);
	pp->PutReal ("PRIOR_DIFFSCALE", sdr);
	pp->PutReal ("PRIOR_PMTHRESHOLD", fT);
    } else {
	pp->PutString ("PRIOR_KAPREFIMG", "NONE");
    }
}

bool GenericSigma::SetLocalScaling (const RVector &scale_all)
{
    kref = SetKref (scale_all);
    KRefIsLocal = true;
    krefimg_name = NULL;
    return true;
}

RVector *GenericSigma::SetKref (const RVector &gkref_all)
{
    int j, sofs, gofs;
    int slen = raster->SLen();
    int glen = raster->GLen();
    int nprm = gkref_all.Dim() / glen;

    RVector *kref_all = new RVector(slen*nprm);
    for (j = sofs = gofs = 0; j < nprm; j++, sofs += slen, gofs += glen) {
        RVector gkref (gkref_all, gofs, glen);
	RVector kref (*kref_all, sofs, slen); // single parameter solution
	raster->Map_GridToSol (gkref, kref);
    }
    return kref_all;
}

RVector *GenericSigma::MakeKref (const RVector &gimgref_all, double sdr,
    double fT)
{
    const IVector &gdim = raster->GDim();

    int i, j, sofs, gofs;
    int dim = raster->Dim();
    int slen = raster->SLen();
    int glen = raster->GLen();
    int nprm = gimgref_all.Dim() / glen;
    xASSERT (gimgref_all.Dim() == glen*nprm, "Invalid image size");
    // create a "reference kappa". Can only do 2D for now.
    RVector simgref_all(slen*nprm);
    RVector *kref_all = new RVector(slen*nprm);
    RVector *kap1d  = new RVector[dim];
    for(int idim = 0; idim < dim; idim++) 
	kap1d[idim].New(gdim[idim]);

    for (j = sofs = gofs = 0; j < nprm; j++, sofs += slen, gofs += glen) {
	RVector kref (*kref_all, sofs, slen); // single parameter solution
	RVector gimgref(gimgref_all, gofs, glen);
	RVector simgref(simgref_all, sofs, slen);
	raster->Map_GridToSol (gimgref, simgref); // initial image
	kref = raster->ImageGradient(simgref,sdr); // arbitrary scale for now

	double max, T;
	max = -1;
	for(i = 0; i < slen; i++)
	    max = (max > kref[i] ? max : kref[i]);
	T = fT*max;
	cout << "ref image gradient threshold : " << T << endl;
	for(i = 0; i < slen; i++)
	    kref[i] = exp( -SQR( kref[i]/T));
    }
    return kref_all;
}

RDenseMatrix *GenericSigma::MakeKrefTens (const RVector &gimgref, double sdr,
    double fT)
{
    int i;
    int slen = raster->SLen();
    int dim = raster->Dim();
    RDenseMatrix *kreftens = new RDenseMatrix [2*slen];

    RVector simgref_all(slen*2);
    RVector simgref(simgref_all,0,slen); // 1st half of solution
    RVector simgref2(simgref_all,slen,slen); // 2nd half of solution
    raster->Map_GridToSol (gimgref, simgref); // initial image
    raster->Map_GridToSol (gimgref, simgref2); // initial image

    int ngim = (dim ==3 ? 10 : 6);
    //	bool *iflag = new bool [dimflag ==3 ? 10 : 6];
    bool *iflag = new bool[ngim]; 
    for (int ng =0; ng <  ngim; ng++) iflag[ng] = false;
    iflag[1] = iflag[2] = true;
    if(dim == 3)
	iflag[3] = true;

    RVector *kref = raster->ImageJet(simgref,sdr,iflag); 
    RVector *kref2 = raster->ImageJet(simgref2,sdr,iflag);
    RVector g1, g2;
    for (i = 0; i < dim ; i++) {
	g1 += kref[i+1]*kref[i+1];
        g2 += kref2[i+1]*kref2[i+1];
    }
    g1 = sqrt(g1);
    g2 = sqrt(g2);
    
    double max, T;
#ifdef LOCAL_TENSOR
    max = -1;
    for(i = 0; i < slen; i++)
	max = (max > g1[i] ? max : g1[i]);
    T = fT*max;
    cout << "ref image gradient threshold : " << T << endl;
    for(i = 0; i < slen; i++) {
	double gam = 1- exp( -SQR( g1[i]/T)); // gamma in [0,1]
	kreftens[i].New(dim,dim);
	kreftens[i].Identity(dim);
	kreftens[i](0,0) -= gam*SQR(kref[1][i]/g1[i]); // xx term
	kreftens[i](1,1) -= gam*SQR(kref[2][i]/g1[i]); // yy term
	kreftens[i](0,1) -= gam*(kref[1][i]*kref[2][i])/SQR(g1[i]); // xy term
	kreftens[i](1,0) -= gam*(kref[1][i]*kref[2][i])/SQR(g1[i]); // xy term
    }

    max = -1;
    for(i = 0; i < slen; i++)
	max = (max > g2[i] ? max : g2[i]);
    T = fT*max;
    cout << "ref image gradient threshold : " << T << endl;
    for(i = 0; i < slen; i++) {
	// 2D ONLY !!
	double gam = 1- exp( -SQR( g2[i]/T)); // gamma in [0,1]
	kreftens[slen+i].New(dim,dim);
	kreftens[slen+i].Identity(dim);
	kreftens[slen+i](0,0) -= gam*SQR(kref2[1][i]/g2[i]); // xx term
	kreftens[slen+i](1,1) -= gam*SQR(kref2[2][i]/g2[i]); // yy term
	kreftens[slen+i](0,1) -= gam*(kref2[1][i]*kref2[2][i])/SQR(g2[i]); // xy term
	kreftens[slen+i](1,0) -= gam*(kref2[1][i]*kref2[2][i])/SQR(g2[i]); // xy term
    }
#else
    // pass the image gradients and do the work in regularisation class
    max = -1;
    for(i = 0; i < slen; i++)
	max = (max > g1[i] ? max : g1[i]);
    T = fT*max;
    cout << "ref image gradient threshold : " << T << endl;
    for(i = 0; i < slen; i++) {
	kreftens[i].New(dim,dim);
	for(int idim = 0; idim < dim; idim++)
	    kreftens[i](idim,idim) = kref[1+idim][i];	 
	kreftens[i](0,1) = T;
    }
       
    max = -1;
    for(i = 0; i < slen; i++)
	max = (max > g2[i] ? max : g2[i]);
    T = fT*max;
    cout << "ref image gradient threshold : " << T << endl;
    for(i = 0; i < slen; i++) {
	kreftens[slen+i].New(dim,dim);
	for(int idim = 0; idim < dim; idim++) {
	    kreftens[slen+i](idim,idim) = kref[1+idim][i];
	}
	kreftens[slen+i](0,1) = T; 
    }
#endif
    delete []kref;
    delete []kref2;
    delete []iflag;
    return kreftens;
}
void GenericSigma::CreateDerivativeOperators ()
{
    // compute structure of sparse matrix operators Dx, Dy, Dz 

    int row, col, i, i0, i1, idim;
    double **valptrb = new double * [dim];
    double **valptrf = new double * [dim];
    for (idim = 0; idim < dim; idim++){
      valptrb[idim] = new double[nzero];
      valptrf[idim] = new double[nzero];
    }

    for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];

	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) {
	        IVector eij = raster->NeighbourShift (NG, row, col);
		for (idim = 0; idim < dim; idim++){
		  valptrb[idim][i] = (eij[idim] < 0 ? -1 : 0);        
		  valptrf[idim][i] = (eij[idim] > 0 ? -1 : 0);        
		}
	    }
	    else
	      for (idim = 0; idim < dim; idim++)
		valptrb[idim][i] = valptrf[idim][i] = 1;    //  central pixel
	}
    }

    for (idim = 0; idim < dim; idim++) {
      Db[idim].New(slen,slen);
      Db[idim].Initialise (rowptr, colidx, valptrb[idim]);
      Db[idim].Shrink();
//      cout << D[idim] << endl;
      delete [] valptrb[idim];
      Df[idim].New(slen,slen);
      Df[idim].Initialise (rowptr, colidx, valptrf[idim]);
      Df[idim].Shrink();
      delete [] valptrf[idim];
    }
    delete [] valptrb;
    delete [] valptrf;
/* the following are almost all the same, apart from boundary effects
   I.e. the product D^T D is almost same as Hessian, which is seen by
   the subtraction of the "exact" Hessian (found in base class)
*/
//   cout << Df[0] *Db[0] + Df[1]*Db[1] - Hs1p;
//   cout << transp(Db[0])*Db[0] + transp(Db[1])*Db[1] - Hs1p;
//   cout << Db[0]*Df[0] + Db[1]*Df[1] -Hs1p;
//   cout << transp(Df[0])*Df[0] + transp(Df[1])*Df[1] - Hs1p;
}
#ifdef NO_EXPLICIT_DIFF_MATRICES
double GenericSigma::GetValue (const RVector &x) const
{
    //LOGOUT("Entering GenericSigma::GetValue");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector dx (x - *x0);
    double result=0.0;

    int row, col, i, i0, i1;

    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	RVector spdx(pdx);
	if(sdf > 0.0) {
	    spdx = raster->SmoothImage(pdx,sdf);
	    cout << "Smooth in GetValue " << sdf << endl;
        }
	//LOGOUT("parameter loop : p = %d",p);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) { 
		    if(KRefTensFlag) {
			RDenseMatrix ktens(dim,dim);
			if(kreftens)
			    ktens =tinterp(kreftens, row,col);
			else
			    ktens.Identity(dim);
			double fij = spdx[row] - spdx[col];
			RVector eij = raster->RNeighbourShift (NG, row, col);
			// cout << kreftens[row] << kreftens[col] << ktens << endl;
			// cout <<raster->NeighbourShift (NG, row, col) << " "<< eij << endl;
			double val = fabs(fij)* sqrt(eij & (ktens * eij));
			result += func(val,npreg,preg);
		    } else {
			double kx;
			kx = (kref ? sinterp(*kref, row+p*pdim ,col+p*pdim): 1);
			// kx = (kref ?  0.5 * ((*kref)[row] + (*kref)[col]) : 1);
			result += kx* func(spdx[row] - spdx[col],npreg,preg) ;
		    }
		}
	    } // end of neighbours of this pixel
	} // end of pixels	
	//LOGOUT_2PRM("row %d result %f",row,result);
    } // end of parameter
    
    //LOGOUT("Leaving GenericSigma::GetValue");
    return tau*result;
}
#else
double GenericSigma::GetValue (const RVector &x) const
{
    //LOGOUT("Entering GenericSigma::GetValue");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector dx (x - *x0);
    double result=0.0;

    int row, col, i, i0, i1, idim;

    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	RVector spdx(pdx);
	if(sdf > 0.0) {
	    spdx = raster->SmoothImage(pdx,sdf);
	    cout << "Smooth in GetValue " << sdf << endl;
        }
	//LOGOUT("parameter loop : p = %d",p);
	if(KRefTensFlag) {
	  RVector *dp = new RVector [dim];  // for the gradients
	  for (idim = 0; idim < dim; idim++)
	    dp[idim] = Db[idim] * spdx; // backward difference
	  RVector g(dim);
	  RDenseMatrix ktens(dim,dim);
	  for( i = 0; i < slen; i++) {
	    if(kreftens)
	      ktens =tinterp(kreftens, i,i);
	    else
	      ktens.Identity(dim);
	    for (idim = 0; idim < dim; idim++)
	      g[idim] = dp[idim][i];
	    double val = sqrt(g & (ktens * g));
	    result += func(val,npreg,preg);
	  }
	  delete [] dp;
	}
	else { // do it the old way !
 	  for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) { 
			double kx;
			kx = (kref ? sinterp(*kref, row+p*pdim ,col+p*pdim): 1);
			// kx = (kref ?  0.5 * ((*kref)[row] + (*kref)[col]) : 1);
			result += kx* func(spdx[row] - spdx[col],npreg,preg) ;
		}
	    } // end of neighbours of this pixel
	} // end of pixels	
	  //LOGOUT_2PRM("row %d result %f",row,result);
      }
    } // end of parameter
    
    //LOGOUT("Leaving GenericSigma::GetValue");
    return tau*result;
}
#endif

#ifdef OLD_GETKAPPA
RVector GenericSigma::GetKappa (const RVector &x) const
{
    //LOGOUT("Entering GenericSigma::GetKappa");
    //cout << "Entering GenericSigma::GetKappa\n";
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  xx(x), ktot(x.Dim());

    int row, col, i, i0, i1;
    
    for (int p = 0; p < nset; p++) { // loop over parameter sets
        cout << "parameter set " << p << endl;
	RVector px (xx, p*pdim, pdim);
	RVector spx(px);
	if(sdx > 0.0) {
	  spx = raster->SmoothImage(spx,sdx);
	  cout << "Smooth in GetKappa " << sdx << endl;
	}
	RVector pkap (ktot, p*pdim, pdim);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    int ncount=0; // neighbour count
	    for (i = i0 ; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) {
		  pkap[row] += kfunc(spx[row] - spx[col],npreg,preg);
		  ncount ++;
		}
	    } // end of neighbours of this pixel
	    dASSERT (ncount >0, "Found a pixel with no neighbours");  
	    pkap[row] *= (kref ? (*kref)[row+p*pdim] : 1.0)/ncount; // average of edges.
	} // end of pixels	
    } // end of parameter

    //cout << "Leaving GenericSigma::GetKappa\n";
    //LOGOUT("Leaving GenericSigma::GetKappa");
    return ktot;
}
#else
#ifdef NO_EXPLICIT_DIFF_MATRICES
RVector GenericSigma::GetKappa (const RVector &x) const
{
    //LOGOUT("Entering GenericSigma::GetKappa");
    //cout << "Entering GenericSigma::GetKappa\n";
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  xx(x), ktot(x.Dim());

    int row, col, i, i0, i1;
    
    for (int p = 0; p < nset; p++) { // loop over parameter sets
      //        cout << "parameter set " << p << endl;
	RVector px (xx, p*pdim, pdim);
	RVector spx(px);
	if(sdx > 0.0) {
	  spx = raster->SmoothImage(spx,sdx);
	  cout << "Smooth in GetKappa " << sdx << endl;
	}
	RVector pkap (ktot, p*pdim, pdim);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    int ncount=0; // neighbour count
	    double *kvals = new double [i1-i0];
	    double kx;
	    for (i = i0 ; i < i1; i++) { // find neighbours
	      col = colidx[i];
	      if (col != row) {
//		  pkap[row] += kfunc(spx[row] - spx[col],npreg,preg);
 	       if(KRefTensFlag) {
		RDenseMatrix ktens(dim,dim);
		if(kreftens)
		      ktens =tinterp(kreftens, row,col);
		else
		      ktens.Identity(dim);
		double fij = spx[row] - spx[col];
		RVector eij = raster->RNeighbourShift (NG, row, col);
		//		cout << eij << endl;
		kx = sqrt (eij & (ktens * eij));
//		cout << eij << " " << ktens << " " << kx << endl;
		kx *= kfunc( fij*kx,npreg,preg);
	       }
	       else {
 //		kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
 		kx = (kref ? sinterp(*kref, row+p*pdim, col+p*pdim) : 1);
		kx *= kfunc(spx[row] - spx[col],npreg,preg);
	       }
	      //		pkap[row] += kx;	      }
	       kvals[ncount] = kx;
	       ncount ++;
	      }
	    } // end of neighbours of this pixel
	    dASSERT (ncount >0, "Found a pixel with no neighbours");  
	    //#define KAPPA_IS_ROW_AVERAGE
#ifdef KAPPA_IS_ROW_AVERAGE
	    for(int ik = 0; ik < ncount; ik++)
	      pkap[row] += kvals[ik]/ncount; // average of edges.
#else
	    double min = 1e6;
	    for(int ik = 0; ik < ncount; ik++)
	       min = (kvals[ik] < min ? kvals[ik] : min);
	    pkap[row] = min;
#endif
//	    cout << "Total " << ncount << "->" << pkap[row] << endl;
	    delete [] kvals;
	} // end of pixels	
    } // end of parameter

    //cout << "Leaving GenericSigma::GetKappa\n";
    //LOGOUT("Leaving GenericSigma::GetKappa");
    return ktot;
}
#else
RVector GenericSigma::GetKappa (const RVector &x) const
{
    //LOGOUT("Entering GenericSigma::GetKappa");
    //cout << "Entering GenericSigma::GetKappa\n";
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  xx(x), ktot(x.Dim());

    int row, col, i, i0, i1, idim;
    
    for (int p = 0; p < nset; p++) { // loop over parameter sets
      //        cout << "parameter set " << p << endl;
	RVector px (xx, p*pdim, pdim);
	RVector spx(px);
	if(sdx > 0.0) {
	  spx = raster->SmoothImage(spx,sdx);
	  cout << "Smooth in GetKappa " << sdx << endl;
	}
	RVector pkap (ktot, p*pdim, pdim);
	if(KRefTensFlag) {

	  RVector *dp = new RVector [dim];  // for the gradients
	  for (idim = 0; idim < dim; idim++)
	    dp[idim] = Db[idim] * spx; // backward difference
	  RVector g(dim);
	  RDenseMatrix ktens(dim,dim);

	  for( i = 0; i < slen; i++) {
	    if(kreftens)
	      ktens =tinterp(kreftens, i,i);
	    else
	      ktens.Identity(dim);
	    for (idim = 0; idim < dim; idim++)
	      g[idim] = dp[idim][i];
	    double val = sqrt(g & (ktens * g));
	    //	    cout << val << " ";
	    pkap[i] = kfunc(val,npreg,preg) * det(ktens); 
	  }
	  delete [] dp;
	}
	else {
	  for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    int ncount=0; // neighbour count
	    double *kvals = new double [i1-i0];
	    double kx;
	    for (i = i0 ; i < i1; i++) { // find neighbours
	      col = colidx[i];
	      if (col != row) {
//		  pkap[row] += kfunc(spx[row] - spx[col],npreg,preg);
 	       if(KRefTensFlag) {
		RDenseMatrix ktens(dim,dim);
		if(kreftens)
		      ktens =tinterp(kreftens, row,col);
		else
		      ktens.Identity(dim);
		double fij = spx[row] - spx[col];
		RVector eij = raster->RNeighbourShift (NG, row, col);
		//		cout << eij << endl;
		kx = sqrt (eij & (ktens * eij));
//		cout << eij << " " << ktens << " " << kx << endl;
		kx *= kfunc( fij*kx,npreg,preg);
	       }
	       else {
 //		kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
		kx = (kref ? sinterp(*kref, row+p*pdim, col+p*pdim) : 1);
		kx *= kfunc(spx[row] - spx[col],npreg,preg);
	       }
	      //		pkap[row] += kx;	      }
	       kvals[ncount] = kx;
	       ncount ++;
	      }
	    } // end of neighbours of this pixel
	    dASSERT (ncount >0, "Found a pixel with no neighbours");  
	    //#define KAPPA_IS_ROW_AVERAGE
#ifdef KAPPA_IS_ROW_AVERAGE
	    for(int ik = 0; ik < ncount; ik++)
	      pkap[row] += kvals[ik]/ncount; // average of edges.
#else
	    double min = 1e6;
	    for(int ik = 0; ik < ncount; ik++)
	       min = (kvals[ik] < min ? kvals[ik] : min);
	    pkap[row] = min;
#endif
//	    cout << "Total " << ncount << "->" << pkap[row] << endl;
	    delete [] kvals;
	  } // end of pixels	
	}
    } // end of parameter

    //cout << "Leaving GenericSigma::GetKappa\n";
    //LOGOUT("Leaving GenericSigma::GetKappa");
    return ktot;
}
#endif
#endif
void GenericSigma::SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap)
{
    //LOGOUT("Entering GenericSigma::SetHess1FromKappa");

    int row, col, i, i0, i1, cnt_idx, j;
    double valsum, *valptr;
    valptr = new double[nzero];

    for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];

	valsum = 0.0;
	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) {
	        j  = row;
		double kx;
	//		kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
      		kx = (kref ? sinterp(*kref,row,col) : 1);
		kx *= 0.5*(kap[row] + kap[col]);
		valptr[i] = -kx;
		valsum += valptr[i];
	    }
	    else
		cnt_idx = i;      // find central pixel
	}
	valptr[cnt_idx] = -valsum;
	//	cout << valptr[cnt_idx] << " ";
    }
    Hess.New(slen,slen);
    Hess.Initialise (rowptr, colidx, valptr);
    //    cout << Hess << endl;
    Hess *= tau;
    delete []valptr;
    //LOGOUT("Leaving GenericSigma::SetHess1FromKappa");
}
#ifdef NO_EXPLICIT_DIFF_MATRICES
void GenericSigma::SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p)
{
    //LOGOUT("Entering GenericSigma::SetHess1");
    RVector dx (x - *x0);
    RVector pdx (dx, p*pdim, pdim);
    RVector spdx(pdx);
    if(sdx > 0.0) {
          spdx = raster->SmoothImage(pdx,sdx);
	  cout << "Smooth in SetHess1 " << sdx << endl;;
    }

    int row, col, i, i0, i1, cnt_idx, cntj, j;
    double valsum, *valptr;
    valptr = new double[nzero];

    for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];
	valsum = 0.0;
	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) {
	      //	        j  = row;
	      if(KRefTensFlag) {
		RDenseMatrix ktens(dim,dim);
		if(kreftens)
		      ktens =tinterp(kreftens, row,col);
		else
		      ktens.Identity(dim);
		double fij = spdx[row] - spdx[col];
		RVector eij = raster->RNeighbourShift (NG, row, col);
		//		cout << eij << endl;
		double kx = sqrt (eij & (ktens * eij));
		kx *= kfunc( fij*kx,npreg,preg);
		valptr[i] = -kx;
		valsum += valptr[i];	      }
	      else {
		double kx;
 //		kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
 		kx = (kref ? sinterp(*kref, row+p*pdim, col+p*pdim) : 1);
		kx *= kfunc(spdx[row] - spdx[col],npreg,preg);
		valptr[i] = -kx;
		valsum += valptr[i];
	      }
	    }
	    else
		cnt_idx = i;      // find central pixel
	}
	valptr[cnt_idx] = -valsum;
	//	cout << valptr[cnt_idx] << " ";
    }
    Hess1.New(slen,slen);
    //    cout << "Calling Hess1.Initialise(() \n";
    Hess1.Initialise (rowptr, colidx, valptr);
    //    cout << Hess << endl;
    Hess1 *= tau;

    delete []valptr;
    //LOGOUT("Leaving GenericSigma::SetHess1");
}
#else
void GenericSigma::SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p)
{
    //LOGOUT("Entering GenericSigma::SetHess1");
    RVector dx (x - *x0);
    RVector pdx (dx, p*pdim, pdim);
    RVector spdx(pdx);
    if(sdx > 0.0) {
          spdx = raster->SmoothImage(pdx,sdx);
	  cout << "Smooth in SetHess1 " << sdx << endl;;
    }
    int row, col, i, i0, i1, cnt_idx, idim, jdim;
    double valsum, *valptr;

    Hess1.New(slen,slen);
#define KDIAG_IS_VECTOR
    if(KRefTensFlag) {

	  RVector *dp = new RVector [dim];  // for the gradients
	  for (idim = 0; idim < dim; idim++)
	    dp[idim] = Db[idim] * spdx; // backward difference
	  RVector g(dim);
	  RDenseMatrix ktens(dim,dim);
#ifdef KDIAG_IS_VECTOR
	  RVector **Kdiags = new  RVector * [dim];
	  //	  cout << "Creating some diagonal matrices\n";
	  for (idim = 0; idim < dim; idim++) {
	    Kdiags[idim] = new  RVector [dim];
	    for (jdim = 0; jdim < dim; jdim++) {
	      Kdiags[idim][jdim].New(slen);
	    }
	  }
#else
	  RCompRowMatrix **Kdiags = new  RCompRowMatrix * [dim];
	  //	  cout << "Creating some diagonal matrices\n";
	  for (idim = 0; idim < dim; idim++) {
	    Kdiags[idim] = new  RCompRowMatrix [dim];
	    for (jdim = 0; jdim < dim; jdim++) {
	      Kdiags[idim][jdim].New(slen,slen);
	      Kdiags[idim][jdim].Identity(slen);
	    }
	  }
#endif
	  //	  cout << "filling the terms\n";
	  for( i = 0; i < slen; i++) {
	    if(kreftens)
	      ktens =tinterp(kreftens, i,i);
	    else
	      ktens.Identity(dim);
	    for (idim = 0; idim < dim; idim++)
	      g[idim] = dp[idim][i];
	    double val = sqrt(g & (ktens * g));
	    val = kfunc(val,npreg,preg);
	    for (idim = 0; idim < dim; idim++){
	      for (jdim = 0; jdim < dim; jdim++) {
#ifdef KDIAG_IS_VECTOR
		Kdiags[idim][jdim][i] = val*ktens(idim,jdim);
#else
		Kdiags[idim][jdim](i,i) = val*ktens(idim,jdim);
#endif
	      }
	    }
	  }
	  cout << "summing terms to create Hessian\n";
	  for (idim = 0; idim < dim; idim++) {
	    cout << idim << ": ";cout.flush();
	    for (jdim = 0; jdim < dim; jdim++) {
	        cout << jdim << " ";cout.flush();
#ifdef KDIAG_IS_VECTOR
		RCompRowMatrix tmp = Db[jdim];
		tmp.RowScale(Kdiags[idim][jdim]);
		tmp = Df[idim] * tmp;
		Hess1 = Hess1 + tmp;
#else
		Hess1 = Hess1 + Df[idim] * Kdiags[idim][jdim] * Db[jdim];
#endif
	    }
	    cout << endl;
	  }
	  for (idim = 0; idim < dim; idim++) {
	    delete [] Kdiags[idim];
	  }
	  delete [] Kdiags;
	  delete [] dp;
    }
    else {
      valptr = new double[nzero];
      for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];
	valsum = 0.0;
	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) {
	      //	        j  = row;
	      if(KRefTensFlag) {
		RDenseMatrix ktens(dim,dim);
		if(kreftens)
		      ktens =tinterp(kreftens, row,col);
		else
		      ktens.Identity(dim);
		double fij = spdx[row] - spdx[col];
		RVector eij = raster->RNeighbourShift (NG, row, col);
		//		cout << eij << endl;
		double kx = sqrt (eij & (ktens * eij));
		kx *= kfunc( fij*kx,npreg,preg);
		valptr[i] = -kx;
		valsum += valptr[i];	      }
	      else {
		double kx;
 //		kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
 		kx = (kref ? sinterp(*kref, row+p*pdim, col+p*pdim) : 1);
		kx *= kfunc(spdx[row] - spdx[col],npreg,preg);
		valptr[i] = -kx;
		valsum += valptr[i];
	      }
	    }
	    else
		cnt_idx = i;      // find central pixel
	}
	valptr[cnt_idx] = -valsum;
	//	cout << valptr[cnt_idx] << " ";
      }
      //    cout << "Calling Hess1.Initialise(() \n";
      Hess1.Initialise (rowptr, colidx, valptr);
      //    cout << Hess << endl;
      delete []valptr;
    }
    Hess1 *= tau;
    //LOGOUT("Leaving GenericSigma::SetHess1");
}
#endif
RVector GenericSigma::GetHess1f(const RVector &x,const RVector &f) const
{
  //LOGOUT("Entering GenericSigma::GetHess1f");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
    dASSERT (f.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
    
    RVector dx (x - *x0), grad(x.Dim());

    int row, col, i, i0, i1;
    
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	RVector spdx(pdx);
	if(sdx > 0.0) {
	  spdx = raster->SmoothImage(pdx,sdx);
	  cout << "Smooth in GetHess1f (x) " << sdx << endl;
        }
	RVector pf (f, p*pdim, pdim);
	RVector spf (pf); // it must be a copy else will corrupt f
	if(sdf > 0.0 ) {
	  spf = raster->SmoothImage(pf,sdf);
	  cout << "Smooth in GetHess1f (f) " << sdf << endl;
        }

	RVector pgd (grad, p*pdim, pdim);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) {
		  if(KRefTensFlag) {
		    RDenseMatrix ktens(dim,dim);
		    if(kreftens)
		      ktens =tinterp(kreftens, row,col);
		    else
		      ktens.Identity(dim);
		    double fij = spdx[row] - spdx[col];
		    RVector eij = raster->RNeighbourShift (NG, row, col);
		    //		    cout << eij << endl;
		    double kx = sqrt (eij & (ktens * eij));
		    kx *= kfunc( fij*kx,npreg,preg);
		    pgd[row] += kx*(spf[row] - spf[col]);
		  }
		  else {
		    double kx;
 //		  kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
		    kx = (kref ? sinterp(*kref, row+p*pdim, col+p*pdim) : 1);
		    kx *= kfunc(spdx[row] - spdx[col],npreg,preg);
		    pgd[row] += kx*(spf[row] - spf[col]); 
		  }
		}
	    } // end of neighbours of this pixel
	} // end of pixels	
    } // end of parameter

    //LOGOUT("Leaving GenericSigma::GetHess1f");
    return grad * tau;
}

void GenericSigma::SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p)
{
    LOGOUT("GenericSigma::SetFullHess not implemented yet");
}
RVector GenericSigma::GetFullHessf (const RVector &x, const RVector &f) const
{
    LOGOUT("GenericSigma::GetFullHessf not implemented yet");
    return RVector();
}

int GenericSigma::GetHessianRow (const RVector &x, int row, idxtype *colidx,
    double *val) const
  // DEFINITELY need to check with Martin !
{
    //  LOGOUT("Entering GenericSigma::GetHessianRow");
    int i,j, col, cnt_idx, cntj, nz, row0, ofs;

    RVector sx(x);

    for(int p = 0; p < nset; p++) { // produced smoothed version of image
      RVector spx(sx, p*pdim, pdim);
      if(sdx>0.0) {
	spx = raster->SmoothImage(spx,sdx); 
	cout << "Hess ";
      }
    }

    row0 = row%pdim;         // reduce row index
    ofs = (row/pdim) * pdim; // parameter offset
    nz = Hs1p.SparseRow (row0, colidx, val); // structure of row
    
    for (i = 0; i < nz; i++) { // find central pixel
	col = colidx[i];
	if (col == row0) {
	    val[i] = 0.0;// reset value of diagonal pixel
	    cnt_idx = i;
	    cntj = col;
	    //cntj = raster->Sol2Grid (col);
	    break;
	}
    } // found central pixel

    for (i = 0; i < nz; i++) {
	col = colidx[i];
	if (col != row0) {
	    j = col;
#ifdef REGUL_FULLHESS
	    double hx = d2func(sx[cntj+ofs] - sx[j+ofs],npreg,preg);
#else
	    double hx = kfunc(sx[cntj+ofs] - sx[j+ofs],npreg,preg);
#endif
	    double kx;
     //	    kx = (kref ? 0.5 * ((*kref)[cntj+ofs] + (*kref)[j+ofs]) : 1);
	    kx = (kref ? sinterp(*kref, cntj+ofs, j+ofs) : 1);
 	    val[i] = -kx*hx ;
	    val[cnt_idx] -= val[i] ;
	}
    }
    /*  output values;
        for (i = 0; i < nz ; i++) cout << val[i] << " ";
        cout << endl;
    */
    // offset coefficient set
    
    for (i = 0; i < nz; i++) colidx[j] += ofs;
    //  LOGOUT("Leaving GenericSigma::GetHessianRow");
    return nz;
}

RVector GenericSigma::GetHessianDiag (const RVector &x) const
{
    //LOGOUT("Entering GenericSigma::GetHessianDiag");
  // DEFINITELY need to check with Martin !
    RVector dx (x - *x0);


    int row, col, i, i0, i1, cntj, j;
    RVector diag(nset*pdim); 
    //    diag =1e-6 ;//set default value so no zero divides

    for (int p = 0; p < nset; p++){// loop over parameter sets
      RVector pdx (dx, p*pdim, pdim);
      if(sdx> 0.0) {
	pdx = raster->SmoothImage(pdx,sdx);
	cout << "HessD " << sdx << endl;
      }
      for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];
	cntj = row;

	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) { 
		j = col;
#ifdef REGUL_FULLHESS
 		double hx = d2func(pdx[row] - pdx[col],npreg,preg);
#else
 		double hx = kfunc(pdx[row] - pdx[col],npreg,preg);
#endif
		double kx;
//		kx = (kref ? 0.5 * ((*kref)[row] + (*kref)[col]) : 1);
		kx = (kref ? sinterp(*kref,row+p*pdim,col+p*pdim) : 1);
               	diag[p*pdim+row] += kx*hx ;
	    }
	} // end of neighbours of this pixel
      } // end of pixels
      //LOGOUT_3PRM("row %d pixel %d result %f",row-1, cntj,diag[p*pdim+row-1])
    } // end of set

    //LOGOUT("Leaving GenericSigma::GetHessianDiag");
    return tau*diag;
}


// ==========================================================================

void TVSigma::ReadParams (ParamParser *pp)
{
    GenericSigma::ReadParams (pp);

    double beta;
    if (!pp->GetReal ("TV_BETA", beta)) {
	cout << "\nSelect TV regularisation parameter beta:\n>> ";
	cin >> beta;
    }
    preg[0] = beta;
}

void TVSigma::WriteParams (ParamParser *pp)
{
    GenericSigma::WriteParams (pp);
    pp->PutReal ("TV_BETA", preg[0]);
}

// ==========================================================================

void HuberSigma::ReadParams (ParamParser *pp)
{
    GenericSigma::ReadParams (pp);

    double eps;
    if (!pp->GetReal ("HUBER_EPS", eps)) {
	cout << "\nSelect Huber regularisation parameter eps:\n>> ";
	cin >> eps;
    }
    preg[0] = eps;
}

void HuberSigma::WriteParams (ParamParser *pp)
{
    GenericSigma::WriteParams (pp);
    pp->PutReal ("HUBER_EPS", preg[0]);
}

// ==========================================================================

GenericScaleSpace::GenericScaleSpace (const Raster *_raster,
    const RVector *_x0, PSIREG _func, PSIREG _kfunc, PSIREG _k1func,
    const int _npreg, double * _preg, double _tau, double _sdx, double _sdf,
    const RDenseMatrix *_kreftens)
    : Regularisation (_raster, _tau, _x0), func(_func), kfunc(_kfunc),
      k1func(_k1func), npreg(_npreg), sdx(_sdx), sdf(_sdf), kreftens(_kreftens)
{
    //LOGOUT("Entering GenericScaleSpace constructor");
    //cout << "Entering GenericScaleSpace constructor\n" ;

    if(npreg > 0) {
      preg = new double[npreg];
      for (int ip = 0; ip < npreg; ip++)
	preg[ip] = _preg[ip];
    }

    if(sdx > 0.0) {
      LOGOUT("Smoothing Kappa ",sdx);
      cout << "Smoothing Kappa "<< sdx << endl;
    }
    if(sdf > 0.0) {
      LOGOUT("Smoothing image argument ",sdf);
      cout << "Smoothing image argument " << sdf << endl;
    }
    // allocate space for diffusion Tensor
    KappaTens = new RDenseMatrix * [nset];
    KappaSet = false;
    //LOGOUT("Leaving GenericScaleSpace constructor");
    //cout << "Leaving GenericScaleSpace constructor\n";
}
GenericScaleSpace::~GenericScaleSpace ()
{
    //LOGOUT("Entering GenericScaleSpace destructor");
    if(npreg > 0)
      delete [] preg;
    for(int p = 0; p < nset; p++) {
      if(KappaTens[p]) {
	delete [] KappaTens[p];
      }
     }
    delete [] KappaTens;
    //LOGOUT("Leaving GenericScaleSpace destructor");
}

double GenericScaleSpace::GetValue (const RVector &x) const
{

    //LOGOUT("Entering GenericSigma::GetValue");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    static int dimflag = raster->Dim();
    static int ngim =  (dimflag == 3 ? 10 : 6);
    bool *iflgs = new bool[ngim];
    for (int ib = 0; ib < ngim; ib++) iflgs[ib] = false;
    iflgs[1] = iflgs[2] = true; // x and y derivative
    if (dimflag ==3)
      iflgs[3] = true;          // z derivative


    RVector dx (x - *x0);
    RVector pgrad(dimflag);  // pixel gradient

    double result=0.0, val;

    int i;
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	RVector *spdx_jet;

	spdx_jet = raster->ImageJet(pdx, sdf, iflgs);
	cout << "Smooth in GetValue " << sdf << endl;
	//LOGOUT("parameter loop : p = %d",p);
	for (i = 0; i < pdx.Dim() ; i++) {
	  for(int k = 0; k < dimflag; k++)
	    pgrad[k] = spdx_jet[k+1][i];  // gradient of image
	  if(kreftens)
	   val = pgrad & (kreftens[ p*pdim + i] * pgrad ); //quadratic inner product
	  else
	   val= pgrad&pgrad; // scalar product
	  result += ( *func)(sqrt(val),npreg,preg);
	} // end of pixels	
	delete [] spdx_jet;
    } // end of parameter
    delete []iflgs;

    return tau*result;
}

/*------ Get Gradient implemented inline in  header */
//RVector GenericScaleSpace::GetGradient (const RVector &x) const
//{
//    return (x - *x0)  * (2.0*tau);
//}

RVector GenericScaleSpace::GetKappa    (const RVector &x) const
{
    //LOGOUT("Entering GenericScaleSpace::GetKappa");
    //cout << "Entering GenericScaleSpace::GetKappa\n";
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  ktot(x.Dim());
    static int dimflag = raster->Dim();
    static int ngim =  (dimflag == 3 ? 10 : 6);
    bool *iflgs = new bool[ngim];
    for(int ib = 0; ib < ngim; ib++) iflgs[ib] = false;
    iflgs[1] = iflgs[2] = true; // x and y derivative
    if (dimflag ==3)
      iflgs[3] = true;          // z derivative


    RVector dx (x - *x0);
    RVector pgrad(dimflag);  // pixel gradient

    double val;

    int i;
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (x, p*pdim, pdim);
	RVector pkap (ktot, p*pdim, pdim);
	RVector *spdx_jet;

	spdx_jet = raster->ImageJet(pdx, sdx, iflgs);
	cout << "Smooth in GetKappa " << sdx << endl;
	//LOGOUT("parameter loop : p = %d",p);
	for (i = 0; i < pdx.Dim() ; i++) {
	  for(int k = 0; k < dimflag; k++)
	    pgrad[k] = spdx_jet[k+1][i];  // gradient of image
	  if(kreftens)
	   val = pgrad & (kreftens[ p*pdim + i] * pgrad ); //quadratic inner product
	  else
	   val= pgrad&pgrad; // scalar product
	  pkap[i] = ( *kfunc)(sqrt(val),npreg,preg);
	} // end of pixels	
	delete [] spdx_jet;
    } // end of parameter
    //cout << "Leaving GenericScaleSpace::GetKappa\n";
    //LOGOUT("Leaving GenericScaleSpace::GetKappa");
    delete []iflgs;
    return ktot;    
}

#ifdef EXPLICIT_IMAGE_JET
RVector GenericScaleSpace::GetHess1f(const RVector &x,const RVector &f) const
{

  // need : order 2 jet of x and f
  // This could be made much more efficient if we know that x is fixed...

    //LOGOUT("Entering GenericScaleSpace::GetHess1f");
    //cout << "Entering GenericScaleSpace::GetHess1f\n";
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  h(x.Dim());
    static int dimflag = raster->Dim();
    static int ngim =  (dimflag == 3 ? 10 : 6);
    bool iflgsx[ngim],iflgsf[ngim]; 
    for(int ib = 0; ib < ngim; ib++) {
      iflgsx[ib] = iflgsf[ib] = true;
    }
    iflgsx[0] = iflgsf[0] = false;

    RVector dx (x - *x0);
    RVector pgradx(dimflag), pgradf(dimflag);  // pixel gradients
    RDenseMatrix phessx(dimflag,dimflag);      // pixel Hessian;

    double result=0.0, valx, valf, val0, val1;

    int i,k;
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (x, p*pdim, pdim);
	RVector pf (f, p*pdim, pdim);
	RVector ph (h, p*pdim, pdim);
	RVector *spdx_jet, *spf_jet;

	spdx_jet = raster->ImageJet(pdx, sdx, iflgsx);
	spf_jet = raster->ImageJet(pf, sdf, iflgsf);
	cout << "Smooth in GetHess1f " << sdx << endl;
	//LOGOUT("parameter loop : p = %d",p);
	for (i = 0; i < pdx.Dim() ; i++) {
	  for(k = 0; k < dimflag; k++) {
	    pgradx[k] = spdx_jet[k+1][i];  // gradient of image
	    pgradf[k] = spf_jet[k+1][i];  // gradient of image
	    phessx(k,k) = spdx_jet[k+1+dimflag][i];
	  }
	  if(dimflag == 3) {
	    phessx(0,1) = phessx(1,0) = spdx_jet[7][i];	    
	    phessx(0,2) = phessx(2,0) = spdx_jet[8][i];	    
	    phessx(1,2) = phessx(2,1) = spdx_jet[9][i];	    
	  }
	  else {
	    phessx(0,1) = phessx(1,0) = spdx_jet[5][i];
	  }
	  if(kreftens) { //quadratic inner product
	   valx = pgradx & (kreftens[ p*pdim + i] * pgradx ); 
	  }
	  else {
	   valx= pgradx&pgradx; // scalar product
	  }
	  val0 = 0;

	  if(kreftens) {
	    for(k = 0; k < dimflag; k++)
		val0 += spf_jet[k+1+dimflag][i]*kreftens[ p*pdim + i](k,k);
	    if(dimflag == 3) {
	      val0 += 2*spf_jet[7][i]*kreftens[ p*pdim + i](0,1);
	      val0 += 2*spf_jet[8][i]*kreftens[ p*pdim + i](0,2);
	      val0 += 2*spf_jet[9][i]*kreftens[ p*pdim + i](1,2);
	    }
	    else {
	      val0 += 2*spf_jet[5][i]*kreftens[ p*pdim + i](0,1);
	    }
	  }
	  else { // Dtensor is identity
	    for(k = 0; k < dimflag; k++) 
	      val0 += spf_jet[k+1+dimflag][i];
	  }
	  val0 *= ( *kfunc)(sqrt(valx),npreg,preg);

	  if(kreftens) { // rotate gradient
	   pgradx = (kreftens[ p*pdim + i] * pgradx); 
	   pgradf = (kreftens[ p*pdim + i] * pgradf); 
	  }
	  val1 =pgradf & (phessx * pgradx); //quadratic inner product
	  val1 *= ( *k1func)(sqrt(valx),npreg,preg);
	  ph[i] = val0 + val1;
	} // end of pixels	
	delete [] spdx_jet;
	delete [] spf_jet;
    } // end of parameter
    //cout << "Leaving GenericScaleSpace::GetHess1f\n";
    //LOGOUT("Leaving GenericScaleSpace::GetHess1f");
    return h;    
}
#else
//RVector GenericScaleSpace::GetHess1f_KappaSet(const RVector &f) const
//{
//  if(KappaSet)
//    return GetHess1f_KappaSet(f);
//  else
//    return GetHess1f_KappaNotSet(f);
//}
RVector GenericScaleSpace::GetHess1f_KappaSet(const RVector &f) const
{

  // This assumes x is fixed and KappaTens is set
  // this version takes div D grad rather than just scale space derivatives.

    //LOGOUT("Entering GenericScaleSpace::GetHess1f");
    //    cout << "Entering GenericScaleSpace::GetHess1f\n";
    //dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  h(f.Dim());
    static int dimflag = raster->Dim();
    static int ngim =  (dimflag == 3 ? 10 : 6);
    bool *iflgsf = new bool[ngim]; 
    iflgsf[0] = false;
    for(int ib = 1; ib < dimflag+1; ib++) {  // only the first derivatives
      iflgsf[ib] = true;
    }
    for(int ib = dimflag+1; ib < ngim; ib++) {
      iflgsf[ib] = false;
    }

    RVector ff (f);
    RVector pgradf(dimflag);  // pixel gradients
    int  local_nset = f.Dim() / pdim; // may not be whole image


    int i,k;
    for (int p = 0; p < local_nset; p++) { // loop over parameter sets
	//LOGOUT("parameter loop : p = %d",p);
	//	cout << "parameter loop : p = " << p << endl;
	RVector pf (ff, p*pdim, pdim);
	RVector ph (h, p*pdim, pdim);
	RVector *spf_jet;

	spf_jet = raster->ImageJet(pf, sdf, iflgsf);
	//	cout << "Smooth in GetHess1f " << sdx << endl;

	for (i = 0; i < pf.Dim() ; i++) {
	  for(k = 0; k < dimflag; k++) {
	    pgradf[k] = spf_jet[k+1][i];  // gradient of image
	  }
	  pgradf = KappaTens[p][i] * pgradf; // rotate
	  for(k = 0; k < dimflag; k++) {  // pjet is rotated and scaled
	    spf_jet[k+1][i] = pgradf[k]; 
	  }	
	} // end of pixels;
	//	cout << "done tensor, now for divergence\n";
	RVector *spfx, *spfy, *spfz;
	iflgsf[2] = false; iflgsf[3] = false;
	spfx =  raster->ImageJet(spf_jet[1], sdf, iflgsf);
	iflgsf[1] = false; iflgsf[2] = true; 
	spfy =  raster->ImageJet(spf_jet[2], sdf, iflgsf);
	ph = spfx[1] + spfy[2];
	//	cout << "ph successful \n" ;
	delete [] spfx;
	//	cout << "deleted spfx \n" ;
	delete [] spfy;
	//	cout << "deleted spfy \n" ;
	if(dimflag==3) {
	  iflgsf[2] = false; iflgsf[3] = true;
	  spfz =  raster->ImageJet(spf_jet[3], sdf, iflgsf);
	  ph += spfz[3];
      	  delete [] spfz;
	  //	  cout << "deleted spfz \n" ;
	}
       	delete [] spf_jet;
	//	cout << "deleted spf_jet \n" ;
	iflgsf[1] = true; iflgsf[2] = true;  // reset
    } // end of parameter
    //    cout << "Leaving GenericScaleSpace::GetHess1f\n";
    //LOGOUT("Leaving GenericScaleSpace::GetHess1f");
    delete []iflgsf;
    return -tau*h;    
}

RVector GenericScaleSpace::GetHess1f_KappaNotSet(const RVector &x,const RVector &f) const
{

  // need : order 2 jet of x and f
  // This could be made much more efficient if we know that x is fixed...
  // this version takes div D grad rather than just scale space derivatives.

    //LOGOUT("Entering GenericScaleSpace::GetHess1f");
    //    cout << "Entering GenericScaleSpace::GetHess1f\n";
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    RVector  h(x.Dim());
    static int dimflag = raster->Dim();
    static int ngim =  (dimflag == 3 ? 10 : 6);
    bool *iflgsx = new bool[ngim];
    bool *iflgsf = new bool[ngim]; 
    iflgsx[0] = iflgsf[0] = false;
    for(int ib = 1; ib < dimflag+1; ib++) {  // only the first derivatives
      iflgsx[ib] = iflgsf[ib] = true;
    }
    for(int ib = dimflag+1; ib < ngim; ib++) {
      iflgsx[ib] = iflgsf[ib] = false;
    }

    RVector dx (x - *x0);
    RVector ff (f);
    RVector pgradx(dimflag), pgradf(dimflag);  // pixel gradients
    RDenseMatrix PsiD(dimflag,dimflag);
    PsiD.Identity();
    int  local_nset = f.Dim() / pdim; // may not be whole image

    double valx;

    int i,k;
    for (int p = 0; p < local_nset; p++) { // loop over parameter sets
	//LOGOUT("parameter loop : p = %d",p);
	//	cout << "parameter loop : p = " << p << endl;
	RVector pdx (x, p*pdim, pdim);
	RVector pf (ff, p*pdim, pdim);
	RVector ph (h, p*pdim, pdim);
	RVector *spdx_jet, *spf_jet;

	spdx_jet = raster->ImageJet(pdx, sdx, iflgsx);
	spf_jet = raster->ImageJet(pf, sdf, iflgsf);
	//	cout << "Smooth in GetHess1f " << sdx << endl;

	for (i = 0; i < pdx.Dim() ; i++) {
	  for(k = 0; k < dimflag; k++) {
	    pgradx[k] = spdx_jet[k+1][i];  // gradient of image
	    pgradf[k] = spf_jet[k+1][i];  // gradient of image
	  }
	  if(kreftens) { //quadratic inner product
	   valx = pgradx & (kreftens[ p*pdim + i] * pgradx ); 
	   PsiD = kreftens[ p*pdim + i] * kreftens[ p*pdim + i]; //CAREFUL!
	  }
	  else {
	   valx= pgradx&pgradx; // scalar product
	  }
	   pgradf = (*kfunc)(sqrt(valx),npreg,preg)*(PsiD * pgradf); // rotate

	  for(k = 0; k < dimflag; k++) {  // pjet is rotated and scaled
	    spf_jet[k+1][i] = pgradf[k]; 
	  }	
	} // end of pixels;
	delete [] spdx_jet;
	//	cout << "done tensor, now for divergence\n";
	RVector *spfx, *spfy, *spfz;
	iflgsf[2] = false; iflgsf[3] = false;
	spfx =  raster->ImageJet(spf_jet[1], sdf, iflgsf);
	iflgsf[1] = false; iflgsf[2] = true; 
	spfy =  raster->ImageJet(spf_jet[2], sdf, iflgsf);
	ph = spfx[1] + spfy[2];
	//	cout << "ph successful \n" ;
	delete [] spfx;
	//	cout << "deleted spfx \n" ;
	delete [] spfy;
	//	cout << "deleted spfy \n" ;
	if(dimflag==3) {
	  iflgsf[2] = false; iflgsf[3] = true;
	  spfz =  raster->ImageJet(spf_jet[3], sdf, iflgsf);
	  ph += spfz[3];
      	  delete [] spfz;
	  //	  cout << "deleted spfz \n" ;
	}
       	delete [] spf_jet;
	//	cout << "deleted spf_jet \n" ;
	iflgsf[1] = true; iflgsf[2] = true;  // reset
    } // end of parameter
    delete []iflgsx;
    delete []iflgsf;

    //    cout << "Leaving GenericScaleSpace::GetHess1f\n";
    //LOGOUT("Leaving GenericScaleSpace::GetHess1f");
    return -tau*h;    
}
#endif
void GenericScaleSpace::SetHess1 (RCompRowMatrix &Hess, const RVector &x, const int p)
{
  //  LOGOUT("SetHess1 not implemented for GenericScaleSpace class");
  /* In ScaleSpace class, Hessian is not created directly.
     Instead, create an array of dim X dim matrices for anisotropic diffusion
     We will use the Hessian Matrix class just for some scalings. This might
     be used in future for preconditioning for example 
  */
    Hess.Identity(Hess.nRows());
    if(KappaTens[p]) // clear out old
      delete [] KappaTens[p];
    KappaTens[p] = new RDenseMatrix [pdim];

  // assume that x is "full vector" from which we select parameter set p
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");  

    static int dimflag = raster->Dim();
    static int ngim =  (dimflag == 3 ? 10 : 6);
    bool *iflgsx = new bool[ngim]; 
    iflgsx[0]  = false;
    for(int ib = 1; ib < dimflag+1; ib++) {  // only the first derivatives
      iflgsx[ib] = true;
    }
    for(int ib = dimflag+1; ib < ngim; ib++) {
      iflgsx[ib]  = false;
    }

    RVector dx (x - *x0);
    RVector pgradx(dimflag);  // pixel gradients
    RDenseMatrix PsiD(dimflag,dimflag);
    PsiD.Identity();

    double valx;

    int i,k;

    //LOGOUT("parameter loop : p = %d",p);
    RVector pdx (x, p*pdim, pdim);
    RVector *spdx_jet;

    spdx_jet = raster->ImageJet(pdx, sdx, iflgsx);

    for (i = 0; i < pdx.Dim() ; i++) {
	  for(k = 0; k < dimflag; k++) {
	    pgradx[k] = spdx_jet[k+1][i];  // gradient of image
	  }
	  if(kreftens) { //quadratic inner product
	   valx = pgradx & (kreftens[ p*pdim + i] * pgradx ); 
	   PsiD = kreftens[ p*pdim + i] * kreftens[ p*pdim + i]; //CAREFUL!
	  }
	  else {
	   valx= pgradx&pgradx; // scalar product
	  }
   	  KappaTens[p][i] = PsiD* (*kfunc)(sqrt(valx),npreg,preg);
    } // end of pixels;
    delete []spdx_jet;
    delete []iflgsx;

    KappaSet = true;
}
void GenericScaleSpace::SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kappa)
{

  LOGOUT("SetHess1FromKappa not implemented for GenericScaleSpace class");
}
void GenericScaleSpace::SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p)
{
    LOGOUT("GenericScaleSpace::SetFullHess not implemented yet");
}
RVector GenericScaleSpace::GetFullHessf (const RVector &x, const RVector &f) const
{
    LOGOUT("GenericScaleSpace::GetFullHessf not implemented yet");
    return RVector();
}


int GenericScaleSpace::GetHessianRow (const RVector &x, int i, idxtype *colidx,
    double *val) const
{
    colidx[0] = i;  // diagonal only
    val[0] = 2.0*tau ;
    return 1;
}

RVector GenericScaleSpace::GetHessianDiag (const RVector &x) const
{
  RVector tmp;
  return tmp;
}

// ==========================================================================

MRF::MRF (double _tau, const RVector *_x0, const Raster *_raster,
    const RVector *_kap)
    : Regularisation (_raster, _tau, _x0)
{
    MRFa = 0;
    MRFb = 0;
    SetKaprefImg (_kap);
}

MRF::~MRF ()
{
    if (MRFa) delete MRFa;
    if (MRFb) delete MRFb;
}

void MRF::SetKaprefImg (const RVector *kap)
{
    if (MRFa) delete MRFa;
    if (MRFb) delete MRFb;

    // Build up the neighbour graph
    int i, j, k, m1, m2, bm1, bm2;
    int *rp, *ci, nz;
    int slen = raster->SLen();
    int blen = raster->BLen();
    int dim = raster->Dim();
    IVector bdim = raster->BDim();
    raster->NeighbourGraph (rp, ci, nz);
    MRFa = new RCompRowMatrix(slen,slen,rp,ci);
    MRFb = new RCompRowMatrix(slen,slen,rp,ci);
    RCompRowMatrix &rMRFa = *MRFa;
    RCompRowMatrix &rMRFb = *MRFb;
    
    const double *muaim = 0, *musim = 0;
    if (kap) {
	muaim = kap->data_buffer();
	musim = (kap->Dim() == blen*2 ? muaim + blen : 0);
    } else {
	double *sol = new double[slen];
	for (i = 0; i < slen; i++) sol[i] = 1.0;
	muaim = new double[raster->BLen()];
	musim = muaim;
	delete []sol;
    }

    double minmrf = 1e-2;
    double mrfdiff = 1.0 - minmrf;
    double diff, sum;
    
    int bx = bdim[0], by = bdim[1], bz = (dim == 3 ? bdim[2]:1);
    for (k = 0; k < bz; k++) {
	for (j = 0; j < by; j++) {
	    for (i = 0; i < bx; i++) {
		bm1 = i + (j + k*by)*bx;
		m1 = raster->Basis2Sol (bm1);
		if (m1 < 0) continue;
		if (i < bx-1) {
		    bm2 = i+1 + (j + k*by)*bx;
		    m2 = raster->Basis2Sol (bm2);
		    if (m2 >= 0) {
			if (muaim) {
			    diff = fabs(muaim[bm1]-muaim[bm2]);
			    sum  = muaim[bm1]+muaim[bm2];
			    rMRFa(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFa(m1,m2) = -minmrf - mrfdiff;
			if (musim) {
			    diff = fabs(musim[bm1]-musim[bm2]);
			    sum  = musim[bm1]+musim[bm2];
			    rMRFb(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFb(m1,m2) = -minmrf - mrfdiff;
		    }
		}
		if (i > 0) {
		    bm2 = i-1 + (j + k*by)*bx;
		    m2 = raster->Basis2Sol (bm2);
		    if (m2 >= 0) {
			if (muaim) {
			    diff = fabs(muaim[bm1]-muaim[bm2]);
			    sum  = muaim[bm1]+muaim[bm2];
			    rMRFa(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFa(m1,m2) = -minmrf - mrfdiff;
			if (musim) {
			    diff = fabs(musim[bm1]-musim[bm2]);
			    sum  = musim[bm1]+musim[bm2];
			    rMRFb(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFb(m1,m2) = -minmrf - mrfdiff;
		    }
		}
		
		if (j < by-1) {
		    bm2 = i + (j+1 + k*by)*bx;
		    m2 = raster->Basis2Sol (bm2);
		    if (m2 >= 0) {
			if (muaim) {
			    diff = fabs(muaim[bm1]-muaim[bm2]);
			    sum  = muaim[bm1]+muaim[bm2];
			    rMRFa(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFa(m1,m2) = -minmrf - mrfdiff;
			if (musim) {
			    diff = fabs(musim[bm1]-musim[bm2]);
			    sum  = musim[bm1]+musim[bm2];
			    rMRFb(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFb(m1,m2) = -minmrf - mrfdiff;
		    }
		}
		if (j > 0) {
		    bm2 = i + (j-1 + k*by)*bx;
		    m2 = raster->Basis2Sol (bm2);
		    if (m2 >= 0) {
			if (muaim) {
			    diff = fabs(muaim[bm1]-muaim[bm2]);
			    sum  = muaim[bm1]+muaim[bm2];
			    rMRFa(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFa(m1,m2) = -minmrf - mrfdiff;
			if (musim) {
			    diff = fabs(musim[bm1]-musim[bm2]);
			    sum  = musim[bm1]+musim[bm2];
			    rMRFb(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFb(m1,m2) = -minmrf - mrfdiff;
		    }
		}
		
		if (k < bz-1) {
		    bm2 = i + (j + (k+1)*by)*bx;
		    m2 = raster->Basis2Sol (bm2);
		    if (m2 >= 0) {
			if (muaim) {
			    diff = fabs(muaim[bm1]-muaim[bm2]);
			    sum  = muaim[bm1]+muaim[bm2];
			    rMRFa(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFa(m1,m2) = -minmrf - mrfdiff;
			if (musim) {
			    diff = fabs(musim[bm1]-musim[bm2]);
			    sum  = musim[bm1]+musim[bm2];
			    rMRFb(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFb(m1,m2) = -minmrf - mrfdiff;
		    }
		}
		if (k > 0) {
		    bm2 = i + (j + (k-1)*by)*bx;
		    m2 = raster->Basis2Sol (bm2);
		    if (m2 >= 0) {
			if (muaim) {
			    diff = fabs(muaim[bm1]-muaim[bm2]);
			    sum  = muaim[bm1]+muaim[bm2];
			    rMRFa(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFa(m1,m2) = -minmrf - mrfdiff;
			if (musim) {
			    diff = fabs(musim[bm1]-musim[bm2]);
			    sum  = musim[bm1]+musim[bm2];
			    rMRFb(m1,m2) = -minmrf - mrfdiff * (diff<0.1*sum ? 1:0);
			} else rMRFb(m1,m2) = -minmrf - mrfdiff;
		    }
		}
	    }
	}
    }

    // Populate the diagonal
    const double *pMRFa = MRFa->ValPtr();
    const double *pMRFb = MRFb->ValPtr();
    for (i = 0; i < slen; i++) {
	double suma = 0.0, sumb = 0.0;
	for (j = rp[i]; j < rp[i+1]; j++) {
	    suma += pMRFa[j];
	    sumb += pMRFb[j];
	}
	rMRFa(i,i) = -suma;
	rMRFb(i,i) = -sumb;
    }

    delete []rp;
    delete []ci;
    if (!kap)
	delete []muaim;
}

double MRF::GetValue (const RVector &x) const
{
    // DEBUG
    ofstream ofs1("dbg.dat");
    MRFa->ExportRCV(ofs1);

    int slen = raster->SLen();
    RVector xa(x, 0, slen);
    RVector xb(x, slen, slen);
    RVector res;
    MRFa->Ax(xa, res);
    double va = res & xa;
    MRFb->Ax(xb, res);
    double vb = res & xb;
    return tau * (va+vb);
}

RVector MRF::GetGradient (const RVector &x) const
{ // TO BE DONE
    return RVector();
}

RVector MRF::GetKappa (const RVector &x)  const
{ // DUMMY
    return RVector();
}

void MRF::SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p)
{ // DUMMY
}

void MRF::SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap)
{ // DUMMY
}

RVector MRF::GetHess1f (const RVector &x, const RVector &f) const
{ // DUMMY
    return RVector();
}

void MRF::SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p)
{ // DUMMY
}

RVector MRF::GetFullHessf (const RVector &x, const RVector &f) const
{ // DUMMY
    return RVector();
}

int MRF::GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)
    const
{ // DUMMY
    return 0;
}

RVector MRF::GetHessianDiag (const RVector &x) const
{ // DUMMY
    return RVector();
}

//#define OLD_REGUL2
#ifdef OLD_REGUL2
// ==========================================================================

double NullRegularisation::GetValue (const RVector &x) const
{
    return 0.0;
}

RVector NullRegularisation::GetGradient (const RVector &x) const
{
    return RVector (x.Dim());
}
RVector NullRegularisation::GetKappa    (const RVector &x) const
{
    return RVector (x.Dim());
}
RVector NullRegularisation::GetHess1f(const RVector &x,const RVector &f) const
{
    return RVector (x.Dim());
}
int NullRegularisation::GetHessianRow (const RVector &x, int i, idxtype *colidx,
    double *val) const
{
    return 0;
}

RVector NullRegularisation::GetHessianDiag (const RVector &x) const
{
    return RVector (x.Dim());
}

// ==========================================================================

Tikhonov0::Tikhonov0 (double _tau, const RVector *_x0, const RVector *_xs)
    : Regularisation (0, _tau, _x0)
{
    if (_x0)
	x0 = new RVector(*_x0);
    else
	x0 = NULL;

    if (_xs)
	xs = new RVector(*_xs);
    else
	xs = x0;
    use_tauvec = false;
    kaprefimg_name = NULL;
}

Tikhonov0::~Tikhonov0()
{
    if (xs && xs != x0) delete xs;
    if (kaprefimg_name) delete []kaprefimg_name;
}

void Tikhonov0::SetTau (double _tau)
{
    tau = _tau;
    use_tauvec = false;
}

void Tikhonov0::SetTau (const RVector &_tau)
{
    xASSERT(_tau.Dim() == x0->Dim(), "Wrong vector size (expected %d, got %d)",
	    x0->Dim(), _tau.Dim());
    xASSERT(vmax(_tau) >= 0.0, "Negative tau element detected");
    tauvec = _tau;
    use_tauvec = true;
}

double Tikhonov0::GetValue (const RVector &x) const
{
    RVector dx ((x - *x0) / *xs);

    if (use_tauvec) return ((dx*tauvec) & dx) * tau;
    else            return (dx & dx) * tau;
}

RVector Tikhonov0::GetGradient (const RVector &x) const
{
  if (use_tauvec) return ((x - *x0)*tauvec) / sqr (*xs) * (2.0*tau);
    else            return (x - *x0) / sqr (*xs) * (2.0*tau);
}

RVector Tikhonov0::GetKappa    (const RVector &x) const
{
  LOGOUT("Diffusivity is not defined for zero order Tikhonov");
  cerr << "Diffusivity is not defined for zero order Tikhonov\n";
    return RVector (x.Dim());
}

RVector Tikhonov0::GetHess1f(const RVector &x,const RVector &f) const
{
    if (use_tauvec) return (f*(2.0*tau)) * (tauvec/sqr(*xs));
    else return f*(2.0*tau)/sqr(*xs);
}

void Tikhonov0::SetHessian (RCompRowMatrix &Hess, const RVector &x)
{
    /*
      RDiagMatrix Ident(x.Dim());
      Hess = Iden;
      for(int i = 0; i < x.Dim(); i++) 
      Hess(i,i) = 1;
    */
    cerr << "SetHessian not implemented for Tikhonov0 class" << endl;
    LOGOUT("SetHessian not implemented for Tikhonov0 class");
}

void Tikhonov0::SetHessianFromKappa (RCompRowMatrix &Hess, const RVector &kappa){
    /*
      RDiagMatrix Ident(x.Dim());
      Hess = Iden;
      for(int i = 0; i < x.Dim(); i++) 
      Hess(i,i) = kappa(i);
    */
    cerr << "SetHessianFromKappa not implemented for Tikhonov0 class" << endl;
    LOGOUT("SetHessianFromKappa not implemented for Tikhonov0 class");
}


int Tikhonov0::GetHessianRow (const RVector &x, int i, idxtype *colidx,
    double *val) const
{
    colidx[0] = i;  // diagonal only
    if (use_tauvec) val[0] = (2.0*tau)*tauvec[i] / ((*xs)[i] * (*xs)[i]);
    else val[0] = 2.0*tau / ((*xs)[i] * (*xs)[i]);
    return 1;
}

RVector Tikhonov0::GetHessianDiag (const RVector &x) const
{
    if (use_tauvec) return (2.0*tau)*(tauvec/sqr(*xs));
    else            return (2.0*tau)/sqr(*xs);
}

void Tikhonov0::ReadParams (ParamParser *pp)
{
    int cmd;
    char cbuf[256];

    Regularisation::ReadParams (pp);

    // === Image for reference kappa ===
    if (!pp->GetString ("PRIOR_KAPREFIMG", cbuf)) {
	for (cmd = -1; cmd <0 || cmd > 1; ) {
	    cout << "\nImage for reference diffusivity:\n";
	    cout << "(0) None\n";
	    cout << "(1) From file\n";
	    cout << "[0|1] >> ";
	    cin >> cmd;
	}
	switch (cmd) {
	case 0:
	    cbuf[0] = '\0';
	    break;
	case 1:
	    cout << "Image data file name:\n>> ";
	    cin >> cbuf;
	    break;
	}
    }
    if (kaprefimg_name) {
	delete []kaprefimg_name;
	kaprefimg_name = NULL;
    }
    if (cbuf[0] && strcasecmp (cbuf, "NONE")) {
	RVector vtau;
	ifstream ifs (cbuf);
	ifs >> vtau;
	SetTau (vtau);
	kaprefimg_name = new char[strlen(cbuf)+1];
	strcpy (kaprefimg_name, cbuf);
    }
}

void Tikhonov0::WriteParams (ParamParser *pp)
{
    Regularisation::WriteParams (pp);

    pp->PutString ("PRIOR_KAPREFIMG", kaprefimg_name ? kaprefimg_name : "NONE");
}

// ==========================================================================

Tikhonov1::Tikhonov1 (double _tau, const RVector *_x0,
    const Raster *_raster, const RVector *_kap)
    : Regularisation (_raster, _tau, _x0)
{
    // Create Laplacian
    Laplacian = new RCompRowMatrix (slen, slen);
    if(_kap==0) {
      cout << "Tikhonov1 constructor default kref\n";
      kappa = new RVector(slen);
      //      *kappa = 1.0; // doesn't seem to work..
      for(int k = 0; k < slen; k++)
	(*kappa)[k] = 1;
    }
    else
      kappa = _kap;
    cout << "Tikhonov1 constructor : calling CreateHessian\n";
    CreateHessian (raster, *kappa, *Laplacian);
    //CreateHessStruct1param (*Laplacian); // MS 080307
    pdim      = Laplacian->nCols();
    nset      = x0->Dim() / pdim;
    xASSERT (nset*pdim == x0->Dim(), "Invalid parameter dimensions");
    if(_kap==0)   {// kappa is not needed if Hessian is explicitly created.
      delete kappa;
      kappa = 0;
    }
}

Tikhonov1::~Tikhonov1 ()
{
    delete Laplacian;
}

void Tikhonov1::CreateHessian (const Raster *raster,
    const RVector &kappa, RCompRowMatrix &Hess)
{
    // compute sparse regularisation matrix as
    // grad kappa grad
    // where input image kappa is defined in BDim basis
    // For Laplacian, kappa=-1 (constant) 

    int slen = raster->SLen();
    idxtype *rowptr, *colidx;
	int nzero;
    int row, col, i, i0, i1, cnt_idx, cntj, j;
    raster->NeighbourGraph (rowptr, colidx, nzero);
    double valsum, *valptr = new double[nzero];
    //IVector cpix(dim), npix(dim);

    for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];
	for (i = i0; i < i1; i++) { // find central pixel
	    col = colidx[i];
	    if (col == row) {
		cnt_idx = i;
		cntj = raster->Sol2Basis (col);
		break;
	    }
	}
	valsum = 0.0;
	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) {
		j = raster->Sol2Basis (col);
		valptr[i] = -0.5 * (kappa[cntj] + kappa[j]);
		valsum += valptr[i];
	    }
	}
	valptr[cnt_idx] = -valsum;
    }
    Hess.Initialise (rowptr, colidx, valptr);
    delete []rowptr;
    delete []colidx;
    delete []valptr;
}

double Tikhonov1::GetValue (const RVector &x) const
{
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector dx (x - *x0), Lx(x.Dim());
    for (int i = 0; i < nset; i++) { // loop over parameter sets
	RVector pdx (dx, i*pdim, pdim);
	RVector pLx (Lx, i*pdim, pdim);
	pLx = *Laplacian * pdx;
    }
    return (dx & Lx) * tau;
}

RVector Tikhonov1::GetGradient (const RVector &x) const
{
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector dx (x - *x0), grad(x.Dim());
    for (int i = 0; i < nset; i++) { // loop over parameter sets
	RVector pdx (dx, i*pdim, pdim);
	RVector pgrad (grad, i*pdim, pdim);
	pgrad = *Laplacian * pdx;
    }
    return grad * tau;
}
RVector Tikhonov1::GetKappa    (const RVector &x) const
{
  RVector kx(x.Dim());
  kx = 1;
  return kx;
}
RVector Tikhonov1::GetHess1f(const RVector &x,const RVector &f) const
{
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
    dASSERT (f.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector  grad(x.Dim());
    for (int i = 0; i < nset; i++) { // loop over parameter sets
	RVector pf (f, i*pdim, pdim);
	RVector pgrad (grad, i*pdim, pdim);
	pgrad = *Laplacian * pf;
    }
    return grad * tau;
}

int Tikhonov1::GetHessianRow (const RVector &x, int i, idxtype *colidx,
    double *val) const
{
    int j, nz = Laplacian->SparseRow (i%pdim, colidx, val);
    for (j = 0; j < nz; j++) val[j] *= tau;

    // offset coefficient set
    int ofs = (i/pdim) * pdim;
    for (j = 0; j < nz; j++) colidx[j] += ofs;

    return nz;
}

void Tikhonov1::SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p)
{
    int i, j, nz;
    idxtype *colidx = new idxtype[pdim];
    double *val = new double[pdim];

    //RVector dx (x - *x0);
    //RVector pdx (dx, p*pdim, pdim);
    RVector h(pdim);

    Hess1.New(pdim,pdim);
    for (i = 0; i < pdim; i++) {
	nz = GetHessianRow (x, i/*+p*pdim*/, colidx, val);
	h.Clear();
	for (j = 0; j < nz; j++) h[colidx[j]] = val[j];
	Hess1.SetRow (i, h);
    }
}

RVector Tikhonov1::GetHessianDiag (const RVector &x) const
{
    RVector diag(nset*pdim);
    RVector Ldiag = Laplacian->Diag();
    for (int i = 0; i < nset; i++)
	for (int j = 0; j < pdim; j++)
	    diag[i*pdim+j] = Ldiag[j] * tau;
    return diag;
}

// ===================== TV with presmoothing ==============================

TVSigmaOLD::TVSigmaOLD (double _tau, double _beta, double _sd,
    const RVector *_x0, const Raster *_raster, bool _SmoothKappaOnly,
    const RVector *_kap)
    : Regularisation (_raster, _tau, _x0)
{
    //LOGOUT("Entering TVSigmaOLD constructor");
    //cout << "Entering TVSigmaOLD constructor\n" ;
    beta      = _beta;
    sd        = _sd;

    SmoothKappaOnly = _SmoothKappaOnly;
    betasq    = beta*beta;
    if(SmoothKappaOnly) {
      LOGOUT("Smoothing Kappa Only");
      cout << "Smoothing Kappa Only\n";
    }
    else {
      LOGOUT("Smoothing Kappa and x");
      cout << "Smoothing Kappa and x\n";
    }
    // Create Hessian Structure
    TVhess = new RCompRowMatrix (slen, slen);


    if(_kap==0) {
      cout << "TV constructor default kref\n";
      delete_kappa = true;
      kappa = new RVector(slen);
      //      *kappa = 1.0; // doesn't seem to work..
      for(int k = 0; k < slen; k++)
	(*kappa)[k] = 1;
    }
    else {
      delete_kappa = false;
      kappa = _kap;
    }
    //    CreateHessian (raster, *kappa, *TVhess); // values are irrelevant

    pdim      = TVhess->nCols();
    nset      = x0->Dim() / pdim;

    xASSERT (nset*pdim == x0->Dim(), "Invalid parameter dimensions");
    //    nset = 1;

    //LOGOUT("Leaving TVSigmaOLD constructor");
    //cout << "Leaving TVSigmaOLD constructor\n";
}

TVSigmaOLD::~TVSigmaOLD ()
{
    //LOGOUT("Entering TVSigmaOLD detructor");
    delete TVhess;
    if(delete_kappa)
      delete kappa;
    //LOGOUT("Leaving TVSigmaOLD constructor");
}
void TVSigmaOLD::OLDSetHessian (const RVector &x) const
  // in TVSigmaOLD, Hessian is not constant : depends on x
  // need to think about this. It will need TWO hessian matrices !
{
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector dx (x - *x0);


    int row, col, i, i0, i1;
 
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	if(sd > 0.0)
	  pdx = raster->SmoothImage(pdx,sd);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) { 
		    double hx = sqrt( SQR(pdx[row] - pdx[col]) + betasq);
		    double kx = 0.5 * ((*kappa)[row] + (*kappa)[col]);
		    (*TVhess)(col,row) = -2*betasq*kx/(CUBE(hx)) ;
		    (*TVhess)(col,col) += 2*betasq*kx/(CUBE(hx)) ;
		}
		//	    else  // this might not be necessary...
		//  	TVhess->operator()(col,col) += 2/sqrt(betasq) ;
	    } // end of neighbours of this pixel
	} // end of pixels	
    } // end of parameter

}


double TVSigmaOLD::GetValue (const RVector &x) const
{
    //LOGOUT("Entering TVSigmaOLD::GetValue");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");

    RVector dx (x - *x0);
    double result=0.0;

    int row, col, i, i0, i1;

    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	if(sd > 0.0 && !SmoothKappaOnly)
	  pdx = raster->SmoothImage(pdx,sd);
	//LOGOUT("parameter loop : p = %d",p);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) { 
		    double kx = 0.5 * ((*kappa)[row] + (*kappa)[col]);
		    result += kx*sqrt( SQR(pdx[row] - pdx[col]) + betasq) ;
		}
	    } // end of neighbours of this pixel
	} // end of pixels	
	//LOGOUT_2PRM("row %d result %f",row,result);
    } // end of parameter
    

    //LOGOUT("Leaving TVSigmaOLD::GetValue");
    return tau*result;
}

RVector TVSigmaOLD::GetGradient (const RVector &x) const
{
    //LOGOUT("Entering TVSigmaOLD::GetGradient");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
    
    RVector dx (x - *x0), grad(x.Dim());

    int row, col, i, i0, i1;
    
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	RVector spdx(pdx);
	if(sd > 0.0)
	  spdx = raster->SmoothImage(pdx,sd);

	RVector pgd (grad, p*pdim, pdim);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) {
		  double kx = 0.5 * ((*kappa)[row] + (*kappa)[col]);
		  if(SmoothKappaOnly)
		    pgd[row] += 2*kx*(pdx[row] - pdx[col]) /
			        sqrt( SQR(spdx[row] - spdx[col]) + betasq) ;
		  else
		    pgd[row] += 2*kx*(spdx[row] - spdx[col]) /
 		                sqrt( SQR(spdx[row] - spdx[col]) + betasq) ;
		}
	    } // end of neighbours of this pixel
	} // end of pixels	
    } // end of parameter

    //LOGOUT("Leaving TVSigmaOLD::GetGradient");
    return grad * tau;
}
RVector TVSigmaOLD::GetKappa    (const RVector &x) const
{
    //LOGOUT("Entering TVSigmaOLD::GetKappa");
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
    
    RVector dx (x - *x0), ktot(x.Dim());

    int row, col, i, i0, i1;
    
    for (int p = 0; p < nset; p++) { // loop over parameter sets
	RVector pdx (dx, p*pdim, pdim);
	RVector spdx(pdx);
	if(sd > 0.0)
	  spdx = raster->SmoothImage(pdx,sd);

	RVector pkap (ktot, p*pdim, pdim);
	for (row = 0; row < slen; row++) {
	    i0 = rowptr[row]; i1 = rowptr[row+1];
	    int ncount = 0;
	    for (i = i0; i < i1; i++) { // find neighbours
		col = colidx[i];
		if (col != row) {
		  ktot[row] += 2/sqrt( SQR(spdx[row] - spdx[col]) + betasq);
		  ncount ++;
		}
	    } // end of neighbours of this pixel
	    dASSERT (ncount >0, "Found a pixel with no neighbours");  
	    pkap[row] *= (*kappa)[row]/ncount; // average of edges.
	} // end of pixels	
    } // end of parameter

    //LOGOUT("Leaving TVSigmaOLD::GetKappa");
    return ktot;
}
RVector TVSigmaOLD::GetHess1f(const RVector &x,const RVector &f) const
{
    dASSERT (x.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
    dASSERT (f.Dim() == x0->Dim(), "Invalid coefficient vector dimension");
  LOGOUT("Hess1f is not implemented for TV");
  cerr << "Hess1f is not implemented for TV\n";
  return RVector (x.Dim());
}


int TVSigmaOLD::GetHessianRow (const RVector &x, int row, idxtype *colidx,
    double *val) const
  // DEFINITELY need to check with Martin !
{
    //  LOGOUT("Entering TVSigmaOLD::GetHessianRow");
    int i,j, col, cnt_idx, cntj, nz, row0, ofs;

    RVector sx(x);

    for(int p = 0; p < nset; p++) { // produced smoothed version of image
      RVector spx(sx, p*pdim, pdim);
      if(sd>0.0)
	spx = raster->SmoothImage(spx,sd); 
    }

    row0 = row%pdim;         // reduce row index
    ofs = (row/pdim) * pdim; // parameter offset
    nz = TVhess->SparseRow (row0, colidx, val); // structure of row
    
    for (i = 0; i < nz; i++) { // find central pixel
	col = colidx[i];
	if (col == row0) {
	    val[i] = 0.0;// reset value of diagonal pixel
	    cnt_idx = i;
	    cntj = col;
	    //cntj = raster->Sol2Grid (col);
	    break;
	}
    } // found central pixel

    for (i = 0; i < nz; i++) {
	col = colidx[i];
	if (col != row0) {
	    j = col;
	    //j = raster->Sol2Grid (col);
	    double hx = sqrt( SQR(sx[cntj+ofs] - sx[j+ofs]) + betasq);
	    double kx = 0.5 * ((*kappa)[cntj+ofs] + (*kappa)[j+ofs]);
	    val[i] = -2*kx*tau*betasq/(CUBE(hx)) ;
	    val[cnt_idx] -= val[i] ;
	}
    }
    /*  output values;
        for (i = 0; i < nz ; i++) cout << val[i] << " ";
        cout << endl;
    */
    // offset coefficient set
    
    for (i = 0; i < nz; i++) colidx[j] += ofs;
    //  LOGOUT("Leaving TVSigmaOLD::GetHessianRow");
    return nz;
}

RVector TVSigmaOLD::GetHessianDiag (const RVector &x) const
{
    //LOGOUT("Entering TVSigmaOLD::GetHessianDiag");
  // DEFINITELY need to check with Martin !
    RVector dx (x - *x0);


    int row, col, i, i0, i1, cntj, j;
    RVector diag(nset*pdim); 
    //    diag =1e-6 ;//set default value so no zero divides

    for (int p = 0; p < nset; p++){// loop over parameter sets
      RVector pdx (dx, p*pdim, pdim);
      if(sd> 0.0)
	pdx = raster->SmoothImage(pdx,sd);
      for (row = 0; row < slen; row++) {
	i0 = rowptr[row]; i1 = rowptr[row+1];
	cntj = row;
	/* no need for this
	for (i = i0; i < i1; i++) { // find central pixel
	    col = colidx[i];
	    if (col == row) {
		cnt_idx = i;
		cntj = col;
		//cntj = raster->Sol2Grid (col);
		break;
	    }
	} // found central pixel
	*/
	for (i = i0; i < i1; i++) { // find neighbours
	    col = colidx[i];
	    if (col != row) { 
		j = col;
		//j = raster->Sol2Grid (col);
		double hx = sqrt( SQR(pdx[cntj] - pdx[j]) + betasq);
	        double kx = 0.5 * ((*kappa)[cntj] + (*kappa)[j]);
       		diag[p*pdim+row] += 2*betasq*kx/(CUBE(hx)) ;
	    }
	    //	    else  // this might not be necessary...
	    //         	diag[p*pdim+cntj] += 2*tau/sqrt(betasq) ;
	} // end of neighbours of this pixel
      } // end of pixels
      //LOGOUT_3PRM("row %d pixel %d result %f",row-1, cntj,diag[p*pdim+row-1])
    } // end of set

    //LOGOUT("Leaving TVSigmaOLD::GetHessianDiag");
    return tau*diag;
}

#endif // !OLD_REGUL2
