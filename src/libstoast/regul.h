// -*-C++-*-
#ifndef __ST_REGUL_H
#define __ST_REGUL_H

#include "stoastlib.h"

class Raster;
class ParamParser;

#ifndef SQR
#define SQR(_X) ((_X)*(_X))
#endif
#ifndef CUBE
#define CUBE(_X) ((_X)*(_X)*(_X))
#endif
typedef double (* PSIREG)(const double, const int np, double *p); //pointer to generic function

typedef enum {
    PRIOR_NONE,
    PRIOR_LAPLACIAN,
    PRIOR_DIAG,
    PRIOR_TK1,
    PRIOR_TV,
    PRIOR_HUBER,
    PRIOR_PM,
    PRIOR_QPM,
    PRIOR_TK1SIGMA,
    PRIOR_TVSIGMA,
    PRIOR_HUBERSIGMA,
    PRIOR_PMSIGMA,
    PRIOR_QPMSIGMA,
    PRIOR_TUKEYSIGMA,
    PRIOR_GGMRF,
    PRIOR_GENERIC,
    PRIOR_GENERICSIGMA,
    PRIOR_GENERICSCALESPACE,
    PRIOR_MRF
} PRIOR;


// ======================   CLASSES    =====================================
class STOASTLIB Regularisation {
public:
    Regularisation (const Raster *_raster = 0, double _tau = 1,
        const RVector *_x0 = 0);

    virtual ~Regularisation ();

    static Regularisation *Create (ParamParser *pp, const RVector *_x0,
        const Raster *_raster, const RVector *_xs = 0);
    // create a regularisation instance from a parameter file

    virtual void ReadParams (ParamParser *pp);
    virtual void WriteParams (ParamParser *pp);

    virtual const char *GetName() const = 0;
    // return a name to identify the regulariser

    inline double GetTau() const
    { return tau; }

    inline int GetNParam() const
    { return nset; }

    virtual double GetValue (const RVector &x) const = 0;
    // value of prior psi for given coefficient vector x

    virtual RVector GetGradient (const RVector &x) const = 0;
    // gradient dpsi/dx_i for given coefficient vector x

    virtual RVector GetKappa (const RVector &x) const = 0;
    // diffusivity function for given coefficient vector x

    virtual void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p) =0;
    // set the first order Hessian for input vector x
    virtual void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) =0;
    virtual RVector GetHess1f (const RVector &x, const RVector &f) const = 0;
    // Effect of 1st order Hessian L on vector f at given coefficient vector x

    virtual void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p) =0;
    // set the total Hessian H for input vector x

    virtual RVector GetFullHessf (const RVector &x, const RVector &f) const = 0;
    // Effect of total Hessian H on vector f at given coefficient vector x

  //    virtual void SetHessian (RCompRowMatrix &Hess, const RVector &x) =0;

    virtual int GetHessianRow (const RVector &x, int i, idxtype *colidx,
        double *val) const = 0;
    // returns single row i of the 2nd derivative matrix
    // d^2 psi / dx_i dx_j in sparse vector format. Return value is
    // number of nonzeros. column index list is returned in colidx,
    // matrix values are returned in val. Both colidx and val must
    // be allocated by the caller with sufficient size

    virtual RVector GetHessianDiag (const RVector &x) const = 0;
    // returns the diagonal of the 2nd derivative matrix

    virtual bool SetLocalScaling (const RVector &scale_all)
    { return false; }
    // set local regularisation scaling factor. scale_all is a concatenation
    // of local scaling factors for all parameters in 'G' basis
    // returns true if regularisation supports local scaling, false otherwise

    virtual PRIOR Type() const = 0; 

protected:
    void CreateHessStruct1param (RCompRowMatrix &Hess) ;
    // set up the Hessian structure for a single parameter.
    // multiple Hessians can be setup by concatenation

    int nset, pdim;
  //    void CreateNonLinearHessian (RCompRowMatrix &Hess, const RVector &x);
    double tau;  // regularisation hyperparameter
    PRIOR ptype;
    const Raster *raster;
    const RVector *x0;
    int slen, glen, dim;
    idxtype *rowptr, *colidx;
	int nzero;
    RCompRowMatrix Hs1p; // will hold sparse matrix structure for 1 Hessian
};
#ifdef REG_1ST_ORDER
// =========================================================================
class Regularisation_1stOrder : public Regularisation {
  // base class for Regularisation using gradient information
    Regularisation_1stOrder (const Raster *_raster = 0,const RVector *_x0 = 0);
    virtual ~Regularisation_1stOrder () {};

protected:
    void CreateHessStruct1param (RCompRowMatrix &Hess);
    // set up the Hessian structure for a single parameter.
    // multiple Hessians can be setup by concatenation

};

// =========================================================================

class STOASTLIB GenericSigma : public Regularisation_1stOrder {
#else
class STOASTLIB GenericSigma : public Regularisation {
#endif
public:
  //    GenericSigma (const Raster *_raster, const RVector *_x0, const int _npreg, double * _preg, double _tau, double _sdx = 0, double _sdf = 0, const RVector *_kref=0);
    GenericSigma (const Raster *_raster, const RVector *_x0, const int _npreg,
        double * _preg, double _tau = 1e-2, double _sdx = 0, double _sdf = 0,
        const void *_kref=0, bool _KRefTensFlag = false);

    GenericSigma (const Raster *_raster, const RVector *_x0, int _npreg,
        double *_preg, double _tau, const RVector &_refimg, double _sdr,
	double _sdx = 0, double _sdf = 0, double _fT = 0,
        bool _KRefTensFlag = false);

    virtual ~GenericSigma ();
    
    const char *GetName() const
    { return "GenericSigma"; }

    void ReadParams (ParamParser *pp);
    void WriteParams (ParamParser *pp);

    RVector *SetKref (const RVector &gkref_all);
    // Set the regularisation scaling for all parameters directly
    // gkref_all: concatenated scaling images for all parameters in
    // 'G' basis

    RVector *MakeKref (const RVector &gimgref, double sdr, double fT);
    // Return scalar field for isotropic diffusivity distribution given
    // image gimgref. sdr: value for reference sigma; fT: fraction of max to
    // use for Perona-Malik threshold [0-1]

    RDenseMatrix *MakeKrefTens (const RVector &gimgref, double sdr, double fT);
    // Return tensor field for anisotropic diffusivity distribution given
    // image gimgref. sdr: value for reference sigma; fT: fraction of max to
    // use for Perona-Malik threshold [0-1]

    PRIOR Type() const {return PRIOR_GENERICSIGMA; } 

  //    void SetHessian (RCompRowMatrix &Hess, const RVector &x) {
  //	CreateNonLinearHessian(raster, x, Hess);
  //    };

    double GetValue (const RVector &x) const;
    RVector GetGradient (const RVector &x) const {
	// MS: modified 100201
	// the second parameter of GetHess1f doesn't get x0 subtracted
	// inside GetHess1f, which can lead to a nonzero gradient even
	// if x = x0

	//return GetHess1f(x,x);
	return GetHess1f (x, x - *x0);
    }
    RVector GetKappa (const RVector &x) const;
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p);
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) ;
    RVector GetHess1f (const RVector &x, const RVector &f) const;
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p);
    RVector GetFullHessf (const RVector &x, const RVector &f) const;
  //    void SetHessian (RCompRowMatrix &Hess, const RVector &x) ;

    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)const;

    RVector GetHessianDiag (const RVector &x) const ;

    inline RDenseMatrix tinterp (const RDenseMatrix *mat, int i, int j) const;

    bool SetLocalScaling (const RVector &scale_all);

protected:
    double sdr, fT;
    const double sdx, sdf;
    virtual double func(const double, const int np, double *p) const = 0;
    virtual double kfunc(const double, const int np, double *p)const = 0;
    virtual double d2func(const double, const int np, double *p)const = 0;
    void CreateDerivativeOperators ();
    const int npreg;
    char *krefimg_name;
    const RVector *kref;
    const RDenseMatrix *kreftens;
    const bool KRefTensFlag;
    bool KRefIsLocal;
    double *preg;
    ICompRowMatrix NG; // used only if tensor flag set
    RCompRowMatrix *Db, *Df;  //will hold derivative oeprators
};

// =====================   Some Utility functions ==========================
inline double sinterp(const RVector & v, const int i, const int j)
{
    return 0.5*(v[i] + v[j]);
}
#ifdef LINTERP
inline RDenseMatrix GenericSigma::tinterp (const RDenseMatrix *mat,
    int i, int j) const
{
    RDenseMatrix tmp(mat[i]);
    for(int r = 0; r < tmp.nRows(); r++)
	for(int c = 0; c < tmp.nCols(); c++)
	    tmp(r,c) = 0.5*(mat[i](r,c) + mat[j](r,c));
    return tmp;
}
#else
#ifdef MINDET
//return matrix with minimum determinant!
inline RDenseMatrix GenericSigma::tinterp (const RDenseMatrix *mat,
    int i, int j) const
{
    double d1 = det(mat[i]);
    double d2 = det(mat[j]);
    return (d1 < d2 ? mat[i] : mat[j]);
}
#else
inline RDenseMatrix GenericSigma::tinterp (const RDenseMatrix *mat,
    int i, int j) const
// calculate tensor from interpolated gradients
{
#ifdef TINTERP2D // old version
    static const int dim = 2;       // only 2D for now.

    double fx = 0.5*( mat[i](0,0) +  mat[j](0,0));
    double fy = 0.5*( mat[i](1,1) +  mat[j](1,1));
    double fxsq = SQR(fx);
    double fysq = SQR(fy);
    double gsq = fxsq + fysq;
    double T  = 0.5*( mat[i](0,1) +  mat[j](0,1));
    double Tsq = SQR(T);
    RDenseMatrix tmp(dim,dim);
    tmp.Identity(dim);
    //  double gam = 1- exp( -gsq/Tsq);
    double gam = 1- exp( -sqrt(gsq)/T);
    tmp(0,0) -= (fx == 0.0 ? 0 : gam*fxsq/gsq); // xx term
    tmp(1,1) -= (fy == 0.0 ? 0 : gam*fysq/gsq); // yy term
    tmp(0,1) -= (gsq == 0.0 ? 0 : gam*fx*fy/gsq);   // xy term
    tmp(1,0) = tmp(0,1); // xy term
    return tmp;
#else // new, but untested
    int m, n;
    double gsq = 0.0;
    RVector f(dim);
    for (m = 0; m < dim; m++) {
	f[m] = 0.5*(mat[i](m,m) + mat[j](m,m));
	gsq += SQR(f[m]);
    }
    double T = 0.5*(mat[i](0,1) + mat[j](0,1));
    RDenseMatrix tmp(dim,dim);
    tmp.Identity (dim);
    double gam = 1.0-exp(-sqrt(gsq)/T);
    for (m = 0; m < dim; m++)
	for (n = 0; n < dim; n++)
	    tmp(m,n) -= (gsq == 0.0 ? 0 : gam*f[m]*f[n]/gsq);
    return tmp;
#endif
}
#endif
#endif

// ========================"Classic" TV========================================
// Scale Space literature refers to this function as "Charbonnier"
//
class STOASTLIB TVSigma: public GenericSigma {
public:
    TVSigma (double _tau, double _beta, double _sd, const RVector *_x0,
        const Raster *_raster, bool _SmoothKappaOnly = true,
	const void *_kref = 0, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1,  (new double [1]), _tau, _sd,
	  (_SmoothKappaOnly ? 0 : _sd), _kref, _KRefTensFlag)
    {preg[0] = _beta;}
    
    TVSigma (double _tau, double _beta, double _sd, const RVector *_x0,
	const Raster *_raster, const RVector &_refimg, double _sdr,
	double _fT, bool _SmoothKappaOnly = true, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1,  (new double [1]), _tau, _refimg,
	  _sdr, _sd, (_SmoothKappaOnly ? 0 : _sd), _fT, _KRefTensFlag)
    {preg[0] = _beta;};

    ~TVSigma ()
    {}

    const char *GetName() const
    { return "TV_SIGMA"; }

    void ReadParams (ParamParser *pp);
    void WriteParams (ParamParser *pp);

protected:
    inline double func(const double t, const int np, double *p) const{
	double betasq = SQR(p[0]);
	return p[0]*sqrt(SQR(t) + betasq) - betasq ;
    }

    inline double kfunc(const double t, const int np, double *p) const{
	double betasq = SQR(p[0]);
	return p[0]/sqrt(SQR(t) + betasq) ;
    }

    inline double d2func(const double t, const int np, double *p) const{
	double betasq = SQR(p[0]);
	//  return CUBE(p[0]/hx);
	return p[0]/sqrt(SQR(t) + betasq);
	// this is "Gauss-Newton" Hessian : convex
    }
};
    
class STOASTLIB TV: public TVSigma {
public:
    TV (double _tau, double _beta, const RVector *_x0, const Raster *_raster,
	const void *_kref=0, bool _KRefTensFlag = false)
	: TVSigma ( _tau, _beta, 0, _x0, _raster,false, _kref, _KRefTensFlag)
    {}

    TV (double _tau, double _beta, const RVector *_x0, const Raster *_raster,
	const RVector &_refimg, double _sdr, double _fT,
	bool _KRefTensFlag = false)
	: TVSigma (_tau, _beta, 0, _x0, _raster, _refimg, _sdr, _fT, false,
	  _KRefTensFlag)
    {}

    ~TV ()
    {}

    const char *GetName() const
    { return "TV"; }
};

// ======================== Laplacian = 1st order Tikhonov ====================
class STOASTLIB TK1Sigma: public GenericSigma {
public:
    TK1Sigma (double _tau, double _sd, const RVector *_x0,
        const Raster *_raster, bool _SmoothKappaOnly = true,
	const void *_kref=0, bool _KRefTensFlag =false)
	: GenericSigma (_raster, _x0, 0, (double*)0, _tau, _sd,
			(_SmoothKappaOnly ? 0 : _sd), _kref, _KRefTensFlag)
    {}

    TK1Sigma (double _tau, double _sd, const RVector *_x0,
        const Raster *_raster, const RVector &_refimg, double _sdr, double _fT,
	bool _SmoothKappaOnly = true, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 0, (double*)0, _tau, _refimg,
	  _sdr, _sd, (_SmoothKappaOnly ? 0 : _sd), _fT, _KRefTensFlag)
    {}

    ~TK1Sigma ()
    {}

    const char *GetName() const
    { return "TK1_SIGMA"; }

protected:
    inline double func(const double t, const int np, double *p) const{
	return 0.5*SQR(t);
    }

    inline double kfunc(const double t, const int np, double *p) const{
	return 1;
    }

    inline double d2func(const double t, const int np, double *p) const{
	return 1;
    }
};

class STOASTLIB TK1: public TK1Sigma {
public:
    TK1 (double _tau, const RVector *_x0, const Raster *_raster,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: TK1Sigma ( _tau, 0, _x0, _raster, false, _kref, _KRefTensFlag)
    {}

    TK1 (double _tau, const RVector *_x0, const Raster *_raster,
	 const RVector &_refimg, double _sdr, double _fT,
	 bool _KRefTenseFlag = false)
	: TK1Sigma (_tau, 0, _x0, _raster, _refimg, _sdr, _fT, false,
	  _KRefTenseFlag)
    {}

    ~TK1 ()
    {}

    const char *GetName() const
    { return "TK1"; }
};

// =========================== Huber Function ===========================
class STOASTLIB HuberSigma: public GenericSigma {
public:
    HuberSigma (double _tau, double _eps, double _sd, const RVector *_x0,
        const Raster *_raster, bool _SmoothKappaOnly = true,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1,  (new double [1]), _tau, _sd,
	  (_SmoothKappaOnly ? 0 : _sd), _kref, _KRefTensFlag)
    {preg[0] = _eps;}

    HuberSigma (double _tau, double _eps, double _sd, const RVector *_x0,
        const Raster *_raster, const RVector _refimg, double _sdr, double _fT,
        bool _SmoothKappaOnly = true, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1, (new double[1]), _tau, _refimg,
	  _sdr, _sd, (_SmoothKappaOnly ? 0 : _sd), _fT, _KRefTensFlag)
    {preg[0] = _eps;}

    ~HuberSigma ()
    {}

    const char *GetName() const
    { return "HUBER_SIGMA"; }

    void ReadParams (ParamParser *pp);
    void WriteParams (ParamParser *pp);

protected:
    inline double func(const double t, const int np, double *p) const{
	double eps = p[0];
	double at = fabs(t);
	return (at <= eps ? 0.5*SQR(t) : eps*(at - 0.5*eps));
    }

    inline double kfunc(const double t, const int np, double *p) const{
	double eps = p[0];
	double at = fabs(t);
	return (at <= eps ? 1 : eps/at);
    }

    inline double d2func(const double t, const int np, double *p) const{
	double eps = p[0];
	double at = fabs(t);
	return (at <= eps ? 1 : CUBE(eps/at)); // not standard!
    }
};

class STOASTLIB Huber: public HuberSigma {
public:
    Huber (double _tau, double _eps, const RVector *_x0, const Raster *_raster,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: HuberSigma (_tau, _eps, 0, _x0, _raster, false, _kref, _KRefTensFlag)
    {}

    Huber (double _tau, double _eps, const RVector *_x0, const Raster *_raster,
        const RVector &_refimg, double _sdr, double _fT,
        bool _KRefTensFlag = false)
	: HuberSigma (_tau, _eps, 0, _x0, _raster, _refimg, _sdr, _fT, false,
	  _KRefTensFlag)
    {}

    ~Huber ()
    {}

    const char *GetName() const
    { return "HUBER"; }
};
// =========================================================================
class STOASTLIB PMSigma: public GenericSigma {
public:
    PMSigma (double _tau, double _T, double _sd, const RVector *_x0,
        const Raster *_raster, bool _SmoothKappaOnly = true,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1,  (new double [1]), _tau, _sd,
	  (_SmoothKappaOnly ? 0 : _sd), _kref, _KRefTensFlag)
    {preg[0] = _T;}

    PMSigma (double _tau, double _T, double _sd, const RVector *_x0,
        const Raster *_raster, const RVector &_refimg, double _sdr, double _fT,
        bool _SmoothKappaOnly = true, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1, (new double[1]), _tau, _refimg,
	  _sdr, _sd, (_SmoothKappaOnly ? 0 : _sd), _fT, _KRefTensFlag)
    {preg[0] = _T;}

    ~PMSigma ()
    {}

protected:
    inline double func(const double t, const int np, double *p) const{
	double T = p[0];
	return 0.5*SQR(T)*(1-exp(-SQR(t/T)));
    }

    inline double kfunc(const double t, const int np, double *p) const{
	double T = p[0];
	return exp(-SQR(t/T));
    }

    inline double d2func(const double t, const int np, double *p) const{
	double T = p[0];
	//  return (1- 2*SQR(t/T))*exp(-SQR(t/T));
	return exp(-SQR(t/T)); // this is "Gauss-Newton" Hessian : convex
    }
};

class STOASTLIB PM: public PMSigma {
public:
    PM (double _tau, double _T, const RVector *_x0, const Raster *_raster,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: PMSigma (_tau, _T, 0, _x0, _raster, false, _kref, _KRefTensFlag)
    {}

    PM (double _tau, double _T, const RVector *_x0, const Raster *_raster,
	const RVector &_refimg, double _sdr, double _fT,
	bool _KRefTensFlag = false)
	: PMSigma (_tau, _T, 0, _x0, _raster, _refimg, _sdr, _fT, false,
	  _KRefTensFlag)
    {}

    ~PM ()
    {}
};
// =========================================================================
class STOASTLIB QPMSigma: public GenericSigma {
public:
    QPMSigma (double _tau, double _T, double _sd, const RVector *_x0,
        const Raster *_raster, bool _SmoothKappaOnly = true,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1,  (new double[1]), _tau, _sd,
	  (_SmoothKappaOnly ? 0 : _sd), _kref, _KRefTensFlag)
    {preg[0] = _T;}

    QPMSigma (double _tau, double _T, double _sd, const RVector *_x0,
        const Raster *_raster, const RVector &_refimg, double _sdr, double _fT,
        bool _SmoothKappaOnly = true, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1, (new double[1]), _tau, _refimg,
	  _sdr, _sd, (_SmoothKappaOnly ? 0 : _sd), _fT, _KRefTensFlag)
    {preg[0] = _T;}

    ~QPMSigma ()
    {}

protected:
    inline double func(const double t, const int np, double *p) const{
	double  T = p[0];
	return 0.5*SQR(T)*log(SQR(t/T) + 1);
    }

    inline double kfunc(const double t, const int np, double *p) const{
	double T = p[0];
	return 1/(SQR(t/T) + 1) ;
    }

    inline double d2func(const double t, const int np, double *p) const{
	double T = p[0];
	double hx = (SQR(t/T) + 1);
	//  return (1-SQR(t/T))/SQR(hx); // this looks a bit dodgy : not sym.pos.def
	return hx;   // this is "Gauss-Newton" Hessian : convex
    }
};

class STOASTLIB QPM: public QPMSigma {
public:
    QPM (double _tau, double _T, const RVector *_x0, const Raster *_raster,
	const void *_kref = 0, bool _KRefTensFlag = false)
	: QPMSigma (_tau, _T, 0, _x0, _raster, false, _kref, _KRefTensFlag)
    {}

    QPM (double _tau, double _T, const RVector *_x0, const Raster *_raster,
	const RVector &_refimg, double _sdr, double _fT,
        bool _KRefTensFlag = false)
	: QPMSigma (_tau, _T, 0, _x0, _raster, _refimg, _sdr, _fT, false,
          _KRefTensFlag)
    {}

    ~QPM ()
    {}
};

// ======================   Tukey BiWeight ===============================
class STOASTLIB TukeySigma: public GenericSigma {
public:
    TukeySigma (double _tau, double _T, double _sd, const RVector *_x0,
        const Raster *_raster, bool _SmoothKappaOnly = true,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1, (new double[1]), _tau, _sd,
	  (_SmoothKappaOnly ? 0 : _sd), _kref, _KRefTensFlag)
    {preg[0] = _T;}

    TukeySigma (double _tau, double _T, double _sd, const RVector *_x0,
        const Raster *_raster, const RVector &_refimg, double _sdr, double _fT,
        bool _SmoothKappaOnly = true, bool _KRefTensFlag = false)
	: GenericSigma (_raster, _x0, 1, (new double[1]), _tau, _refimg,
	  _sdr, _sd, (_SmoothKappaOnly ? 0 : _sd), _fT, _KRefTensFlag)
    {preg[0] = _T;}

    ~TukeySigma ()
    {}

protected:
    inline double func(const double t, const int np, double *p) const{
	double  T = p[0];
	return 0.5*SQR(t)*(1 - SQR(t/T) + SQR(SQR(t/T))/3.0 );
    }

    inline double kfunc(const double t, const int np, double *p) const{
	double T = p[0];
	double at = fabs(t);
	return (at <= T ? SQR(1 - (SQR(t/T))) : 0);
    }

    inline double d2func(const double t, const int np, double *p) const{
	double T = p[0];
	double at = fabs(t);
	return (at <= T ? SQR(1 - (SQR(t/T))) : 0);
    }
};

class STOASTLIB Tukey: public TukeySigma {
public:
    Tukey (double _tau, double _T, const RVector *_x0, const Raster *_raster,
        const void *_kref = 0, bool _KRefTensFlag = false)
	: TukeySigma (_tau, _T, 0, _x0, _raster, false, _kref, _KRefTensFlag)
    {}

    Tukey (double _tau, double _T, const RVector *_x0, const Raster *_raster,
        const RVector &_refimg, double _sdr, double _fT,
        bool _KRefTensFlag = false)
	: TukeySigma (_tau, _T, 0, _x0, _raster, _refimg, _sdr, _fT, false,
	  _KRefTensFlag)
    {}

    ~Tukey ()
    {}
};

// =========================================================================

class STOASTLIB GenericScaleSpace : public Regularisation {
public:
    GenericScaleSpace (const Raster *_raster, const RVector *_x0, PSIREG _func, PSIREG _kfunc, PSIREG _k1func, const int _npreg, double * _preg, double _tau, double _sdx = 0, double _sdf = 0, const RDenseMatrix *_kreftens=0);

    ~GenericScaleSpace ();

    PRIOR Type() const {return PRIOR_GENERICSCALESPACE; } 

  //    void SetHessian (RCompRowMatrix &Hess, const RVector &x) {
  //	CreateNonLinearHessian(raster, x, Hess);
  //    };

    double GetValue (const RVector &x) const;
    RVector GetGradient (const RVector &x) const {
	// MS: modified 100201
	// the second parameter of GetHess1f doesn't get x0 subtracted
	// inside GetHess1f, which can lead to a nonzero gradient even
	// if x = x0

	//return GetHess1f(x,x);
	return GetHess1f (x, x - *x0);
	//return GetHess1f(x,x);
    }
    RVector GetKappa (const RVector &x) const;
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p);
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) ;
    RVector GetHess1f (const RVector &x, const RVector &f) const {
    if(KappaSet)
      return GetHess1f_KappaSet(f);
    else
      return GetHess1f_KappaNotSet(x,f);
    }
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p);
    RVector GetFullHessf (const RVector &x, const RVector &f) const;
  //    void SetHessian (RCompRowMatrix &Hess, const RVector &x) ;

    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val) const;

    RVector GetHessianDiag (const RVector &x) const ;

protected:
    RVector GetHess1f_KappaSet (const RVector &f) const;
    RVector GetHess1f_KappaNotSet (const RVector &x, const RVector &f) const;

    const double sdx, sdf ;
    const PSIREG func, kfunc, k1func;
    const int npreg;
    const RDenseMatrix *kreftens;
    RDenseMatrix **KappaTens;
    bool  KappaSet;
    double *preg;
};
// =====================Pretty Much Redundant==================================

class STOASTLIB Generic : public GenericSigma {// zero scale version of GenericSigma

public:
  Generic (const Raster *_raster, const RVector *_x0,  const int _npreg, double * _preg,  double _tau, const RVector *_kap=0) : GenericSigma (_raster, _x0, _npreg,  _preg, _tau, 0, 0,_kap) {};
};

#define OLD_REGUL2
#ifdef OLD_REGUL2
// =========================================================================

class STOASTLIB NullRegularisation: public Regularisation {
public:
    NullRegularisation (): Regularisation() {}
    PRIOR Type() const { return PRIOR_NONE; } 
    const char *GetName() const { return "NONE"; }
    double GetValue     (const RVector &x) const;
    RVector GetGradient (const RVector &x) const;
    RVector GetKappa    (const RVector &x) const;
    RVector GetHess1f   (const RVector &x, const RVector &f) const;

    void SetHessian (RCompRowMatrix &Hess, const RVector &x) {}
    void SetHessianFromKappa(RCompRowMatrix &Hess, const RVector &kappa) {}
    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)
        const;
    RVector GetHessianDiag (const RVector &x) const;

  /* undefined functions */
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p) {}
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) {}
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p) {}
    RVector GetFullHessf (const RVector &x, const RVector &f) const { return RVector(); }

};
// =========================================================================

class STOASTLIB Tikhonov0: public Regularisation {
public:
    Tikhonov0 (double _tau, const RVector *_x0, const RVector *_xs);
    ~Tikhonov0();
    PRIOR Type() const {return PRIOR_DIAG; } 

    const char *GetName() const { return "TK0"; }
    double GetValue (const RVector &x) const;
    RVector GetGradient (const RVector &x) const;
    RVector GetKappa    (const RVector &x) const;
    RVector GetHess1f   (const RVector &x, const RVector &f) const;
    void SetHessian (RCompRowMatrix &Hess, const RVector &x);
    void SetHessianFromKappa(RCompRowMatrix &Hess, const RVector &kappa);
    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)
        const;
    RVector GetHessianDiag (const RVector &x) const;
    void SetTau (double _tau);
    void SetTau (const RVector &_tau);

  /* undefined functions */
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p) {}
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) {}
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p) {}
    RVector GetFullHessf (const RVector &x, const RVector &f) const { return RVector(); }

    void ReadParams (ParamParser *pp);
    void WriteParams (ParamParser *pp);

private:
    RVector tauvec;
    bool use_tauvec;
    const  RVector *xs;
    char *kaprefimg_name;
};

// =========================================================================

class STOASTLIB Tikhonov1: public Regularisation {
public:
    Tikhonov1 (double _tau, const RVector *_x0, const Raster *raster, const RVector *_kap=0);
    ~Tikhonov1 ();
    PRIOR Type() const { return PRIOR_LAPLACIAN; } 
    const char *GetName() const { return "Tikhonov1"; }
    double GetValue (const RVector &x) const;
    RVector GetGradient (const RVector &x) const;
    RVector GetKappa    (const RVector &x) const;
    RVector GetHess1f   (const RVector &x, const RVector &f) const;
    void SetHessian (RCompRowMatrix &Hess, const RVector &x) {
      Hess= *Laplacian;
    };
    void SetHessianFromKappa(RCompRowMatrix &Hess, const RVector &kappa) {
      //      CreateHessian(raster,kappa,Hess);
    };
    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)
        const;
    RVector GetHessianDiag (const RVector &x) const;
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p);

  /* undefined functions */
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) {}
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p) {}
    RVector GetFullHessf (const RVector &x, const RVector &f) const { return RVector(); }

private:
    void CreateHessian (const Raster *raster, const RVector &kappa,
        RCompRowMatrix &Hess);

    const RVector *kappa;
    RCompRowMatrix *Laplacian;
  //    int nset, pdim;
};

// =========================================================================

class STOASTLIB TVSigmaOLD: public Regularisation {
public:
    TVSigmaOLD (double _tau, double _beta, double _sd, const RVector *_x0, const Raster *_raster, bool _SmoothKappaOnly = true, const RVector *_kap=0);
    ~TVSigmaOLD ();

    PRIOR Type() const {return PRIOR_TVSIGMA; } 
    double GetValue (const RVector &x) const;
    RVector GetGradient (const RVector &x) const;
    RVector GetKappa    (const RVector &x) const;
    RVector GetHess1f   (const RVector &x, const RVector &f) const;
    void SetHessian (RCompRowMatrix &Hess, const RVector &x) {
      //      CreateNonLinearHessian(raster, x, Hess) ;
    };
    void SetHessianFromKappa(RCompRowMatrix &Hess, const RVector &_kap) {
      //      CreateHessian(raster,_kap,Hess);
    };
    void OLDSetHessian (const RVector &x) const;
    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)
        const;
    RVector GetHessianDiag (const RVector &x) const;
  /* undefined functions */
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p) {}
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap) {}
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p) {}
    RVector GetFullHessf (const RVector &x, const RVector &f) const { return RVector(); }

protected:
    double tau, beta, betasq, sd;
    bool SmoothKappaOnly ;
    const RVector *kappa;
    bool  delete_kappa;   // this isn't at all pretty...
    RCompRowMatrix *TVhess;
  //    int nset, pdim;
};

// =========================================================================

class STOASTLIB TVOLD: public TVSigmaOLD { // zero scale version of TVSigma
public:
  TVOLD (double _tau, double _beta, const RVector *_x0, const Raster *_raster,  const RVector *_kap=0): TVSigmaOLD(_tau,_beta,0,_x0, _raster, false, _kap) {};
    PRIOR Type() const {return PRIOR_TV; } 

};

// =========================================================================

class STOASTLIB MRF: public Regularisation {
public:
    MRF (double _tau, const RVector *_x0, const Raster *raster, const RVector *_kap = 0);
    ~MRF ();
    PRIOR Type() const { return PRIOR_MRF; }
    const char *GetName() const { return "MRF"; }

    void SetKaprefImg (const RVector *kap = 0);

    double GetValue (const RVector &x) const;
    RVector GetGradient (const RVector &x) const;

    // not implemented (yet)
    RVector GetKappa (const RVector &x)  const;
    void SetHess1 (RCompRowMatrix &Hess1, const RVector &x, const int p);
    void SetHess1FromKappa (RCompRowMatrix &Hess, const RVector &kap);
    RVector GetHess1f (const RVector &x, const RVector &f) const;
    void SetFullHess (RCompRowMatrix &Hess, const RVector &x, const int p);
    RVector GetFullHessf (const RVector &x, const RVector &f) const;
    int GetHessianRow (const RVector &x, int i, idxtype *colidx, double *val)
	const;
    RVector GetHessianDiag (const RVector &x) const;

private:
    RCompRowMatrix *MRFa, *MRFb;
};

#endif // !OLD_REGUL2
#endif // !__ST_REGUL_H
