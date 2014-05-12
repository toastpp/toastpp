#include "PN_incl.h"
#include "mathlib.h"
/** factorial
**/
double factorial(int n);

/** double factorial
**/
double doublefactorial(int n);

/** Computes C^{n} where C is a 'complex' number and 'n' is a positive integer
**/
toast::complex cpowi(toast::complex &c, const int m);

/** Computes (-1)^{m}
**/
int sign(int m);

/** Returns 0, 1, and 2 respectively based on whether m=0, m<0 and m>0.
**/
int signum(int m);

/* evaluates \int_{S^{2}} Y_{\ell_{1}, m_{1}} Y_{\ell_{2}, m_{2}} where both Y's are complex spherical harmonics
*/
toast::complex kronD(const IVector &a, const IVector &b);

/** Computes associated Legendre polynomials of a given order on a set of points 
*	l -> Maximum order
*	numpts -> number of points
*	pt -> the three dimensional point set
*	LT[] -> Legendre polynomials in a table form
**/
void LegendreTable(const int l, const int numpts, const RDenseMatrix &pt, RDenseMatrix* &LT);

/**Real spherical harmonics computed on a set of points
* order -> Maximum order of spherical harmonics
* numpts ->  number of points on which to evaluate
* pt -> 3D point set
* sphHarm[l] -> spherical harmonics over point set of order 'l'
**/
void sphericalHarmonics(const int order, const int numpts, const RDenseMatrix& pt, RDenseMatrix* &Ylm);

/** Adds aB to its appropriate place in the system matrix
* spatrow -> spatial row where 'a' is drawn from
* spatcol -> spatial column where 'a' is drawn from
* node_angN -> number of angular degrees of freedom for all the spatial nodes
* offset -> starting location in the system matrix for each spatial node
* a_ij -> 'a'
* B -> B
* C -> output (System matrix) 
**/
void kronplus(const int spatrow, const int spatcol, const IVector& node_angN, const IVector& offset, const double Aval, const RCompRowMatrix &B, RCompRowMatrix& C);

/*Computes the integral on the sphere of the form 
	 * \int_{S^{n-1}} (s.n)_{plusminus} \psi_{i} \psi_{j}
	 * Inputs:
	 * 	size1 -> number of rows
	 *	size2 -> number of columns
	 *	sphOrder1 -> order of spherical harmonics along rows
         * 	sphOrder2 -> order of spherical harmonics along columns
	 *      pts -> quadrature points
	 *      wts -> quadrature weights
	 *      bnormal -> outward pointing normal to the boundary
	 *      Ylm -> precomputed table of spherical harmonics over the quadrature points
	 * Output:
	 * 	bintplus -> boundary integral on the outward facing half sphere
	 *      bintminus -> boundary integral on the inward facing half sphere
*/
void BIntUnitSphere(const int, const int, const int, const int, const RDenseMatrix&, const RVector&, const RVector&, RDenseMatrix* &, RCompRowMatrix&, RCompRowMatrix&);

/*Computes the integral on the sphere of the form 
	 * \int_{S^{n-1}} (s.n)_{plus} \psi_{i} \psi_{j}, \int_{S^{n-1}} (s.n)_{minus} \psi_{i} \psi_{j} and \int_{S^{n-1}} (s.n)_{minus} R(s) \psi_{i}(Hs) \psi_{j} (s) for reflection boundary condition
	 * Inputs:
	 *      ref_in -> refractive index inside the domain
	 *      ref_out -> refractive index of the medium outside the domain
	 * 	size1 -> number of rows
	 *	size2 -> number of columns
	 *	sphOrder1 -> order of spherical harmonics along rows
         * 	sphOrder2 -> order of spherical harmonics along columns
	 *      pts -> quadrature points
	 *      wts -> quadrature weights
	 *      bnormal -> outward pointing normal to the boundary
	 *      Ylm -> precomputed table of spherical harmonics over the quadrature points
	 * Output:
	 * 	bintplus -> boundary integral on the outward facing half sphere
	 *      bintminus -> boundary integral on the inward facing half sphere
*/
void BRIntUnitSphere(const double, const double, const int, const int, const int, const int, const RDenseMatrix&, const RVector&, const RVector&, RDenseMatrix* &, RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&);


/** Thresholds and shrinks a real dense matrix to give a real sparse matrix
* NOTE!! The threshold is set to 1e-15
**/
RCompRowMatrix shrink(const RDenseMatrix &dnsmat);

/**Computes the maximum spherical harmonic order used in a given element
**/
void findMaxLocalSphOrder(const Mesh &mesh, const IVector& sphOrder, const IVector& node_angN, const int el, int &maxSphOrder, int &maxAngN);

/* Sin(\theta)Cos(\phi)Y_{l, m}^{R} = Sin(\theta)Cos(\phi)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m+1}, Y_{l+1, m+1}, Y_{l-1, m-1}, Y_{l+1, m-1}, Y_{l-1, -m+1}, 
Y_{l+1, -m+1}, Y_{l-1, -m-1} and Y_{l+1, -m-1}.
*/
void sincosY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c);

/* Sin(\theta)Sin(\phi)Y_{l, m}^{R} = Sin(\theta)Sin(\phi)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m+1}, Y_{l+1, m+1}, Y_{l-1, m-1}, Y_{l+1, m-1}, Y_{l-1, -m+1}, 
Y_{l+1, -m+1}, Y_{l-1, -m-1} and Y_{l+1, -m-1}. 
*/
void sinsinY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c);

/* Cos(\theta)Y_{l, m}^{R} = Cos(\theta)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m}, Y_{l+1, m}, Y_{l-1, -m}, Y_{l+1, -m}. 
*/
void cosY(const int l, const int m, CVector& e, CVector& f, IDenseMatrix& e1c, IDenseMatrix& f1c);

/* Y_{l, m}^{R} = (a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l, m}, Y_{l, -m}. 
*/
void sphY(const int l, const int m, CVector& p, IDenseMatrix& p1c);

/*
Computes integral of the form
\int (a_{0}Y_{l1, m1} + a_{1}Y_{l2, m2})(b_{0}Y_{l3, m3} + b_{1}Y_{l4, m4})
*/
double Integrate2(CVector &a, CVector &b, IDenseMatrix &a1c, IDenseMatrix &b1c);


double intSinCosY(const int l, const int m);

double intSinSinY(const int l, const int m);

double intCosY(const int l, const int m);


/**
Computes all the angular integrals required by variable order PN approximation
**/
void genmat_angint(const int sphOrder, const int angN, RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Aintscsc, RCompRowMatrix& Aintscss, RCompRowMatrix& Aintscc, RCompRowMatrix& Aintssss, RCompRowMatrix& Aintssc, RCompRowMatrix& Aintcc);

/**
Phase function discretization
NOTE!! Only Henyey-Greenstein phase function has been implemented currently
**/
RVector phaseFuncDisc(int sphOrder, toast::complex (*phaseFunc)(const double g, const double costheta), const double g);

/** Phase integrals 
**/
void genmat_apu(toast::complex (*phaseFunc)(const double g, const double costheta), const double g, const int angN, const int sphOrder, RCompRowMatrix& apu1, RCompRowMatrix& apu1sc, RCompRowMatrix& apu1ss, RCompRowMatrix& apu1c);

