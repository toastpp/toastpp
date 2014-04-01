#include "mathlib.h"
#include "felib.h"
#include "PN_incl.h"
#define MIN(A,B) ( (A) < (B) ? (A) : (B))
#define MAX(A,B) ( (A) > (B) ? (A) : (B))

/**Computes the maximum spherical harmonic order used in a given element
**/
void findMaxLocalSphOrder(const Mesh &mesh, const IVector& sphOrder, const IVector& node_angN, const int el, int &maxSphOrder, int &maxAngN);

/** Adds aB to its appropriate place in the system matrix
* spatrow -> spatial row where 'a' is drawn from
* spatcol -> spatial column where 'a' is drawn from
* node_angN -> number of angular degrees of freedom for all the spatial nodes
* offset -> starting location in the system matrix for each spatial node
* a_ij -> 'a'
* B -> B
* C -> output (System matrix) 
**/
void kronplus(const int spatrow, const int spatcol, const IVector& node_angN, const IVector& offset, const double a_ij, const RCompRowMatrix &B, RCompRowMatrix& C);

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
void BIntUnitSphere(const int size1, const int size2, const int sphOrder1, const int sphOrder2,  const RDenseMatrix& pts, const RVector& wts, const RVector& bnormal, RDenseMatrix* &Ylm, RCompRowMatrix& bintplus, RCompRowMatrix& bintminus);

void BRIntUnitSphere(const double nin, const double nout, const int size1, const int size2, const int sphOrder1, const int sphOrder2,  const RDenseMatrix& pts, const RVector& wts, const RVector& bnormal, RDenseMatrix* &Ylm, RCompRowMatrix& bintplus, RCompRowMatrix& bintminus, RCompRowMatrix& brintminus);

/** Preallocating memory for boundary integral terms.
The matrix is considered to be spatially sparse and dense angularly(since the angular sparsity pattern is not known)
which is eventually shrunk after computation.
**/
void initialiseA2b1(const Mesh &mesh, const IVector &node_angN, const IVector &offset, RCompRowMatrix& A2, RCompRowMatrix& b1);

/**Compute the boundary integral terms using quadrature
**/
void genmat_boundint(const Mesh& mesh,  const RVector &ref, const double ref_out, const char srctp, const IVector& sphOrder, const IVector& node_angN, const IVector& offset, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm, RCompRowMatrix& A2, RCompRowMatrix& b1);


