#ifndef __SOURCE_H
#define __SOURCE_H

#include "mathlib.h"
#include "felib.h"
#include "sphericalHarmonic_algebra.h"
#include "pparse.h"

typedef enum {
    SRCMODE_NEUMANN,
    SRCMODE_ISOTROPIC
} SourceMode;

void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth, SourceMode &srctp);

RVector QVec_Gaussian (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode);

RVector QVec_Cosine (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode);

RVector QVec_Point (const Mesh &mesh, const Point &cnt, SourceMode mode);

/* Computes the source vectors for all the boundary sources when a QM file has been specified
*/
void genmat_toastsource(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh, const RCompRowMatrix qvec, const int ns, RVector* &dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm);

/** Computes source vectors for point sources on the boundary 
**/
void genmat_source(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh,  const IVector& Nsource, const int ns, RVector* &dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm);

/** Computes source vectors for point sources in the interior 
**/
void genmat_intsource(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh,  const IVector& Nsource, const int ns, RVector* &dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm, const RVector &delta);

/** When used in conjunction with genmat_source function with point sources, this function simulates an uncollided source 
**/
void genmat_intsourceuncollided(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* &Source, Mesh& mesh, const IVector& Nsource, const int ns, RVector* &dirVec, const RDenseMatrix& pts, const RVector& wts, const RVector &mua, const RVector &mus, const RVector &delta, const double g);


//
#endif // !__SOURCE_H
