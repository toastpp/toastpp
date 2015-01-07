#ifndef __DETECTOR_H
#define __DETECTOR_H

#include "mathlib.h"
#include "felib.h"
#include "pparse.h"

void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);

void genmat_toastdetector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Detector, const Mesh& mesh, const RCompRowMatrix mvec, const int nd, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);

void genmat_toastdetectorvalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Dvec, const Mesh& mesh, const RCompRowMatrix mvec, const int iq, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);

void genmat_detector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Detector, const Mesh& mesh,  const IVector& Ndetector, const int nd, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);

void genmat_detectorvalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Dvec, const Mesh& mesh, const int Ndetector, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);



#endif // !__DETECTOR_H
