// -*-C++-*-
// =========================================================================
// jacobian.h (libstoast)
// Interface for Jacobian calculation
//
// Jacobian of forward operator: derivative of measurables with respect to
// parameters (frequency domain case)
// The Jacobian is a real matrix is composed of 4 blocks:
//
//     |  dlnamp  |  dlnamp  |
//     |  ------  |  ------  |
//     |   dmua   |  dkappa  |
//     |----------+----------|
//     |  dphase  |  dphase  |
//     |  ------  |  ------  |
//     |   dmua   |  dkappa  |
//
// Each row corresponds to a measurement (source-detector pair), each column
// corresponds to a voxel in solution space.
// The returned Jacobian does not contain any data or parameter scalings.
// =========================================================================

#ifndef __JACOBIAN_H
#define __JACOBIAN_H

#include "mathlib.h"
#include "felib.h"
#include "toasttype.h"

// Flags for computing Jacobian on nodal/element basis rather than
// basis provided by a raster object
#define INT_NDBASIS 0
#define INT_ELBASIS -1
#define RASTER_NDBASIS ((Raster*)INT_NDBASIS)
#define RASTER_ELBASIS ((Raster*)INT_ELBASIS)

/**
 * \defgroup jacobian Jacobian matrix calculation
 *
 * This section contains functions that calculate the Jacobian matrix for
 * the frequency-domain diffusion equation, in terms of absorption and
 * diffusion coefficient, and either for real and imaginary parts of the
 * complex data, or for log amplitude and phase.
 */
/// @{
/**
 * \brief Generate Jacobian matrix of forward operator (complex case)
 * \param raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param mesh FEM forward mesh
 * \param mvec measurement vectors (rows of sparse matrix)
 * \param dphi array of forward fields
 * \param aphi array of adjoint fields
 * \param dscale data scaling flag (lin/log)
 * \param J Jacobian
 * \param elbasis (only used if raster==0) If true, returned Jacobian is
 *   generated for piecewise constant element basis rather than shape
 *   function basis
 * \note This version generates the boundary projection values on the fly
 *   from mvec, if required (i.e. if using log data).
 */
STOASTLIB void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrix &J);

/**
 * \brief Generate Jacobian matrix of forward operator (complex case)
 * \param raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param mesh FEM forward mesh
 * \param dphi array of forward fields
 * \param aphi array of adjoint fields
 * \param dscale data scaling flag (lin/log)
 * \param J Jacobian
 * \param elbasis (only used if raster==0) If true, returned Jacobian is
 *   generated for piecewise constant element basis rather than shape
 *   function basis
 * \note This version passes the boundary projection values as input arguments.
 */
STOASTLIB void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrix &J);

/**
 * \brief Generate raw Jacobian matrix for CW intensity data (real case)
 * \param [in] raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param [in] mesh FEM forward mesh
 * \param [in] dphi array of direct fields
 * \param [in] aphi array of adjoint fields
 * \param [out] Jmua pointer to Jacobian matrix (absorption component), or
 *   NULL if not required
 * \param [out] Jkap pointer to Jacobian matrix (diffusion component), or
 *   NULL if not required
 * \note This version performs no data or parameter scaling or mapping to
 *   logarithmic values. It calculates dy/dx, where y are boundary cw intensity
 *   data, and x is mua (absorption) and/or kappa (diffusion).
 * \note If Jmua == NULL or Jkap == NULL, then the corresponding component
 *   is not calculated.
 */
STOASTLIB void GenerateJacobian_cw (const Raster *raster, const QMMesh *mesh,
    const RVector *dphi, const RVector *aphi,
    RDenseMatrix *Jmua = NULL, RDenseMatrix *Jkap = NULL);

/**
 * \brief Generate Jacobian matrix for CW intensity data (real case).
 * \param [in] raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param [in] mesh FEM forward mesh
 * \param [in] mvec measurement vectors (rows of sparse matrix)
 * \param [in] dphi array of forward fields
 * \param [in] aphi array of adjoint fields
 * \param [in] dscale data scaling flag (lin/log)
 * \param [out] Jmua Jacobian matrix (absorption component)
 * \param [out] Jkap Jacobian matrix (diffusion component)
 * \note This version generates the boundary projection values on the fly
 *   from mvec, if required (i.e. if using log data).
 * \note If Jmua == NULL or Jkap == NULL, then the corresponding component
 *   is not calculated.
 * \note If dscale == DATA_LOG, then the Jacobian is calculated in terms of
 *   logarithmic amplitude (d log y)/dx.
 */
STOASTLIB void GenerateJacobian_cw (const Raster *raster,
    const QMMesh *mesh, const RCompRowMatrix &mvec,
    const RVector *dphi, const RVector *aphi,
    DataScale dscale, RDenseMatrix *Jmua = NULL, RDenseMatrix *Jkap = NULL);

/**
 * \brief Generate mua Jacobian matrix for CW intensity data (real case).
 * \param [in] raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param [in] mesh FEM forward mesh
 * \param [in] dphi array of forward fields
 * \param [in] aphi array of adjoint fields
 * \param [in] dscale data scaling flag (lin/log)
 * \param [out] Jmua Jacobian matrix (absorption component)
 * \param [out] Jkap Jacobian matrix (diffusion component)
 * \note This version passes the boundary projection values as input arguments.
 * \note If Jmua == NULL or Jkap == NULL, then the corresponding component
 *   is not calculated.
 * \note If dscale == DATA_LOG, then the Jacobian is calculated in terms of
 *   logarithmic amplitude (d log y)/dx.
 */
STOASTLIB void GenerateJacobian_cw (const Raster *raster,
    const QMMesh *mesh, const RVector *dphi, const RVector *aphi,
    const RVector *proj, DataScale dscale,
    RDenseMatrix *Jmua = NULL, RDenseMatrix *Jkap = NULL);

/// @}

template<class T>
STOASTLIB void ImageGradient (const IVector &dim, const RVector &size,
    const TVector<T> &im, TVector<T> *grad, const int *mask = 0);

#endif // !__JACOBIAN_H
