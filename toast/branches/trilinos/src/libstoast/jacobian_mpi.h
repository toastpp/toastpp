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

#ifndef __JACOBIANMPI_H
#define __JACOBIANMPI_H

#include "mathlib.h"
#include "felib.h"
#include "toasttype.h"

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
 * \param Jmod Jacobian (log amplitude block)
 * \param Jarg Jacobian (phase block)
 * \note This version generates the boundary projection values on the fly
 *   from mvec, if required (i.e. if using log data).
 */
STOASTLIB void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CCompRowMatrix &mvec, const CVector *dphi, const CVector *aphi,
    DataScale dscale, RDenseMatrixMPI &Jmod, RDenseMatrixMPI &Jarg);

/**
 * \brief Generate Jacobian matrix of forward operator (complex case)
 * \param raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param mesh FEM forward mesh
 * \param mvec measurement vectors (rows of sparse matrix)
 * \param dphi array of forward fields
 * \param aphi array of adjoint fields
 * \param dscale data scaling flag (lin/log)
 * \param Jmod Jacobian (log amplitude block)
 * \param Jarg Jacobian (phase block)
 * \note This version passes the boundary projection values as input arguments.
 */
STOASTLIB void GenerateJacobian (const Raster *raster, const QMMesh *mesh,
    const CVector *dphi, const CVector *aphi, const CVector *proj,
    DataScale dscale, RDenseMatrixMPI &Jmod, RDenseMatrixMPI &Jarg);

/**
 * \brief Generate Jacobian matrix for CW intensity data (real case).
 * \param raster solution basis mapper, or NULL to generate the Jacobian
 *   directly on the forward mesh
 * \param mesh FEM forward mesh
 * \param mvec measurement vectors (rows of sparse matrix)
 * \param dphi array of forward fields
 * \param aphi array of adjoint fields
 * \param dscale data scaling flag (lin/log)
 * \param J Jacobian
 * \note This version generates the boundary projection values on the fly
 *   from mvec, if required (i.e. if using log data).
 */
STOASTLIB void GenerateJacobian_cw_mua (const Raster *raster,
    const QMMesh *mesh, const RCompRowMatrix &mvec,
    const RVector *dphi, const RVector *aphi,
    DataScale dscale, RDenseMatrixMPI &J);

/// @}

// ============================================================================
// PMDF calculations

// absorption PMDF
template<class T>
inline TVector<T> PMDF_mua (const TVector<T> &dphi, const TVector<T> &aphi)
{
    return -dphi*aphi;
}

// diffusion PMDF
template<class T>
inline TVector<T> PMDF_kap (const TVector<T> *Gdphi, const TVector<T> *Gaphi,
    int dim)
{
    TVector<T> pmdf(Gdphi[0].Dim());
    for (int j = 0; j < dim; j++)
	pmdf -= Gdphi[j]*Gaphi[j];
    return pmdf;
}

// convert to log data PMDF for measurement m
template<class T>
inline TVector<T> PMDF_log (const TVector<T> &pmdf, T &m)
{
    return pmdf/m;
}


template<class T>
void ImageGradient (const IVector &dim, const RVector &size,
    const TVector<T> &im, TVector<T> *grad, const int *mask = 0);

#endif // !__JACOBIANMPI_H
