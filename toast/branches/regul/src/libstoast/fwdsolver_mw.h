// -*-C++-*-
// ==========================================================================
// Forward model: finite element method (multi-wavelength version)
// ==========================================================================

#ifndef __FWDSOLVER_MW_H
#define __FWDSOLVER_MW_H

#include "fwdsolver.h"
#ifdef TOAST_MPI
#include "toast_mpi.h"
#endif

// =========================================================================

/**
 * \brief Templated forward solver class for multi-wavelength problems
 *
 * %TFwdSolverMW is an extension to TFwdSolver which allows to calculate
 * the boundary projections for multiple wavelengths with a single call
 * to \ref ProjectAll_wavel.
 */

template<class T> class TFwdSolverMW: public TFwdSolver<T> {
public:
    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param mesh pointer to associated QMMesh object
     * \param linsolver linear solver type (see \ref lsolver)
     * \param tol linear solver tolerance (only used for iterative methods)
     */
    TFwdSolverMW (const QMMesh *mesh, LSOLVER linsolver, double tol = 1e-10);

    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param mesh pointer to associated QMMesh object
     * \param solver linear solver type, provided as a string (see notes)
     * \param tol linear solver tolerance (only used for iterative methods)
     * \note For valid strings to define the linear solver, see
     * \ref SetLinSolver.
     */
    TFwdSolverMW (const QMMesh *mesh, char *solver, double tol = 1e-10);

    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param mesh pointer to associated QMMesh object
     * \param pp parser to read solver options from.
     * \note For recognised items in the configuration file used by the parser,
         see the \ref ReadParams method.
    */
    TFwdSolverMW (const QMMesh *mesh, ParamParser &pp);

    /**
     * \brief Destructor. Destroys the forward solver instance.
     */
    ~TFwdSolverMW ();

    /**
     * \brief Return boundary data for all sources and all detectors at all
     *   wavelengths.
     * \param qvec array of source column vectors
     * \param mvec array of measurement column vectors
     * \param sol multi-wavelength solution instance
     * \param omega modulation frequency (cycles/ps)
     * \param scl output data scaling: linear or logarithmic
     * \return data vector (size nw*nQM, with nw=number of wavelengths)
     * \note This method sequentially resets the system matrix for all
     *   wavelengths according to the optical parameters stored in sol,
     *   then calculates the fields and projections.
     * \note The returned vector is arranged in blocks of data for each
     *   wavelength. Each block is structured in the same way as the vector
     *   returned by \ref TFwdSolver::ProjectAll.
     */
    TVector<T> ProjectAll_wavel (const TCompRowMatrix<T> &qvec,
        const TCompRowMatrix<T> &mvec, const MWsolution &sol, double omega,
        DataScale scl = DATA_DEFAULT);

    /**
     * \brief Return boundary data for the complex case in a real vector
     * \param qvec array of source column vectors
     * \param mvec array of measurement column vectors
     * \param sol multi-wavelength solution instance
     * \param omega modulation frequency (cycles/ps)
     * \param scl output data scaling: linear or logarithmic
     * \return data vector (size nw*nQM*2, with nw=number of wavelengths)
     * \note This method sequentially resets the system matrix for all
     *   wavelengths according to the optical parameters stored in sol,
     *   then calculates the fields and projections.
     * \note The returned vector consists of a real and imaginary block.
     *   Both blocks consist of sub-blocks for each wavelength, and each
     *   sub-block contains further blocks for each source, consisting of
     *   the measurements for that source.
     */
    RVector ProjectAll_wavel_real (const TCompRowMatrix<T> &qvec,
        const TCompRowMatrix<T> &mvec, const MWsolution &sol, double omega,
        DataScale scl = DATA_DEFAULT);

protected:
    /**
     * \brief Forward solver initialisation routines
     */
    void Setup();

    /**
     * \brief Clean-up routines
     */
    void Cleanup();

private:
#ifdef TOAST_MPI
    MPI_Datatype mpitp;
    int sze, rnk;
    int *projall_count;
    int *projall_ofs;
    int *qidx;
#endif
};

// ==========================================================================
// template typedefs

typedef TFwdSolverMW<double>         RFwdSolverMW;
typedef TFwdSolverMW<std::complex<double> > CFwdSolverMW;

// ==========================================================================
// extern declarations of FwdSolverMW (only required for VS)

#ifndef __FWDSOLVER_MW_CC
extern template class STOASTLIB TFwdSolverMW<double>;
extern template class STOASTLIB TFwdSolverMW<std::complex<double> >;
#endif // !__FWDSOLVER_MW_CC

#endif // __FWDSOLVER_MW_H
