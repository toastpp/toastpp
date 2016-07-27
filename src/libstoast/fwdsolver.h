// -*-C++-*-
// ==========================================================================
// Forward model: finite element method
// ==========================================================================

#ifndef __FWDSOLVER_H
#define __FWDSOLVER_H

//#include "supermatrix.h"
#include "toasttype.h"
#include "mathlib.h"

class Solution;
class MWsolution;

/// \defgroup lsolver Linear solver methods
//@{
enum LSOLVER{
    LSOLVER_UNDEFINED,      ///<   undefined
    LSOLVER_DIRECT,         ///<   direct solver (LU)
    LSOLVER_ITERATIVE       ///<   iterative solver
};
//@}

/// \defgroup dscale Scaling methods for projections
//@{
enum DataScale {
    DATA_DEFAULT,           ///<   default method
    DATA_LIN,               ///<   linear scaling
    DATA_LOG                ///<   logarithmic scaling
};
//@}

// =========================================================================
// Nonmember declarations

template<class T> class TFwdSolver;

template<class T>
STOASTLIB TVector<T> ProjectSingle (const QMMesh *mesh, int q,
    const TCompRowMatrix<T> &mvec, const TVector<T> &phi,
    DataScale dscale = DATA_LIN);

template<class T>
STOASTLIB TVector<T> ProjectAll (const QMMesh *mesh,
    const TCompRowMatrix<T> &mvec, const TVector<T> *phi,
    DataScale dscale = DATA_LIN);

// =========================================================================

/**
 * \brief Templated forward solver class
 *
 * %TFwdSolver is a template class which encapsulates the FEM diffusion forward
 * solver. It can be instantiated either as TFwdSolver<double> (or RFwdSolver
 * for short) to describe real-valued problems (continuous-wave measurements),
 * or as TFwdSolver<complex> (or CFwdSolver) for complex-valued problems
 * (frequency domain problems).
 *
 * The %TFwdSolver class has an associated linear solver (defined in the
 * constructor) for solving the linear FEM system. This can be either a
 * direct solver (Cholesky for TFwdSolver<double> and LU for
 * TFwdSolver<complex>), or an iterative method (CG, BiCG, BiCGSTAB, GMRES).
 * For iterative methods, a tolerance limit is also required.
 *
 * Each %TFwdSolver instance has an FEM mesh associated with it. The mesh is
 * defined by a call to \ref Allocate, which also constructs the sparse system
 * matrix structure based on the mesh connectivity graph.
 *
 * The system matrix is built with a call to \ref Reset for a given set of
 * nodal optical coefficients. Subsequently, the photon density field (given
 * a source vector) can be evaluated with \ref CalcField and \ref CalcFields,
 * and the measurement values on the boundary can be obtained with
 * \ref Project and \ref ProjectAll.
 */

template<class T> class TFwdSolver {
public:
    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param mesh pointer to associated QMMesh object
     * \param linsolver linear solver type (see \ref lsolver)
     * \param tol linear solver tolerance (only used for iterative methods)
     */
    TFwdSolver (const QMMesh *mesh, LSOLVER linsolver, double tol = 1e-10);

    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param mesh pointer to associated QMMesh object
     * \param solver linear solver type, provided as a string (see notes)
     * \param tol linear solver tolerance (only used for iterative methods)
     * \note For valid strings to define the linear solver, see
     * \ref SetLinSolver.
     */
    TFwdSolver (const QMMesh *mesh, const char *solver, double tol = 1e-10, int nth=1);

    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param mesh pointer to associated QMMesh object
     * \param pp parser to read solver options from.
     * \note For recognised items in the configuration file used by the parser,
         see the \ref ReadParams method.
    */
    TFwdSolver (const QMMesh *mesh, ParamParser &pp);

    /**
     * \brief Destructor. Destroys the forward solver instance.
     */
    ~TFwdSolver ();

    /**
     * \brief Set the default scaling for data returned by the projection
     *   methods.
     * \param scl scaling method (DATA_LIN or DATA_LOG)
     * \param The default scaling defined by this method is used by projection
     *   methods whenever they use a scaling parameter of DATA_DEFAULT.
     * \note Before the first call to SetDataScaling, the default setting
     *   is DATA_LIN (linear data scaling)
     * \sa GetDataScaling, ProjectAll, ProjectAll_real
     */
    void SetDataScaling (DataScale scl);

    /**
     * \brief Returns the current default setting for data scaling.
     * \return current data scaling method: DATA_LIN (linear scaling) or
     *   DATA_LOG (logarithmic scaling).
     * \sa SetDataScaling, ProjectAll, ProjectAll_real
     */
    DataScale GetDataScaling () const;

    /**
     * \brief Set the preconditioning method.
     * \param type Preconditioner type.
     * \note The default preconditioner is None.
     * \note If the preconditioner is constructed from a ParamParser object,
     *   the preconditioner may also be set via the LINSOLVER_PRECON tag in
     *   the parameter file.
     * \sa ReadParams
     */
    void SetPrecon (PreconType type);

    void SetPhaseUnwrap (bool unwrap)
    { unwrap_phase = unwrap; }

    /**
     * \brief Read solver parameters from a parameter file.
     * \param pp parser instance
     * \note This method recognises the following items in a parameter file:
     * <table col="3">
     * <tr><td>LINSOLVER</td><td><i>string</i></td><td>Linear solver method.
     *   The following choices are available: DIRECT (LU solver) CG, (conjugate
     *   gradient solver) BICG (biconjugate gradient solver), BICGSTAB
     *   (biconjugate gradient stabilised solver), GMRES (generalised minimum
     *   residual solver), GAUSSSEIDEL (Gauss-Seidel solver)</td></tr>
     * <tr><td>LINSOLVER_TOL</td><td><i>real</i></td>
     *   <td>Tolerance value for iterative solvers</td></tr>
     * </table>
     * \note The LINSOLVER_TOL item is not required if LINSOLVER = DIRECT.
     * \note Any required items not found in the file are queried interactively
     *   from the user.
     */
    void ReadParams (ParamParser &pp);

    /**
     * \brief Write the current solver settings to a parameter file.
     * \param pp parser instance
     * \note This method writes the LINSOLVER and (for iterative methods only)
     *    the LINSOLVER_TOL items to the parameter file. See \ref ReadParams
     *    for details on the items.
     */
    void WriteParams (ParamParser &pp);

    /**
     * \brief Set the solver type and tolerance.
     * \param solver linear solver type, provided as a string (see notes)
     * \param tol  linear solver tolerance (only used for iterative methods)
     * \note Valid options for the linear solver are:
     * <table col="2">
     * <tr><td>DIRECT</td><td>Direct solver (LU)</td></tr>
     * <tr><td>CG</td><td>Conjugate gradient solver</td></tr>
     * <tr><td>BICG</td><td>Bi-conjugate gradient solver</td></tr>
     * <tr><td>BICGSTAB</td>
     *   <td>Bi-conjugate gradient stabilised solver</td></tr>
     * <tr><td>GMRES</td><td>Generalised minimum residual solver</td></tr>
     * </table>
     */
    void SetLinSolver (const char *solver, double tol = 1e-10);

    /**
     * \brief Returns the current solver type.
     * \return Solver type enumeration value (see \ref lsolver).
     */
    inline LSOLVER LinSolver() const { return solvertp; }

    /**
     * \brief Returns the solver tolerance.
     * \return Current solver tolerance value.
     * \note The tolerance is only relevant for iterative solvers.
     */
    inline double GetLinSolverTol() const { return iterative_tol; }

    /**
     * \brief Set max iteration count for iterative solver
     * \param maxit max number of iterations (0 for auto)
     * \note If the max iteration count is set to 0, it will be
     *   set to the dimension of the linear problem.
     */
    inline void SetLinSolverMaxit (int maxit) { iterative_maxit = maxit; }

    /**
     * \brief Returns a pointer to the FEM mesh associated with the forward
     *   solver instance.
     * \return Mesh pointer
     */
    inline const QMMesh *MeshPtr() const { return meshptr; }

    /**
     * \brief Evaluates fill structure of system matrices and allocates
     *   dynamic memory for them.
     * \note This function evaluates the fill structure of the sparse
     *   system matrix from the element connectivity graph of the associated
     *   mesh and allocates storage for it.
     * \note If a direct solver was selected, memory for the factorised matrix
     *   is also allocated, based on a symbolic Cholesky factorisation of the
     *   system matrix.
     * \note For iterative solvers, a simple diagonal precoditioner is
     *   generated.
     */
    void Allocate ();

    /**
     * \brief Construct the FEM system matrix from a set of parameter
     *   distributions for a given modulation frequency.
     * \param sol Solution instance containing parameter vectors
     * \param omega modulation frequency (cycles/ps)
     * \note For template type \e double, parameter omega must be 0.
     * \note A call to \ref Allocate must have preceeded the first call to
     *   AssembleSystemMatrix.
     * \sa Reset
     */
    void AssembleSystemMatrix (const Solution &sol, double omega = 0,
        bool elbasis = false);

    void AssembleSystemMatrixComponent (const RVector &prm, int type);

    /**
     * \brief Construct the FEM mass matrix for time-dependent problems.
     * \param mesh mesh pointer
     * \note If the provided mesh pointer is NULL, then the mesh from the
     *   preceeding call to \ref Allocate is used instead. If no mesh has been
     *   assigned, the method fails.
     * \sa Allocate, AssembleSystemMatrix
     */
    void AssembleMassMatrix (const Mesh *mesh = 0);

    /**
     * \brief Reset the forward solver by re-building the system matrices
     *   from a new set of parameters.
     * \param sol Solution instance containing the parameter vectors
     * \param omega modulation frequency (cycles/ps)
	 * \param elbasis parameter vectors in sol are provided in piecewise
	 *   constant element basis rather than nodal basis
     * \note For template type \e double, parameter omega must be 0.
     * \note The system matrices must have been allocated (once) with a
     *   previous call to \ref Allocate.
     * \note This function rebuilds the system matrix with a call to
     *   \ref AssembleSystemMatrix. If a direct solver was selected, it also
     *   performs a Cholesky (double) or LU factorisation (complex).
     * \note If the mass matrix exists, it is rebuilt as well
     * \bug Resetting the mass matrix shouldn't be necessary, since it doesn't
     *   depend on the parameters.
     * \sa AssembleSystemMatrix, AssembleMassMatrix
     */
    void Reset (const Solution &sol, double omega = 0, bool elbasis = false);

    /**
     * \brief Calculate photon density field for a given source vector.
     * \param [in] qvec source vector (h-basis)
     * \param [out] phi photon density field (h-basis)
     * \param [out] res optional; if used, and if an iterative solver is
     *   used, the structure is filled with iteration count and relative
     *   error.
     * \sa CalcFields
     */
    void CalcField (const TVector<T> &qvec, TVector<T> &phi,
        IterativeSolverResult *res=0, int th=0) const;

    /**
     * \brief Calculate photon density fields for all sources.
     * \param [in] qvec array of column source vectors
     * \param [out] phi array of photon density fields
     * \param [out] res optional; if used, and if an iterative solver is
     *   used, the structure is filled with the maximum iteration count
     *   and the maximum relative error for any source
     * \note On call, phi must be an array of nq vectors, where nq is the
     *   number of rows in qvec.
     * \sa CalcField
     */
    void CalcFields (const TCompRowMatrix<T> &qvec, TVector<T> *phi,
        IterativeSolverResult *res = 0) const;

    /**
     * \brief Calculate photon density fields for a range of sources.
     * \param q0 min source index
     * \param q1 max source index + 1
     * \param qvec complete array of column source vectors
     * \param [out] phi array of photon density fields
     * \param [out] res optional; if used, and if an iterative solver is
     *   used, the structure is filled with the maximum iteration count
     *   and the maximum relative error for any source
     * \note This function calculates the photon density fields of a sub-range
     *   q of sources (q0 <= q < q1)
     * \note The returned array phi contains q1-q0 fields.
     */
    void CalcFields (int q0, int q1, const TCompRowMatrix<T> &qvec,
        TVector<T> *phi, IterativeSolverResult *res = 0) const;

    /**
     * \brief Return boundary data for a single source, given the
     *   corresponding photon density field.
     * \param q source index (>= 0)
     * \param mvec array of measurement column vectors
     * \param phi photon density field for source q (h-basis representation)
     * \return data vector (size nQMref[q])
     */
    TVector<T> ProjectSingle (int q, const TCompRowMatrix<T> &mvec,
        const TVector<T> &phi, DataScale scl = DATA_DEFAULT) const;

    /**
     * \brief Return boundary data for all sources and all
     *   detectors, given photon density fields for all sources.
     * \param mvec array of measurement column vectors
     * \param dphi array of photon density fields for all sources (h-basis
     *   representation)
     * \param scl output data scaling: linear or logarithmic
     * \return data vector (size nQM)
     */
    TVector<T> ProjectAll (const TCompRowMatrix<T> &mvec,
        const TVector<T> *dphi, DataScale scl = DATA_DEFAULT);

    /**
     * \brief Return boundary data for all sources and all detectors,
     *   given parameter distributions.
     * \param qvec array of source column vectors
     * \param mvec array of measurement column vectors
     * \param sol solution containing optical parameters in h-basis
     * \param omega modulation frequency (cycles/ps)
     * \param scl output data scaling (linear or logarithmic)
     */
    TVector<T> ProjectAll (const TCompRowMatrix<T> &qvec,
        const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
        DataScale scl = DATA_DEFAULT);

    /**
     * \brief Return boundary data for the complex case in a real vector
     * \param mvec array of measurement column vectors
     * \param dphi array of photon density fields for all sources (h-basis
     *   representation)
     * \param scl output data scaling: linear or logarithmic
     * \return data vector (size nQM*2)
     * \note This function splits the complex data vector into real and
     *   imaginary parts, and returns them concatenated in a real vector,
     *   to be compliant with the standard TOAST data vector format.
     */
    RVector ProjectAll_real (const TCompRowMatrix<T> &mvec,
        const TVector<T> *dphi, DataScale scl = DATA_DEFAULT);
    FVector ProjectAll_singlereal (const TCompRowMatrix<T> &mvec,
        const TVector<T> *dphi, DataScale scl = DATA_DEFAULT);

    /**
     * \brief Return boundary data for all sources and all detectors, given
     *   parameters defined via a Solution instance.
     * \param qvec array of source column vectors
     * \param mvec array of measurement column vectors
     * \param sol solution instance (h-basis representation)
     * \param omega modulation frequency (cycles/ps)
     * \param scl output data scaling: linear or logarithmic
     * \return data vector (size nQM or nQM*2)
     * \note This version always returns a real vector. For complex data, it
     *   returns the real and imaginary parts in sequential blocks of the
     *   vector.
     * \note This function constructs the system matrix from the solution,
     *   then calculates the fields for all sources, and obtains the
     *   projections by applying the measurement vectors.
     */
    RVector ProjectAll_real (const TCompRowMatrix<T> &qvec,
        const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
        DataScale scl = DATA_DEFAULT);
    FVector ProjectAll_singlereal (const TCompRowMatrix<T> &qvec,
        const TCompRowMatrix<T> &mvec, const Solution &sol, double omega,
        DataScale scl = DATA_DEFAULT);

    const QMMesh *meshptr;  ///< pointer to the associated FEM mesh
    LSOLVER solvertp;       ///< linear solver type
    IterativeMethod method; ///< iterative solver method, if applicable
#ifdef MPI_FWDSOLVER
    TCompRowMatrixMPI<T> *F;  ///< Distributed FEM system matrix
#else
    TCompRowMatrix<T> *F;   ///< FEM system matrix
#endif
    TCompRowMatrix<T> *FL;  ///< lower triangle of system matrix decomposition
    TVector<T> *Fd;         ///< diagonal of Cholesky factorisation
    TPreconditioner<T> *precon; ///< preconditioner instance
    TCompRowMatrix<T> *B;   ///< mass matrix; only used for time-domain problems
    //mutable SuperLU_data<T> lu_data; ///< parameters for LU solver
    void *SuperLU;          ///< SuperLU solver engine
    double iterative_tol;   ///< iterative solver tolerance
    int iterative_maxit;    ///< iterative solver max iterations (0 for auto)
    PreconType precontp;    ///< preconditioner

protected:
    void Setup (int nth=1);
    void SetupType (int nth=1);
    void DeleteType ();

    /**
     * \brief Unfold real and imaginary parts of a complex vector.
     * \param vec Input vector argument of template type
     * \return Unfolded real vector.
     * \note For complex template types, this creates a real output vector of
     *   twice the size of the input vector, containing the real and imaginary
     *   parts of the vector.
     * \note For real template types, this simply returns the input vector.
     */
    RVector UnfoldComplex (const TVector<T> &vec) const;
    FVector UnfoldSComplex (const TVector<T> &vec) const;

#ifdef UNDEF
    void UnfoldComplex (const TVector<T> &vec, RVector &res) const;
    void UnfoldComplex (const TVector<T> &vec, FVector &res) const;
#endif

    DataScale dscale;       ///< default data scaling: DATA_LIN or DATA_LOG
    TVector<T> *pphi;       ///< work buffer for field calculation
    bool unwrap_phase;      ///< use phase unwrapping?

    // ===============================================================
    // MPI-specific functions and data members

#ifdef MPI_FWDSOLVER

public:

    /**
     * \brief Set MPI parameters.
     * \note Sets up the MPI-specific parameters (stores MPI size and
     *   rank, defines the initial distribution of sources over processes)
     */
    void Setup_MPI();

    /**
     * \brief Deallocate MPI data structures.
     */
    void Cleanup_MPI();

    /**
     * \brief Define the distribution of nq sources over the available
     *   processes.
     * \param nq number of sources
     * \note An even distribution is assumed. Modify this to implement
     *   a load-balancing strategy.
     */
    void DistributeSources_MPI (int nq) const;

    /**
     * \brief Distributed field calculation (non-blocking)
     * \param [in] qvec array of source column vectors
     * \param [out] phi array of field vectors
     * \note No synchronisation is performed at the end of the computation.
     *   Each process only computes the field vectors for the sources it
     *   is responsible for. Other field vectors remain unchanged.
     *   The calling function is responsible to either continue the distributed
     *   computation or perform a synchronisation.
     * \note On call, phi must be an array of nq vectors, where nq is the
     *   number of rows in qvec.
     */
    void CalcFields_proc (const TCompRowMatrix<T> &qvec,
        TVector<T> *phi) const;

    /**
     * \brief Distributed boundary projection calculation (non-blocking)
     * \param [in] mvec array of measurement column vectors
     * \param [in] phi array of field vectors
     * \param [in] scl data scaling
     * \param [out] proj projection vector
     * \note Each process calculates the projections for the sources it is
     *   responsible for, and only updates the corresponding section of
     *   \e proj. No synchronisation is performed. Each process only
     *   requires those fields phi which correspond to processed sources, so
     *   this method can be combined with \ref CalcFields_proc.
     */
    void ProjectAll_proc (const TCompRowMatrix<T> &mvec,
        const TVector<T> *phi, DataScale scl,  TVector<T> &proj);

protected:

    MPI_Datatype mpitp;   ///< MPI type corresponding to template type
    int rnk;              ///< MPI process number (>= 0)
    int sze;              ///< number of MPI processes (>= 1)
    mutable int nQ;       ///< current number of sources for load balancing
    mutable int *Q0, *Q1; ///< range of sources (q0<=q<q1) for all processes

#endif // MPI_FWDSOLVER

#if THREAD_LEVEL==2
    int nthread;
#endif
};


// ==========================================================================
// Some macros

#ifdef MPI_FWDSOLVER
#define SETUP_MPI() Setup_MPI()
#define CLEANUP_MPI() Cleanup_MPI()
#else
#define SETUP_MPI()
#define CLEANUP_MPI()
#endif

// ==========================================================================
// template typedefs

typedef TFwdSolver<float>          FFwdSolver;
typedef TFwdSolver<double>         RFwdSolver;
typedef TFwdSolver<std::complex<float> > SCFwdSolver;
typedef TFwdSolver<std::complex<double> > CFwdSolver;

// ==========================================================================
// extern declarations of FwdSolver (only required for VS)

#ifndef __FWDSOLVER_CC
extern template class STOASTLIB TFwdSolver<float>;
extern template class STOASTLIB TFwdSolver<double>;
extern template class STOASTLIB TFwdSolver<std::complex<float> >;
extern template class STOASTLIB TFwdSolver<std::complex<double> >;
#endif // !__FWDSOLVER_CC

#endif // !__FWDSOLVER_H
