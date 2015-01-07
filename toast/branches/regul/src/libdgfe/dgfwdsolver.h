// -*-C++-*-
// ==========================================================================
// Forward model: finite element method
// ==========================================================================
#ifndef __DGFWDSOLVER_H
#define __DGFWDSOLVER_H
//#include "supermatrix.h"
//#include "zsp_defs.h"
//#include "toasttype.h"
#include "stoastlib.h"
#include "nonconformingMesh.h"
#include <iostream>
using namespace std;
// =========================================================================

/**
 * \brief Templated forward solver class
 *
 * %TDGFwdSolver is a template class which encapsulates the FEM diffusion forward
 * solver. It can be instantiated either as TDGFwdSolver<double> (or RDGFwdSolver
 * for short) to describe real-valued problems (continuous-wave measurements),
 * or as TDGFwdSolver<complex> (or CDGFwdSolver) for complex-valued problems
 * (frequency domain problems).
 *
 * The %TDGFwdSolver class has an associated linear solver (defined in the
 * constructor) for solving the linear FEM system. This can be either a
 * direct solver (Cholesky for TDGFwdSolver<double> and LU for
 * TDGFwdSolver<complex>), or an iterative method (CG, BiCG, BiCGSTAB, GMRES).
 * For iterative methods, a tolerance limit is also required.
 *
 * Each %TDGFwdSolver instance has an FEM mesh associated with it. The mesh is
 * defined by a call to \ref Allocate, which also constructs the sparse system
 * matrix structure based on the mesh connectivity graph.
 *
 * The system matrix is built with a call to \ref Reset for a given set of
 * nodal optical coefficients. Subsequently, the photon density field (given
 * a source vector) can be evaluated with \ref CalcField and \ref CalcFields,
 * and the measurement values on the boundary can be obtained with
 * \ref Project and \ref ProjectAll.
 */
template <class T>
class TDGFwdSolver{
public:
    TDGFwdSolver (LSOLVER linsolver, double tol = 1e-10);

    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param solver linear solver type, provided as a string (see notes)
     * \param tol linear solver tolerance (only used for iterative methods)
     * \note For valid strings to define the linear solver, see
     * \ref SetLinSolver.
     */
    TDGFwdSolver (char *solver, double tol = 1e-10);

    /**
     * \brief Constructor. Creates a forward solver instance.
     * \param pp parser to read solver options from.
     * \note For recognised items in the configuration file used by the parser,
         see the \ref ReadParams method.
    */
    TDGFwdSolver (ParamParser &pp);

    /**
     * \brief Destructor. Destroys the forward solver instance.
     */
    ~TDGFwdSolver ();
    
     inline NonconformingMesh *MeshPtr() const { return meshptr; }

    /**
     * \brief Evaluates fill structure of system matrices and allocates
     *   dynamic memory for them.
     * \param mesh FEM mesh to be associated with the forward solver
     * \note This function associates the specified mesh with the forward
     *   solver, evaluates the fill structure of the resulting sparse
     *   system matrix and allocates storage for it.
     * \note If a direct solver was selected, memory for the factorised matrix
     *   is also allocated, based on a symbolic Cholesky factorisation of the
     *   system matrix.
     * \note For iterative solvers, a simple diagonal precoditioner is
     *   generated.
     */
    void Allocate (NonconformingMesh &mesh);

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
    void AssembleSystemMatrix (const Solution &sol, double omega = 0, bool elbasis = false);
    void ReadParams(ParamParser &pp);
    void SetLinSolver(char *sol, double tol);
    void CalcField (const TVector<T> &qvec, TVector<T> &phi) const;



    /**
     * \brief Construct the FEM mass matrix for time-dependent problems.
     * \param mesh mesh pointer
     * \note If the provided mesh pointer is NULL, then the mesh from the
     *   preceeding call to \ref Allocate is used instead. If no mesh has been
     *   assigned, the method fails.
     * \sa Allocate, AssembleSystemMatrix
     */
   
    NonconformingMesh *meshptr;  ///< pointer to the associated FEM mesh
    LSOLVER solvertp;       ///< linear solver type
    IterativeMethod method; ///< iterative solver method, if applicable
    TCompRowMatrix<T> *F;   ///< FEM system matrix
    TCompRowMatrix<T> *FL;  ///< lower triangle of system matrix decomposition
    TVector<T> *Fd;         ///< diagonal of Cholesky factorisation
    TPreconditioner<T> *precon; ///< preconditioner instance
    RCompRowMatrix *B;      ///< mass matrix; only used for time-domain problems
#ifdef ENABLE_DIRECTSOLVER
    mutable SuperLU_data lu_data; ///< parameters for LU solver
#endif

    double iterative_tol;   
};

// ==========================================================================
// template typedefs

typedef TDGFwdSolver<double>         RDGFwdSolver;
typedef TDGFwdSolver<std::complex<double> > CDGFwdSolver;

#ifndef __DGFWDSOLVER_CC
extern template class STOASTLIB TDGFwdSolver<double>;
extern template class STOASTLIB TDGFwdSolver<std::complex<double> >;
#endif // !__DGFWDSOLVER_CC


#endif
