/**
 * \page opt2 opt2: Mesh optimiser
 *
 * opt2 is a filter that re-sorts the node order of a mesh to optimise its
 * performance with the linear solvers of the finite element package. opt2
 * reads a mesh from standard input, and writes the optimised mesh to
 * standard output.
 *
 * \section syntax Syntax:
 * \code
 * opt2 -H
 * opt2 < inmesh > outmesh
 * opt2 -m < inmesh > outmesh
 * opt2 -b < inmesh > outmesh
 * opt2 -t < inmesh > outmesh
 * \endcode
 *
 * \section descr Description
 * opt2 re-arranges the mesh nodes such that the amount of
 * fill-in of the decomposition of the system matrix is minimised.
 *
 * Fill-in occurs because the decomposition is less sparse than the original
 * system matrix. It can increase the memory demand for the decomposed matrix
 * significantly.
 *
 * Mesh optimisation is most critical when using direct solver methods
 * (Cholesy and LU), but can also affect the performance of iterative
 * solvers, in particular in connection with preconditioners (ICH, ILU).
 *
 * opt2 supports three different optimisation strategies:
 *
 * \subsection sub1 Minimum degree (command line option -m)
 * Graph-based optimisation. This is the default method, so the -m flag is
 * optional.
 *
 * \subsection sub2 Tinney-2 optimisation (command line option -t)
 * Another graph-based method, usually somewhat less effective than -m.
 *
 * \subsection sub3 Bandwidth minimisation (command line option -b)
 * Arranges nonzero elements close to the diagonal. Decomposition does not
 * increase the matrix bandwidth, so the smaller the bandwidth of the
 * original matrix, the less fill-in occurs in the decomposition.
 *
 * opt2 -H lists the supported command line options.
 */
