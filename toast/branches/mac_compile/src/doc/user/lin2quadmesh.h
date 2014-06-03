/**
 * \page lin2quadmesh lin2quadmesh: Convert to quadratic elements
 *
 * lin2quadmesh is a filter to convert a toast mesh from linear to quadratic
 * element types.
 * The mesh is read from standard input and written to standard output.
 *
 * \section syntax Syntax
 * \code
 * lin2quadmesh < inmesh > outmesh
 * \endcode
 *
 * \section descr Description
 * The following element types are converted:
 * - 3-noded triangles -> 6-noded triangles
 * - 4-noded tetrahedra -> 10-noded tetrahedra
 *
 * The total number of elements is retained, while the total number of nodes
 * increases.
 */
