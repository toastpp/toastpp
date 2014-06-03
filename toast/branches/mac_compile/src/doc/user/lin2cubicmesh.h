/**
 * \page lin2cubicmesh lin2cubicmesh: Convert to cubic elements
 *
 * lin2cubicmesh is a filter to convert a toast mesh from linear to cubic
 * element types.
 * The mesh is read from standard input and written to standard output.
 *
 * \section syntax Syntax
 * \code
 * lin2cubicmesh < inmesh > outmesh
 * \endcode
 *
 * \section descr Description
 * Currently, triangluar elements are the only element types in toast that
 * support cubic shape function:
 * - 3-noded triangles -> 10-noded triangles
 *
 * The total number of elements is retained, while the total number of nodes
 * increases.
 */
