/**
 * \page quad2linmesh quad2linmesh: Convert to linear elements
 *
 * quad2linmesh is a filter to convert a toast mesh from quadratic to linear
 * element types.
 * The mesh is read from standard input and written to standard output.
 *
 * \section syntax Syntax
 * \code
 * quad2linmesh < inmesh > outmesh
 * quad2linmesh -r < inmesh > outmesh
 * \endcode
 *
 * \section descr Description
 * The following element types are converted:
 * - 6-noded triangles -> 3-noded triangles
 * - 10-noded tetrahedra -> 4-noded tetrahedra
 *
 * There are two modes of operation:
 * - By default, the number of elements is retained, and all non-vertex nodes
 *   are removed.
 * - With command line parameter -r ('refine'), all nodes are retained, and
 *   each quadratic element is subdivided into multiple linear elements.
 *
 * The second method is a convenient way to subdivide meshes:
 * \code
 * lin2quadmesh < coarse_lin.msh | quad2linmesh -r > fine_lin.msh
 * \endcode
 */
