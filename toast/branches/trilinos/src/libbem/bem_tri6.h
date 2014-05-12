#ifndef __BEM_TRIANGLE6_H
#define __BEM_TRIANGLE6_H

#include "bemlib.h"

/**
 * \brief BEM surface element: 6-noded triangle.
 */
class BEMLIB BEM_Triangle6: public BEM_Element {
public:
	BEM_Triangle6 (BEM_Surface *s, IVector &ndidx);
	~BEM_Triangle6();

	int nNode() const { return 6; }

	/**
	 * \brief Singular integration over element for specified point
	 * \param nodep node index of loading point
	 * \param point reference point
	 * \param invert flag for inverting element normals	
	 */
	CVector Integrate_Singular (BEM_Kernel *kernel, int nodep, const Point3D &point, bool invert = false);

	/**
	 * \brief Nonsingular integral of a kernel over the element
	 * \param kernel kernel instance
	 * \param load load point (global coordinates)
	 * \param invert ???
	 */
	CVector Integrate_Nonsingular (BEM_Kernel *kernel, const Point3D &load, bool invert = false);

	/**
	 * \brief Returns the shape functions for all element nodes at
	 *   point 'loc', given in local (2D) element coordinates.
	 * \param loc local point coordinates
	 */
	RVector ShapeF (Point2D &loc) const;

	/**
	 * \brief Returns the shape function derivatives for all element nodes at
	 *  point 'loc', given in local (2D) element coordinates
	 */
	RVector ShapeD (Point2D &loc) const;

	/**
	 * \brief Returns the coordinates of quadrature point i in local coordinates.
	 * \param i quadrature point index
	 */
	Point2D &QuadraturePoint (int i) const;
	/**
	 * \brief Returns matrix of dimension m x n of m quadrature points for
	 *   n nodes.
	 * \param i quadrature point index
	 */
	RVector &QuadratureShapeF (int i) const;

	/**
	 * \brief Returns matrix of dimension m x 2n of m quadrature points for
	 *   n nodes.
	 * \param i quadrature point index
	 */
	RVector &QuadratureShapeD (int i) const;

	/**
	 * \brief Returns the Jacobi matrix for quadrature point i
	 */
	RVector &QuadratureJacobi (int i) const;

private:
	mutable RVector *qj; // stores the Jacobi quadrature points
};

#endif // !__BEM_TRIANGLE6_H
