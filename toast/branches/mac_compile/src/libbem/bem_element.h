// -*-C++-*-

#ifndef __BEM_ELEMENT_H
#define __BEM_ELEMENT_H

class BEM_Surface;
class BEM_Kernel;

// ==========================================================================
// Module libbem
// File bem_element.h
// Base class for BEM element types
// (Should later be derived from libfe/element.h)
// ==========================================================================

/**
 * \brief Base class for 3D surface elements for the BEM solver.
 */
class BEMLIB BEM_Element {
public:
	/**
	 * \brief Constructs an element as a member of surface 'surf'.
	 * \param surf pointer to surface
	 */
	BEM_Element (BEM_Surface *s, IVector &ndidx);

	/**
	 * \brief Singular integration over element for specified point
	 * \param kernel kernel to use for a particular Green's function
	 * \param nodep node index of loading point
	 * \param point reference point
	 * \param invert flag for inverting element normals	
	 */
	virtual CVector Integrate_Singular (BEM_Kernel *kernel, int nodep, const Point3D &point, bool invert = false) = 0;

	/**
	 * \brief Nonsingular integration over element for specified point
	 * \param kernel kernel to use for a particular Green's function
	 * \param point reference point
	 * \param invert flag for inverting element normals	
	 */
	virtual CVector Integrate_Nonsingular (BEM_Kernel *kernel, const Point3D &point, bool invert = false) = 0;

	/**
	 * \brief Returns the global node number of the i-th node of the element.
	 */
	inline int NodeIndex (int i) const { return node[i]; }

	/**
	 * \brief Returns the number of nodes in the element.
	 */
	virtual int nNode() const = 0;

	/**
	 * \brief Returns the surface the element belongs to
	 */
	const BEM_Surface *Surface() const { return surf; }

	/**
	 * \brief Returns the shape functions for all element nodes at
	 *   point 'loc', given in local (2D) element coordinates.
	 * \param loc local point coordinates
	 */
	virtual RVector ShapeF (Point2D &loc) const = 0;

	/**
	 * \brief Returns matrix of dimension m x n of m quadrature points for
	 *   n nodes.
	 */
	virtual RVector &QuadratureShapeF (int i) const = 0;

	/**
	 * \brief Returns the shape function derivatives for all element nodes at
	 *  point 'loc', given in local (2D) element coordinates
	 */
	virtual RVector ShapeD (Point2D &loc) const = 0;

	/**
	 * \brief Calculates the Jacobian for the element at point 'loc'
	 */
	RVector Jacobi(Point2D &loc) const;

protected:
	virtual void Initialise (IVector &ndidx);

	/**
     * \brief List of element node indices (>= 0) relative to surface
	 */
	int *node;

	/**
	 * \brief Pointer to surface
	 */
	const BEM_Surface *surf;

private:
};

#endif // !__BEM_ELEMENT_H
