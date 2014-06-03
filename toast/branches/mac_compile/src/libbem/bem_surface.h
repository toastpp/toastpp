#ifndef __BEM_SURFACE_H
#define __BEM_SURFACE_H

class BEM_Element;
/**
 * \brief Class representation of a single closed BEM surface.
 */
class BEMLIB BEM_Surface {
public:
	/**
	 * \brief Constructs a surface from a matrix N of node
	 *   coordinates, and a matrix E of element indices.
	 * \param N node coordinate matrix (n x d) where n is
	 *   number of nodes, and d is dimension (probably always 3)
	 * \param E element index matrix (m x s) where m is number
	 *   of elements, and s is nodes per element (assumed the
	 *   same for all elements). Index entries are zero-based.
	 */
	BEM_Surface (RDenseMatrix &N, IDenseMatrix &E); 

	/**
	 * \brief Construct a BEM surface from a FEM volume mesh.
	 * \param mesh FEM mesh (currently must be 10-noded tetrahedra)
	 */
	BEM_Surface (Mesh &mesh);

	/**
	 * \brief Constructs a surface from a stream.
	 * \param is input stream instance
	 * \note The stream is assumed to contain a mesh
	 *   in toast format.
	 */
	BEM_Surface (std::istream &is);

	/**
	 * \brief Returns the node coordinates as a matrix
	 */
	RDenseMatrix CoordinateMatrix() const;

	/**
	 * \brief Returns the element index array as a matrix
	 */
	IDenseMatrix IndexMatrix() const;

	/**
	 * \brief Integration over element 'el' for an arbitrary
	 *   reference point outside the element.
	 * \param point reference point
	 * \param el surface element index
	 * \param invert flag for inverting element normals
	 */
	CVector Integrate_Nonsingular (BEM_Kernel *kernel, const Point3D &point, int el, bool invert = false);

	/**
	 * \brief Integration over element 'el' for a reference
	 *   point given by node 'nd' of the surface. This selects
	 *   singular or nonsingular integration as required.
	 * \param nd Surface node index
	 * \param el Surface element index
	 * \param invert flag for inverting element normals
	 */
	CVector Integrate (BEM_Kernel *kernel, int nd, int el, bool invert = false);

	/**
	 * \brief Returns the node list for the surface
	 * \return List of pointers to nodes.
	 */
	inline Point3D *NodeList() const { return Nd; }

	/**
	 * \brief Returns the number of nodes in the surface
	 */
	inline int nNodes() const { return nNd; }

	/**
	 * \brief Returns pointer to an element
	 * \param i element index (>= 0)
	 */
	inline BEM_Element *Element (int i) const { return El[i]; }

	/**
	 * \brief Returns the number of elements in the surface
	 */
	inline int nElements() const { return nEl; }

    // I/O
	friend BEMLIB std::ostream& operator<< (std::ostream& i, BEM_Surface& surf);

protected:
	/**
	 * \brief List of surface nodes
	 */
	Point3D *Nd;

	/**
	 * \brief Number of surface nodes
	 */
	int nNd;

	/**
	 * \brief List of surface elements
	 */
	BEM_Element **El;

	/**
	 * \brief Number of surface elements
	 */
	int nEl;
};

#endif // !__BEM_SURFACE_H
