#ifndef __BEM_MESH_H
#define __BEM_MESH_H

class BEM_Surface;
class BEM_Region;

/**
 * \brief BEM mesh class.
 *
 * A mesh is a collection of surfaces.
 */
class BEMLIB BEM_Mesh {
public:
	BEM_Mesh();

	int AddSurface (BEM_Surface *surf);

	void AddRegion (BEM_Region *region, double mua, double mus, double ref, BEM_Region *parent = NULL);

protected:
	/**
	 * \brief List of regions
	 */
	BEM_Region **Region;

	/**
	 * \brief Number of regions
	 */
	int nRegion;

	/**
	 * \brief List of surfaces
	 */
	BEM_Surface **Surf;

	/**
	 * \brief Number of surfaces
	 */
	int nSurf;

	/**
	 * \brief Modulation frequency [rad/s]
	 */
	double omega;
};

#endif // !__BEM_MESH_H
