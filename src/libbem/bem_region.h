#ifndef __BEM_REGION_H
#define __BEM_REGION_H

class BEM_Surface;
class BEM_Kernel;

/**
 * \brief Class representation of a single BEM region.
 *
 * Regions are enclosed by surfaces (see \ref BEM_Surface) and define an
 * area of homogeneous optical parameters.
 */
class BEMLIB BEM_Region {
public:
	/**
	 * \brief Region constructor.
	 * Creates a region with an outer surface.
	 * \param outer pointer to outer surface
	 * \param parent pointer to parent region
	 * \note If no parent is specified, the region is assumed to
	 *   be the outermost region of the mesh.
	 */
	BEM_Region (BEM_Surface *outer, BEM_Region *parent = NULL);

	/**
	 * \brief Region destructor.
	 */
	~BEM_Region ();

	/**
	 * \brief Replace the outer surface of the region.
	 * \param outer Pointer to new outer surface
	 */
	void ResetOuterSurface (BEM_Surface *outer);

	/**
	 * \brief reset absorption parameter for the region.
	 * \param newmua new absorption coefficient [1/mm]
	 */
	void SetMua (double newmua) { mua = newmua; }

	/**
	 * \brief reset scattering parameter for the region.
	 * \param newmus new scattering coefficient [1/mm]
	 */
	void SetMus (double newmus) { mus = newmus; }

	/**
	 * \brief reset refractive index for the region.
	 * \param newref new refractive index
	 */
	void SetRef (double newref) { ref = newref; }

	/**
	 * \brief reset modulation frequency for the region
	 * \param newfreq new modulation frequency [Hz]
	 */
	void SetFreq (double newfreq) { freq = newfreq; }

	/**
	 * \brief reset kernel to use for integrations
	 * \param newkernel pointer to new kernel
	 */
	void SetKernel (BEM_Kernel *newkernel) { kernel = newkernel; }

	/**
	 * \brief Add a new region as a child of this region.
	 * \param child pointer to new child region.
	 * \return child index (>= 0)
	 */
	int AddChild (BEM_Region *child);

	/**
	 * \brief Remove a child entry from the list.
	 * \param idx child index (>= 0)
	 */
	void DeleteChild (int idx);

	/**
	 * \brief Notification of child surface change.
	 * \param child pointer to child sending the notification
	 * \note This function is called by a child if it
	 *   changes its outer surface (corresponding to an
	 *   inner surface of the parent).
	 */
	void NotifyChildSurface (BEM_Region *child);

	/**
	 * \brief Perform integration over region surfaces and
	 *   return as region matrices A and B.
	 * \param [out] A region matrix A
	 * \param [out] B region matrix B
	 */
	void ConstructRegionMatrix (CDenseMatrix &A, CDenseMatrix &B);

	/**
	 * \brief Returns an inner surface of the region (given
	 *   by the outer surface of child 'idx').
	 * \param idx child index
	 * \return surface pointer
	 */
	BEM_Surface *GetInnerSurface (int idx);

	std::complex<double> WaveNumberMinus() const;

	std::complex<double> WaveNumber() const;

protected:
	/**
	 * \brief Pointer to outer surface of region
	 * This is only a reference to an externally allocated Surface
	 * object, not a local copy. The Region will not de-allocate
	 * surfaces, and surfaces must remain valid during the lifetime
	 * of the region.
	 */
	BEM_Surface *outerSurf;

	/**
	 * \brief Pointer to parent region, enclosing this region.
	 *   The parent region must contain our outer surface as an
	 *   inner surface. If Parent == NULL, the region has no
	 *   parent (outermost region).
	 */
	BEM_Region *Parent;

	/**
	 * \brief List of pointers to child regions, enclosed directly
	 *   by this region.
	 */
	 BEM_Region **Child;

	 /**
	  * \brief Number of children
	  */
	 int nChildren;

    /**
	 * \Brief Kernel associated with the region
	 */
	BEM_Kernel *kernel;

	/**
	 * \brief Absorption coefficient [1/mm]
	 */
	double mua;

	/**
	 * \brief Scattering coefficient [1/mm]
	 */
	double mus;

	/**
	 * \brief Refractive index.
	 */
	double ref;

	/**
	 * \brief frequency [Hz]
	 */
	double freq;
};

#endif // !__BEM_REGION_H
