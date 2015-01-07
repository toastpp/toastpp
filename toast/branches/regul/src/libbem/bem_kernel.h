#ifndef __BEM_KERNEL
#define __BEM_KERNEL

#define BEMLIB_IMPLEMENTATION

#include "bemlib.h"

/**
 * \brief Base class for BEM kernels.
 */
class BEMLIB BEM_Kernel {
public:
	BEM_Kernel () {}

	/**
	 * \brief reset the wave number for the following calculations
	 * \param newk new wave number
	 */
	void SetWavenumber (std::complex<double> newk) { k = newk; }

	/**
	 * \brief Calculates the Greens function for the element at
	 *   local point 'loc' from loading point 'load'.
	 * \param loc local point on the element
	 * \param load loading point
	 * \param k wave number
	 */
	virtual CVector Calculate (BEM_Element *el, Point2D &loc, const Point3D &load) = 0;

protected:
	std::complex<double> k; ///< wave number
};

/**
 * \brief BEM Helmholtz kernel.
 */
class BEMLIB BEM_Kernel_Helmholtz: public BEM_Kernel {
public:
	BEM_Kernel_Helmholtz (): BEM_Kernel () {}

	CVector Calculate (BEM_Element *el, Point2D &loc, const Point3D &load);
};


/**
 * \brief BEM Helmholtz derivative kernel dG/domega (used for moments calculation)
 */
class BEMLIB BEM_Kernel_Domega: public BEM_Kernel {
public:
    BEM_Kernel_Domega (): BEM_Kernel () {}
    CVector Calculate (BEM_Element *el, Point2D &loc, const Point3D &load);
};


/**
 * \brief BEM Helmholtz derivative kernel dG/dmua (used for derivatives calculation)
 */
class BEMLIB BEM_Kernel_Dk: public BEM_Kernel {
public:
    BEM_Kernel_Dk (): BEM_Kernel () {}
    CVector Calculate (BEM_Element *el, Point2D &loc, const Point3D &load);
};
#endif // !__BEM_KERNEL
