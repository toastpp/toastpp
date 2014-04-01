// Multi-Wavelength Solution Class 17 Aug 2006

#ifndef __MWSOLUTION_H
#define __MWSOLUTION_H
#include "mathlib.h"
#include "solution.h"

class MWsolution: public Solution {

public:
  MWsolution (int nparam, int _meshlength, int _nofwavel,
	      RDenseMatrix _extcoef, RVector _wlength);
  // nparam is the total number of absorption chromophores +
  // 2 parameters, scattering prefactor A, and power b

  MWsolution (const MWsolution &mwsol);
  ~MWsolution();
  
  void RegisterChange();
  // convert chromophores, Scatter Params and N to cmua, ckappa, n

  RVector GetJacobianCoeff_A (int wavelind) const;
  // del(cKappa)/del(A)
  
  RVector GetJacobianCoeff_b (int wavelind) const;
  // del(cKappa)/del(b)

  RVector GetC2A () const;
  // refractive index parameter c/(2A)

  Solution **swsol; // single wavelength mesh basis solutions

  RDenseMatrix extcoef;
  int nofwavel;
  int meshlength;
  int nmuaChromo;
  RVector wlength;
};
  
#endif // !__MWSOLUTION_H
