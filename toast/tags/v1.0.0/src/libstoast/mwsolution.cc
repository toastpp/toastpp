// Multi-Wavelength Solution Class 17 Aug 2006
// Alper Corlu
//
// In addition to chromphores and Scattering parameters
// this class holds the OT params for each wavelength
// in solutions swsol.
// 
// MWsolution SIZE must match MESH number of nodes!
// 
// The parameters are listed as:
// param[nparam-1] = N (index of Refraction) = param[nmuaChromo+2]
// param[nparam-2] = b (scatter power) = param[nmuaChromo+1]
// param[nparam-3] = A (scatter prefactor) = param[nmuaChromo]
// param[nparam-4] = (nparam-3)rd chromophore
// param[nparam-5] = (nparam-2)nd chromophore
//        :        =      :           :
//        :        =      :           :
//        :        =      :           :

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

#include "mwsolution.h"
#include "solution.h"
#include "param.h"
//#ifdef TOAST_MPI
//#include <mpi.h>
//#endif

using namespace std;

MWsolution::MWsolution (int nparam, int _meshlength, 
			int _nofwavel, RDenseMatrix _extcoef, 
			RVector _wlength)
  :Solution (nparam, _meshlength)
{
  int i;
  nofwavel = _nofwavel;
  meshlength = _meshlength;
  extcoef = _extcoef;
  wlength = _wlength;
  nmuaChromo = extcoef.nCols();

  swsol = new Solution *[nofwavel];
  for (i = 0; i < nofwavel; i++)
    swsol[i] = new Solution (OT_NPARAM, meshlength);
  
}

MWsolution::~MWsolution ()
{
  for (int i = 0; i < nofwavel; i++) 
    delete  swsol[i];   
  delete []swsol;

}

MWsolution::MWsolution (const MWsolution &mwsol)
  : Solution (mwsol)
{
  int i;
  nofwavel = mwsol.nofwavel;
  meshlength = mwsol.meshlength;
  extcoef = mwsol.extcoef;
  wlength = mwsol.wlength;
  nmuaChromo = mwsol.nmuaChromo;

  swsol = new Solution *[nofwavel];
  for (i = 0; i < nofwavel; i++)
    swsol[i] = new Solution (OT_NPARAM, meshlength);
  
}

void MWsolution::RegisterChange()
{
    const double c0 = 0.3;
    RVector covern(meshlength);
    RVector tmp_c2A = GetC2A();

    covern = c0/GetParam(nmuaChromo+2);
    
    for (int i = 0; i < nofwavel; i++) {
	RVector tmp_cMua(meshlength);
	RVector tmp_ckappa(meshlength);
	RVector tmp_Musp(meshlength);
	tmp_cMua = GetParam(nmuaChromo+3+i); // background mua
	for (int j = 0; j < nmuaChromo; j++)
	    tmp_cMua += extcoef(i, j) * GetParam(j);

	for (int k = 0; k < meshlength; k++) {
	    tmp_Musp[k] = param[nmuaChromo][k] * 
		pow(wlength[i],-param[nmuaChromo+1][k]);
	}

	tmp_ckappa = covern/(3.0*(tmp_cMua + tmp_Musp));
	tmp_cMua *= covern;

	swsol[i]->SetParam(OT_CMUA, tmp_cMua);
	swsol[i]->SetParam(OT_CKAPPA, tmp_ckappa);
	swsol[i]->SetParam(OT_C2A, tmp_c2A);
    }
}

RVector MWsolution::GetJacobianCoeff_A (int wavelind) const
{
  const double c0 = 0.3;
  RVector jcoeff(meshlength);
  RVector covern(meshlength);
  covern = c0/GetParam(nmuaChromo+2);
  
  RVector tmp_ckappa(meshlength);
  tmp_ckappa = swsol[wavelind]->GetParam(OT_CKAPPA);
  
  for (int k = 0; k < meshlength; k++) 
    jcoeff[k] = -3.0/covern[k] * tmp_ckappa[k] * tmp_ckappa[k]
      * pow(wlength[wavelind],-param[nmuaChromo+1][k]);

  return jcoeff;
}

RVector MWsolution::GetJacobianCoeff_b (int wavelind) const
{
  
  RVector jcoeff(meshlength);
  jcoeff = -1.0 * log(wlength[wavelind]) * GetParam(nmuaChromo) 
                * GetJacobianCoeff_A (wavelind);

  return jcoeff; 
}

RVector MWsolution::GetC2A () const
{
    RVector c2a(meshlength);
    for (int k = 0; k < meshlength; k++) {
        double n = param[nmuaChromo+2][k];
	c2a[k] = Parameter::C2A(REFLECTION_KEIJZER, n);
    }
    return c2a;
}
