#define BEMLIB_IMPLEMENTATION

#include "bemlib.h"
#include "bem_element.h"
#include "bem_surface.h"

BEM_Element::BEM_Element (BEM_Surface *s, IVector &ndidx)
{
	surf = s;
}

void BEM_Element::Initialise (IVector &ndidx)
{
	int i;

	node = new int[nNode()];
	for (i = 0; i < nNode(); i++)
		node[i] = ndidx[i];
}

RVector BEM_Element::Jacobi(Point2D &loc) const
{
 
  //==================================================================
  //  purpose: Calculate the Jacobian matrix entries
  //  
  //==================================================================

  int icd, iic;
  int dim = nNode(), dim2 = dim*2;
  double dx_d_ksi1, dy_d_ksi1, dz_d_ksi1, dx_d_ksi2, dy_d_ksi2, dz_d_ksi2, det_jacob, nx, ny, nz;
  Point *p = surf->NodeList();
  static RVector jac_tab(dim+1), shapd(dim2);

  // local coordinate system
  dx_d_ksi1 = 0.; dy_d_ksi1 = 0.; dz_d_ksi1 = 0.; 
  dx_d_ksi2 = 0.; dy_d_ksi2 = 0.; dz_d_ksi2 = 0.;

  shapd=ShapeD(loc);
  
  for(int ic=0; ic<dim;++ic)
    {
      iic=node[ic]; icd=ic+dim;
      dx_d_ksi1 += shapd[ic]*p[iic][0];
      dy_d_ksi1 += shapd[ic]*p[iic][1];
      dz_d_ksi1 += shapd[ic]*p[iic][2];
      
      dx_d_ksi2 += shapd[icd]*p[iic][0];
      dy_d_ksi2 += shapd[icd]*p[iic][1];
      dz_d_ksi2 += shapd[icd]*p[iic][2];
    }   
  
  // calculate the unit normal components
  nx = dy_d_ksi1*dz_d_ksi2-dz_d_ksi1*dy_d_ksi2;
  ny = dz_d_ksi1*dx_d_ksi2-dx_d_ksi1*dz_d_ksi2;
  nz = dx_d_ksi1*dy_d_ksi2-dy_d_ksi1*dx_d_ksi2;
   
  // calculate the determinant of Jacobian as a function of point "ksi1, ksi2"
  det_jacob = 1./sqrt(nx*nx + ny*ny + nz*nz);
  
  // jac_tab table collect the Jacobian and directional components 
  // of the outward normal vector to the boundary (numbering system-anticlocwise) 
  jac_tab[0] = 1./det_jacob;
  jac_tab[1] = -nx*det_jacob; 
  jac_tab[2] = -ny*det_jacob; 
  jac_tab[3] = -nz*det_jacob;

  return jac_tab;
}
