#include "PN_spatint.h"
/*Generating all the spatial matrices required by variable order PN method*/
void gen_spatint(const QMMesh& mesh, const RVector& muabs, const RVector& muscat, const RVector& ref, const RVector& delta, double w, double c, CCompRowMatrix& Sint, CCompRowMatrix& Sdx, CCompRowMatrix& Sdy, CCompRowMatrix& Sdz, RCompRowMatrix& Sx, RCompRowMatrix& Sy, RCompRowMatrix& Sz, RCompRowMatrix& Sdxx, RCompRowMatrix& Sdxy, RCompRowMatrix& Sdyx, RCompRowMatrix& Sdyy, RCompRowMatrix& Sdxz, RCompRowMatrix& Sdzx, RCompRowMatrix& Sdyz, RCompRowMatrix& Sdzy, RCompRowMatrix& Sdzz, RCompRowMatrix& spatA3_rte, RCompRowMatrix& spatA3_sdmx, RCompRowMatrix& spatA3_sdmy, RCompRowMatrix& spatA3_sdmz, RCompRowMatrix& SPS, RCompRowMatrix& SPSdx, RCompRowMatrix& SPSdy, RCompRowMatrix& SPSdz)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

   int *rowptr, *colidx, nzero;
   mesh.SparseRowStructure (rowptr, colidx, nzero);
   
   Sint.New (sysdim, sysdim);
   Sint.Initialise (rowptr, colidx);
   Sx.New (sysdim, sysdim);
   Sx.Initialise (rowptr, colidx);
   Sy.New (sysdim, sysdim);
   Sy.Initialise (rowptr, colidx);
   Sz.New (sysdim, sysdim);
   Sz.Initialise (rowptr, colidx);
   SPS.New (sysdim, sysdim);
   SPS.Initialise (rowptr, colidx); 
   spatA3_rte.New (sysdim, sysdim);
   spatA3_rte.Initialise (rowptr, colidx);
   Sdx.New (sysdim, sysdim);
   Sdx.Initialise (rowptr, colidx);
   Sdy.New (sysdim, sysdim);
   Sdy.Initialise (rowptr, colidx);
   Sdz.New (sysdim, sysdim);
   Sdz.Initialise (rowptr, colidx);
   Sdxx.New (sysdim, sysdim);
   Sdxx.Initialise (rowptr, colidx);
   Sdxy.New(sysdim, sysdim);
   Sdxy.Initialise (rowptr, colidx);
   Sdyx.New (sysdim, sysdim);
   Sdyx.Initialise (rowptr, colidx);
   Sdyy.New (sysdim, sysdim);
   Sdyy.Initialise (rowptr, colidx);
   Sdxz.New (sysdim, sysdim);
   Sdxz.Initialise (rowptr, colidx);
   Sdzx.New(sysdim, sysdim);
   Sdzx.Initialise(rowptr, colidx);
   Sdyz.New (sysdim, sysdim);
   Sdyz.Initialise (rowptr, colidx);
   Sdzy.New(sysdim, sysdim);
   Sdzy.Initialise(rowptr, colidx);
   Sdzz.New (sysdim, sysdim);
   Sdzz.Initialise (rowptr, colidx);
   SPSdx.New (sysdim, sysdim);
   SPSdx.Initialise (rowptr, colidx); 
   spatA3_sdmx.New (sysdim, sysdim);
   spatA3_sdmx.Initialise (rowptr, colidx);
   SPSdy.New (sysdim, sysdim);
   SPSdy.Initialise (rowptr, colidx); 
   spatA3_sdmy.New (sysdim, sysdim);
   spatA3_sdmy.Initialise (rowptr, colidx);
   SPSdz.New (sysdim, sysdim);
   SPSdz.Initialise (rowptr, colidx); 
   spatA3_sdmz.New (sysdim, sysdim);
   spatA3_sdmz.Initialise (rowptr, colidx);

   
   double elk_ij, elb_ij, elsx_ij, elsy_ij, elsz_ij;
   int el, nodel, i, j, k,is, js;


   
   for (el = 0; el < mesh.elen(); el++) {
	nodel = mesh.elist[el]->nNode();
	
	double dss = delta[el]; //streamline diffusion value for this element.
	double sigmatot = muabs[el] + muscat[el];
	
	RSymMatrix eldd = mesh.elist[el]->Intdd(); // "all at once!""
	for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		
		elb_ij = mesh.elist[el]->IntFF (i, j);
		Sint(is,js) += toast::complex(0, w*ref[el]*elb_ij/c); 
		SPS(is, js) += elb_ij*muscat[el];
		spatA3_rte(is, js) += elb_ij*sigmatot;


		elsx_ij = mesh.elist[el]->IntFd (j,i,0);
		Sx(is,js) += elsx_ij;
		Sdx(is,js) += toast::complex(0, w*ref[el]*dss*elsx_ij/c);
		SPSdx(is, js) += dss*elsx_ij*muscat[el];
		spatA3_sdmx(is, js) += dss*elsx_ij*sigmatot;

		elsy_ij = mesh.elist[el]->IntFd (j,i,1);
		Sy(is,js) += elsy_ij;
		Sdy(is,js) += toast::complex(0, w*ref[el]*dss*elsy_ij/c);
		SPSdy(is, js) += dss*elsy_ij*muscat[el];
		spatA3_sdmy(is, js) += dss*elsy_ij*sigmatot;
	
		if(mesh.elist[el]->Dimension() == 3)
		{
			elsz_ij = mesh.elist[el]->IntFd (j,i,2);
			Sz(is,js) += elsz_ij;
  			Sdz(is,js) += toast::complex(0, w*ref[el]*dss*elsz_ij/c);
			SPSdz(is, js) += dss*elsz_ij*muscat[el];
			spatA3_sdmz(is, js) += dss*elsz_ij*sigmatot;
		}
		int dim = mesh.elist[el]->Dimension();
		Sdxx(is,js) += dss * eldd(i*dim,j*dim);
	       	Sdxy(is,js) += dss * eldd(i*dim,j*dim+1);
     		Sdyx(is,js) += dss * eldd(i*dim+1,j*dim);
       		Sdyy(is,js) += dss * eldd(i*dim+1,j*dim+1);
		if(mesh.elist[el]->Dimension() == 3)
		{
	       		Sdxz(is,js) += dss * eldd(i*dim,j*dim+2);
     			Sdzx(is,js) += dss * eldd(i*dim+2,j*dim);
	       		Sdyz(is,js) += dss * eldd(i*dim+1,j*dim+2);
     			Sdzy(is,js) += dss * eldd(i*dim+2,j*dim+1);
       			Sdzz(is,js) += dss * eldd(i*dim+2,j*dim+2);	
		}
	    }
	}
   }

  delete []rowptr;
  delete []colidx;
}

