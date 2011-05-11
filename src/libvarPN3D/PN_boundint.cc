#include "PN_boundint.h"
#include "toast.h"
#include "sphericalHarmonic_algebra.h"
#include <stdio.h>
/**Computes the maximum spherical harmonic order used in a given element
**/
void findMaxLocalSphOrder(const Mesh &mesh, const IVector& sphOrder, const IVector& node_angN, const int el, int &maxSphOrder, int &maxAngN)
{
  int nodel = mesh.elist[el]->nNode();
  maxAngN = 0; maxSphOrder =0; 
  for(int i=0; i<nodel; i++)
  {
	int is = mesh.elist[el]->Node[i];
	maxSphOrder = MAX(sphOrder[is], maxSphOrder);
	maxAngN = MAX(node_angN[is], maxAngN);

  }

}

/** Adds aB to its appropriate place in the system matrix
* spatrow -> spatial row where 'a' is drawn from
* spatcol -> spatial column where 'a' is drawn from
* node_angN -> number of angular degrees of freedom for all the spatial nodes
* offset -> starting location in the system matrix for each spatial node
* a_ij -> 'a'
* B -> B
* C -> output (System matrix) 
**/
void kronplus(const int spatrow, const int spatcol, const IVector& node_angN, const IVector& offset, const double a_ij, const RCompRowMatrix &B, RCompRowMatrix& C)
{
    double *Cval = C.ValPtr();
    int row_offset = 0, col_offset = 0;
    const int *rowptr, *colidx;
    B.GetSparseStructure(&rowptr, &colidx);
    const double *bval = B.ValPtr(); 
    int nr = B.nRows();
    int nc = B.nCols();
 
    row_offset = offset[spatrow];
	
    col_offset = offset[spatcol];
  
    for(int ib=0; ib < nr; ib++)
    {
	for(int j=rowptr[ib]; j < rowptr[ib+1]; j++)
	{
		int jb = colidx[j];
		if(ib<node_angN[spatrow] && jb <node_angN[spatcol]){
		C(row_offset+ib , col_offset+jb) = C.Get(row_offset+ib , col_offset+jb) + a_ij*bval[j]; 
		}
	}
    }
    
}

/*Computes the integral on the sphere of the form 
	 * \int_{S^{n-1}} (s.n)_{plusminus} \psi_{i} \psi_{j}
	 * Inputs:
	 * 	size1 -> number of rows
	 *	size2 -> number of columns
	 *	sphOrder1 -> order of spherical harmonics along rows
         * 	sphOrder2 -> order of spherical harmonics along columns
	 *      pts -> quadrature points
	 *      wts -> quadrature weights
	 *      bnormal -> outward pointing normal to the boundary
	 *      Ylm -> precomputed table of spherical harmonics over the quadrature points
	 * Output:
	 * 	bintplus -> boundary integral on the outward facing half sphere
	 *      bintminus -> boundary integral on the inward facing half sphere
*/
void BIntUnitSphere(const int size1, const int size2, const int sphOrder1, const int sphOrder2,  const RDenseMatrix& pts, const RVector& wts, const RVector& bnormal, RDenseMatrix* &Ylm, RCompRowMatrix& bintplus, RCompRowMatrix& bintminus)
{
	
   int numpts = pts.nRows();

   /*Finding points which are on the outward facing or inward facing sides*/ 
   RDenseMatrix dnssdotnplus(1, numpts), dnssdotnminus(1, numpts);
   for(int i=0; i < numpts; i++){
	   double tmp = bnormal[0]*pts.Get(i, 0) + bnormal[1]*pts.Get(i, 1) + bnormal[2]*pts.Get(i, 2);
	   if(tmp>0) dnssdotnplus(0, i) = tmp; 
	   else if(tmp<0) dnssdotnminus(0, i) = -tmp;
     }

    // making them sparse for computational efficiency
    RCompRowMatrix sdotnplus = shrink(dnssdotnplus);
    RCompRowMatrix sdotnminus = shrink(dnssdotnminus);
    const int *rowptr1, *colidx1, *rowptr2, *colidx2;
    sdotnplus.GetSparseStructure(&rowptr1, &colidx1);
    sdotnminus.GetSparseStructure(&rowptr2, &colidx2);
       
     //computing integrals by quadrature
     int is, js, indl1, indm1, indl2, indm2;
     RDenseMatrix dnsbintplus(size1, size2), dnsbintminus(size1, size2);
     for(int l1 = 0; l1 <= sphOrder1; l1++){
	indl1 = l1*l1; 
	for(int m1 = -1*l1; m1 <= l1; m1++){
		indm1 = l1 + m1;
		is = indl1 + indm1;
	        for(int l2 = 0; l2 <= sphOrder2; l2++){
			indl2 = l2*l2;
			for(int m2 = -l2; m2 <= l2; m2++){
				indm2 = l2 + m2;
				js = indl2 + indm2;
				for(int i=0; i<sdotnplus.nRows(); i++)
				{
					for(int j=rowptr1[i]; j < rowptr1[i+1]; j++)
						dnsbintplus(is, js) += sdotnplus(i, colidx1[j])*Ylm[l1](indm1, colidx1[j])*Ylm[l2](indm2, colidx1[j])*wts[colidx1[j]]*4*M_PI;
				}
				for(int i=0; i<sdotnminus.nRows(); i++)
				{
					for(int j=rowptr2[i]; j < rowptr2[i+1]; j++)
						dnsbintminus(is, js) += sdotnminus(i, colidx2[j])*Ylm[l1](indm1, colidx2[j])*Ylm[l2](indm2, colidx2[j])*wts[colidx2[j]]*4*M_PI;
				}
			
			}
		}
	
	     }
	}
	
         //shrinking the matrices finally		
         bintplus = shrink(dnsbintplus);
	 bintminus = shrink(dnsbintminus);	

} 

void BRIntUnitSphere(const double nin, const double nout, const int size1, const int size2, const int sphOrder1, const int sphOrder2,  const RDenseMatrix& pts, const RVector& wts, const RVector& bnormal, RDenseMatrix* &Ylm, RCompRowMatrix& bintplus, RCompRowMatrix& bintminus, RCompRowMatrix& brintminus)
{
	
   int numpts = pts.nRows();

   /*Finding points which are on the outward facing or inward facing sides*/ 
   RDenseMatrix dnssdotnplus(1, numpts), dnssdotnminus(1, numpts), dnsrdotnminus(1, numpts);
   for(int i=0; i < numpts; i++)
   {
	   double costheta = bnormal[0]*pts.Get(i, 0) + bnormal[1]*pts.Get(i, 1) + bnormal[2]*pts.Get(i, 2);
	   if(costheta>0) dnssdotnplus(0, i) = costheta; 
	   else if(costheta<0) 
	   {
		dnssdotnminus(0, i) = -costheta;
		double sintheta_sq = 1 - costheta*costheta;
		double term_sq = 1 - (nin*nin*sintheta_sq)/(nout*nout);
		if(term_sq >0)
		{
			double eps = std::numeric_limits<double>::epsilon();

			double nipnt = nin/nout;
			double costht = sqrt( 1.0 - pow(nipnt, 2.0) * (1.0 - pow(costheta, 2.0)) );
			double thi;
				
			if(costheta > 0.0) thi = acos(costheta);
			else thi = acos(-costheta);
			double tht = acos(costht);
			double R;
			if( !(sin(thi + tht) > eps) ) R = pow( (nipnt - 1.0) / (nipnt + 1.0), 2.0 );
			else R = 0.5 * ( pow(sin(thi - tht) / sin(thi + tht), 2.0) + pow(tan(thi - tht) / tan(thi + tht), 2.0) );

				/*double term = -sqrt(term);
				//std::cout<<"nin: "<<nin<<" nout: "<<nout<<" costheta: "<<costheta<<"term in sqrt: "<< (1 - (nin*nin*sintheta_sq)/(nout*nout))<<" term: "<<term<<std::endl; 
				double R_s = (nin*costheta - nout*term)/(nin*costheta + nout*term);
				//std::cout<<"R_s: "<<R_s<<std::endl;
				R_s *= R_s;
				//std::cout<<"After squaring R_s: "<<R_s<<std::endl;
				double R_p = (nin*term - nout*costheta)/(nin*term + nout*costheta);
				//std::cout<<"R_p: "<<R_p<<std::endl;
				R_p *= R_p;
				//std::cout<<"After squaring R_p: "<<R_p<<std::endl;
				//cout<<"term in sqrt: "<< (1 - (nin*nin*sintheta_sq)/(nout*nout))<<" R_s: "<<R_s<<" R_p: "<<R_p<<endl;
				dnsrdotnminus(0, i) = 0.5*costheta*(R_s + R_p);
				std::cout<<"if block "<<0.5*(R_s + R_p)<<std::endl;*/
			dnsrdotnminus(0, i) = costheta*R;
			//std::cout<<"if block "<<R<<std::endl;
		}
		else
		{
		     dnsrdotnminus(0, i) = costheta;
		     //std::cout<<"else block 1"<<std::endl;
		}
	   }
     }
    
    RDenseMatrix H(3, 3);
    for(int i=0; i < 3; i++)
    {
	for(int j=0; j < 3; j++)
	{
		if(i == j) H(i, i) = -2*bnormal[i]*bnormal[i] + 1;
		else H(i, j) = -2*bnormal[i]*bnormal[j];
	}
    }
    // making them sparse for computational efficiency
    RCompRowMatrix sdotnplus = shrink(dnssdotnplus);
    RCompRowMatrix sdotnminus = shrink(dnssdotnminus);
    RCompRowMatrix rsdotnminus = shrink(dnsrdotnminus);
    const int *rowptr1, *colidx1, *rowptr2, *colidx2, *rowptr3, *colidx3;
    sdotnplus.GetSparseStructure(&rowptr1, &colidx1);
    sdotnminus.GetSparseStructure(&rowptr2, &colidx2);
    rsdotnminus.GetSparseStructure(&rowptr3, &colidx3);
 
    RDenseMatrix subPts(rowptr3[1], 3);
    int idx = 0;
    for(int i=0; i < rowptr3[1]; i++)
    {
	subPts(idx, 0) = H(0, 0)*pts.Get(colidx3[i], 0) + H(0, 1)*pts.Get(colidx3[i], 1) + H(0, 2)*pts.Get(colidx3[i], 2);
	subPts(idx, 1) = H(1, 0)*pts.Get(colidx3[i], 0) + H(1, 1)*pts.Get(colidx3[i], 1) + H(1, 2)*pts.Get(colidx3[i], 2);	
	subPts(idx, 2) = H(2, 0)*pts.Get(colidx3[i], 0) + H(2, 1)*pts.Get(colidx3[i], 1) + H(2, 2)*pts.Get(colidx3[i], 2);	
	idx++;	
    } 
    
    RDenseMatrix *YHlm;
    YHlm = new RDenseMatrix[sphOrder1+1];
    for(int l=0; l<=sphOrder1; l++)
		YHlm[l].New(2*l+1, subPts.nRows());
    sphericalHarmonics(sphOrder1, subPts.nRows(), subPts, YHlm);
 
     //computing integrals by quadrature
     int is, js, indl1, indm1, indl2, indm2;
     RDenseMatrix dnsbintplus(size1, size2), dnsbintminus(size1, size2), dnsbrintminus(size1, size2);
     for(int l1 = 0; l1 <= sphOrder1; l1++){
	indl1 = l1*l1; 
	for(int m1 = -1*l1; m1 <= l1; m1++){
		indm1 = l1 + m1;
		is = indl1 + indm1;
	        for(int l2 = 0; l2 <= sphOrder2; l2++){
			indl2 = l2*l2;
			for(int m2 = -l2; m2 <= l2; m2++){
				indm2 = l2 + m2;
				js = indl2 + indm2;
				for(int i=0; i<sdotnplus.nRows(); i++)
				{
					for(int j=rowptr1[i]; j < rowptr1[i+1]; j++)
						dnsbintplus(is, js) += sdotnplus(i, colidx1[j])*Ylm[l1](indm1, colidx1[j])*Ylm[l2](indm2, colidx1[j])*wts[colidx1[j]]*4*M_PI;
				}
				for(int i=0; i<sdotnminus.nRows(); i++)
				{
					for(int j=rowptr2[i]; j < rowptr2[i+1]; j++)
						dnsbintminus(is, js) += sdotnminus(i, colidx2[j])*Ylm[l1](indm1, colidx2[j])*Ylm[l2](indm2, colidx2[j])*wts[colidx2[j]]*4*M_PI;
				}
				for(int i=0; i<rsdotnminus.nRows(); i++)
				{
					int idx =0;
					for(int j=rowptr3[i]; j < rowptr3[i+1]; j++)
					{
						dnsbrintminus(is, js) += rsdotnminus(i, colidx3[j])*YHlm[l1](indm1, idx)*Ylm[l2](indm2, colidx3[j])*wts[colidx3[j]]*4*M_PI;
						idx++;
					}
				}
			
			}
		}
	
	     }
	}
	delete []YHlm;	
         //shrinking the matrices finally		
         bintplus = shrink(dnsbintplus);
	 bintminus = shrink(dnsbintminus);
	 brintminus = shrink(dnsbrintminus);
	
} 


/** Preallocating memory for boundary integral terms.
The matrix is considered to be spatially sparse and dense angularly(since the angular sparsity pattern is not known)
which is eventually shrunk after computation.
**/
void initialiseA2b1(const Mesh &mesh, const IVector &node_angN, const IVector &offset, RCompRowMatrix& A2, RCompRowMatrix& b1)
{
   int el, nodel, i, j, k,is, js;
   int *crrowptr, *crcolidx, nzero;
   int sysdim = mesh.nlen();
   mesh.SparseRowStructure(crrowptr, crcolidx, nzero);
   
   int *status = new int[crrowptr[sysdim]];
   for(i=0; i<nzero; i++)
	status[i] = 0; // 1 implies nonzero 0 denotes zero;

   /*Computing the spatial sparsity pattern corresponding to the boundary*/
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->HasBoundarySide ()) continue;
	nodel = mesh.elist[el]->nNode();
	// now determine the element integrals
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++)  {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
	  for (int i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (int j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		for (int rp = crrowptr[is]; rp < crrowptr[is+1]; rp++){
        		if (crcolidx[rp] == js) status[rp] = 1;
	    }
	  } 
	 }
	}
    }
   
   /*Computing the new spatial sparsity pattern corresponding to the boundary based on 'status' flag*/ 
    int tnzero=0;
    for(i=0; i < nzero; i++)
	if(status[i]) tnzero++; 
   int *spatrowptr, *spatcolidx;
   spatrowptr = new int[sysdim + 1];
   spatcolidx = new int[tnzero];
   spatrowptr[0] = 0;
   j=0;
   for(i = 0; i < sysdim; i++)
   {
	int rp1 = crrowptr[i];
	int rp2 = crrowptr[i+1];
        k=0;
	for(int rp = rp1; rp < rp2; rp++)
	{
		if(status[rp]){ 
			k++;
			spatcolidx[j] = crcolidx[rp];
			j++;			
		} 
	}
	spatrowptr[i+1] = spatrowptr[i] + k;
   } 
   delete []status;

    /*Computing the overall sparsity pattern where angular degrees of freedom are considered dense*/
    int ia, ib, ka, kb, ja, jb, idx;
    int va_i, vb_i, v_i;
    int na = sysdim, ma = sysdim, va = spatrowptr[sysdim], col;
    int *rowptr = new int[sum(node_angN)+1];
    for(i = 0; i < sum(node_angN) + 1; i++)
    	rowptr[i] = 0;
    
    k=1; 
    for (ia = 0; ia < na; ia++)
    {
	for(j = 0; j < node_angN[ia]; j++)
	{
		for(i = spatrowptr[ia]; i < spatrowptr[ia+1]; i++)
		{
			col = spatcolidx[i];
			rowptr[k] += node_angN[col];
		 }

		k++;
        } 
    }
  
    for(i = 1; i < sum(node_angN)+1; i++)
	rowptr[i] += rowptr[i-1];
	
    
   int *colidx = new int[rowptr[sum(node_angN)]];
   int col_offset;
   k=0;
   for(ia = 0; ia < na; ia++)
   {
	for(j = 0; j < node_angN[ia]; j++)
	{	
		for(i = spatrowptr[ia]; i < spatrowptr[ia+1]; i++)
		{
			col = spatcolidx[i];

			col_offset = offset[col];

			/*if(ia == 1669)
				cout<<col<< "  "<<node_angN[col]<<"  "<<node_angN[1669]<<endl;*/

			for(int l = 0; l < node_angN[col]; l++)
			{
				colidx[k] = col_offset + l;
				k++;
			}
		}
	}
   }
   A2.Initialise(rowptr, colidx);
   b1.Initialise(rowptr, colidx);
   A2.Zero(); 
   b1.Zero();

   delete []rowptr;
   delete []colidx;
   delete []spatrowptr;
   delete []spatcolidx;

}

/**Compute the boundary integral terms using quadrature
**/
void genmat_boundint(const Mesh& mesh,  const RVector &ref, const double ref_out, const char bctype, const IVector& sphOrder, const IVector& node_angN, const IVector& offset, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm, RCompRowMatrix& A2, RCompRowMatrix& b1)
{
  
   const int sysdim = mesh.nlen();       // dimensions are size of nodes.
   const int fullsysdim = sum(node_angN);     // full size of angles X space nodes
   
   A2.New(fullsysdim, fullsysdim);
   b1.New(fullsysdim, fullsysdim);
   initialiseA2b1(mesh, node_angN, offset, A2, b1);
   
   double ela_ij;
   int el, nodel, is, js;
   for (el = 0; el < mesh.elen(); el++) {
	if(!(el*100%mesh.elen()))
		std::cout<<el*100/mesh.elen() <<"% elements assembled"<<std::endl;
        if(!mesh.elist[el]->HasBoundarySide ()) continue;

	nodel = mesh.elist[el]->nNode();
       
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++)  {
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
            
	  RVector nhat(3);
	  RVector temp = mesh.ElDirectionCosine(el,sd);
	  if(mesh.Dimension() == 3){
		 nhat[0] = temp[0]; nhat[1] = temp[1]; nhat[2] = temp[2];
	  }
	  else
	  {
		nhat[0] = temp[0]; nhat[1] = temp[1]; nhat[2] = 0;
	  }
	  	  
	  RCompRowMatrix  Angbintplus, Angbintminus, Angbrintminus;
	  int lmaxAngN, lmaxSphOrder;
          findMaxLocalSphOrder(mesh, sphOrder, node_angN, el, lmaxSphOrder, lmaxAngN);
	  BRIntUnitSphere(ref[el], ref_out, lmaxAngN, lmaxAngN, lmaxSphOrder, lmaxSphOrder, pts, wts, nhat, Ylm, Angbintplus, Angbintminus, Angbrintminus);

	  for (int i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (int j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		ela_ij = mesh.elist[el]->BndIntFFSide (i, j,sd);
		if(bctype == 'r'){
			 kronplus(is, js, node_angN, offset, ela_ij, Angbintplus + Angbrintminus, A2);
		}
		else{
			 kronplus(is, js, node_angN, offset, ela_ij, Angbintplus, A2);
		}
	 	kronplus(is, js, node_angN, offset, ela_ij, Angbintminus, b1);	
	    }
	  } 
	
        } // end loop on sides
   	
    
   } // end loop on elements

   //Shrinking the matrices 
   A2.Shrink();
   b1.Shrink();
}

