#include "detector.h"

void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth)
{
    char cbuf[256];
    mtype = -1;
    if (pp.GetString ("MEASUREMENTPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    mtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    mtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    mtype = 2;
	}
    }
    while (mtype < 0) {
	cout << "\nMeasurement profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> mtype;
	cout << mtype<<endl;
	xASSERT(mtype == 1 || mtype == 2 || mtype == 3, "Unidentified measurement profile.");
	mtype -= 1;
    }
    if (mtype > 0 && !pp.GetReal ("MEASUREMENTWIDTH", mwidth)) {
	switch (mtype) {
	case 1:
	    cout << "\nMeasurement 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nMeasurement support radius [mm]:\n>> ";
	    break;
	}
	cin >> mwidth;
	cout << mwidth << endl;
    }
    switch (mtype) {
    case 0:
	pp.PutString ("MEASUREMENTPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("MEASUREMENTPROFILE", "GAUSSIAN");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    case 2:
	pp.PutString ("MEASUREMENTPROFILE", "COSINE");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    }
}

void genmat_toastdetector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Detector, const Mesh& mesh, const RCompRowMatrix mvec, const int nd, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
{
   int sysdim = mesh.nlen();       
   int fullsysdim = sum(node_angN);   
   for (int i = 0; i < nd; i++) {
     Detector[i].New(fullsysdim,1);
     genmat_toastdetectorvalvector(sphOrder, node_angN, offset, Detector[i], mesh, mvec,i, pts, wts, Ylm);
   }
}

/* Computes the source vector per a boundary source when a QM file has been specified
*/
void genmat_toastdetectorvalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Dvec, const Mesh& mesh, const RCompRowMatrix mvec, const int iq, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       
   int fullsysdim = sum(node_angN);      
   int *rowptr, *colidx;
   rowptr = new int[fullsysdim+1];
   colidx = new int[fullsysdim];
   rowptr[0] = 0;
   for(int i=1;  i <= fullsysdim; i++)
	rowptr[i] = i;
   for(int i=0; i<fullsysdim; i++)
	colidx[i] = 0;
   Dvec.New (fullsysdim,1);
   Dvec.Initialise(rowptr, colidx);

   int row_offset = 0;
   RDenseMatrix Ystarlm;
   RDenseMatrix dirMat(1, 3);

   for (int jq = mvec.rowptr[iq]; jq < mvec.rowptr[iq+1]; jq++) {
   	int Ndetector = mvec.colidx[jq];
   	double sweight = norm(mvec.Get(iq,Ndetector)); 
   	for (el = 0; el < mesh.elen(); el++) {
        	if(!mesh.elist[el]->IsNode(Ndetector)) continue; // source not in this el
        	if(!mesh.elist[el]->HasBoundarySide()) continue;
		for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  		// if sd is not a boundary side. skip 
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
	  		for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    			is = mesh.elist[el]->SideNode(sd,nd);
	    			js = mesh.elist[el]->Node[is];
				double ela_i =  mesh.elist[el]->IntF(is);
				RDenseMatrix Angsvec(node_angN[js], 1);
				for(int l=0; l<= sphOrder[js]; l++){
					int indl = l*l;
					for(int m=-l; m<=l; m++){
						int indm = l + m;
						for(int i=0; i < pts.nRows(); i++)
						{
							double sdotn = nhat[0]*pts.Get(i, 0) + nhat[1]*pts.Get(i, 1) + nhat[2]*pts.Get(i, 2);
							if(sdotn>0)
								Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(sdotn)*wts[i]*4*M_PI;
						}
					}
				}
				for(int i = 0; i < node_angN[js]; i++)
					Dvec(offset[js] + i, 0) += Angsvec(i, 0)*sweight*ela_i;

			}

		} // end loop on element sides

   	} // end loop on elements
  }
   delete []rowptr;
   delete []colidx;
}  

/** Computes source vectors for point sources on the boundary 
**/
void genmat_detector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Detector, const Mesh& mesh,  const IVector& Ndetector, const int nd, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
{
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);       // full size of angles X space nodes

   //RCompRowMatrix Dvec;
   for (int i = 0; i < nd; i++) {
     Detector[i].New(fullsysdim,1);
     genmat_detectorvalvector(sphOrder, node_angN, offset, Detector[i], mesh, Ndetector[i], pts, wts, Ylm);
     //b1.AB(Dvec, Detector[i]);
   }
}

/**Computes source vector for a single point source on the boundary
**/
void genmat_detectorvalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Dvec, const Mesh& mesh, const int Ndetector, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   
   int *rowptr, *colidx;
   rowptr = new int[fullsysdim+1];
   colidx = new int[fullsysdim];
   rowptr[0] = 0;
   for(int i=1;  i <= fullsysdim; i++)
	rowptr[i] = i;
   for(int i=0; i<fullsysdim; i++)
	colidx[i] = 0;
   Dvec.New (fullsysdim,1);
   Dvec.Initialise(rowptr, colidx);

   int row_offset = 0;
   RDenseMatrix dirMat(1, 3);
   RDenseMatrix Ystarlm;

   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Ndetector)) continue; // source not in this el
        if(!mesh.elist[el]->HasBoundarySide()) continue;
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
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
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    double ela_i =  mesh.elist[el]->IntF(is);
	    RDenseMatrix Angsvec(node_angN[js], 1);
	    for(int l=0; l<= sphOrder[js]; l++){
		int indl = l*l;
		for(int m=-l; m<=l; m++){
			int indm = l + m;
			for(int i=0; i < pts.nRows(); i++)
			{
				double sdotn = nhat[0]*pts.Get(i, 0) + nhat[1]*pts.Get(i, 1) + nhat[2]*pts.Get(i, 2);
				if(sdotn>0)
					Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(sdotn)*wts[i]*4*M_PI;
			}
		}
	    }
	   for(int i = 0; i < node_angN[js]; i++)
		Dvec(offset[js] + i, 0) = Angsvec(i, 0)*ela_i;
			
	}
    } 


   } // end loop on elements
   delete []rowptr;
   delete []colidx;
}


