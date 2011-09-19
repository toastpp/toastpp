#include "source.h"
#include "PN_incl.h"

using namespace std;
void genmat_toastsourcevalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const RCompRowMatrix qvec, const int iq, RVector& dirVec, const int srctp,  const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);
void genmat_sourcevalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const int Nsource, RVector& dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);
void genmat_intsourcevalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const int Nsource, RVector& dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm, const RVector &delta);
void genmat_intsourceuncollidedvalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, Mesh& mesh, const int Nsource, RVector& dirVec, const RDenseMatrix& pts, const RVector& wts, const RVector &mua, const RVector &mus, const RVector &delta, const double g);

void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth, SourceMode &srctp)
{
    char cbuf[256];
    int cmd;

    bool typeok = false;
    if (pp.GetString ("SOURCETYPE", cbuf)) {
	if (!strcasecmp (cbuf, "NEUMANN")) {
	    srctp = SRCMODE_NEUMANN;
	    typeok = true;
	} else if (!strcasecmp (cbuf, "ISOTROPIC")) {
	    srctp = SRCMODE_ISOTROPIC;
	    typeok = true;
	}
    }
    while (!typeok) {
	cout << "\nSource type:\n";
	cout << "(1) Neumann boundary source\n";
	cout << "(2) Isotropic point source\n";
	cout << "[1|2] >> ";
	cin  >> cmd;
	cout<<cmd<<endl;
	xASSERT(cmd == 1 || cmd == 2, "Unidentified source type.");
	switch (cmd) {
	    case 1: srctp = SRCMODE_NEUMANN;   typeok = true; break;
	    case 2: srctp = SRCMODE_ISOTROPIC; typeok = true; break;
	}
    }
    pp.PutString ("SOURCETYPE",
        srctp == SRCMODE_NEUMANN ? "NEUMANN" : "ISOTROPIC");

    qtype = -1;
    if (pp.GetString ("SOURCEPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    qtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    qtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    qtype = 2;
	}
    }
    while (qtype < 0) {
	cout << "\nSource profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> qtype;
	cout<<qtype<<endl;
	xASSERT(qtype == 1 || qtype == 2 || qtype == 3, "Unidentified source profile.");
	qtype -= 1;
    }
    if (qtype > 0 && !pp.GetReal ("SOURCEWIDTH", qwidth)) {
	switch (qtype) {
	case 1:
	    cout << "\nSource 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nSource support radius [mm]:\n>> ";
	    break;
	}
	cin >> qwidth;
	cout<<qwidth<<endl;
    }
    switch (qtype) {
    case 0:
	pp.PutString ("SOURCEPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("SOURCEPROFILE", "GAUSSIAN");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    case 2:
	pp.PutString ("SOURCEPROFILE", "COSINE");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    }
}

RVector QVec_Gaussian (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode)
{
    int n = mesh.nlen();
    RVector qvec(n);
    Element *pel;
    int i, j, is, js, el, nnode, *node;
    double d, q, w2 = w*w;
    double fac1 = 1.0 / sqrt (2.0*M_PI * w2);
    double fac2 = -0.5/w2;

    // assemble source vector
    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    if ((is = node[i]) >= n) continue;
	    if (mode == SRCMODE_NEUMANN && !mesh.nlist[is].isBnd()) continue;

	    // source strength at node is
	    d = cnt.Dist (mesh.nlist[is]);
	    q = exp (d*d*fac2) * fac1;

	    for (j = 0; j < nnode; j++) {
		if ((js = node[j]) >= n) continue;
		if (mode == SRCMODE_NEUMANN) {
		    if (mesh.nlist[is].isBnd())
		        qvec[js] += q * pel->BndIntFF (i,j);
		} else
		    qvec[js] += q * pel->IntFF (i,j);
	    }
	}
    }
    return qvec;
}


// ============================================================================

RVector QVec_Cosine (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode)
{
    double scale = 0.5*M_PI/w;
    int n = mesh.nlen();
    int i, j, el, sd, nnode, *node, nabsc = 0;
    double *wght, val, tval = 0.0, d;
    Point *absc;
    RVector *F;
    RDenseMatrix *D;
    Element *pel;
    RVector qvec (n);
    double fac1 = 1.0 / sqrt (2.0*M_PI * w*w);
    double fac2 = -0.5/(w*w);

    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	if (mode == SRCMODE_NEUMANN) {
	    for (sd = 0; sd < pel->nSide(); sd++) {
		if (!pel->IsBoundarySide (sd)) continue;
//#define USE_SUBSAMPLING
#ifdef USE_SUBSAMPLING
		nabsc = pel->GetBndSubsampleFD (sd, nabsc, wght, absc, F, D,
						mesh.nlist);
		for (i = 0; i < nabsc; i++) {
		    d = cnt.Dist (absc[i]);
		    val = (d < w ? cos (scale*d) : 0.0);
		    if (val) {
			val *= wght[i];
			// scale with interval width
			val *= M_PI/(4.0*w);
			// normalise with integral over cosine
			for (j = 0; j < nnode; j++)
			    if (mesh.nlist[node[j]].isBnd())
				qvec[node[j]] += val * F[i][j];
		    }
		}
#else
		RVector ucos = pel->BndIntFCos (sd, cnt, w, mesh.nlist);
		ucos *= M_PI/(4.0*w);
		// normalise with integral over cosine
		for (i = 0; i < pel->nSideNode (sd); i++)
		    qvec[node[pel->SideNode (sd, i)]] += ucos[i];
#endif
	    }
	} else {
	    nabsc = pel->GetSubsampleFD (nabsc, wght, absc, F, D, mesh.nlist);
	    for (i = 0; i < nabsc; i++) {
		d = cnt.Dist (absc[i]);
		val = (d < w ? cos (scale*d) : 0.0);
		if (val) {
		    val *= wght[i];
		    tval += val;
		    for (j = 0; j < nnode; j++)
			qvec[node[j]] += val * F[i][j];
		}
	    }
	}
    }
    if (mode == SRCMODE_ISOTROPIC) qvec /= tval;
    return qvec;
}

// ============================================================================

RVector QVec_Point (const Mesh &mesh, const Point &cnt, SourceMode mode)
{
    int i, nd, n = mesh.nlen();
    const NodeList &nlist = mesh.nlist;
    const ElementList &elist = mesh.elist;

    int el = mesh.ElFind (cnt);
    if (el < 0) { // point outside mesh
	double d, d2, dmin = 1e10, dmin2 = 1e10;
	const double eps = 1e-8;
	int j, k;
	for (i = 0; i < mesh.elen(); i++) {
	    for (j = 0; j < elist[i]->nNode(); j++) {
		nd = elist[i]->Node[j];
		d = nlist[nd].Dist(cnt);
		if (d < dmin+eps) {
		    for (k = 0; k < elist[i]->nNode(); k++) {
			if (k == j) continue;
			d2 = nlist[elist[i]->Node[k]].Dist(cnt);
			if (d2 < dmin2) {
			    dmin2 = d2;
			    dmin = d;
			    el = i;
			}
		    }
		}
	    }
	}
    }

    RVector qvec (n);
    Point loc = elist[el]->Local (mesh.nlist, cnt);
    RVector fun = elist[el]->LocalShapeF (loc);
    for (i = 0; i < elist[el]->nNode(); i++) {
	nd = elist[el]->Node[i];
	if (nlist[nd].isBnd() || mode != SRCMODE_NEUMANN)
	    qvec[nd] = fun[i];
    }
    return qvec;
}

/* Computes the source vectors for all the boundary sources when a QM file has been specified
*/
void genmat_toastsource(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh, const RCompRowMatrix qvec, const int ns, RVector* &dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm)
{
   int sysdim = mesh.nlen();       
   int fullsysdim = sum(node_angN);  
   if(srctp ==2)
	cerr<<" Not implemented for uncollided sources. Try using point sources instead"<<endl;
   RCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     Source[i].New(fullsysdim,1);
     genmat_toastsourcevalvector(sphOrder, node_angN, offset, Svec, mesh, qvec,i, dirVec[i], srctp,  pts, wts, Ylm);
     b1.AB(Svec, Source[i]);
   }
}

/* Computes the source vector per a boundary source when a QM file has been specified
*/
void genmat_toastsourcevalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const RCompRowMatrix qvec, const int iq, RVector& dirVec, const int srctp,  const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
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
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);

   int row_offset = 0;
   RDenseMatrix Ystarlm;
   RDenseMatrix dirMat(1, 3);

   for (int jq = qvec.rowptr[iq]; jq < qvec.rowptr[iq+1]; jq++) {
   	int Nsource = qvec.colidx[jq];
   	double sweight = norm(qvec.Get(iq,Nsource)); 
   	for (el = 0; el < mesh.elen(); el++) {
        	if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        	if(!mesh.elist[el]->HasBoundarySide()) continue;
		for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  		// if sd is not a boundary side. skip 
	  		if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
	  		RVector nhat = mesh.ElDirectionCosine(el,sd);
			if(length(dirVec) == 0){
				for(int i=0; i < mesh.Dimension(); i++)
				{
					dirMat(0, i) = -nhat[i];
					dirVec[i] = -nhat[i];
				}
			}
			else{
			    for(int i=0; i < mesh.Dimension(); i++)
					dirMat(0, i) = dirVec[i];
			}
	  		for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    			is = mesh.elist[el]->SideNode(sd,nd);
	    			js = mesh.elist[el]->Node[is];
				double ela_i =  mesh.elist[el]->IntF(is);
	    			if(srctp == 1)
				{
					RDenseMatrix Angsvec(node_angN[js], 1);
					for(int l=0; l<= sphOrder[js]; l++){
						int indl = l*l;
						for(int m=-l; m<=l; m++){
							int indm = l + m;
							for(int i=0; i < pts.nRows(); i++)
							{
								double sdotn = 0;
								for(int j=0; j < mesh.Dimension(); j++)
									sdotn += nhat[j]*pts.Get(i, j);
								if(sdotn<0)
									Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(-sdotn)*wts[i]*4*M_PI;
							}
						}
					}
					for(int i = 0; i < node_angN[js]; i++)
						Svec(offset[js] + i, 0) += Angsvec(i, 0)*sweight*ela_i;

				}
				else if(srctp == 0){
					RDenseMatrix Angsvec(node_angN[js], 1);
					RDenseMatrix *Ystarlm;
					Ystarlm = new RDenseMatrix[sphOrder[js]+1];
					for(int l=0; l<= sphOrder[js]; l++)
						Ystarlm[l].New(2*l+1, 1);
					sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
					for(int l=0; l<= sphOrder[js]; l++){
						int indl = l*l;
						for(int m=-l; m<=l; m++){
							int indm = l + m;
							Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0)*(4*M_PI)/(2*l+1);
						}
					}
					for(int i = 0; i < node_angN[js]; i++)
						Svec(offset[js] + i, 0) += Angsvec(i, 0)*sweight*ela_i;

					delete []Ystarlm;
				}
				else{}
	    		}

		} // end loop on element sides

   	} // end loop on elements
  }
   delete []rowptr;
   delete []colidx;
}  

/** Computes source vectors for point sources on the boundary 
**/
void genmat_source(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh,  const IVector& Nsource, const int ns, RVector* &dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm)
{
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);       // full size of angles X space nodes

   RCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     Source[i].New(fullsysdim,1);
     genmat_sourcevalvector(sphOrder, node_angN, offset, Svec, mesh, Nsource[i], dirVec[i], srctp, pts, wts, Ylm);
     b1.AB(Svec, Source[i]);
   }
}

/**Computes source vector for a single point source on the boundary
**/
void genmat_sourcevalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const int Nsource, RVector& dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
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
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);

   int row_offset = 0;
   RDenseMatrix dirMat(1, 3);
   RDenseMatrix Ystarlm;

   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        if(!mesh.elist[el]->HasBoundarySide()) continue;
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  if(length(dirVec) == 0){
		for(int i=0; i < mesh.Dimension(); i++)
		{
			dirMat(0, i) = -nhat[i];
			dirVec[i] = -nhat[i];
		}
	  }
	  else{
		for(int i=0; i < mesh.Dimension(); i++)
			dirMat(0, i) = dirVec[i];
	  }
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    double ela_i =  mesh.elist[el]->IntF(is);
	    if(srctp == 1)
	    {
			RDenseMatrix Angsvec(node_angN[js], 1);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = l*l;
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					for(int i=0; i < pts.nRows(); i++)
					{
						double sdotn = 0;
						for(int j=0; j < mesh.Dimension(); j++)
							sdotn += nhat[j]*pts.Get(i, j);
						if(sdotn<0)
							Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(-sdotn)*wts[i]*4*M_PI;
					}
				}
			}
			for(int i = 0; i < node_angN[js]; i++)
				Svec(offset[js] + i, 0) = Angsvec(i, 0)*ela_i;
			
		}
	    	else{
			RDenseMatrix Angsvec(node_angN[js], 1);
			RDenseMatrix *Ystarlm;
			Ystarlm = new RDenseMatrix[sphOrder[js]+1];
			for(int l=0; l<= sphOrder[js]; l++)
			 Ystarlm[l].New(2*l+1, 1);
			sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = l*l;
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0)*(4*M_PI)/(2*l+1);
				}
			}
			for(int i = 0; i < node_angN[js]; i++)
				Svec(offset[js] + i, 0) = Angsvec(i, 0)*ela_i;
			delete []Ystarlm;
		}
	}
    } 


   } // end loop on elements
   delete []rowptr;
   delete []colidx;
}

void genmat_intsource(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh,  const IVector& Nsource, const int ns, RVector* &dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm, const RVector &delta)
{
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);       // full size of angles X space nodes

   for (int i = 0; i < ns; i++) {
     Source[i].New(fullsysdim,1);
     genmat_intsourcevalvector(sphOrder, node_angN, offset, Source[i], mesh, Nsource[i], dirVec[i], srctp, pts, wts, Ylm, delta);
   }
}

/**Computes source vector for a single point source in the interior of the domain
**/
void genmat_intsourcevalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const int Nsource, RVector& dirVec, const int srctp, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm, const RVector &delta)
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
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);
 
   int row_offset = 0;
   RDenseMatrix dirMat(1, 3);
   RDenseMatrix Ystarlm;

   // now create vector, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  if(length(dirVec) == 0){
		for(int i=0; i < mesh.Dimension(); i++)
		{
			dirMat(0, i) = -nhat[i];
			dirVec[i] = -nhat[i];
		}
	  }
	  else{
		for(int i=0; i < mesh.Dimension(); i++)
			dirMat(0, i) = dirVec[i];
	  }
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    if(mesh.nlist[js].isBnd() || js != Nsource) continue;
	    double ela_i =  mesh.elist[el]->IntF(is);
	    RVector p(mesh.nlen());
	    for(int i=0; i < mesh.nlen(); i++) p[i] = 1.0;

	    if(srctp == 1) {
		// Adding the CG components to the source vector
		Svec(offset[js], 0) += ela_i*sqrt(4*M_PI);
		//Adding SDM components to the source vector
		for(int l=0; l <= sphOrder[js]; l++){
			int indl = l*l;
			for(int m = -l; m <= l; m++){
				int indm = l + m;
				Svec(offset[js] + indl + indm, 0) += (intSinCosY(l, m)*mesh.elist[el]->IntPd(p, is, 0) + intSinSinY(l, m)*mesh.elist[el]->IntPd(p, is, 1))*delta[el];
				if(mesh.Dimension() == 3)
					Svec(offset[js] + indl + indm, 0) += intCosY(l, m)*mesh.elist[el]->IntPd(p, is, 2)*delta[el];

				 
			}
		}
	    }
	    else{
		RDenseMatrix Angsvec(node_angN[js], 1);
		RDenseMatrix *Ystarlm;
		Ystarlm = new RDenseMatrix[sphOrder[js]+1];
		for(int l=0; l<= sphOrder[js]; l++)
		 	Ystarlm[l].New(2*l+1, 1);
		sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
		for(int l=0; l<= sphOrder[js]; l++){
			int indl = l*l;
			for(int m=-l; m<=l; m++){
				int indm = l + m;
				Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0)*(4*M_PI)/(2*l+1);
			}
		}
		//Adding CG components to the source vector
		for(int i = 0; i < node_angN[js]; i++)
			Svec(offset[js] + i, 0) += Angsvec(i, 0)*ela_i;
		//Adding SDM components to the source vector
		for(int l=0; l <= sphOrder[js]; l++){
			int indl = l*l;
			for(int m = -l; m <= l; m++){
				int indm = l + m;
				Svec(offset[js] + indl + indm, 0) += delta[el]*Angsvec(indl+indm, 0)*(dirMat(0, 0)*mesh.elist[el]->IntPd(p, is, 0) + dirMat(0, 1)*mesh.elist[el]->IntPd(p, is, 1));
				if(mesh.Dimension() == 3)
					Svec(offset[js] + indl + indm, 0) += delta[el]*Angsvec(indl+indm, 0)*dirMat(0, 2)*mesh.elist[el]->IntPd(p, is, 2); 	
			}
		}

		delete []Ystarlm;
		}
	    }
	}
   } // end loop on elements
   delete []rowptr;
   delete []colidx;

}
void genmat_intsourceuncollided(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* &Source, Mesh& mesh, const IVector& Nsource, const int ns, RVector* &dirVec, const RDenseMatrix& pts, const RVector& wts, const RVector &mua, const RVector &mus, const RVector &delta, const double g)
{
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);       // full size of angles X space nodes

   for (int i = 0; i < ns; i++) {
     Source[i].New(fullsysdim,1);
     genmat_intsourceuncollidedvalvector(sphOrder, node_angN, offset, Source[i], mesh, Nsource[i], dirVec[i], pts, wts, mua, mus, delta, g);
   }

}
void genmat_intsourceuncollidedvalvector(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, Mesh& mesh, const int Nsource, RVector& dirVec, const RDenseMatrix& pts, const RVector& wts, const RVector &mua, const RVector &mus, const RVector &delta, const double g)
{
   xASSERT(length(dirVec) != 0, "dirVec should be defined for uncollided sources");  
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   NodeList &nlist=mesh.nlist;
   int *rowptr, *colidx;
   rowptr = new int[fullsysdim+1];
   colidx = new int[fullsysdim];
   rowptr[0] = 0;
   for(int i=1;  i <= fullsysdim; i++)
	rowptr[i] = i;
   for(int i=0; i<fullsysdim; i++)
	colidx[i] = 0;
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);
 
   int row_offset = 0;
   double t = 100*mesh.Size();
   cout<<"length of longest bounding box size: "<<t<<endl;
   int dim = mesh.Dimension();
   RDenseMatrix dirMat(1, 3);
   Point p1(dim), p2(dim);
   p1 = nlist[Nsource];
   for(int i=0; i < dim; i++)
   {
	p2[i] = p1[i] + t*dirVec[i];
	dirMat(0, i) = dirVec[i];
   }
   Point p3 = mesh.BndIntersect(p1, p2);
   cout<<"The intersection point on the mesh boundary along this direction is: "<<p2<<endl;

   double len = p1.Dist(p3);
   int num_intervals = 100;
   RDenseMatrix Angsvec(vmax(node_angN), 1);
   RDenseMatrix *Ystarlm;
   Ystarlm = new RDenseMatrix[vmax(sphOrder)+1];
   for(int l=0; l<= vmax(sphOrder); l++)
   	Ystarlm[l].New(2*l+1, 1);
   sphericalHarmonics(vmax(sphOrder), 1, dirMat, Ystarlm);
   for(int l=0; l<= vmax(sphOrder); l++){
	int indl = l*l;
	for(int m=-l; m<=l; m++){
		int indm = l + m;
		Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0)*(4*M_PI)/(2*l+1);
	}
   }
   for(int i=0 ; i< num_intervals-1; i++)
   {
	Point p(dim);
	p = p1 + i*len*dirVec/num_intervals;
	for(int el=0; el < mesh.elen(); el++)
	{
	    if(!mesh.elist[el]->GContains(p, nlist)) continue;
	    RVector fn(mesh.nlen());
	    for(int j=0; j < mesh.nlen(); j++) fn[j] = 1.0;
 
	    for(int is=0; is < mesh.elist[el]->nNode(); is++)
	    {
		int nd = mesh.elist[el]->Node[is];
		double ela_i =  mesh.elist[el]->IntF(is);
		//adding CG components to the source vector
		for(int js = 0; js < node_angN[nd]; js++)
			Svec(offset[nd] + js, 0) += mus[el]*pow(g, floor(sqrt(js)))*Angsvec(js, 0)*exp(-i*len*(mua[el]+mus[el])/num_intervals)*ela_i;

		//adding SDM components to the source vector
		for(int l=0; l <= sphOrder[nd]; l++){
			int indl = l*l;
			for(int m = -l; m <= l; m++){
				int indm = l + m;
				Svec(offset[nd] + indl + indm, 0) += mus[el]*delta[el]*exp(-i*len*(mua[el]+mus[el])/num_intervals)*pow(g,(double)l)*Angsvec(indl+indm, 0)*(intSinCosY(l, m)*mesh.elist[el]->IntPd(fn, is, 0) + intSinSinY(l, m)*mesh.elist[el]->IntPd(fn, is, 1));
				if(mesh.Dimension() == 3)
					Svec(offset[nd] + indl + indm, 0) += mus[el]*delta[el]*exp(-i*len*(mua[el]+mus[el])/num_intervals)*pow(g, (double)l)*Angsvec(indl+indm, 0)*intCosY(l, m)*mesh.elist[el]->IntPd(fn, is, 2); 	
			}
		}

	    }	
		
	}
   } 
   delete []Ystarlm;
   delete []rowptr;
   delete []colidx;
     
}


