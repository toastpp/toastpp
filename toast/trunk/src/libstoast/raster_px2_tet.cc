// Implementation of those routines in Raster_Pixel2 specific to tetrahedral
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tet_qr.h"
#include "tetsplit.h"

// ==========================================================================
// Create the mixed-basis mass matrix Buv by subdividing tetrahedral elements
// at voxel boundaries

RCompRowMatrix *Raster_Pixel2::CreateMixedMassmat_tet4 () const
{
    int i, j, k, r, m, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, jj, idx_i, idx_j;
    int imin, imax, jmin, jmax, kmin, kmax;
    RVector fun;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double dz = zrange/(bdim[2]-1.0);
    double xmin, xmax, ymin, ymax, zmin, zmax, djac, vb, v;
    int stride_i = bdim[0], stride_j = stride_i*bdim[1];

    // quadrature rule for local tetrahedron
    const double *wght;
    const Point *absc;
    int np = QRule_tet_4_14 (&wght, &absc);

    // pass 1: determine matrix fill structure
    int *nimin = new int[n];
    int *nimax = new int[n];
    int *njmin = new int[n];
    int *njmax = new int[n];
    int *nkmin = new int[n];
    int *nkmax = new int[n];
    for (i = 0; i < n; i++) {
	nimin[i] = bdim[0];
	njmin[i] = bdim[1];
	nkmin[i] = bdim[2];
	nimax[i] = njmax[i] = nkmax[i] = -1;
    }
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	// element bounding box
	double exmin = meshptr->nlist[pel->Node[0]][0];
	double exmax = meshptr->nlist[pel->Node[0]][0];
	double eymin = meshptr->nlist[pel->Node[0]][1];
	double eymax = meshptr->nlist[pel->Node[0]][1];
	double ezmin = meshptr->nlist[pel->Node[0]][2];
	double ezmax = meshptr->nlist[pel->Node[0]][2];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	    ezmin = min (ezmin, meshptr->nlist[pel->Node[j]][2]);
	    ezmax = max (ezmax, meshptr->nlist[pel->Node[j]][2]);
	}
	// determine which pixels overlap the element
	imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));
	kmin = max (0, (int)floor((bdim[2]-1) * (ezmin-bbmin[2])/zrange));
	kmax = min (bdim[2]-2, (int)floor((bdim[2]-1) * (ezmax-bbmin[2])/zrange));
	for (i = 0; i < pel->nNode(); i++) {
	    int nidx = pel->Node[i];
	    if (imin < nimin[nidx]) nimin[nidx] = imin;
	    if (imax > nimax[nidx]) nimax[nidx] = imax;
	    if (jmin < njmin[nidx]) njmin[nidx] = jmin;
	    if (jmax > njmax[nidx]) njmax[nidx] = jmax;
	    if (kmin < nkmin[nidx]) nkmin[nidx] = kmin;
	    if (kmax > nkmax[nidx]) nkmax[nidx] = kmax;
	}
    }

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (r = 0; r < n; r++) {
	int nentry = (nimax[r]-nimin[r]+2)*(njmax[r]-njmin[r]+2)*(nkmax[r]-nkmin[r]+2);
	rowptr[r+1] = rowptr[r]+nentry;
    }
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (r = m = 0; r < n; r++) {
	for (k = nkmin[r]; k <= nkmax[r]+1; k++) {
	    for (j = njmin[r]; j <= njmax[r]+1; j++) {
		for (i = nimin[r]; i <= nimax[r]+1; i++) {
		    colidx[m++] = i + j*bdim[0] + k*bdim[0]*bdim[1];
		}
	    }
	}
    }

    RCompRowMatrix *Buv = new RCompRowMatrix (n, blen, rowptr, colidx);
    delete []rowptr;
    delete []colidx;
    delete []nimin;
    delete []nimax;
    delete []njmin;
    delete []njmax;
    delete []nkmin;
    delete []nkmax;

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    double t_split = 0.0, t_integrate = 0.0;

    // pass 2: fill the matrix
    for (el = 0; el < nel; el++) {
	std::cout << "Buv: processing el " << el << " of " << nel << std::endl;
	Element *pel = meshptr->elist[el];
	xASSERT(pel->Type() == ELID_TET4, "Currently only implemented for 4-noded tetrahedra");

	double orig_size = pel->Size();

	// element bounding box
	double exmin = meshptr->nlist[pel->Node[0]][0];
	double exmax = meshptr->nlist[pel->Node[0]][0];
	double eymin = meshptr->nlist[pel->Node[0]][1];
	double eymax = meshptr->nlist[pel->Node[0]][1];
	double ezmin = meshptr->nlist[pel->Node[0]][2];
	double ezmax = meshptr->nlist[pel->Node[0]][2];
	for (j = 1; j < pel->nNode(); j++) {
	    exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	    exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	    eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	    eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	    ezmin = min (ezmin, meshptr->nlist[pel->Node[j]][2]);
	    ezmax = max (ezmax, meshptr->nlist[pel->Node[j]][2]);
	}

	// determine which pixels overlap the element
	imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));
	kmin = max (0, (int)floor((bdim[2]-1) * (ezmin-bbmin[2])/zrange));
	kmax = min (bdim[2]-2, (int)floor((bdim[2]-1) * (ezmax-bbmin[2])/zrange));

	tic();

	// Create a mesh containing this single element
	for (i = 0; i < 4; i++) {
	    submesh.elist[0]->Node[i] = i;
	    submesh.nlist[i] = meshptr->nlist[pel->Node[i]];
	}
	submesh.elist[0]->SetRegion(0);
	submesh.elen_used = 1;
	submesh.nlen_used = 4;

	// DEBUG
	//submesh.SubSetup();
	//ofstream ofs1("dbg1.msh");
	//ofs1 << submesh << std::endl;
	//double orig_size = submesh.CalcFullSize();
	//std::cout << "Initial element size: " << orig_size << std::endl;

	// perform subdivisions along x-cutting planes
	for (i = imin; i < imax; i++) {
	    int elen = submesh.elen_used;
	    double cut_pos = bbmin[0] + (i+1)*dx;
	    for (m = 0; m < elen; m++) {
		int subreg = submesh.elist[m]->Region() & 0x3FF;
		if (subreg >= i-imin)
		    Tetsplit (&submesh, m, 0, cut_pos);
	    }
	}
	    
#ifdef UNDEF
	// DEBUG: compare integral over original element with sum
	// over sub elements
	submesh.SubSetup();
	double sum = 0.0, sum_orig = 0.0;
	for (m = 0; m < np; m++) {
	    djac = pel->DetJ(absc[m], &meshptr->nlist);
	    fun = pel->LocalShapeF (absc[m]);
	    v = wght[m] * djac;
	    for (ii = 0; ii < fun.Dim(); ii++) {
		sum_orig += v*fun[ii];
	    }
	}

	for (i = 0; i < submesh.elen_used; i++) {
	    Element *psel = submesh.elist[i];
	    for (m = 0; m < np; m++) {
		djac = psel->DetJ(absc[m], &submesh.nlist);
		Point glob = psel->Global (submesh.nlist, absc[m]);
		Point loc = pel->Local(meshptr->nlist, glob);
		fun = pel->LocalShapeF (loc);
		v = wght[m] * djac;
		for (ii = 0; ii < fun.Dim(); ii++) {
		    sum += v*fun[ii];
		}
	    }
	}
#endif

	// DEBUG
	//submesh.SubSetup();
	//ofstream ofs2("dbg2.msh");
	//ofs2 << submesh << std::endl;
	//double sz = 0.0;
	//std::cout << "Subdiv element sizes: " << std::endl;
	//for (i = 0; i < submesh.elen_used; i++) {
	//    double s = submesh.elist[i]->Size();
	//    std::cout << s << std::endl;
	//    sz += s;
	//    if (s <= 0)
	//	std::cerr << "**** x-split: El size(" << el << ',' << i
	//		  << ")=" << s << std::endl;
	//}
	//std::cout << "Total: " << sz << std::endl;

	// perform subdivisions along y-cutting planes
	for (j = jmin; j < jmax; j++) {
	    int elen = submesh.elen_used;
	    double cut_pos = bbmin[1] + (j+1)*dy;
	    for (m = 0; m < elen; m++) {
		int subreg = (submesh.elist[m]->Region() >> 10) & 0x3FF;
		if (subreg >= j-jmin)
		    Tetsplit (&submesh, m, 1, cut_pos);
	    }
	}
	    
	// DEBUG
	//submesh.SubSetup();
	//ofstream ofs3("dbg3.msh");
	//ofs3 << submesh << std::endl;
	//sz = 0.0;
	//std::cout << "Subdiv element sizes: " << std::endl;
	//for (i = 0; i < submesh.elen_used; i++) {
	//    double s = submesh.elist[i]->Size();
	//    std::cout << s << std::endl;
	//    sz += s;
	//    if (s <= 0)
	//	std::cerr << "**** y-split: El size(" << el << ',' << i
	//		  << ")=" << s << std::endl;
	//}
	//std::cout << "Total: " << sz << std::endl;

	// perform subdivisions along z-cutting planes
	for (k = kmin; k < kmax; k++) {
	    int elen = submesh.elen_used;
	    double cut_pos = bbmin[2] + (k+1)*dz;
	    for (m = 0; m < elen; m++) {
		int subreg = (submesh.elist[m]->Region() >> 20) & 0x3FF;
		if (subreg >= k-kmin)
		    Tetsplit (&submesh, m, 2, cut_pos);
	    }
	}

	// DEBUG
	//submesh.SubSetup();
	//ofstream ofs4("dbg4.msh");
	//ofs4 << submesh << std::endl;
	//sz = 0.0;
	//std::cout << "Subdiv element sizes: " << std::endl;
	//for (i = 0; i < submesh.elen_used; i++) {
	//    double s = submesh.elist[i]->Size();
	//    std::cout << s << std::endl;
	//    sz += s;
	//    if (s <= 0)
	//	std::cerr << "**** z-split: El size(" << el << ',' << i
	//		  << ")=" << s << std::endl;
	//}
	//std::cout << "Total: " << sz << std::endl;

	t_split += toc();
	tic();

	// for any voxel cells completely inside the element,
	// replace the subdivided elements by a simple 6-tetra
	// subdivision to reduce the number of elements for the
	// integral
	for (k = kmin; k <= kmax; k++) {
	    zmax = (zmin = bbmin[2] + dz*k) + dz;
	    for (j = jmin; j <= jmax; j++) {
		ymax = (ymin = bbmin[1] + dy*j) + dy;
		for (i = imin; i <= imax; i++) {
		    xmax = (xmin = bbmin[0] + dx*i) + dx;
		    if (pel->GContains(Point3D(xmin,ymin,zmin),false) &&
			pel->GContains(Point3D(xmax,ymin,zmin),false) &&
			pel->GContains(Point3D(xmin,ymax,zmin),false) &&
			pel->GContains(Point3D(xmax,ymax,zmin),false) &&
			pel->GContains(Point3D(xmin,ymin,zmax),false) &&
			pel->GContains(Point3D(xmax,ymin,zmax),false) &&
			pel->GContains(Point3D(xmin,ymax,zmax),false) &&
			pel->GContains(Point3D(xmax,ymax,zmax),false)) {
			int reg = (i-imin) + (j-jmin) << 10 + (k-kmin) << 20;
			for (m = 0; m < submesh.elen_used; m++)
			    if (submesh.elist[m]->Region() == reg)
				submesh.elist[m]->SetRegion(-1); // mark disabled
			Tetsplit_cube (&submesh, xmin, xmax, ymin, ymax,
				       zmin, zmax, reg);
		    }
		}
	    }
	}


	// now perform quadrature on all sub-elements
	submesh.SubSetup();
	//double sz = 0.0; // DEBUG
	for (i = 0; i < submesh.elen_used; i++) {
	    Element *psel = submesh.elist[i];

	    // find the voxel cell containing the sub-tet
	    int reg = psel->Region();
	    if (reg < 0) continue; // disabled
	    int ci = imin + reg & 0x3FF;
	    int cj = jmin + (reg >> 10) & 0x3FF;
	    int ck = kmin + (reg >> 20) & 0x3FF;
	    int cellidx = ci + (cj + ck*bdim[1])*bdim[0];
	    // map quadrature points into global frame
	    for (m = 0; m < np; m++) {
		djac = psel->DetJ(absc[m], &submesh.nlist);
		Point glob = psel->Global (submesh.nlist, absc[m]);
		Point loc = pel->Local(meshptr->nlist, glob);
		fun = pel->LocalShapeF (loc);
		v = wght[m] * djac;
		// loop over the vertices of the voxel cell
		for (jj = 0; jj < 8; jj++) {
		    idx_j = cellidx + (jj&1) + ((jj>>1)&1)*stride_i
		    	+ (jj>>2)*stride_j;
		    //idx_j = cellidx + jj%2 + ((jj/2)%2)*bdim[0] +
		    //	(jj/4)*bdim[0]*bdim[1];
		    vb = v*Value_nomask (glob, idx_j, false);
		    for (ii = 0; ii < fun.Dim(); ii++) {
			idx_i = pel->Node[ii];
			(*Buv)(idx_i, idx_j) += vb*fun[ii];
		    }
		}
	    }
	    //sz += psel->Size();
	}
	t_integrate += toc();
    }

    std::cout << "#### t_split = " << t_split << std::endl;
    std::cout << "#### t_integrate = " << t_integrate << std::endl;

    Buv->Shrink();
    return Buv;
}
// ==========================================================================
// Assemble single-element contribution for element "el" into global
// system matrix M, where coefficients (where applicable) are given in pixel
// basis. "mode" is the integration type.

void Raster_Pixel2::AddToElMatrix_tet (int el, RGenericSparseMatrix &M,
    const RVector *pxcoeff, int mode) const
{
    // All "parameter-free" integrals can be done in the standard way
    if (mode != ASSEMBLE_PFF && mode != ASSEMBLE_PDD) {
	::AddToElMatrix (*meshptr, el, M, pxcoeff, mode);
	return;
    }

    dASSERT(pxcoeff, "AddToElMatrix: requires a parameter pointer");

    Element *pel = meshptr->elist[el];
    int nnode = pel->nNode();
    int i, j, k, m, ii, jj, kk, idx_i, idx_j, idx_k, nv;
    int imin, imax, jmin, jmax, kmin, kmax;
    double xmin, xmax, ymin, ymax, zmin, zmax, b, v;
    RVector fun;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double dz = zrange/(bdim[2]-1.0);
    int stride_i = bdim[0], stride_j = stride_i*bdim[1];

    // quadrature rule for local tetrahedron
    const double *wght;
    const Point *absc;
    int np = QRule_tet_4_14 (&wght, &absc);

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    // SPLIT ELEMENT

    // compute bounding box
    double exmin = meshptr->nlist[pel->Node[0]][0];
    double exmax = meshptr->nlist[pel->Node[0]][0];
    double eymin = meshptr->nlist[pel->Node[0]][1];
    double eymax = meshptr->nlist[pel->Node[0]][1];
    double ezmin = meshptr->nlist[pel->Node[0]][2];
    double ezmax = meshptr->nlist[pel->Node[0]][2];
    for (j = 1; j < pel->nNode(); j++) {
	exmin = min (exmin, meshptr->nlist[pel->Node[j]][0]);
	exmax = max (exmax, meshptr->nlist[pel->Node[j]][0]);
	eymin = min (eymin, meshptr->nlist[pel->Node[j]][1]);
	eymax = max (eymax, meshptr->nlist[pel->Node[j]][1]);
	ezmin = min (ezmin, meshptr->nlist[pel->Node[j]][2]);
	ezmax = max (ezmax, meshptr->nlist[pel->Node[j]][2]);
    }
    // determine which pixels overlap the element
    imin = max (0, (int)floor((bdim[0]-1) * (exmin-bbmin[0])/xrange));
    imax = min (bdim[0]-2, (int)floor((bdim[0]-1) * (exmax-bbmin[0])/xrange));
    jmin = max (0, (int)floor((bdim[1]-1) * (eymin-bbmin[1])/yrange));
    jmax = min (bdim[1]-2, (int)floor((bdim[1]-1) * (eymax-bbmin[1])/yrange));
    kmin = max (0, (int)floor((bdim[2]-1) * (ezmin-bbmin[2])/zrange));
    kmax = min (bdim[2]-2, (int)floor((bdim[2]-1) * (ezmax-bbmin[2])/zrange));
    
    // Create a mesh containing this single element
    for (i = 0; i < 4; i++) {
	submesh.elist[0]->Node[i] = i;
	submesh.nlist[i] = meshptr->nlist[pel->Node[i]];
    }
    submesh.elist[0]->SetRegion(0);
    submesh.elen_used = 1;
    submesh.nlen_used = 4;
    
    // perform subdivisions along x-cutting planes
    for (i = imin; i < imax; i++) {
	int elen = submesh.elen_used;
	double cut_pos = bbmin[0] + (i+1)*dx;
	for (m = 0; m < elen; m++) {
	    int subreg = submesh.elist[m]->Region() & 0x3FF;
	    if (subreg >= i-imin)
		Tetsplit (&submesh, m, 0, cut_pos);
	}
    }
    // perform subdivisions along y-cutting planes
    for (j = jmin; j < jmax; j++) {
	int elen = submesh.elen_used;
	double cut_pos = bbmin[1] + (j+1)*dy;
	for (m = 0; m < elen; m++) {
	    int subreg = (submesh.elist[m]->Region() >> 10) & 0x3FF;
	    if (subreg >= j-jmin)
		Tetsplit (&submesh, m, 1, cut_pos);
	}
    }
    // perform subdivisions along z-cutting planes
    for (k = kmin; k < kmax; k++) {
	int elen = submesh.elen_used;
	double cut_pos = bbmin[2] + (k+1)*dz;
	for (m = 0; m < elen; m++) {
	    int subreg = (submesh.elist[m]->Region() >> 20) & 0x3FF;
	    if (subreg >= k-kmin)
		Tetsplit (&submesh, m, 2, cut_pos);
	}
    }

#ifdef UNDEF
    // for any voxel cells completely inside the element,
    // replace the subdivided elements by a simple 6-tetra
    // subdivision to reduce the number of elements for the
    // integral
    for (k = kmin; k <= kmax; k++) {
	zmax = (zmin = bbmin[2] + dz*k) + dz;
	for (j = jmin; j <= jmax; j++) {
	    ymax = (ymin = bbmin[1] + dy*j) + dy;
	    for (i = imin; i <= imax; i++) {
		xmax = (xmin = bbmin[0] + dx*i) + dx;
		if (pel->GContains(Point3D(xmin,ymin,zmin),false) &&
		    pel->GContains(Point3D(xmax,ymin,zmin),false) &&
		    pel->GContains(Point3D(xmin,ymax,zmin),false) &&
		    pel->GContains(Point3D(xmax,ymax,zmin),false) &&
		    pel->GContains(Point3D(xmin,ymin,zmax),false) &&
		    pel->GContains(Point3D(xmax,ymin,zmax),false) &&
		    pel->GContains(Point3D(xmin,ymax,zmax),false) &&
		    pel->GContains(Point3D(xmax,ymax,zmax),false)) {
		    int reg = (i-imin) + (j-jmin) << 10 + (k-kmin) << 20;
		    for (m = 0; m < submesh.elen_used; m++)
			if (submesh.elist[m]->Region() == reg)
			    submesh.elist[m]->SetRegion(-1); // mark disabled
		    Tetsplit_cube (&submesh, xmin, xmax, ymin, ymax,
				   zmin, zmax, reg);
		}
	    }
	}
    }
#endif

    // now perform quadrature on all sub-elements
    RDenseMatrix M0(pel->nNode(), pel->nNode());
    submesh.SubSetup();
    for (i = 0; i < submesh.elen_used; i++) {
	Element *psel = submesh.elist[i];
	RDenseMatrix M1(pel->nNode(), pel->nNode());

	// find the voxel cell containing the sub-tet
	int reg = psel->Region();
	if (reg < 0) continue; // disabled
	int ci = imin + reg & 0x3FF;
	int cj = jmin + (reg >> 10) & 0x3FF;
	int ck = kmin + (reg >> 20) & 0x3FF;
	int cellidx = ci + (cj + ck*bdim[1])*bdim[0];

	switch (mode) {
	case ASSEMBLE_PFF:
	    for (m = 0; m < np; m++) {
		v = psel->DetJ(absc[m], &submesh.nlist) * wght[m];
		Point glob = psel->Global (submesh.nlist, absc[m]);
		Point loc = pel->Local(meshptr->nlist, glob);
		fun = pel->LocalShapeF (loc);
		for (kk = 0; kk < 8; kk++) {
		    idx_k = cellidx + (kk&1) + ((kk>>1)&1)*stride_i
			+ (kk>>2)*stride_j;
		    b = v * (*pxcoeff)[idx_k] *
			Value_nomask (glob, idx_k, false);
		    for (ii = 0; ii < fun.Dim(); ii++) {
			//idx_i = pel->Node[ii];
			for (jj = 0; jj < fun.Dim(); jj++) {
			    //idx_j = pel->Node[jj];
			    M1(ii,jj) += b * fun[ii] * fun[jj];
			    //M(idx_i,idx_j) += b * fun[ii] * fun[jj];
			}
		    }
		}
	    }
	    break;
	case ASSEMBLE_PDD: {
	    RDenseMatrix IJ(3,3), der(3,4);
	    double dd;
	    // assumes that the derivative is constant over the element
	    RVector vtxfun[4];
	    for (ii = 0; ii < 4; ii++)
		vtxfun[ii] = pel->GlobalShapeF(meshptr->nlist,
					       submesh.nlist[psel->Node[ii]]);
	    for (ii = 0; ii < 4; ii++) {
		der(0,ii) = vtxfun[1][ii] - vtxfun[0][ii];
		der(1,ii) = vtxfun[2][ii] - vtxfun[0][ii];
		der(2,ii) = vtxfun[3][ii] - vtxfun[0][ii];
	    }

	    for (m = 0; m < np; m++) {
		Point glob = psel->Global (submesh.nlist, absc[m]);
		Point loc = pel->Local(meshptr->nlist, glob);
		fun = pel->LocalShapeF (loc);
		v = psel->IJacobian (absc[m], &submesh.nlist, IJ) * wght[m];

#ifdef UNDEF
		// DEBUG
		static bool logout=true;
		if (logout && !v) {
		    std::cerr << "Problem in element " << el << ", subel "
			      << i << std::endl;
		    for (kk=0; kk < 4; kk++) {
			std::cerr << psel->Node[kk] << ": " 
				  << submesh.nlist[psel->Node[kk]] << std::endl;
		    }
		    logout = false;
		}
#endif

		for (kk = 0; kk < 8; kk++) {
		    idx_k = cellidx + (kk&1) + ((kk>>1)&1)*stride_i
			+ (kk>>2)*stride_j;
		    b = v * (*pxcoeff)[idx_k] *
			Value_nomask (glob, idx_k, false);
		    for (ii = 0; ii < fun.Dim(); ii++) {
			//idx_i = pel->Node[ii];
			for (jj = 0; jj < fun.Dim(); jj++) {
			    //idx_j = pel->Node[jj];
			    dd = 
				((IJ(0,0)*der(0,ii) + IJ(0,1)*der(1,ii) + IJ(0,2)*der(2,ii)) *
				 (IJ(0,0)*der(0,jj) + IJ(0,1)*der(1,jj) + IJ(0,2)*der(2,jj)) +
				 (IJ(1,0)*der(0,ii) + IJ(1,1)*der(1,ii) + IJ(1,2)*der(2,ii)) *
				 (IJ(1,0)*der(0,jj) + IJ(1,1)*der(1,jj) + IJ(1,2)*der(2,jj)) +
				 (IJ(2,0)*der(0,ii) + IJ(2,1)*der(1,ii) + IJ(2,2)*der(2,ii)) *
				 (IJ(2,0)*der(0,jj) + IJ(2,1)*der(1,jj) + IJ(2,2)*der(2,jj)));

			    M1(ii,jj) += b * dd;
			    //M(idx_i,idx_j) += b * dd;
			}
		    }
		}
	    }
	    } break;
	default:
	    xERROR("Integration type not supported");
	}
	for (ii = 0; ii < M0.nRows(); ii++) {
	    for (jj = 0; jj < M0.nCols(); jj++) {
	        M0(ii,jj) += M1(ii,jj);
	    }
        }
	if (!((i+1)%20) || i == submesh.elen_used-1) {
	    for (ii = 0; ii < M0.nRows(); ii++) {
		idx_i = pel->Node[ii];
		for (jj = 0; jj < M0.nCols(); jj++) {
		    idx_j = pel->Node[jj];
		    M(idx_i, idx_j) += M0(ii,jj);
		}
	    }
	    M0.Zero();
	}
    }
}
