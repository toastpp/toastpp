// Implementation of those routines in Raster_CPixel specific to tetrahedral
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tet_qr.h"
#include "tetsplit.h"

RCompRowMatrix *Raster_CPixel::CreateBuv_tet4 () const
{
    int i, j, k, r, m, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, idx_i, idx_j;
    int imin, imax, jmin, jmax, kmin, kmax;
    RVector fun;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/bdim[0];
    double dy = yrange/bdim[1];
    double dz = zrange/bdim[2];
    double djac, v;
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
        std::cerr << "pass 1, " << el << "/" << nel << std::endl;
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
	imin = max (0, (int)floor(bdim[0] * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-1, (int)floor(bdim[0] * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor(bdim[1] * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-1, (int)floor(bdim[1] * (eymax-bbmin[1])/yrange));
	kmin = max (0, (int)floor(bdim[2] * (ezmin-bbmin[2])/zrange));
	kmax = min (bdim[2]-1, (int)floor(bdim[2] * (ezmax-bbmin[2])/zrange));
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
	int nentry = (nimax[r]-nimin[r]+1)*(njmax[r]-njmin[r]+1)*
	    (nkmax[r]-nkmin[r]+1);
	rowptr[r+1] = rowptr[r]+nentry;
    }
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (r = m = 0; r < n; r++) {
	for (k = nkmin[r]; k <= nkmax[r]; k++) {
	    for (j = njmin[r]; j <= njmax[r]; j++) {
		for (i = nimin[r]; i <= nimax[r]; i++) {
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

    // pass 2: fill the matrix
    for (el = 0; el < nel; el++) {
        std::cerr << "pass 2, " << el << "/" << nel << std::endl;
	Element *pel = meshptr->elist[el];
	xASSERT(pel->Type() == ELID_TET4,
		"Currently only implemented for 4-noded tetrahedra");

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
	imin = max (0, (int)floor(bdim[0] * (exmin-bbmin[0])/xrange));
	imax = min (bdim[0]-1, (int)floor(bdim[0] * (exmax-bbmin[0])/xrange));
	jmin = max (0, (int)floor(bdim[1] * (eymin-bbmin[1])/yrange));
	jmax = min (bdim[1]-1, (int)floor(bdim[1] * (eymax-bbmin[1])/yrange));
	kmin = max (0, (int)floor(bdim[2] * (ezmin-bbmin[2])/zrange));
	kmax = min (bdim[2]-1, (int)floor(bdim[2] * (ezmax-bbmin[2])/zrange));

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
				submesh.elist[m]->SetRegion(-1);
			        // mark disabled
			Tetsplit_cube (&submesh, xmin, xmax, ymin, ymax,
				       zmin, zmax, reg);
		    }
		}
	    }
	}
#endif

	// now perform quadrature on all sub-elements
	submesh.SubSetup();
	for (i = 0; i < submesh.elen_used; i++) {
	    Element *psel = submesh.elist[i];

	    // find the voxel cell containing the sub-tet
	    int reg = psel->Region();
	    if (reg < 0) continue; // disabled
	    int ci = imin + reg & 0x3FF;
	    int cj = jmin + (reg >> 10) & 0x3FF;
	    int ck = kmin + (reg >> 20) & 0x3FF;
	    idx_j = ci + cj*stride_i + ck*stride_j;
	    // map quadrature points into global frame
	    for (m = 0; m < np; m++) {
		djac = psel->DetJ(absc[m], &submesh.nlist);
		Point glob = psel->Global (submesh.nlist, absc[m]);
		Point loc = pel->Local(meshptr->nlist, glob);
		fun = pel->LocalShapeF (loc);
		v = wght[m] * djac;
		for (ii = 0; ii < fun.Dim(); ii++) {
		    idx_i = pel->Node[ii];
		    (*Buv)(idx_i, idx_j) += v*fun[ii];
		}
	    }
	}
    }

    Buv->Shrink();
    return Buv;
}
