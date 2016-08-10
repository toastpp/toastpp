// Implementation of those routines in Raster_Blob2 specific to tetrahedral
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_blob2.h"
#include "tet_qr.h"
#include "tetsplit.h"
#include "timing.h"

// ==========================================================================
// Local prototypes

void PlaneCut (BufMesh &mesh, const Point *pt, const RDenseMatrix &R);
void IcosaCut (BufMesh &mesh, double scale, const Point &cnt);

// ==========================================================================

#ifdef UNDEF
RCompRowMatrix *Raster_Blob2::CreateBvv_tet4 () const
{
    int i, i0, j0, k0, i1, j1, k1, tmp, nd, s, m, idx_i, idx_j;
    int el, nel = meshptr->elen();
    bool intersect0, intersect1;
    double px0, py0, pz0, px1, py1, pz1, rx, ry, rz, djac, v;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double dz = zrange/(bdim[2]-1.0);
    double radlimit2 = sup*sup;
    int stride1 = bdim[0];
    int stride2 = bdim[0]*bdim[1];

    // quadrature rule for local tetrahedron
    const double *wght;
    const Point *absc;
    int np = QRule_tet_6_24 (&wght, &absc);

    int *npx = new int[blen];
    for (i = 0; i < blen; i++) npx[i] = 0;

    // pass 1: determine memory requirements
    for (idx_i = 0; idx_i < blen; idx_i++) {
	std::cerr << "pass 1, " << idx_i << "/" << blen << std::endl;
        k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	for (idx_j = 0; idx_j < blen; idx_j++) {
	    if (idx_j == idx_i) {
	        npx[idx_i]++;
		continue;
	    }
	    k1 = idx_j / stride2;
	    tmp = idx_j - k1*stride2;
	    j1 = tmp / stride1;
	    i1 = tmp % stride1;
	    px1 = bbmin[0] + i1*dx;
	    py1 = bbmin[1] + j1*dy;
	    pz1 = bbmin[2] + k1*dz;
	    rx = px0-px1;
	    ry = py0-py1;
	    rz = pz0-pz1;
	    if (rx*rx + ry*ry + rz*rz >= radlimit2*4.0)
	        continue; // blobs don't overlap
	    for (el = 0; el < nel; el++) {
	        Element *pel = meshptr->elist[el];
		intersect0 = intersect1 = false;
		for (s = 0; s < pel->nNode() && !(intersect0 && intersect1); s++) {
		    nd = pel->Node[s];
		    if (!intersect0) {
		        rx = px0-meshptr->nlist[nd][0];
			ry = py0-meshptr->nlist[nd][1];
			rz = pz0-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2)
			  intersect0 = true;
		    }
		    if (!intersect1) {
		        rx = px1-meshptr->nlist[nd][0];
			ry = py1-meshptr->nlist[nd][1];
			rz = pz1-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2)
			  intersect1 = true;
		    }
		}
		if (intersect0 && intersect1) {
		    npx[idx_i]++;
		    break;
		}
	    }
	}
    }

    int *rowptr = new int[blen+1];
    rowptr[0] = 0;
    for (i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen];
    int *colidx = new int[nz];

    // pass 2: create the matrix
    for (idx_i = 0; idx_i < blen; idx_i++) {
	std::cerr << "pass 2, " << idx_i << "/" << blen << std::endl;
        k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	for (idx_j = 0; idx_j < blen; idx_j++) {
	    if (idx_j == idx_i) {
		colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
		continue;
	    }
	    k1 = idx_j / stride2;
	    tmp = idx_j - k1*stride2;
	    j1 = tmp / stride1;
	    i1 = tmp % stride1;
	    px1 = bbmin[0] + i1*dx;
	    py1 = bbmin[1] + j1*dy;
	    pz1 = bbmin[2] + k1*dz;
	    rx = px0-px1;
	    ry = py0-py1;
	    rz = pz0-pz1;
	    if (rx*rx + ry*ry + rz*rz >= radlimit2*4.0)
	        continue; // blobs don't overlap
	    for (el = 0; el < nel; el++) {
	        Element *pel = meshptr->elist[el];
		intersect0 = intersect1 = false;
		for (s = 0; s < pel->nNode() && !(intersect0 && intersect1); s++) {
		    nd = pel->Node[s];
		    if (!intersect0) {
		        rx = px0-meshptr->nlist[nd][0];
			ry = py0-meshptr->nlist[nd][1];
			rz = pz0-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2)
			  intersect0 = true;
		    }
		    if (!intersect1) {
		        rx = px1-meshptr->nlist[nd][0];
			ry = py1-meshptr->nlist[nd][1];
			rz = pz1-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2)
			  intersect1 = true;
		    }
		}
		if (intersect0 && intersect1) {
		    colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
		    break;
		}
	    }
	}
    }

    delete []npx;
    RCompRowMatrix *bvv = new RCompRowMatrix (blen, blen, rowptr, colidx);

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    for (el = 0; el < nel; el++) {
	std::cerr << "pass 3, " << el << "/" << nel << std::endl;
        Element *pel = meshptr->elist[el];
	for (idx_i = 0; idx_i < blen; idx_i++) {
	    k0 = idx_i / stride2;
	    tmp = idx_i - k0*stride2;
	    j0 = tmp / stride1;
	    i0 = tmp % stride1;
	    px0 = bbmin[0] + i0*dx;
	    py0 = bbmin[1] + j0*dy;
	    pz0 = bbmin[2] + k0*dz;
	    intersect0 = false;
	    for (s = 0; s < pel->nNode(); s++) {
	        nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx + ry*ry + rz*rz < radlimit2) {
		    intersect0 = true;
		    break;
		}
	    }
	    if (!intersect0) continue;
	    for (s = 0; s < 4; s++) {
	        submesh.elist[0]->Node[s] = s;
		submesh.nlist[s] = meshptr->nlist[pel->Node[s]];
	    }
	    for (idx_j = 0; idx_j <= idx_i; idx_j++) {
	        k1 = idx_j / stride2;
		tmp = idx_j - k1*stride2;
		j1 = tmp / stride1;
		i1 = tmp % stride1;
		px1 = bbmin[0] + i1*dx;
		py1 = bbmin[1] + j1*dy;
		pz1 = bbmin[2] + k1*dz;
		rx = px1-px0;
		ry = py1-py0;
		rz = pz1-pz0;
		if (rx*rx + ry*ry + rz*rz >= radlimit2*4.0)
		    continue;
		intersect1 = false;
		for (s = 0; s < pel->nNode(); s++) {
		    nd = pel->Node[s];
		    rx = px1-meshptr->nlist[nd][0];
		    ry = py1-meshptr->nlist[nd][1];
		    rz = pz1-meshptr->nlist[nd][2];
		    if (rx*rx + ry*ry + rz*rz < radlimit2) {
		        intersect1 = true;
			break;
		    }
		}
		if (!intersect1) continue;
		submesh.elen_used = 1;
		submesh.nlen_used = 4;
		submesh.elist[0]->SetRegion(0);
		Point cnt(3);
		cnt[0] = px0, cnt[1] = py0, cnt[2] = pz0;
		IcosaCut (submesh, sup, cnt);
		cnt[0] = px1, cnt[1] = py1, cnt[2] = pz1;
		IcosaCut (submesh, sup, cnt);
		submesh.SubSetup();

		for (s = 0; s < submesh.elen_used; s++) {
		    Element *psel = submesh.elist[s];
		    if (psel->Region()) continue;
		    for (m = 0; m < np; m++) {
		        djac = psel->DetJ(absc[m], &submesh.nlist);
			Point glob = psel->Global (submesh.nlist, absc[m]);
			v = wght[m] * djac *
			    Value_nomask(glob,idx_i,false) *
			    Value_nomask(glob,idx_j,false);
			(*bvv)(idx_i,idx_j) += v;
			if (idx_i != idx_j)
			    (*bvv)(idx_j,idx_i) += v;
		    }
		}
	    }
	}
    }

    delete []rowptr;
    delete []colidx;

    // diagonal conditioning
    RVector d = bvv->Diag();
    double dmean = mean(d);
    for (i = 0; i < blen; i++)
	(*bvv)(i,i) += dmean * dgscale;

    bvv->Shrink();
    std::cerr << "exit BasisMatrix" << std::endl;
    return bvv;
}
#endif

// ==========================================================================

// currently active
RCompRowMatrix *Raster_Blob2::CreateBvv_tet4 () const
{
    int i, j, i0, j0, k0, i1, j1, k1, s, m, el, nd, idx_i, idx_j, tmp;
    int nel = meshptr->elen();
    bool intersect, intersect0, intersect1;
    double px0, py0, pz0, px1, py1, pz1, rx, ry, rz, dstx, dsty, dstz, djac, v;
    Point glob(3);

    double dx = grid[0];
    double dy = grid[1];
    double dz = grid[2];
    double radlimit2 = sup*sup;
    int stride1 = bdim[0];
    int stride2 = bdim[0]*bdim[1];

    // quadrature rule for local tet
    const double *wght;
    const Point *absc;
    int np = QRule_tet_6_24 (&wght, &absc);

    // pre-compute local shape functions at quadrature points
    // for faster access from the inner loops
    Tetrahedron4 tet;
    RVector *fun = new RVector[np];
    for (i = 0; i < np; i++)
        fun[i] = tet.LocalShapeF(absc[i]);

    int *npx = new int[blen];
    for (i = 0; i < blen; i++) npx[i] = 0;

    // pass 1: determine storage requirements
    for (idx_i = 0; idx_i < blen; idx_i++) {
	std::cerr << "Basismat pass 1, " << idx_i << "/" << blen << std::endl;
	k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	intersect0 = false;
	for (el = 0; el < nel; el++) {
	    Element *pel = meshptr->elist[el];
	    intersect = false;
	    for (s = 0; s < pel->nNode(); s++) {
	        nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx+ry*ry+rz*rz < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (intersect) {
	        intersect0 = true;
		pel->SetRegion(0);
	    } else
	        pel->SetRegion(-1);
	}
	for (idx_j = 0; idx_j < blen; idx_j++) {
	    if (idx_i == idx_j) { // diagonal elements always exist
	        npx[idx_i]++;
		continue;
	    }
	    if (!intersect0) continue;
	    k1 = idx_j / stride2;
	    tmp = idx_j - k1*stride2;
	    j1 = tmp / stride1;
	    i1 = tmp % stride1;
	    px1 = bbmin[0] + i1*dx;
	    py1 = bbmin[1] + j1*dy;
	    pz1 = bbmin[2] + k1*dz;
	    dstx = px1-px0;
	    dsty = py1-py0;
	    dstz = pz1-pz0;
	    if (dstx*dstx + dsty*dsty + dstz*dstz >= radlimit2*4.0)
		continue;  // no overlap between blobs
	    for (el = 0, intersect1 = false; el < nel && !intersect1; el++) {
		Element *pel = meshptr->elist[el];
		if (pel->Region()) continue; // el doesn't intersect first blob
		for (s = 0; s < 4; s++) {
		    Node &node = meshptr->nlist[pel->Node[s]];
		    rx = px1-node[0];
		    ry = py1-node[1];
		    rz = pz1-node[2];
		    if (rx*rx + ry*ry + rz*rz < radlimit2) {
		        intersect1 = true;
			break;
		    }
		}
	    }
	    if (intersect1)
		npx[idx_i]++; // overlap with mesh domain
	}
    }

    int *rowptr = new int[blen+1];
    rowptr[0] = 0;
    for (i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen];
    int *colidx = new int[nz];

    // pass 2: determine matrix sparsity pattern
    for (idx_i = 0; idx_i < blen; idx_i++) {
	std::cerr << "Basismat pass 2, " << idx_i << "/" << blen << std::endl;
	k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	intersect0 = false;
	for (el = 0; el < nel; el++) {
	    Element *pel = meshptr->elist[el];
	    intersect = false;
	    for (s = 0; s < pel->nNode(); s++) {
	        nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx+ry*ry+rz*rz < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (intersect) {
	        intersect0 = true;
		pel->SetRegion(0);
	    } else
	        pel->SetRegion(-1);
	}
	for (idx_j = 0; idx_j < blen; idx_j++) {
	    if (idx_i == idx_j) {
		colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
		continue;
	    }
	    if (!intersect0) continue;
	    k1 = idx_j / stride2;
	    tmp = idx_j - k1*stride2;
	    j1 = tmp / stride1;
	    i1 = tmp % stride1;
	    px1 = bbmin[0] + i1*dx;
	    py1 = bbmin[1] + j1*dy;
	    pz1 = bbmin[2] + k1*dz;
	    dstx = px1-px0;
	    dsty = py1-py0;
	    dstz = pz1-pz0;
	    if (dstx*dstx + dsty*dsty + dstz*dstz >= radlimit2*4.0)
		continue;  // no overlap between blobs
	    for (el = 0, intersect1 = false; el < nel && !intersect1; el++) {
		Element *pel = meshptr->elist[el];
		if (pel->Region()) continue; // element doesn't intersect first blob
		for (s = 0; s < 4; s++) {
		    Node &node = meshptr->nlist[pel->Node[s]];
		    rx = px1-node[0];
		    ry = py1-node[1];
		    rz = pz1-node[2];
		    if (rx*rx+ry*ry+rz*rz < radlimit2) {
		        intersect1 = true;
			break;
		    }
		}
	    }
	    if (intersect1)
		colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
	}
    }

    delete []npx;
    RCompRowMatrix *bvv = new RCompRowMatrix (blen, blen, rowptr, colidx);

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    tic();
    for (el = 0; el < nel; el++) {
	std::cerr << "Basismat pass 3, " << el << "/" << nel << std::endl;
        Element *pel = meshptr->elist[el];
	for (idx_i = 0; idx_i < blen; idx_i++) {
	    k0 = idx_i / stride2;
	    tmp = idx_i - k0*stride2;
	    j0 = tmp / stride1;
	    i0 = tmp % stride1;
	    px0 = bbmin[0] + i0*dx;
	    py0 = bbmin[1] + j0*dy;
	    pz0 = bbmin[2] + k0*dz;
	    intersect0 = false;
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx + ry*ry + rz*rz < radlimit2) {
		    intersect0 = true;
		    break;
		}
	    }
	    if (!intersect0)
  	        continue; // 1st blob doesn't intersect element
	    
	    for (idx_j = 0; idx_j <= idx_i; idx_j++) {
	        if (idx_i != idx_j) {
		    k1 = idx_j / stride2;
		    tmp = idx_j - k1*stride2;
		    j1 = tmp / stride1;
		    i1 = tmp % stride1;
		    px1 = bbmin[0] + i1*dx;
		    py1 = bbmin[1] + j1*dy;
		    pz1 = bbmin[2] + k1*dz;
		    
		    rx = px1-px0;
		    ry = py1-py0;
		    rz = pz1-pz0;
		    if (rx*rx + ry*ry + rz*rz >= radlimit2*4.0)
		        continue;  // blobs don't intersect each other
		    
		    intersect1 = false;
		    for (s = 0; s < pel->nNode(); s++) {
		        nd = pel->Node[s];
			rx = px1-meshptr->nlist[nd][0];
			ry = py1-meshptr->nlist[nd][1];
			rz = pz1-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect1 = true;
			    break;
			}
		    }
		    if (!intersect1)
		        continue; // 2nd blob doesn't intersect element
		}
		for (s = 0; s < 4; s++) {
		    submesh.elist[0]->Node[s] = s;
		    submesh.nlist[s] = meshptr->nlist[pel->Node[s]];
		}
		submesh.elen_used = 1;
		submesh.nlen_used = 4;
		submesh.elist[0]->SetRegion(0);
		Point cnt(3);
		cnt[0] = px0, cnt[1] = py0, cnt[2] = pz0;
		IcosaCut (submesh, sup, cnt);
		if (idx_i != idx_j) {
		    cnt[0] = px1, cnt[1] = py1, cnt[2] = pz1;
		    IcosaCut (submesh, sup, cnt);
		}
		submesh.SubSetup();
		v = 0.0;
		for (s = 0; s < submesh.elen_used; s++) {
		    Element *psel = submesh.elist[s];
		    if (psel->Region()) continue;
		    for (m = 0; m < np; m++) {
		        djac = psel->DetJ(absc[m], &submesh.nlist);
			// quadrature point in global frame
			for (i = 0; i < 3; i++)
			    for (glob[i] = 0.0, j = 0; j < 4; j++)
			        glob[i] += fun[m][j] * submesh.nlist[psel->Node[j]][i];
			//Point glob = psel->Global (submesh.nlist, absc[m]);
			v += wght[m] * djac *
			    Value_nomask(glob,idx_i,false) *
			    Value_nomask(glob,idx_j,false);
		    }
		}
		if (v) {
		    (*bvv)(idx_i,idx_j) += v;
		    if (idx_i != idx_j)
		        (*bvv)(idx_j,idx_i) += v;
		}
	    }
	}
    }
    std::cerr << "time: " << toc() << std::endl;

    delete []rowptr;
    delete []colidx;
	delete []fun;

    // diagonal conditioning
    RVector d = bvv->Diag();
    double dmean = mean(d);
    for (i = 0; i < blen; i++)
	(*bvv)(i,i) += dmean * dgscale;

    bvv->Shrink();
    std::cerr << "exit BasisMatrix" << std::endl;
    return bvv;
}






#ifdef UNDEF
// SIMPLE DIAGONAL VERSION
RCompRowMatrix *Raster_Blob2::CreateBvv_tet4 () const
{
    int i, i0, j0, k0, i1, j1, k1, s, m, el, nd, idx_i, idx_j, tmp;
    int nel = meshptr->elen();
    bool intersect, intersect0, intersect1;
    double px0, py0, pz0, px1, py1, pz1, rx, ry, rz, dstx, dsty, dstz, djac, v;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double dz = zrange/(bdim[2]-1.0);
    double radlimit2 = sup*sup;
    int stride1 = bdim[0];
    int stride2 = bdim[0]*bdim[1];

    // quadrature rule for local tet
    const double *wght;
    const Point *absc;
    int np = QRule_tet_6_24 (&wght, &absc);

    int *npx = new int[blen];
    for (i = 0; i < blen; i++) npx[i] = 0;

    // pass 1: determine storage requirements
    for (idx_i = 0; idx_i < blen; idx_i++) {
        npx[idx_i] = 1;
    }

    int *rowptr = new int[blen+1];
    rowptr[0] = 0;
    for (i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen];
    int *colidx = new int[nz];

    // pass 2: determine matrix sparsity pattern
    for (idx_i = 0; idx_i < blen; idx_i++) {
        colidx[idx_i] = idx_i;
    }

    delete []npx;
    RCompRowMatrix *bvv = new RCompRowMatrix (blen, blen, rowptr, colidx);

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    for (el = 0; el < nel; el++) {
        Element *pel = meshptr->elist[el];
	for (idx_i = 0; idx_i < blen; idx_i++) {
	    k0 = idx_i / stride2;
	    tmp = idx_i - k0*stride2;
	    j0 = tmp / stride1;
	    i0 = tmp % stride1;
	    px0 = bbmin[0] + i0*dx;
	    py0 = bbmin[1] + j0*dy;
	    pz0 = bbmin[2] + k0*dz;
	    intersect0 = false;
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx + ry*ry + rz*rz < radlimit2) {
		    intersect0 = true;
		    break;
		}
	    }
	    if (!intersect0) continue;
	    for (s = 0; s < 4; s++) {
	        submesh.elist[0]->Node[s] = s;
		submesh.nlist[s] = meshptr->nlist[pel->Node[s]];
	    }
	    submesh.elen_used = 1;
	    submesh.nlen_used = 4;
	    submesh.elist[0]->SetRegion(0);
	    Point cnt(3);
	    cnt[0] = px0, cnt[1] = py0, cnt[2] = pz0;
	    IcosaCut (submesh, sup, cnt);
	    submesh.SubSetup();

	    for (s = 0; s < submesh.elen_used; s++) {
		Element *psel = submesh.elist[s];
		if (psel->Region()) continue;
		for (m = 0; m < np; m++) {
		    djac = psel->DetJ(absc[m], &submesh.nlist);
		    Point glob = psel->Global (submesh.nlist, absc[m]);
		    v = wght[m] * djac *
			Value_nomask(glob,idx_i,false) *
			Value_nomask(glob,idx_i,false);
		    (*bvv)(idx_i,idx_i) += v;
		}
	    }
	}
    }

#ifdef UNDEF
    // pass 3: fill the matrix
    BufMesh submesh1, submesh2, submesh(*meshptr);
    double t_copy = 0.0, t_cut = 0.0, t_setup = 0.0, t_comp = 0.0;
    for (idx_i = 0; idx_i < blen; idx_i++) {
	std::cerr << "pass 3, " << idx_i << "/" << blen
		  << ", t_copy=" << t_copy << ", t_cut=" << t_cut
		  << ", t_setup=" << t_setup << ", t_comp=" << t_comp
		  << std::endl;
	k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;

	intersect = false;
	for (el = 0; el < nel; el++) {
	    Element *pel = meshptr->elist[el];
	    intersect0 = false;
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx + ry*ry + rz*rz < radlimit2) {
		    intersect0 = true;
		    break;
		}
	    }
	    if (intersect0) intersect = true;
	    submesh.elist[el]->SetRegion (intersect0 ? 0 : -1);
	}
	if (!intersect) continue;
	Point cnt(3);
	cnt[0] = px0, cnt[1] = py0, cnt[2] = pz0;
	submesh2.Copy (submesh);
	//submesh1.Shrink();
	IcosaCut (submesh2, sup, cnt); // cut mesh with first blob
	if (!submesh2.elen_used) continue; // blob has no mesh support

	for (idx_j = 0; idx_j <= idx_i; idx_j++) {
	    if (idx_j != idx_i) continue;
	    tic();
	    submesh2.SubSetup();
	    t_setup += toc();
	    tic();
	    for (s = 0; s < submesh2.elen_used; s++) {
		Element *psel = submesh2.elist[s];
		if (psel->Region()) continue;
		for (m = 0; m < np; m++) {
		    djac = psel->DetJ(absc[m], &submesh2.nlist);
		    Point glob = psel->Global (submesh2.nlist, absc[m]);
		    v = wght[m] * djac *
			Value_nomask(glob,idx_i,false) *
			Value_nomask(glob,idx_j,false);
		    (*bvv)(idx_i,idx_j) += v;
		    if (idx_i != idx_j)
			(*bvv)(idx_j,idx_i) += v;
		}
	    }
	    t_comp += toc();
	}
    }
#endif
    delete []rowptr;
    delete []colidx;

    // diagonal conditioning
    RVector d = bvv->Diag();
    double dmean = mean(d);
    for (i = 0; i < blen; i++)
	(*bvv)(i,i) += dmean * dgscale;

    bvv->Shrink();
    std::cerr << "exit BasisMatrix" << std::endl;
    return bvv;
}
#endif


// ==========================================================================

RCompRowMatrix *Raster_Blob2::CreateBuv_tet4 () const
{
    int i, j, k, s, m, nd, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, idx_i, idx_j;
    bool intersect;
    RVector fun;
    int *npx = new int[n];
    for (i = 0; i < n; i++) npx[i] = 0;

    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double zrange = bbmax[2]-bbmin[2];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double dz = zrange/(bdim[2]-1.0);
    double djac, v;
    double rx, ry, rz;
    int stride1 = bdim[0], stride2 = bdim[0]*bdim[1];
    double radlimit2 = sup*sup;

    // quadrature rule for local tetrahedron
    const double *wght;
    const Point *absc;
    int np = QRule_tet_6_24 (&wght, &absc);

    // pass 1: determine matrix fill structure
    double px, py, pz;
    for (el = 0; el < nel; el++) {
	std::cerr << "pass1, el=" << el << "/" << nel << std::endl;
	Element *pel = meshptr->elist[el];
	for (k = 0; k < bdim[2]; k++) {
	    pz = bbmin[2] + k*dz;
	    for (j = 0; j < bdim[1]; j++) {
		py = bbmin[1] + j*dy;
		for (i = 0; i < bdim[0]; i++) {
		    px = bbmin[0] + i*dx;
		    intersect = false;
		    for (s = 0; s < pel->nNode(); s++) {
			nd = pel->Node[s];
			rx = px-meshptr->nlist[nd][0];
			ry = py-meshptr->nlist[nd][1];
			rz = pz-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    break;
			}
		    }
		    if (intersect) {
			for (s = 0; s < pel->nNode(); s++)
			    npx[pel->Node[s]]++;
		    }
		}
	    }
	}
    }

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (i = 0; i < n; i++)
	rowptr[i+1] = rowptr[i]+npx[i];
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (i = 0; i < n; i++)
	npx[i] = 0;
    for (el = 0; el < nel; el++) {
	std::cerr << "pass2, el=" << el << "/" << nel << std::endl;
	Element *pel = meshptr->elist[el];
	for (k = 0; k < bdim[2]; k++) {
	    pz = bbmin[2] + k*dz;
	    for (j = 0; j < bdim[1]; j++) {
		py = bbmin[1] + j*dy;
		for (i = 0; i < bdim[0]; i++) {
		    px = bbmin[0] + i*dx;
		    intersect = false;
		    for (s = 0; s < pel->nNode(); s++) {
			nd = pel->Node[s];
			rx = px-meshptr->nlist[nd][0];
			ry = py-meshptr->nlist[nd][1];
			rz = pz-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    break;
			}
		    }
		    if (intersect) {
			for (s = 0; s < pel->nNode(); s++) {
			    nd = pel->Node[s];
			    colidx[rowptr[nd]+npx[nd]++] =
				i + j*stride1 + k*stride2;
			}
		    }
		}
	    }
	}
    }

    RCompRowMatrix *buv = new RCompRowMatrix (n, blen, rowptr, colidx);
    //RVector scale(blen);

    BufMesh submesh;
    submesh.elist.Append(new Tetrahedron4);
    submesh.nlist.New(4);
    submesh.nlen_used = 0;
    submesh.elen_used = 0;

    for (el = 0; el < nel; el++) {
	std::cerr << "pass3, el=" << el << "/" << nel << std::endl;
	Element *pel = meshptr->elist[el];

	for (k = 0; k < bdim[2]; k++) {
	    pz = bbmin[2] + k*dz;
	    for (j = 0; j < bdim[1]; j++) {
		py = bbmin[1] + j*dy;
		for (i = 0; i < bdim[0]; i++) {
		    px = bbmin[0] + i*dx;
		    intersect = false;
		    for (s = 0; s < pel->nNode(); s++) {
			nd = pel->Node[s];
			rx = px-meshptr->nlist[nd][0];
			ry = py-meshptr->nlist[nd][1];
			rz = pz-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    break;
			}
		    }
		    if (intersect) {

			// create a mesh containing this single element
			for (s = 0; s < 4; s++) {
			    submesh.elist[0]->Node[s] = s;
			    submesh.nlist[s] = meshptr->nlist[pel->Node[s]];
			}
			submesh.elen_used = 1;
			submesh.nlen_used = 4;
			submesh.elist[0]->SetRegion(0);
			Point cnt(3);
			cnt[0] = px, cnt[1] = py, cnt[2] = pz;
			IcosaCut (submesh, sup, cnt);
			submesh.SubSetup();
			idx_j = i + j*stride1 + k*stride2;
			for (s = 0; s < submesh.elen_used; s++) {
			    Element *psel = submesh.elist[s];
			    if (psel->Region()) continue;
			    // map quadrature points into global frame
			    for (m = 0; m < np; m++) {
				djac = psel->DetJ(absc[m], &submesh.nlist);
				Point glob = psel->Global (submesh.nlist, absc[m]);
				Point loc = pel->Local(meshptr->nlist, glob);
				fun = pel->LocalShapeF (loc);
				v = wght[m] * djac * Value_nomask(glob,idx_j,false);
				for (ii = 0; ii < fun.Dim(); ii++) {
				    idx_i = pel->Node[ii];
				    (*buv)(idx_i, idx_j) += v*fun[ii];
				}
				//scale[idx_j] += v;
			    }
			}
		    }
		}
	    }
	}
    }

    delete []rowptr;
    delete []colidx;
    delete []npx;

    buv->Shrink();
    return buv;
}

// ==========================================================================
// Map directly between basis and a regular voxel image with piecewise
// constant basis functions

RCompRowMatrix *Raster_Blob2::CreateBvw_tet4 (const IVector &pxdim) const
{
    int i, j, k, idx_i, idx_j, i0, j0, k0, i1, j1, k1, jj, tmp;
    double px0, py0, pz0, px1, py1, pz1, x, y, z, rx, ry, rz, v;
    int stride1 = bdim[0];
    int stride2 = bdim[0]*bdim[1];
    int pstride1 = pxdim[0];
    int pstride2 = pxdim[0]*pxdim[1];
    int plen = pxdim[0]*pxdim[1]*pxdim[2];
    int *npx = new int[blen];
    for (i = 0; i < blen; i++) npx[i] = 0;
    double dx = grid[0];
    double dy = grid[1];
    double dz = grid[2];
    double pdx = bbsize[0]/pxdim[0];
    double pdy = bbsize[1]/pxdim[1];
    double pdz = bbsize[2]/pxdim[2];
    double radlimit2 = sup*sup;
    bool intersect;
    int subdiv = 10;
    double subwght = 1.0/(subdiv*subdiv*subdiv);
    Point p(3);

    // pass 1: determine storage requirements
    for (idx_i = 0; idx_i < blen; idx_i++) {
 	std::cerr << "BasisPixelmat pass 1, " << idx_i << "/" << blen << std::endl;
        k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	for (idx_j = 0; idx_j < plen; idx_j++) {
	    k1 = idx_j / pstride2;
	    tmp = idx_j - k1*pstride2;
	    j1 = tmp / pstride1;
	    i1 = tmp % pstride1;
	    px1 = bbmin[0] + i1*pdx;
	    py1 = bbmin[1] + j1*pdy;
	    pz1 = bbmin[2] + k1*pdz;
	    intersect = false;
	    for (k = 0; k < 2; k++) {
	        for (j = 0; j < 2; j++) {
		    for (i = 0; i < 2; i++) {
		        rx = px0 - (px1+i*pdx);
			ry = py0 - (py1+j*pdy);
			rz = pz0 - (pz1+k*pdz);
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    i = j = k = 2; // break
			}
		    }
		}
	    }
	    if (intersect)
	        npx[idx_i]++;
	}
    }

    int *rowptr = new int[blen+1];
    rowptr[0] = 0;
    for (i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen];
    int *colidx = new int[nz];

    // pass 2: determine matrix sparsity pattern
    for (idx_i = 0; idx_i < blen; idx_i++) {
  	std::cerr << "BasisPixelmat pass 2, " << idx_i << "/" << blen << std::endl;
	k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	for (idx_j = 0; idx_j < plen; idx_j++) {
	    k1 = idx_j / pstride2;
	    tmp = idx_j - k1*pstride2;
	    j1 = tmp / pstride1;
	    i1 = tmp % pstride1;
	    px1 = bbmin[0] + i1*pdx;
	    py1 = bbmin[1] + j1*pdy;
	    pz1 = bbmin[2] + k1*pdz;
	    intersect = false;
	    for (k = 0; k < 2; k++) {
	        for (j = 0; j < 2; j++) {
		    for (i = 0; i < 2; i++) {
		        rx = px0 - (px1+i*pdx);
			ry = py0 - (py1+j*pdy);
			rz = pz0 - (pz1+k*pdz);
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    i = j = k = 2; // break
			}
		    }
		}
	    }
	    if (intersect)
		colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
	}      
    }

    delete []npx;
    RCompRowMatrix *bvw = new RCompRowMatrix (blen, plen, rowptr, colidx);

    for (idx_i = 0; idx_i < blen; idx_i++) {
 	std::cerr << "BasisPixelmat pass 3, " << idx_i << "/" << blen << std::endl;
        k0 = idx_i / stride2;
	tmp = idx_i - k0*stride2;
	j0 = tmp / stride1;
	i0 = tmp % stride1;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	pz0 = bbmin[2] + k0*dz;
	for (jj = rowptr[idx_i]; jj < rowptr[idx_i+1]; jj++) {
	    idx_j = colidx[jj];
	    k1 = idx_j / pstride2;
	    tmp = idx_j - k1*pstride2;
	    j1 = tmp / pstride1;
	    i1 = tmp % pstride1;
	    px1 = bbmin[0] + i1*pdx;
	    py1 = bbmin[1] + j1*pdy;
	    pz1 = bbmin[2] + k1*pdz;
	    v = 0.0;
	    for (k = 0; k < subdiv; k++) {
	        z = pz1 + (k+0.5)/subdiv*pdz;
		for (j = 0; j < subdiv; j++) {
		    y = py1 + (j+0.5)/subdiv*pdy;
		    for (i = 0; i < subdiv; i++) {
		        x = px1 + (i+0.5)/subdiv*pdx; 
			rx = x-px0;
			ry = y-py0;
			rz = z-pz0;
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    p[0] = x; p[1] = y; p[2] = z;
			    v += Value_nomask(p,idx_i,false);
			}
		    }
		}
	    }
	    (*bvw)(idx_i,idx_j) = v*subwght;
	}
    }

    delete []rowptr;
    delete []colidx;

    bvw->Shrink();
    return bvw;
}

// ==========================================================================
// Split the tetrahedra in 'mesh' along a plane defined by 3 points

void PlaneCut (BufMesh &mesh, const Point *pt, const RDenseMatrix &R)
{
    // rotate the mesh with R
    mesh.Rotate (R);

    // perform plane cutting
    double x = R(0,0)*pt[0][0] + R(0,1)*pt[0][1] + R(0,2)*pt[0][2];
    int i, n = mesh.elen_used;
    for (i = 0; i < n; i++)
	if (!mesh.elist[i]->Region())
	    Tetsplit (&mesh, i, 0, x);

    // rotate back
    mesh.Rotate (transpose(R));
}


// ==========================================================================
// Split the tetrahedra in 'mesh' along the surface of a tetrahedron

void IcosaCut (BufMesh &mesh, double scale, const Point &cnt)
{
    const int nvtx = 12;
    const int nidx = 20;
    const int dim = 3;
    const int trind = 3;

    double ico_vtx[nvtx][dim] = {
	{                   0,                   0,   1.000000000000000},
	{   0.894427190999916,                   0,   0.447213595499958},
	{  -0.894427190999916,                   0,  -0.447213595499958},
	{                   0,                   0,  -1.000000000000000},
	{  -0.723606797749979,   0.525731112119134,   0.447213595499958},
	{  -0.723606797749979,  -0.525731112119134,   0.447213595499958},
	{   0.723606797749979,   0.525731112119134,  -0.447213595499958},
	{   0.723606797749979,  -0.525731112119134,  -0.447213595499958},
	{   0.276393202250021,   0.850650808352040,   0.447213595499958},
	{   0.276393202250021,  -0.850650808352040,   0.447213595499958},
	{  -0.276393202250021,   0.850650808352040,  -0.447213595499958},
	{  -0.276393202250021,  -0.850650808352040,  -0.447213595499958}
    };
    int idx[nidx][trind] = {
	{     0,     1,     8},
	{     0,     8,     4},
	{     0,     4,     5},
	{     0,     5,     9},
	{     0,     9,     1},
	{     1,     6,     8},
	{     8,     6,    10},
	{     8,    10,     4},
	{     4,    10,     2},
	{     4,     2,     5},
	{     5,     2,    11},
	{     5,    11,     9},
	{     9,    11,     7},
	{     9,     7,     1},
	{     1,     7,     6},
	{     3,     6,     7},
	{     3,     7,    11},
	{     3,    11,     2},
	{     3,     2,    10},
	{     3,    10,     6}

    };

    int i, j, k;
    Point pt[3];
    for (j = 0; j < trind; j++)
	pt[j].New(3);

    // compute the rotation matrices to rotate each of the faces of the
    // icosahedron so it has outward normal [1,0,0];
    static bool need_setup = true;
    static RDenseMatrix R[nidx];
    if (need_setup) {
	for (i = 0; i < nidx; i++) {
	    R[i].New(3,3);
	    RDenseMatrix &Ri = R[i];
	    for (j = 0; j < trind; j++) {
		for (k = 0; k < dim; k++) {
		    pt[j][k] = ico_vtx[idx[i][j]][k];
		}
	    }
	    RVector nml(3);
	    nml[0] = pt[0][1]*(pt[1][2]-pt[2][2]) + pt[1][1]*(pt[2][2]-pt[0][2])
		+ pt[2][1]*(pt[0][2]-pt[1][2]);
	    nml[1] = pt[0][0]*(pt[2][2]-pt[1][2]) + pt[1][0]*(pt[0][2]-pt[2][2])
		+ pt[2][0]*(pt[1][2]-pt[0][2]);
	    nml[2] = pt[0][0]*(pt[1][1]-pt[2][1]) + pt[1][0]*(pt[2][1]-pt[0][1])
		+ pt[2][0]*(pt[0][1]-pt[1][1]);
	    nml /= length(nml);
	    
	    double c = nml[0], v2 = nml[2], v3 = -nml[1];
	    double s = hypot(v2,v3);
	    double scale = (1.0-c)/(s*s);
	    Ri(0,0) = 1.0 + (-v2*v2-v3*v3)*scale;
	    Ri(0,1) = -v3;
	    Ri(0,2) = v2;
	    Ri(1,0) = v3;
	    Ri(1,1) = 1.0 - v3*v3*scale;
	    Ri(1,2) = v2*v3*scale;
	    Ri(2,0) = -v2;
	    Ri(2,1) = v2*v3*scale;
	    Ri(2,2) = 1.0 - v2*v2*scale;
	}
	need_setup = false;
    }

    for (i = 0; i < nidx; i++) {
	for (j = 0; j < trind; j++) {
	    for (k = 0; k < dim; k++) {
		pt[j][k] = cnt[k] + ico_vtx[idx[i][j]][k]*scale;
	    }
	}
	PlaneCut (mesh, pt, R[i]);
    }
    for (i = 0; i < mesh.elen_used; i++)
    	if (mesh.elist[i]->Region())
    	    mesh.elist[i]->SetRegion (-1);
    mesh.Shrink();
}
