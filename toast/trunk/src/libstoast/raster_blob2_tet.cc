// Implementation of those routines in Raster_Blob2 specific to tetrahedral
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_blob2.h"
#include "tet_qr.h"
#include "tetsplit.h"

// ==========================================================================
// Local prototypes

void PlaneCut (BufMesh &mesh, const Point *pt, const RDenseMatrix &R);
void IcosaCut (BufMesh &mesh, double scale, const Point &cnt);

// ==========================================================================

RCompRowMatrix *Raster_Blob2::CreateBasisMassmat_tet4 () const
{
    int i, i0, j0, k0, i1, j1, k1, s, el, nd, idx_i, idx_j, tmp;
    int nel = meshptr->elen();
    bool intersect;
    double px0, py0, pz0, px1, py1, pz1, rx, ry, rz, dstx, dsty, dstz;

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

    bool *have_diag = new bool[blen];
    for (i = 0; i < blen; i++) have_diag[i] = false;

    // pass 1: determine matrix fill structure
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (k0 = 0; k0 < bdim[2]; k0++) {
	    pz0 = bbmin[2] + k0+dz;
	    for (j0 = 0; j0 < bdim[1]; j0++) {
		py0 = bbmin[1] + j0*dy;
		for (i0 = 0; i0 < bdim[0]; i0++) {
		    px0 = bbmin[0] + i0*dx;
		    intersect = false;
		    for (s = 0; s < pel->nNode(); s++) {
			nd = pel->Node[s];
			rx = px0-meshptr->nlist[nd][0];
			ry = py0-meshptr->nlist[nd][1];
			rz = pz0-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    break;
			}
		    }
		    if (intersect) {
			for (k1 = 0; k1 < bdim[2]; k1++) {
			    pz1 = bbmin[2] + k1*dz;
			    for (j1 = 0; j1 < bdim[1]; j1++) {
				py1 = bbmin[1] + j1*dy;
				for (i1 = 0; i1 < bdim[0]; i1++) {
				    px1 = bbmin[0] + i1*dx;
				    dstx = px1-px0;
				    dsty = py1-py0;
				    dstz = pz1-pz0;
				    if (dstx*dstx + dsty*dsty + dstz*dstz
					>= radlimit2) continue;
				    intersect = false;
				    for (s = 0; s < pel->nNode(); s++) {
					nd = pel->Node[s];
					rx = px1-meshptr->nlist[nd][0];
					ry = py1-meshptr->nlist[nd][1];
					rz = pz1-meshptr->nlist[nd][2];
					if (rx*rx + ry*ry + rz*rz < radlimit2) {
					    intersect = true;
					    break;
					}
				    }
				    if (!intersect) continue;
				    if (i0 + j0*stride1 + k0*stride2 ==
					i1 + j1*stride1 + k1*stride2)
					have_diag[i0+j0*stride1+k0*stride2] =
					    true;
				    npx[i0 + j0*stride1 + k0*stride2]++;
				}
			    }
			}
		    }
		}
	    }
	}
    }

    for (i = 0; i < blen; i++) {
	if (!have_diag[i]) // reserve space for diag element
	    npx[i]++;
	have_diag[i] = false;
    }

    int *rowptr = new int[blen+1];
    rowptr[0] = 0;
    for (i = 0; i < blen; i++)
	rowptr[i+1] = rowptr[i]+npx[i];
    int nz = rowptr[blen];
    int *colidx = new int[nz];
    for (i = 0; i < blen; i++)
	npx[i] = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (k0 = 0; k0 < bdim[2]; k0++) {
	    pz0 = bbmin[2] + k0*dz;
	    for (j0 = 0; j0 < bdim[1]; j0++) {
		py0 = bbmin[1] + j0*dy;
		for (i0 = 0; i0 < bdim[0]; i0++) {
		    px0 = bbmin[0] + i0*dx;
		    intersect = false;
		    for (s = 0; s < pel->nNode(); s++) {
			nd = pel->Node[s];
			rx = px0-meshptr->nlist[nd][0];
			ry = py0-meshptr->nlist[nd][1];
			rz = pz0-meshptr->nlist[nd][2];
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    break;
			}
		    }
		    idx_i = i0 + j0*stride1 + k0*stride2;
		    if (intersect) {
			for (k1 = 0; k1 < bdim[2]; k1++) {
			    pz1 = bbmin[2] + k1*dz;
			    for (j1 = 0; j1 < bdim[1]; j1++) {
				py1 = bbmin[1] + j1*dy;
				for (i1 = 0; i1 < bdim[0]; i1++) {
				    px1 = bbmin[0] + i1*dx;
				    dstx = px1-px0;
				    dsty = py1-py0;
				    dstz = pz1-pz0;
				    if (dstx*dstx + dsty*dsty + dstz*dstz
					>= radlimit2) continue;
				    intersect = false;
				    for (s = 0; s < pel->nNode(); s++) {
					nd = pel->Node[s];
					rx = px1-meshptr->nlist[nd][0];
					ry = py1-meshptr->nlist[nd][1];
					rz = pz1-meshptr->nlist[nd][2];
					if (rx*rx + ry*ry + rz*rz < radlimit2) {
					    intersect = true;
					    break;
					}
				    }
				    if (!intersect) continue;
				    idx_j = i1 + j1*stride1 + k1*stride2;
				    if (idx_i == idx_j)
					have_diag[idx_i] = true;
				    colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    for (i = 0; i < blen; i++) {
	if (!have_diag[i])
	    colidx[rowptr[i] + npx[i]++] = i;
    }
    delete []have_diag;

    RCompRowMatrix *bvv = new RCompRowMatrix (blen, blen, rowptr, colidx);

    BufMesh submesh0, submesh1;
    submesh0.elist.Append(new Tetrahedron4);
    submesh0.nlist.New(4);
    submesh0.nlen_used = 0;
    submesh0.elen_used = 0;
    submesh1.elist.Append(new Tetrahedron4);
    submesh1.nlist.New(4);
    submesh1.nlen_used = 0;
    submesh1.elen_used = 0;

    // pass 2: fill the matrix
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
	    intersect = false;
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		rz = pz0-meshptr->nlist[nd][2];
		if (rx*rx + ry*ry + rz*rz < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (!intersect) continue;
	    // create a mesh contining this single element
	    for (s = 0; s < 4; s++) {
		submesh0.elist[0]->Node[s] = s;
		submesh0.nlist[s] = meshptr->nlist[pel->Node[s]];
	    }
	    submesh0.elen_used = 1;
	    submesh0.nlen_used = 4;
	    submesh0.elist[0]->SetRegion(0);
	    Point cnt(3);
	    cnt[0] = px0, cnt[1] = py0, cnt[2] = pz0;
	    IcosaCut (submesh0, sup, cnt);
	    for (idx_j = 0; idx_j <= idx_i; idx_j++) {
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
		if (dstx*dstx + dsty*dsty + dstz*dstz >= radlimit2)
		    continue;
		intersect = false;
		for (s = 0; s < pel->nNode(); s++) {
		    nd = pel->Node[s];
		    rx = px1-meshptr->nlist[nd][0];
		    ry = py1-meshptr->nlist[nd][1];
		    rz = pz1-meshptr->nlist[nd][2];
		    if (rx*rx + ry*ry + rz*rz < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (!intersect) continue;
		for (s = 0; s < 4; s++) {
		    submesh1.elist[0]->Node[s] = s;
		    submesh1.nlist[s] = meshptr->nlist[pel->Node[s]];
		}
		submesh1.elen_used = 1;
		submesh1.nlen_used = 4;
		submesh1.elist[0]->SetRegion(0);
		Point cnt(3);
		cnt[0] = px1, cnt[1] = py1, cnt[2] = pz1;
		IcosaCut (submesh1, sup, cnt);
		// what now? need intersection of the two submeshes
	    }
	}
    }

    delete []rowptr;
    delete []colidx;
    delete []npx;

    bvv->Shrink();
    return bvv;
}

// ==========================================================================

RCompRowMatrix *Raster_Blob2::CreateMixedMassmat_tet4 () const
{
    int i, j, k, r, s, m, nd, el, nel = meshptr->elen(), n = meshptr->nlen();
    int ii, jj, idx_i, idx_j;
    int imin, imax, jmin, jmax, kmin, kmax;
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
    double xmin, xmax, ymin, ymax, zmin, zmax, djac, vb, v;
    double rx, ry, rz;
    int stride_i = bdim[0], stride_j = stride_i*bdim[1];
    double radlimit2 = sup*sup;

    // quadrature rule for local tetrahedron
    const double *wght;
    const Point *absc;
    int np = QRule_tet_4_14 (&wght, &absc);

    // pass 1: determine matrix fill structure
    double px, py, pz;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (k = 0; k < bdim[2]; k++) {
	    pz = bbmin[2] + k+dz;
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
			if (rx*rx + ry*ry + rz*rz < radlimit2) {
			    intersect = true;
			    break;
			}
		    }
		    if (intersect) {
			for (s = 0; s < pel->nNode(); s++) {
			    nd = pel->Node[s];
			    colidx[rowptr[nd]+npx[nd]++] =
				i + j*bdim[0] + k*bdim[0]*bdim[1];
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
			idx_j = i + (j + k*bdim[1])*bdim[0];
			for (s = 0; s < submesh.elen_used; i++) {
			    Element *psel = submesh.elist[i];
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
