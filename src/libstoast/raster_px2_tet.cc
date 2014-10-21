// Implementation of those routines in Raster_Pixel2 specific to tetrahedral
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tet_qr.h"

// ==========================================================================
// Mesh class with pre-allocated element and node lists, so we don't
// need to dynamically grow the lists too often
// nlen_used and elen_used are the length of the lists actually in use

class BufMesh: public Mesh
{
public:
    BufMesh(): Mesh()
    { nlen_used = elen_used = 0; }

    void SubSetup()
    {
	for (int el=0; el < elen_used; el++)
	    if (elist[el]->Region() >= 0)
		elist[el]->Initialise(nlist);
	for (int el=0; el < elen_used; el++)
	    if (elist[el]->Region() >= 0)
		elist[el]->PostInitialisation(nlist);
    }

    friend std::ostream& operator<< (std::ostream& o, BufMesh &mesh)
    {
	o << "MeshData 5.0\n\n";
	o << "NodeList " << mesh.nlen_used << " 1" << std::endl;
	o.precision(10);
	for (int i = 0; i < mesh.nlen_used; i++)
	    o << "N[" << mesh.nlist[i][0] << " " << mesh.nlist[i][1] << " "
	      << mesh.nlist[i][2] << "]R0" << std::endl;
	o << "\nElementList " << mesh.elen_used << std::endl;
	for (int i = 0; i < mesh.elen_used; i++)
	    o << "c " << mesh.elist[i]->Node[0]+1 << " " << mesh.elist[i]->Node[1]+1
	      << " " << mesh.elist[i]->Node[2]+1 << " " << mesh.elist[i]->Node[3]+1
	      << " R" << mesh.elist[i]->Region() << std::endl;
	o << "\n[ParameterList]N" << std::endl;
	o << "Size " << mesh.nlen_used << std::endl;
	o << "Param1 MUA" << std::endl;
	o << "Param2 KAPPA" << std::endl;
	o << "Param3 N" << std::endl;
	o << "Data" << std::endl;
	for (int i = 0; i < mesh.nlen_used; i++)
	    o << "0.01 0.33 1.4" << std::endl;
	return o;
    }

    int nlen_used, elen_used;
};

// ==========================================================================
// Create the mixed-basis mass matrix Buv by subdividing tetrahedral elements
// at voxel boundaries

RCompRowMatrix *Raster_Pixel2::CreateMixedMassmat_tet4 () const
{
    void Tetsplit (BufMesh *mesh, int el, int cut_orient, double cut_pos);
    void Tetsplit_cube (BufMesh *mesh, double xmin, double xmax,
			double ymin, double ymax, double zmin, double zmax,
			int reg);

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

// =========================================================================
// =========================================================================
// Tetrahedra splitting code. This is preliminary and should eventually
// be shifted elsewhere (e.g. to the individual element classes)
// =========================================================================
// =========================================================================

// Calculate the intersection point pi of a line defined by two points p0
// and p1, and a cutting plane parallel to two coordinate axes
// cut_orient defines the orientation of the cutting plane
// (0: x=const, 1: y=const, 2: z=const)
// cut_pos is the position of the intersection of the cutting plane with
// the corresponding axis
// It is assumed that there is an intersection, i.e. the line is not
// parallel to the plane
inline void Intersect_line_plane (int cut_orient, double cut_pos,
			   const Point &p0, const Point &p1, Point &pi)
{
    double s = (cut_pos-p0[cut_orient])/(p1[cut_orient]-p0[cut_orient]);
    for (int i = 0; i < 3; i++)
	pi[i] = (i == cut_orient ? cut_pos : p0[i] + (p1[i]-p0[i])*s);
}


// Check if a point is left of a cutting plane
inline bool Point_is_left (const Point &p, int cut_orient, double cut_pos)
{
    const double eps = 1e-8;
    return (p[cut_orient] < cut_pos-eps);
}


// Check if a point is right of a cutting plane
inline bool Point_is_right (const Point &p, int cut_orient, double cut_pos)
{
    const double eps = 1e-8;
    return (p[cut_orient] > cut_pos+eps);
}


// Check if a tetrahedron has a portion left of a cutting plane
bool Tet_extends_left (Mesh *mesh, int el, int cut_orient, double cut_pos)
{
    Element *pel = mesh->elist[el];
    for (int i = 0; i < 4; i++)
	if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
	    return true;
    return false;
}


// Check if a tetrahedron has a portion right of a cutting plane
bool Tet_extends_right (Mesh *mesh, int el, int cut_orient, double cut_pos)
{
    Element *pel = mesh->elist[el];
    for (int i = 0; i < 4; i++)
	if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
	    return true;
    return false;
}


// Calculate the intersection points of the edges of tetrahedron el
// with a cutting plane parallel to two axes.
// cut_orient defines the orientation of the cutting plane
// (0: x=const, 1: y=const, 2: z=const)
// cut_pos is the position of the intersection of the cutting plane with
// the corresponding axis
// Return value is the number of intersection points found (0, 3, or 4)
// The points are written to the array pointed to by isect (must be allocated to
// at least 4 points of dimension 3)
int CalcIntersections_tet(Mesh *mesh, int el, int cut_orient, double cut_pos, Point *isect)
{
    int i, j, k, nleft, nright, nsingle;
    Element *pel = mesh->elist[el];

    for (i = nleft = nright = 0; i < 4; i++) {
	Point &pt = mesh->nlist[pel->Node[i]];
	if      (Point_is_left  (pt, cut_orient, cut_pos))  nleft++;
	else if (Point_is_right  (pt, cut_orient, cut_pos)) nright++;
    }

    if (!nleft || !nright) return 0; // no intersection

    if (nleft == 1 || nright == 1) { // triangular intersection

	// find the single node
	if (nleft == 1) {
	    for (i = 0; i < 4; i++)
		if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		    nsingle = i; break;
		}
	} else {
	    for (i = 0; i < 4; i++)
		if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		    nsingle = i; break;
		}
	}
	Point &psingle = mesh->nlist[pel->Node[nsingle]];
	double vsingle = psingle[cut_orient];

	// find the intersection points on the 3 edges connected to psingle
	for (i = j = 0; i < 4; i++) {
	    if (i == nsingle) continue;
	    Point &p = mesh->nlist[pel->Node[i]];
	    Intersect_line_plane (cut_orient, cut_pos, psingle, p, isect[j++]);
	}

	return 3;

    } else { // quadrilateral intersection

	int nleft[2], nright[2];
	for (i = j = k = 0; i < 4; i++) {
	    if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
		nleft[j++] = i;
	    else
		nright[k++] = i;
	}

	// re-arrange a bit for simpler computation
	if (nleft[0] != 0 && nleft[1] != 0) {
	    nleft[0] = 0;
	    nleft[1] = (nright[0] == 0 ? nright[1]:nright[0]);
	    nright[0] = (nleft[1] == 1 ? 2:1);
	    nright[1] = (nleft[1] == 3 ? 2:3);
	} else {
	    if (nleft[0] != 0) {
		nleft[1] = nleft[0];
		nleft[0] = 0;
	    }
	    if (nright[0] > nright[1]) {
		int tmp = nright[0]; nright[0] = nright[1]; nright[1] = tmp;
	    }
	}
	

	static const int tet_edge[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
	// the 6 edges of the tetrahedron, defined by the node indices

	for (i = j = 0; i < 6; i++) {
	    if ((tet_edge[i][0] == nleft[0] && tet_edge[i][1] == nleft[1]) ||
		(tet_edge[i][0] == nright[0] && tet_edge[i][1] == nright[1]))
		continue; // we skip the 2 edges with nodes on the same side of the cutting plane
	    Intersect_line_plane (cut_orient, cut_pos,
				  mesh->nlist[pel->Node[tet_edge[i][0]]],
				  mesh->nlist[pel->Node[tet_edge[i][1]]],
				  isect[j++]);
	}
	return 4;

    }
}

// Split a tetrahedron that has one node on one side of the cutting plane,
// the other 3 on the other side, into 4 subtets.
// 'inode' is the index of the single node (0..3)
// isect is a list of 3 global points defining the intersections of the
// tet edges with the cutting plane
// On exit, the returned mesh has modified element el, added the 3 additional
// tetrahedra, and appended the additional nodes to its node list
// If mod_el_idx is set, it should point to an integer array of length >= 4,
// and will receive the indices of the 4 modified/added tetrahedra
// Return value is the number of modified/added elements (4)
int Tetsplit_1_3 (BufMesh *mesh, int el, int inode, const Point *isect,
		   int *mod_el_idx=0)
{
    int i, j;

    Element *pel[4];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4, "Currently only works with 4-noded tetrahedra");

    NodeList &nlist = mesh->nlist;

    if (mod_el_idx) {
	mod_el_idx[0] = el;
	for (i = 0; i < 3; i++)
	    mod_el_idx[i+1] = mesh->elen_used+i;
    }

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+3 > ebuf)
	for (i = 0; i < mesh->elen_used+3-ebuf; i++)
	    mesh->elist.Append(new Tetrahedron4);
    
    for (i = 1; i < 4; i++)
	pel[i] = mesh->elist[mesh->elen_used++];

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+3 > nbuf) {
	nlist.Append(mesh->nlen_used+3-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 3;

    int nidx[7];
    for (i = 0; i < 4; i++)
	nidx[i] = pel[0]->Node[i];
    for (i = 0; i < 3; i++)
	nidx[i+4] = nlen+i;

    for (i = 0; i < 3; i++)
	nlist[nlen+i] = isect[i];

    switch (inode) {
    case 0:
	pel[0]->Node[0] = nidx[0];
	pel[0]->Node[1] = nidx[4];
	pel[0]->Node[2] = nidx[5];
	pel[0]->Node[3] = nidx[6];

	pel[1]->Node[0] = nidx[4];
	pel[1]->Node[1] = nidx[1];
	pel[1]->Node[2] = nidx[5];
	pel[1]->Node[3] = nidx[6];
	
	pel[2]->Node[0] = nidx[1];
	pel[2]->Node[1] = nidx[2];
	pel[2]->Node[2] = nidx[5];
	pel[2]->Node[3] = nidx[3];
	
	pel[3]->Node[0] = nidx[6];
	pel[3]->Node[1] = nidx[1];
	pel[3]->Node[2] = nidx[5];
	pel[3]->Node[3] = nidx[3];
	break;
    case 1:
	pel[0]->Node[0] = nidx[1];
	pel[0]->Node[1] = nidx[5];
	pel[0]->Node[2] = nidx[4];
	pel[0]->Node[3] = nidx[6];

	pel[1]->Node[0] = nidx[5];
	pel[1]->Node[1] = nidx[2];
	pel[1]->Node[2] = nidx[4];
	pel[1]->Node[3] = nidx[6];
	
	pel[2]->Node[0] = nidx[2];
	pel[2]->Node[1] = nidx[0];
	pel[2]->Node[2] = nidx[4];
	pel[2]->Node[3] = nidx[3];
	
	pel[3]->Node[0] = nidx[6];
	pel[3]->Node[1] = nidx[2];
	pel[3]->Node[2] = nidx[4];
	pel[3]->Node[3] = nidx[3];
	break;
    case 2:
	pel[0]->Node[0] = nidx[2];
	pel[0]->Node[1] = nidx[4];
	pel[0]->Node[2] = nidx[5];
	pel[0]->Node[3] = nidx[6];

	pel[1]->Node[0] = nidx[4];
	pel[1]->Node[1] = nidx[0];
	pel[1]->Node[2] = nidx[5];
	pel[1]->Node[3] = nidx[6];
	
	pel[2]->Node[0] = nidx[0];
	pel[2]->Node[1] = nidx[1];
	pel[2]->Node[2] = nidx[5];
	pel[2]->Node[3] = nidx[3];
	
	pel[3]->Node[0] = nidx[6];
	pel[3]->Node[1] = nidx[0];
	pel[3]->Node[2] = nidx[5];
	pel[3]->Node[3] = nidx[3];
	break;
    case 3:
	pel[0]->Node[0] = nidx[3];
	pel[0]->Node[1] = nidx[5];
	pel[0]->Node[2] = nidx[4];
	pel[0]->Node[3] = nidx[6];

	pel[1]->Node[0] = nidx[5];
	pel[1]->Node[1] = nidx[1];
	pel[1]->Node[2] = nidx[4];
	pel[1]->Node[3] = nidx[6];
	
	pel[2]->Node[0] = nidx[1];
	pel[2]->Node[1] = nidx[0];
	pel[2]->Node[2] = nidx[4];
	pel[2]->Node[3] = nidx[2];
	
	pel[3]->Node[0] = nidx[6];
	pel[3]->Node[1] = nidx[1];
	pel[3]->Node[2] = nidx[4];
	pel[3]->Node[3] = nidx[2];
	break;
    }
    return 4;
}

// Split a tetrahedron taht has one node on one side, one node on 
// the cutting plane, and two nodes on the other side of the
// cutting plane, into 3 subtets.
int Tetsplit_1_2 (BufMesh *mesh, int el, int inode, int skipnode,
		  const Point *isect, int *mod_el_idx=0)
{
    int i;
    Element *pel[3];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4, "Currently only works with 4-noded tetrahedra");

    NodeList &nlist = mesh->nlist;

    if (mod_el_idx) {
	mod_el_idx[0] = el;
	for (i = 0; i < 2; i++)
	    mod_el_idx[i+1] = mesh->elen_used+i;
    }

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+2 > ebuf)
	for (i = 0; i < mesh->elen_used+2-ebuf; i++)
	    mesh->elist.Append(new Tetrahedron4);
    
    for (i = 1; i < 3; i++)
	pel[i] = mesh->elist[mesh->elen_used++];

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+3 > nbuf) {
	nlist.Append(mesh->nlen_used+3-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 3;

    int nidx[7];
    for (i = 0; i < 4; i++)
	nidx[i] = pel[0]->Node[i];
    for (i = 0; i < 3; i++)
	nidx[i+4] = nlen+i;

    for (i = 0; i < 3; i++)
	nlist[nlen+i] = isect[i];

    switch (inode) {
    case 0:
	pel[0]->Node[0] = nidx[0];
	pel[0]->Node[1] = nidx[4];
	pel[0]->Node[2] = nidx[5];
	pel[0]->Node[3] = nidx[6];
	switch (skipnode) {
	case 1:
	    pel[1]->Node[0] = nidx[5];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    pel[2]->Node[0] = nidx[6];
	    pel[2]->Node[1] = nidx[4];
	    pel[2]->Node[2] = nidx[5];
	    pel[2]->Node[3] = nidx[3];
	    break;
	case 2:
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[5];
	    pel[1]->Node[3] = nidx[3];
	    pel[2]->Node[0] = nidx[6];
	    pel[2]->Node[1] = nidx[4];
	    pel[2]->Node[2] = nidx[5];
	    pel[2]->Node[3] = nidx[3];
	    break;
	case 3:
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[6];
	    pel[2]->Node[0] = nidx[5];
	    pel[2]->Node[1] = nidx[4];
	    pel[2]->Node[2] = nidx[2];
	    pel[2]->Node[3] = nidx[6];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    case 1:
	pel[0]->Node[0] = nidx[1];
	pel[0]->Node[1] = nidx[5];
	pel[0]->Node[2] = nidx[4];
	pel[0]->Node[3] = nidx[6];
	switch (skipnode) {
	case 0:
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[5];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[6];
	    pel[2]->Node[0] = nidx[4];
	    pel[2]->Node[1] = nidx[6];
	    pel[2]->Node[2] = nidx[2];
	    pel[2]->Node[3] = nidx[3];
	    break;
	case 2:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[5];
	    pel[1]->Node[3] = nidx[3];
	    pel[2]->Node[0] = nidx[6];
	    pel[2]->Node[1] = nidx[5];
	    pel[2]->Node[2] = nidx[4];
	    pel[2]->Node[3] = nidx[3];
	    break;
	case 3:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[6];
	    pel[2]->Node[0] = nidx[4];
	    pel[2]->Node[1] = nidx[5];
	    pel[2]->Node[2] = nidx[2];
	    pel[2]->Node[3] = nidx[6];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    case 2:
	pel[0]->Node[0] = nidx[2];
	pel[0]->Node[1] = nidx[4];
	pel[0]->Node[2] = nidx[5];
	pel[0]->Node[3] = nidx[6];
	switch (skipnode) {
	case 0:
	    pel[1]->Node[0] = nidx[1];
	    pel[1]->Node[1] = nidx[5];
	    pel[1]->Node[2] = nidx[4];
	    pel[1]->Node[3] = nidx[6];
	    pel[2]->Node[0] = nidx[1];
	    pel[2]->Node[1] = nidx[6];
	    pel[2]->Node[2] = nidx[4];
	    pel[2]->Node[3] = nidx[3];
	    break;
	case 1:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[5];
	    pel[1]->Node[2] = nidx[4];
	    pel[1]->Node[3] = nidx[3];
	    pel[2]->Node[0] = nidx[5];
	    pel[2]->Node[1] = nidx[6];
	    pel[2]->Node[2] = nidx[4];
	    pel[2]->Node[3] = nidx[3];
	    break;
	case 3:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[5];
	    pel[1]->Node[3] = nidx[6];
	    pel[2]->Node[0] = nidx[0];
	    pel[2]->Node[1] = nidx[5];
	    pel[2]->Node[2] = nidx[4];
	    pel[2]->Node[3] = nidx[3];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    case 3:
	pel[0]->Node[0] = nidx[3];
	pel[0]->Node[1] = nidx[5];
	pel[0]->Node[2] = nidx[4];
	pel[0]->Node[3] = nidx[6];
	switch (skipnode) {
	case 0:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[5];
	    pel[2]->Node[0] = nidx[0];
	    pel[2]->Node[1] = nidx[5];
	    pel[2]->Node[2] = nidx[2];
	    pel[2]->Node[3] = nidx[6];
	    break;
	case 1:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[5];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[4];
	    pel[2]->Node[0] = nidx[4];
	    pel[2]->Node[1] = nidx[5];
	    pel[2]->Node[2] = nidx[2];
	    pel[2]->Node[3] = nidx[6];
	    break;
	case 2:
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[6];
	    pel[1]->Node[3] = nidx[4];
	    pel[2]->Node[0] = nidx[4];
	    pel[2]->Node[1] = nidx[1];
	    pel[2]->Node[2] = nidx[6];
	    pel[2]->Node[3] = nidx[5];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    }
    return 3;
}

int Tetsplit_1_1 (BufMesh *mesh, int el, int inode, int onode,
		  const Point *isect, int *mod_el_idx=0)
{
    int i;
    Element *pel[3];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4, "Currently only works with 4-noded tetrahedra");

    NodeList &nlist = mesh->nlist;

    if (mod_el_idx) {
	mod_el_idx[0] = el;
	mod_el_idx[1] = mesh->elen_used;
    }

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+1 > ebuf)
	for (i = 0; i < mesh->elen_used+1-ebuf; i++)
	    mesh->elist.Append(new Tetrahedron4);
    
    for (i = 1; i < 2; i++)
	pel[i] = mesh->elist[mesh->elen_used++];

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+1 > nbuf) {
	nlist.Append(mesh->nlen_used+1-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 1;

    int nidx[5];
    for (i = 0; i < 4; i++)
	nidx[i] = pel[0]->Node[i];
    nidx[4] = nlen;

    switch (inode) {
    case 0:
	switch (onode) {
	case 1:
	    nlist[nlen] = isect[0];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[4];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 2:
	    nlist[nlen] = isect[1];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[4];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 3:
	    nlist[nlen] = isect[2];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[4];
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    case 1:
	switch (onode) {
	case 0:
	    nlist[nlen] = isect[0];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[4];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 2:
	    nlist[nlen] = isect[1];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[4];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 3:
	    nlist[nlen] = isect[2];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[4];
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    case 2:
	switch (onode) {
	case 0:
	    nlist[nlen] = isect[0];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[4];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 1:
	    nlist[nlen] = isect[1];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[4];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 3:
	    nlist[nlen] = isect[2];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[4];
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[4];
	    pel[1]->Node[3] = nidx[3];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    case 3:
	switch (onode) {
	case 0:
	    nlist[nlen] = isect[0];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[4];
	    pel[1]->Node[0] = nidx[4];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 1:
	    nlist[nlen] = isect[1];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[3];
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[4];
	    pel[1]->Node[2] = nidx[2];
	    pel[1]->Node[3] = nidx[3];
	    break;
	case 2:
	    nlist[nlen] = isect[2];
	    pel[0]->Node[0] = nidx[0];
	    pel[0]->Node[1] = nidx[1];
	    pel[0]->Node[2] = nidx[2];
	    pel[0]->Node[3] = nidx[4];
	    pel[1]->Node[0] = nidx[0];
	    pel[1]->Node[1] = nidx[1];
	    pel[1]->Node[2] = nidx[4];
	    pel[1]->Node[3] = nidx[3];
	    break;
	default:
	    xERROR("Inconsistent node numbers");
	    break;
	}
	break;
    }
    return 2;
}


// Split a tetrahedron that has two nodes on either side of the cutting plane,
// into 6 subtets.
// 'inode' is the list of two node indices on the same side of the cutting
// plane.
// isect is a list of 4 global points defining the intersections of the
// tet edges with the cutting plane
// On exit, the returned mesh has modified element el, added the 5 additional
// tetrahedra, and appended the additional nodes to its node list
// If mod_el_idx is set, it should point to an integer array of length >= 6,
// and will receive the indices of the 6 modified/added tetrahedra
// Return value is the number of modified/added elements (6)

int Tetsplit_2_2 (BufMesh *mesh, int el, int *inode, const Point *isect,
		   int *mod_el_idx=0)
{
    int i, j, k;
    int onode[2];
    for (i = j = 0; i < 4; i++) {
	for (k = 0; k < 2; k++)
	    if (inode[k] == i) break;
	if (k == 2)
	    onode[j++] = i;
    }
	
    Element *pel[6];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4, "Currently only works with 4-noded tetrahedra");
    dASSERT(inode[0] == 0, "Inconsistent node order");

    NodeList &nlist = mesh->nlist;

    if (mod_el_idx) {
	mod_el_idx[0] = el;
	for (i = 0; i < 5; i++)
	    mod_el_idx[i+1] = mesh->elen_used+i;
    }

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+5 > ebuf)
	for (i = 0; i < mesh->elen_used+7-ebuf; i++)
	    mesh->elist.Append (new Tetrahedron4);

    for (i = 1; i < 6; i++)
	pel[i] = mesh->elist[mesh->elen_used++];

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+4 > nbuf) {
	nlist.Append(mesh->nlen_used+4-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 4;

    int nidx[8];
    for (i = 0; i < 4; i++)
	nidx[i] = pel[0]->Node[i];
    for (i = 0; i < 4; i++)
	nidx[i+4] = nlen+i;

    for (i = 0; i < 4; i++)
	nlist[nlen+i] = isect[i];

    int local_nidx_a[6][4] = { // split 01-23
	{5,7,4,3},{4,7,6,2},{3,7,4,2},
	{5,7,0,4},{4,7,1,6},{7,1,0,4}};
    int local_nidx_b[6][4] = { // split 02-13
	{4,6,5,1},{5,6,7,3},{1,6,5,3},
	{4,6,0,5},{5,6,2,7},{6,2,0,5}};
    int local_nidx_c[6][4] = { // split 03-12
	{5,7,4,2},{4,7,6,1},{2,7,4,1},
	{5,7,0,4},{4,7,3,6},{7,3,0,4}};
    int (*local_nidx)[4] = (inode[1] == 1 ? local_nidx_a :
			       inode[1] == 2 ? local_nidx_b :
			       local_nidx_c);

#ifdef UNDEF
    int local_nidx_a[6][4] = {
	{inode[0],inode[1],6,7},{5,inode[0],4,7},{4,inode[0],6,7},
	{onode[1],5,4,7},{onode[1],4,onode[0],7},{4,6,onode[0],7}
    };
    int local_nidx_b[6][4] = {
	{inode[0],inode[1],7,6},{4,inode[0],5,6},{5,inode[0],7,6},
	{onode[1],4,5,6},{onode[1],6,onode[0],4},{5,7,onode[1],6}
    };

    bool case_a = (inode[0]!=0 || inode[1]!=3) && (inode[0]!=1 || inode[1]!=3);
    int (&local_nidx)[6][4] = (case_a ? local_nidx_a : local_nidx_b);
#endif

    for (i = 0; i < 6; i++)
	for (j = 0; j < 4; j++)
	    pel[i]->Node[j] = nidx[local_nidx[i][j]];

    return 6;
}


// split a tetrahedron in the mesh along a cutting plane.
// cut_orient: orientation of the cutting plane (0: x=const, 1: y=const, 2: z=const)
// cut_pos: position of the cutting plane along the corresponding orientation
// The voxel cell for each of the resulting sub-tetrahedra is encoded in their
// region label. Tetrahedra left of the cutting plane retain their region value.
// Tetrahedra right of the cutting plane get their region index along the
// corresponding slice orientation incremented by one.
// The i,j,k cell grid indices are encoded in the region index as follows:
// reg = i + j<<10 + k<<20.
// i,j,k can thus have values between 0 and 1023.
// i,j,k are indices for the sub-cell grid bounding the given element, rather
// than absolute values

void Tetsplit (BufMesh *mesh, int el, int cut_orient, double cut_pos)
{
    Element *pel = mesh->elist[el];
    int regidx = pel->Region();

    // compute the intersection points
    int i, nisect, nelidx;
    int elidx[8];
    static Point isect[4];
    if (!isect[0].Dim())
	for (i = 0; i < 4; i++) isect[i].New(3);
    nisect = CalcIntersections_tet (mesh, el, cut_orient, cut_pos, isect);
    
    if (!nisect) { // no intersections: tet is either fully left or fully right

	if (Tet_extends_right (mesh, el, cut_orient, cut_pos)) // increment slice index
	    pel->SetRegion(regidx + (1 << (cut_orient*10)));

    } else if (nisect == 3) {

	// find the single node to be cut
	int singlend, leftnd, rightnd, nleft, nright, nflat, flatnd;
	for (i = nleft = nright = nflat = 0; i < 4; i++)
	    if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		nleft++; leftnd = i;
	    } else if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos)) {
		nright++; rightnd = i;
	    } else {
		nflat++; flatnd = i;
	    }
	singlend = (nleft == 1 ? leftnd : rightnd);

	if (nflat == 0) {
	    nelidx = Tetsplit_1_3 (mesh, el, singlend, isect, elidx);
	} else if (nflat == 1) {
	    nelidx = Tetsplit_1_2 (mesh, el, singlend, flatnd, isect, elidx);
	} else if (nflat == 2) {
	    int oppnd = (nleft == 1 ? rightnd : leftnd);
	    nelidx = Tetsplit_1_1 (mesh, el, singlend, oppnd, isect, elidx);
	} else {
	    xERROR("Case not supported yet.");
	}

	// modify slice indices for all sub-tetrahedra that ended up right of the cutting plane
	for (i = 0; i < nelidx; i++)
	    if (Tet_extends_right (mesh, elidx[i], cut_orient, cut_pos)) {
		mesh->elist[elidx[i]]->SetRegion (regidx + (1 << (cut_orient*10)));
	    } else {
		mesh->elist[elidx[i]]->SetRegion (regidx);
	    }
	
    } else if (nisect == 4) {

	// find the two nodes each left and right of the cutting plane
	int leftnd[2], rightnd[2], nleft, nright;
	for (i = nleft = nright = 0; i < 4; i++) {
	    if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
		leftnd[nleft++] = i;
	    else
		rightnd[nright++] = i;
	}

	int *grp0 = (leftnd[0] == 0 || leftnd[1] == 0 ? leftnd : rightnd);
	nelidx = Tetsplit_2_2 (mesh, el, grp0, isect, elidx);

	// modify slice indices for all sub-tetrahedra that ended up right of the cutting plane
	for (i = 0; i < nelidx; i++)
	    if (Tet_extends_right (mesh, elidx[i], cut_orient, cut_pos)) {
		mesh->elist[elidx[i]]->SetRegion (regidx + (1 << (cut_orient*10)));
	    } else {
		mesh->elist[elidx[i]]->SetRegion (regidx);
	    }
    }
}


// Add 6 tetrahedra to the mesh subdividing a voxel cell given by its
// extents
void Tetsplit_cube (BufMesh *mesh, double xmin, double xmax,
		    double ymin, double ymax, double zmin, double zmax,
		    int reg)
{
    int i;

    // grow element list if required
    int ebuf = mesh->elen();
    if (mesh->elen_used+6 > ebuf)
	for (i = 0; i < mesh->elen_used+6-ebuf; i++)
	    mesh->elist.Append (new Tetrahedron4);

    Element *pel[6];
    for (i = 0; i < 6; i++) {
	pel[i] = mesh->elist[mesh->elen_used++];
	pel[i]->SetRegion (reg);
    }

    // grow node list if required
    int nbuf = mesh->nlen();
    if (mesh->nlen_used+8 > nbuf) {
	mesh->nlist.Append(mesh->nlen_used+8-nbuf);
	for (i = mesh->nlen_used; i < mesh->nlen(); i++)
	    mesh->nlist[i].New(3);
    }
    int nlen = mesh->nlen_used;
    mesh->nlen_used += 6;

    for (i = 0; i < 8; i++) {
	Node &nd = mesh->nlist[nlen+i];
	nd[0] = (i%2 ? xmax:xmin);
	nd[1] = ((i/2)%2 ? ymax:ymin);
	nd[2] = (i/4 ? zmax:zmin);
    }

    pel[0]->Node[0] = nlen+0;
    pel[0]->Node[1] = nlen+1;
    pel[0]->Node[2] = nlen+2;
    pel[0]->Node[3] = nlen+5;

    pel[1]->Node[0] = nlen+4;
    pel[1]->Node[1] = nlen+6;
    pel[1]->Node[2] = nlen+5;
    pel[1]->Node[3] = nlen+2;

    pel[2]->Node[0] = nlen+0;
    pel[2]->Node[1] = nlen+4;
    pel[2]->Node[2] = nlen+5;
    pel[2]->Node[3] = nlen+2;

    pel[3]->Node[0] = nlen+3;
    pel[3]->Node[1] = nlen+2;
    pel[3]->Node[2] = nlen+1;
    pel[3]->Node[3] = nlen+6;

    pel[4]->Node[0] = nlen+7;
    pel[4]->Node[1] = nlen+5;
    pel[4]->Node[2] = nlen+6;
    pel[4]->Node[3] = nlen+1;

    pel[5]->Node[0] = nlen+3;
    pel[5]->Node[1] = nlen+7;
    pel[5]->Node[2] = nlen+6;
    pel[5]->Node[3] = nlen+1;
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
