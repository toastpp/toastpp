// -*-C++-*-

#ifndef __TETSPLIT_H
#define __TETSPLIT_H

// Mesh class with pre-allocated element and node lists, so we don't
// need to dynamically grow the lists too often
// nlen_used and elen_used are the length of the lists actually in use

class BufMesh: public Mesh
{
public:
    BufMesh ();
    BufMesh (const Mesh &mesh);
    void SubSetup();
    void Copy (const BufMesh &mesh);
    void Shrink();
    void Rotate (const RDenseMatrix &R);
    friend std::ostream& operator<< (std::ostream& o, BufMesh &mesh);

    int nlen_used, elen_used;
};


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
inline bool Tet_extends_left (Mesh *mesh, int el, int cut_orient,
    double cut_pos)
{
    Element *pel = mesh->elist[el];
    for (int i = 0; i < 4; i++)
	if (Point_is_left (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
	    return true;
    return false;
}


// Check if a tetrahedron has a portion right of a cutting plane
inline bool Tet_extends_right (Mesh *mesh, int el, int cut_orient,
    double cut_pos)
{
    Element *pel = mesh->elist[el];
    for (int i = 0; i < 4; i++)
	if (Point_is_right (mesh->nlist[pel->Node[i]], cut_orient, cut_pos))
	    return true;
    return false;
}


// Calculate the intersection point pi of a line defined by two points p0
// and p1, and a cutting plane parallel to two coordinate axes
// cut_orient defines the orientation of the cutting plane
// (0: x=const, 1: y=const, 2: z=const)
// Calculate the intersection points of the edges of tetrahedron el
// with a cutting plane parallel to two axes.
// cut_orient defines the orientation of the cutting plane
// (0: x=const, 1: y=const, 2: z=const)
// cut_pos is the position of the intersection of the cutting plane with
// the corresponding axis
// Return value is the number of intersection points found (0, 3, or 4)
// The points are written to the array pointed to by isect (must be allocated to
// at least 4 points of dimension 3)

int CalcIntersections_tet(Mesh *mesh, int el, int cut_orient, double cut_pos,
    Point *isect);


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
    int *mod_el_idx=0);


// Split a tetrahedron that has one node on one side, one node on 
// the cutting plane, and two nodes on the other side of the
// cutting plane, into 3 subtets.

int Tetsplit_1_2 (BufMesh *mesh, int el, int inode, int skipnode,
    const Point *isect, int *mod_el_idx=0);


// Split a tetrahedron that has two nodes on the cutting plane, and one node
// on each side of the cutting plane, into 2 subtets

int Tetsplit_1_1 (BufMesh *mesh, int el, int inode, int onode,
    const Point *isect, int *mod_el_idx=0);


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
    int *mod_el_idx=0);


// split a tetrahedron in the mesh along a cutting plane.
// cut_orient: orientation of the cutting plane (0: x=const, 1: y=const,
// 2: z=const)
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

void Tetsplit (BufMesh *mesh, int el, int cut_orient, double cut_pos);


// Add 6 tetrahedra to the mesh subdividing a voxel cell given by its
// extents

void Tetsplit_cube (BufMesh *mesh, double xmin, double xmax,
    double ymin, double ymax, double zmin, double zmax, int reg);

#endif // !__TETSPLIT_H
