
// =========================================================================
// =========================================================================
// Tetrahedra splitting code. This is preliminary and should eventually
// be shifted elsewhere (e.g. to the individual element classes)
// =========================================================================
// =========================================================================

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "tetsplit.h"

// ==========================================================================

BufMesh::BufMesh(): Mesh()
{
    nlen_used = elen_used = 0;
}

BufMesh::BufMesh (const Mesh &mesh): Mesh(mesh)
{
    nlen_used = nlen();
    elen_used = elen();
    for (int i = 0; i < elen_used; i++)
	elist[i]->SetRegion(0);  // mark all elements active
}

void BufMesh::SubSetup()
{
    for (int el=0; el < elen_used; el++)
	if (elist[el]->Region() >= 0)
	    elist[el]->Initialise(nlist);
    for (int el=0; el < elen_used; el++)
	if (elist[el]->Region() >= 0)
	    elist[el]->PostInitialisation(nlist);
}

void BufMesh::Copy (const BufMesh &mesh)
{
    nlist = mesh.nlist;
    elist = mesh.elist;
    nlen_used = mesh.nlen_used;
    elen_used = mesh.elen_used;
    for (int i = 0; i < elen_used; i++)
	elist[i]->SetRegion (mesh.elist[i]->Region());
    SubSetup();
}

void BufMesh::Shrink()
{
    int el, i, j, elen_used_new;

    for (el = 0; el < elen_used; el++) {
	if (elist[el]->Region()) {
	    for (i = el+1; i < elen_used; i++) {
		if (elist[i]->Region() == 0) {
		    for (j = 0; j < elist[el]->nNode(); j++)
			elist[el]->Node[j] = elist[i]->Node[j];
		    elist[el]->SetRegion(0);
		    elist[i]->SetRegion(-1);
		    elist[el]->Initialise (nlist);
		    break;
		}
	    }
	}
    }
    for (elen_used_new = el = 0; el < elen_used; el++)
	if (elist[el]->Region() == 0) elen_used_new++;
    elen_used = elen_used_new;
    // also need to shrink the node list
}

void BufMesh::Rotate (const RDenseMatrix &R)
{
    for (int i = 0; i < nlen_used; i++) {
	nlist[i] = R * nlist[i];
    }
}

std::ostream& operator<< (std::ostream& o, BufMesh &mesh)
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


// =========================================================================

int CalcIntersections_tet(Mesh *mesh, int el, int cut_orient, double cut_pos,
    Point *isect)
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

// =========================================================================

int Tetsplit_1_3 (BufMesh *mesh, int el, int inode, const Point *isect,
    int *mod_el_idx)
{
    int i;

    Element *pel[4];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4,
	    "Currently only works with 4-noded tetrahedra");

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

// =========================================================================

int Tetsplit_1_2 (BufMesh *mesh, int el, int inode, int skipnode,
    const Point *isect, int *mod_el_idx)
{
    int i;
    Element *pel[3];
    pel[0] = mesh->elist[el];

    xASSERT(pel[0]->Type() == ELID_TET4,
	    "Currently only works with 4-noded tetrahedra");

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

// =========================================================================

int Tetsplit_1_1 (BufMesh *mesh, int el, int inode, int onode,
    const Point *isect, int *mod_el_idx)
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

// =========================================================================

int Tetsplit_2_2 (BufMesh *mesh, int el, int *inode, const Point *isect,
    int *mod_el_idx)
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

    xASSERT(pel[0]->Type() == ELID_TET4,
	    "Currently only works with 4-noded tetrahedra");
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

// =========================================================================

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

	if (Tet_extends_right (mesh, el, cut_orient, cut_pos))
	    // increment slice index
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

// =========================================================================

void Tetsplit_cube (BufMesh *mesh, double xmin, double xmax,
    double ymin, double ymax, double zmin, double zmax, int reg)
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
