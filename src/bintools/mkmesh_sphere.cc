#include<fstream>
#include<mathlib.h>
#include<felib.h>

using namespace std;

static const RDenseMatrix octa_nd = RDenseMatrix (7, 3,
   " 0  0  0 \
     1  0  0 \
     0  1  0 \
    -1  0  0 \
     0 -1  0 \
     0  0  1 \
     0  0 -1");

static const IDenseMatrix octa_el = IDenseMatrix (8, 4,
   "0 1 2 5 \
    0 2 3 5 \
    0 3 4 5 \
    0 4 1 5 \
    0 1 4 6 \
    0 4 3 6 \
    0 3 2 6 \
    0 2 1 6");

static const RDenseMatrix icosa_nd = RDenseMatrix (13, 3,
   "0 0 0 \
    0 0 1 \
    0.27639357928 0.85065117347 0.44721369208 \
    0.89442568286 0 0.44721369208 \
    -0.72360727137 0.52573134735 0.44721369208 \
    -0.72360727137 -0.52573134735 0.44721369208 \
    0.27639357928 -0.85065117347 0.44721369208 \
    0.72360727137 0.52573134735 -0.44721369208 \
    -0.27639357928 0.85065117347 -0.44721369208 \
    -0.89442568286 0 -0.44721369208 \
    -0.27639357928 -0.85065117347 -0.44721369208 \
    0.72360727137 -0.52573134735 -0.44721369208 \
    0 0 -1");

static const IDenseMatrix icosa_el = IDenseMatrix (20, 4,
 "0 1 3 2 \
  0 1 2 4 \
  0 1 4 5 \
  0 1 5 6 \
  0 1 6 3 \
  0 3 7 2 \
  0 2 8 4 \
  0 4 9 5 \
  0 5 10 6 \
  0 6 11 3 \
  0 7 8 2 \
  0 8 9 4 \
  0 9 10 5 \
  0 10 11 6 \
  0 11 7 3 \
  0 7 12 8 \
  0 8 12 9 \
  0 9 12 10 \
  0 10 12 11 \
  0 11 12 7");

static const int nb[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

int Exists (const Node &nd, const NodeList &nlist, int len)
{
    static const double eps = 1e-4;
    for (int i = 0; i < len; i++)
        if (fabs (Dist (nd, nlist[i])) < eps) return i;
    return -1;
}

int main (void)
{
    const double eps = 1e-4;

    double hmua, hmus, href;
    double rad, r1, r2, scale;
    int res;
    int i, j, el, hreg, r;
    int elen, nlen, nn;
    char fname[256];
    Mesh mesh;
    Tetrahedron4 *tet;
    Node newn(3);

    cout << "mkmesh_sphere: Create spherical 3D mesh using tetrahedral elements"
	 << endl << endl;

    cout << "Mesh name: ";
    cin  >> fname;
    cout << "Sphere radius: ";
    cin  >> rad;

    cout << "Resolution (>= 0): ";
    cin  >> res;

    cout << "Background parameters:" << endl;
    cout << "mua mus refind region: " ;
    cin  >> hmua >> hmus >> href >> hreg;

#ifdef UNDEF
    // build initial octagon
    mesh.nlist.New (7);
    for (i = 0; i < 7; i++) {
        mesh.nlist[i].New (3);
	for (j = 0; j < 3; j++) mesh.nlist[i][j] = octa_nd(i,j)*rad;
	mesh.nlist[i].SetBndTp (i ? BND_DIRICHLET : BND_NONE);
	mesh.nlist[i].SetRegion (hreg);
    }
    mesh.elist.Clear();
    Element **list = new Element*[8];
    for (i = 0; i < 8; i++) {
        list[i] = new Tetrahedron4;
	for (j = 0; j < 4; j++) list[i]->Node[j] = octa_el(i,j);
    }
    mesh.elist.AppendList (8, list);
#endif

    // build initial icosahedron
    mesh.nlist.New (13);
    for (i = 0; i < 13; i++) {
        mesh.nlist[i].New (3);
	for (j = 0; j < 3; j++) mesh.nlist[i][j] = icosa_nd(i,j)*rad;
	mesh.nlist[i].SetBndTp (i ? BND_DIRICHLET : BND_NONE);
	mesh.nlist[i].SetRegion (hreg);
    }
    mesh.elist.Clear();
    Element **list = new Element*[20];
    for (i = 0; i < 20; i++) {
        list[i] = new Tetrahedron4;
	for (j = 0; j < 4; j++) list[i]->Node[j] = icosa_el(i,j);
    }
    mesh.elist.AppendList (20, list);


    delete []list;
    mesh.Setup();

    // now refine elements
    for (r = 0; r < res; r++) {

        cout << "Starting refinement level " << r+1 << " of " << res << endl;
        nlen = mesh.nlen();
        elen = mesh.elen();
        NodeList newnl (elen*6); // max number of new nodes
	int newnl_len = 0;

	for (el = 0; el < elen; el++) {

	    int nno[10];
	    for (i = 0; i < 4; i++) nno[i] = mesh.elist[el]->Node[i];

	    for (nn = 0; nn < 6; nn++) {
	        Node &n1 = mesh.nlist[mesh.elist[el]->Node[nb[nn][0]]];
		Node &n2 = mesh.nlist[mesh.elist[el]->Node[nb[nn][1]]];
		r1 = length(n1);
		r2 = length(n2);
		for (i = 0; i < 3; i++) newn[i] = 0.5 * (n1[i]+n2[i]);
		if (fabs(r1-r2) < eps) { // push to sphere surface
		    scale = r1/length(newn);
		    for (i = 0; i < 3; i++) newn[i] *= scale;
		}
		i = Exists (newn, newnl, newnl_len);
		if (i >= 0) { // got this node already
		    nno[nn+4] = nlen+i;
		} else {
		    nno[nn+4] = nlen+newnl_len;
		    if (n1.BndTp() == n2.BndTp()) newn.SetBndTp (n1.BndTp());
		    else                          newn.SetBndTp (BND_NONE);
		    newnl[newnl_len++].Copy (newn);
		}
	    }

	    mesh.elist[el]->Node[0] = nno[0];
	    mesh.elist[el]->Node[1] = nno[4];
	    mesh.elist[el]->Node[2] = nno[5];
	    mesh.elist[el]->Node[3] = nno[6];

	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[4];
	    tet->Node[1] = nno[1];
	    tet->Node[2] = nno[7];
	    tet->Node[3] = nno[8];
	    mesh.elist.Append (tet);

	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[5];
	    tet->Node[1] = nno[7];
	    tet->Node[2] = nno[2];
	    tet->Node[3] = nno[9];
	    mesh.elist.Append (tet);

	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[6];
	    tet->Node[1] = nno[8];
	    tet->Node[2] = nno[9];
	    tet->Node[3] = nno[3];
	    mesh.elist.Append (tet);

	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[5];
	    tet->Node[1] = nno[4];
	    tet->Node[2] = nno[7];
	    tet->Node[3] = nno[6];
	    mesh.elist.Append (tet);

	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[5];
	    tet->Node[1] = nno[7];
	    tet->Node[2] = nno[9];
	    tet->Node[3] = nno[6];
	    mesh.elist.Append (tet);

	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[7];
	    tet->Node[1] = nno[4];
	    tet->Node[2] = nno[8];
	    tet->Node[3] = nno[6];
	    mesh.elist.Append (tet);
	    
	    tet = new Tetrahedron4;
	    tet->Node[0] = nno[9];
	    tet->Node[1] = nno[6];
	    tet->Node[2] = nno[7];
	    tet->Node[3] = nno[8];
	    mesh.elist.Append (tet);
	}

	mesh.nlist.Append (newnl_len);
	for (i = 0; i < newnl_len; i++)
	    mesh.nlist[nlen+i].Copy (newnl[i]);
	cout << "Mesh now has " << mesh.nlen() << " nodes, " << mesh.elen()
	     << " elements" << endl;
    }
    cout << "Finished refinement" << endl;

    mesh.Setup();
    cout << "Mesh initialised" << endl;

    // write mesh
    ofstream ofs(fname);
    ofs << mesh;
    cout << "Mesh written to " << fname << endl;
    cout << "Bye." << endl;

    return 0;
}
