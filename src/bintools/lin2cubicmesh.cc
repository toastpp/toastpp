// =======================================================================
// Converts TRI3 -> TRI10IP
// and      TET4 -> (not supported yet)
// Other element types are not supported
// =======================================================================

#include <mathlib.h>
#include <felib.h>

using namespace std;

void Convert2DMesh (Mesh &mesh);
void Convert3DMesh (Mesh &mesh);

int Exists (const Node &nd, const NodeList &nlist, int nn, double rad = 1e-8)
{
    // we assume without test that nn <= nlist.Dim()

    for (int i = 0; i < nn; i++) {
        if (Dist (nd, nlist[i]) < rad) return i;
    }
    return -1;
}

int main (void)
{
    CHECK_EXPIRED();

    Mesh mesh;

    cin >> mesh;
    mesh.Setup();
    int dimension = mesh.nlist[0].Dim();
    switch (dimension) {
    case 2:  Convert2DMesh (mesh); break;
    case 3:  Convert3DMesh (mesh); break;
    default: xERROR ("Invalid mesh dimension"); break;
    }

    mesh.Setup();
    cout << mesh;
    return 0;
}

void Convert2DMesh (Mesh &mesh)
{
    int i, j, el, sd, nlen, elen;
    int n0, n1, s, inew;
    int nd[10];

    elen = mesh.elen();
    nlen = mesh.nlen();

    NodeList nnlist (7*elen); // max number of new nodes
    int nnlen = 0;

    for (el = 0; el < elen; el++) {
        if (!(el % 1000)) cerr << "Processing el " << el << endl;
        Element *pel = mesh.elist[el];
	xASSERT (pel->Type() == ELID_TRI3OLD, "Unsupported element type");

	// check node orientation
	double x0 = mesh.nlist[pel->Node[0]][0];
	double y0 = mesh.nlist[pel->Node[0]][1];
	double x1 = mesh.nlist[pel->Node[1]][0];
	double y1 = mesh.nlist[pel->Node[1]][1];
	double x2 = mesh.nlist[pel->Node[2]][0];
	double y2 = mesh.nlist[pel->Node[2]][1];
	double size2 = x1*y2-x2*y1-x0*y2+x2*y0+x0*y1-x1*y0;
	if (size2 < 0.0) { // wrong node orientation
	    int tmp = pel->Node[1];
	    pel->Node[1] = pel->Node[2];
	    pel->Node[2] = tmp;
	}

	Element *qel = new Triangle10;
	inew = 0;
	for (sd = 0; sd < pel->nSide(); sd++) {
	    n0 = pel->Node[pel->SideNode (sd, 0)];
	    n1 = pel->Node[pel->SideNode (sd, 1)];
	    for (s = 0; s < 2; s++) {
	        Node newn(2);
		for (j = 0; j < 2; j++)
		    newn[j] = (2.0-s)/3.0 * mesh.nlist[n0][j] +
		              (1.0+s)/3.0 * mesh.nlist[n1][j];
		if ((nd[inew] = Exists (newn, nnlist, nnlen)) < 0) {
		    nd[inew] = nlen + nnlen;
		    nnlist[nnlen].New (2);
		    for (j = 0; j < 2; j++)
		        nnlist[nnlen][j] = newn[j];
		    if (mesh.nlist[n0].BndTp() == mesh.nlist[n1].BndTp())
		        nnlist[nnlen].SetBndTp (mesh.nlist[n0].BndTp());
		    else
		        nnlist[nnlen].SetBndTp (BND_NONE);
		    nnlen++;
		} else {
		  nd[inew] += nlen;
		}
		inew++;
	    }
	}

	// generate central node
	nnlist[nnlen].New(2);
	for (j = 0; j < 2; j++) {
	    nnlist[nnlen][j] = 1.0/3.0 * (mesh.nlist[pel->Node[0]][j] +
					  mesh.nlist[pel->Node[1]][j] +
					  mesh.nlist[pel->Node[2]][j]);
	}
	nd[inew++] = nlen + nnlen;
	nnlen++;
	
	qel->Node[0] = pel->Node[0];
	qel->Node[1] = pel->Node[1];
	qel->Node[2] = pel->Node[2];
	qel->Node[3] = nd[0];
	qel->Node[4] = nd[1];
	qel->Node[5] = nd[2];
	qel->Node[6] = nd[3];
	qel->Node[7] = nd[4];
	qel->Node[8] = nd[5];
	qel->Node[9] = nd[6];
	delete mesh.elist[el];
	mesh.elist[el] = qel;
    }
    if (nnlen) {
        mesh.nlist.Append (nnlen);
	for (i = 0; i < nnlen; i++) {
	    mesh.nlist[nlen+i].Copy (nnlist[i]);
	}
    }
}

void Convert3DMesh (Mesh &mesh)
{
    ERROR_UNDEF;
#ifdef UNDEF
    const int edgend[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    int el, elen, nlen, edge, i, j;
    int n0, n1, n2, nd[6];

    elen = mesh.elen();
    nlen = mesh.nlen();

    NodeList nnlist (6*elen); // max number of new nodes
    int nnlen = 0;

    for (el = 0; el < elen; el++) {
        if (!(el % 1000)) cerr << "Processing el " << el << endl;
        Element *pel = mesh.elist[el];
	xASSERT (pel->Type() == ELID_TET4, "Unsupported element type");

	Element *qel = new Tetrahedron10;
	for (edge = 0; edge < 6; edge++) {
	    n0 = pel->Node[edgend[edge][0]];
	    n1 = pel->Node[edgend[edge][1]];
	    Node newn(3);
	    for (j = 0; j < 3; j++)
	        newn[j] = 0.5 * (mesh.nlist[n0][j] + mesh.nlist[n1][j]);
	    if ((nd[edge] = Exists (newn, nnlist, nnlen)) < 0) {
	        nd[edge] = nlen + nnlen;
		nnlist[nnlen].New (3);
		for (j = 0; j < 3; j++)
		    nnlist[nnlen][j] = newn[j];
		if (mesh.nlist[n0].BndTp() == mesh.nlist[n1].BndTp())
		    nnlist[nnlen].SetBndTp (mesh.nlist[n0].BndTp());
		else
		    nnlist[nnlen].SetBndTp (BND_NONE);
		nnlen++;
	    } else {
	        nd[edge] += nlen;
	    }
	}
	for (i = 0; i < 4; i++)
	    qel->Node[i] = pel->Node[i];
	for (i = 0; i < 6; i++)
	    qel->Node[4+i] = nd[i];
	delete mesh.elist[el];
	mesh.elist[el] = qel;
    }
    if (nnlen) {
        mesh.nlist.Append (nnlen);
	for (i = 0; i < nnlen; i++) {
	    mesh.nlist[nlen+i].Copy (nnlist[i]);
	}
    }
#endif
}
