// =======================================================================
// Converts  TET10 -> TET4
// By replacing each 10-noded tetrahedron with 1 or 8 4-noded ones
// Other element types are not supported
// =======================================================================

#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

void Quad2Lin (Mesh &mesh);
void Quad2LinRefine (Mesh &mesh);
void Usage();

int main (int argc, char *argv[])
{
  //CHECK_EXPIRED();
    
    int i;
    bool do_stdin = true;
    bool do_stdout = true;
    bool do_refine = false;
    char imeshname[256], omeshname[256];

    Mesh mesh;

    for (i = 1; i < argc; i++) {
        xASSERT (argv[i][0] == '-', "Error parsing command line");
	switch (argv[i][1]) {
	case 'i':
	    do_stdin = false;
	    strcpy (imeshname, argv[++i]);
	    break;
	case 'o':
	    do_stdout = false;
	    strcpy (omeshname, argv[++i]);
	    break;
	case 'r':
	    do_refine = true;
	    break;
	case '?':
	case 'h':
	case 'H':
	    Usage();
	    exit(0);
	default:
	    xERROR("Error parsing command line");
	    break;
	}
    }

    if (do_stdin) {
        cin >> mesh;
    } else {
        ifstream ifs (imeshname);
	ifs >> mesh;
    }
    mesh.Setup();

    if (do_refine)
        Quad2LinRefine (mesh);
    else
        Quad2Lin (mesh);

    if (do_stdout) {
        cout << mesh;
    } else {
        ofstream ofs (omeshname);
	ofs << mesh;
    }

    return 0;
}

void Quad2Lin (Mesh &mesh)
{
    int i, j, nn2 = 0, nn = mesh.nlen(), ne = mesh.elen();
    
    // pass 1: find out how many nodes are required
    bool *nused = new bool[nn];
    for (i = 0; i < nn; i++) nused[i] = false;

    for (i = 0; i < ne; i++) {
        Element *pel = mesh.elist[i];
	switch (pel->Type()) {
	case ELID_TET10:
	    for (j = 0; j < 4; j++)
	        nused[pel->Node[j]] = true;
	    break;
	case ELID_TRI6:
	case ELID_TRI10:
	case ELID_TRI3D3:
	    for (j = 0; j < 3; j++)
	        nused[pel->Node[j]] = true;
	    break;
	default:
	    cerr << "Element type not supported" << endl;
	    exit (1);
	}
    }
    for (i = 0; i < nn; i++)
        if (nused[i]) nn2++;
    cerr << "New mesh has " << nn2 << " nodes" << endl;

    int *nmap = new int[nn];

    // construct new node and element lists
    Node *nlist2 = new Node[nn2];
    Element **elist2 = new Element*[ne];

    for (i = j = 0; i < nn; i++) {
        if (nused[i]) {
	    nlist2[j].New(mesh.nlist[i].Dim());
	    nlist2[j] = mesh.nlist[i];
	    nmap[i] = j++;
	} else nmap[i] = -1;
    }
    for (i = 0; i < ne; i++) {
        switch (mesh.elist[i]->Type()) {
	case ELID_TET10:
	    elist2[i] = new Tetrahedron4;
	    for (j = 0; j < 4; j++)
	        elist2[i]->Node[j] = nmap[mesh.elist[i]->Node[j]];
	    break;
	case ELID_TRI6:
	case ELID_TRI10:
	    //elist2[i] = new Triangle3old;
	    //for (j = 0; j < 3; j++)
	    //    elist2[i]->Node[j] = nmap[mesh.elist[i]->Node[j]];
	    elist2[i] = new Triangle3;
	    for (j = 0; j < 3; j++)
	        elist2[i]->Node[j] = nmap[mesh.elist[i]->Node[j]];
	    break;
	case ELID_TRI3D6:
	    elist2[i] = new Triangle3D3;
	    for (j = 0; j < 3; j++)
	        elist2[i]->Node[j] = nmap[mesh.elist[i]->Node[j]];
	    break;
	default:
	    cerr << "Element type not supported" << endl;
	    exit (1);
	}
    }
    mesh.elist.SetList (ne, elist2);
    mesh.nlist.SetList (nn2, nlist2);
    delete []nused;
    delete []nmap;
}

void Quad2LinRefine (Mesh &mesh)
{
    int i, els = 0;
    int nd0, nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8, nd9;

    // pass 1: find out how many new elements we need
    for (i = 0; i < mesh.elen(); i++) {
        switch (mesh.elist[i]->Type()) {
	case ELID_TRI6:
	case ELID_TRI3D6:
  	    els += 4; // split quadratic triangles into 4 linear ones
	    break;
	case ELID_TET10:
	    els += 8; // split quadratic tetrahedra into 8 linear ones
	    break;
	default:
	    cerr << "Element type not supported" << endl;
	    exit(1);  // no other types are recognised
	}
    }
    Element **pel = new Element*[els];

    // pass 2: construct new elements
    els = 0;
    for (i = 0; i < mesh.elen(); i++) {

        switch (mesh.elist[i]->Type()) {

	case ELID_TET10:
	    nd0 = mesh.elist[i]->Node[0];
	    nd1 = mesh.elist[i]->Node[1];
	    nd2 = mesh.elist[i]->Node[2];
	    nd3 = mesh.elist[i]->Node[3];
	    nd4 = mesh.elist[i]->Node[4];
	    nd5 = mesh.elist[i]->Node[5];
	    nd6 = mesh.elist[i]->Node[6];
	    nd7 = mesh.elist[i]->Node[7];
	    nd8 = mesh.elist[i]->Node[8];
	    nd9 = mesh.elist[i]->Node[9];

	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd0;
	    pel[els]->Node[1] = nd4;
	    pel[els]->Node[2] = nd5;
	    pel[els]->Node[3] = nd6;
	    els++;

	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd4;
	    pel[els]->Node[1] = nd1;
	    pel[els]->Node[2] = nd7;
	    pel[els]->Node[3] = nd8;
	    els++;
	    
	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd5;
	    pel[els]->Node[1] = nd7;
	    pel[els]->Node[2] = nd2;
	    pel[els]->Node[3] = nd9;
	    els++;
	    
	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd6;
	    pel[els]->Node[1] = nd8;
	    pel[els]->Node[2] = nd9;
	    pel[els]->Node[3] = nd3;
	    els++;

	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd6;
	    pel[els]->Node[1] = nd4;
	    pel[els]->Node[2] = nd5;
	    pel[els]->Node[3] = nd8;
	    els++;

	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd4;
	    pel[els]->Node[1] = nd7;
	    pel[els]->Node[2] = nd5;
	    pel[els]->Node[3] = nd8;
	    els++;

	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd5;
	    pel[els]->Node[1] = nd7;
	    pel[els]->Node[2] = nd9;
	    pel[els]->Node[3] = nd8;
	    els++;

	    pel[els] = new Tetrahedron4;
	    pel[els]->Node[0] = nd5;
	    pel[els]->Node[1] = nd9;
	    pel[els]->Node[2] = nd6;
	    pel[els]->Node[3] = nd8;
	    els++;
	    break;

	case ELID_TRI6:
	    nd0 = mesh.elist[i]->Node[0];
	    nd1 = mesh.elist[i]->Node[1];
	    nd2 = mesh.elist[i]->Node[2];
	    nd3 = mesh.elist[i]->Node[3];
	    nd4 = mesh.elist[i]->Node[4];
	    nd5 = mesh.elist[i]->Node[5];

	    pel[els] = new Triangle3;
	    pel[els]->Node[0] = nd0;
	    pel[els]->Node[1] = nd3;
	    pel[els]->Node[2] = nd5;
	    els++;

	    pel[els] = new Triangle3;
	    pel[els]->Node[0] = nd3;
	    pel[els]->Node[1] = nd1;
	    pel[els]->Node[2] = nd4;
	    els++;

	    pel[els] = new Triangle3;
	    pel[els]->Node[0] = nd3;
	    pel[els]->Node[1] = nd4;
	    pel[els]->Node[2] = nd5;
	    els++;

	    pel[els] = new Triangle3;
	    pel[els]->Node[0] = nd5;
	    pel[els]->Node[1] = nd4;
	    pel[els]->Node[2] = nd2;
	    els++;
	    break;

	case ELID_TRI3D6:
	    nd0 = mesh.elist[i]->Node[0];
	    nd1 = mesh.elist[i]->Node[1];
	    nd2 = mesh.elist[i]->Node[2];
	    nd3 = mesh.elist[i]->Node[3];
	    nd4 = mesh.elist[i]->Node[4];
	    nd5 = mesh.elist[i]->Node[5];

	    pel[els] = new Triangle3D3;
	    pel[els]->Node[0] = nd0;
	    pel[els]->Node[1] = nd3;
	    pel[els]->Node[2] = nd5;
	    els++;

	    pel[els] = new Triangle3D3;
	    pel[els]->Node[0] = nd3;
	    pel[els]->Node[1] = nd1;
	    pel[els]->Node[2] = nd4;
	    els++;

	    pel[els] = new Triangle3D3;
	    pel[els]->Node[0] = nd3;
	    pel[els]->Node[1] = nd4;
	    pel[els]->Node[2] = nd5;
	    els++;

	    pel[els] = new Triangle3D3;
	    pel[els]->Node[0] = nd5;
	    pel[els]->Node[1] = nd4;
	    pel[els]->Node[2] = nd2;
	    els++;

	    break;
	}
    }
    mesh.elist.SetList (els, pel);
    //delete []pel;
    //mesh.Setup();
}

void Usage ()
{
    cout << endl << "quad2linmesh" << endl;
    cout << "Convert 2nd order to 1st order meshes" << endl;
    cout << endl << "Parameters:" << endl;
    cout << "-i <meshname>  input mesh name (default is stdin)" << endl;
    cout << "-o <meshname>  output mesh name (default is stdout)" << endl;
    cout << "-r             refine (keep all nodes, split elements)" << endl;
    cout << "               default is to remove mid nodes" << endl;
}
