#include<fstream>
#include<mathlib.h>
#include<felib.h>

using namespace std;

int main (void)
{
    double hmua, hmus, href;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    int dimx, dimy, dimz, yofs, zofs, n0, n1, n2, n3, n4, n5, n6, n7;
    int i, j, k, n, el, nnode, nnode_slice, nels, reg, hreg, tmp;
    double x, y, z;
    char fname[256];
    Mesh mesh;

    cout << "mkmesh_brick: Create simple 3D mesh using tetrahedral elements"
	 << endl << endl;

    cout << "Mesh name: ";
    cin  >> fname;
    cout << "Physical dimensions:" << endl;
    cout << "xmin ymin zmin: ";
    cin  >> xmin >> ymin >> zmin;
    cout << "xmax ymax zmax: ";
    cin  >> xmax >> ymax >> zmax;

    cout << "Number of element cells in each dimension (each element" << endl;
    cout << "cell holds 6 tetrahedra):" << endl;
    cout << "x y z: ";
    cin  >> dimx >> dimy >> dimz;

    cout << "Background parameters:" << endl;
    cout << "mua mus refind region: " ;
    cin  >> hmua >> hmus >> href >> hreg;

    nnode_slice = (dimx+1) * (dimy+1);     // nodes per z-slice
    nnode       = nnode_slice * (dimz+1);  // total number of nodes
    nels        = dimx * dimy * dimz * 5;

    // create node list
    mesh.nlist.New (nnode);
    n = 0;
    for (k = 0; k <= dimz; k++) {
        z = (zmax-zmin)*(double)k/(double)dimz + zmin;
        for (j = 0; j <= dimy; j++) {
	    y = (ymax-ymin)*(double)j/(double)dimy + ymin;
	    for (i = 0; i <= dimx; i++) {
	        x = (xmax-xmin)*(double)i/(double)dimx + xmin;
		mesh.nlist[n].New (3);
	        mesh.nlist[n][0] = x;
		mesh.nlist[n][1] = y;
		mesh.nlist[n][2] = z;
		if (i == 0 || i == dimx || j == 0 || j == dimy ||
		    k == 0 || k == dimz)
		    mesh.nlist[n].SetBndTp (BND_DIRICHLET);
		else
		    mesh.nlist[n].SetBndTp (BND_NONE);
		n++;
	    }
	}
    }
    cout << "Generated node list (" << nnode << " nodes)" << endl;

    // create element list
    mesh.elist.Clear();
    Element **list = new Element*[nels];
    for (el = 0; el < nels; el++)
        list[el] = new Tetrahedron4;
    mesh.elist.AppendList (nels, list);
    delete []list;
    //for (el = 0; el < nels; el++)
    //  mesh.elist.Append (new Tetrahedron4);
    el = 0;
    for (k = 0; k < dimz; k++) {
        zofs = nnode_slice * k;
	for (j = 0; j < dimy; j++) {
	    yofs = (dimx+1) * j;
	    for (i = 0; i < dimx; i++) {
	        n0 = zofs + yofs + i;
		n1 = n0 + 1;
		n2 = n0 + (dimx+1);
		n3 = n2 + 1;
		n4 = n0 + nnode_slice;
		n5 = n4 + 1;
		n6 = n4 + (dimx+1);
		n7 = n6 + 1;
	        if ((i+j+k) & 1) {
		    // rotate element arrangement in every other cube 
		    // to make element edges align across cubes
		    tmp = n0, n0 = n1, n1 = n3, n3 = n2, n2 = tmp;
		    tmp = n4, n4 = n5, n5 = n7, n7 = n6, n6 = tmp;
		}
		mesh.elist[el]->Node[0] = n1;
		mesh.elist[el]->Node[1] = n3;
		mesh.elist[el]->Node[2] = n0;
		mesh.elist[el]->Node[3] = n5;
		el++;
		mesh.elist[el]->Node[0] = n2;
		mesh.elist[el]->Node[1] = n0;
		mesh.elist[el]->Node[2] = n3;
		mesh.elist[el]->Node[3] = n6;
		el++;
		mesh.elist[el]->Node[0] = n4;
		mesh.elist[el]->Node[1] = n5;
		mesh.elist[el]->Node[2] = n0;
		mesh.elist[el]->Node[3] = n6;
		el++;
		mesh.elist[el]->Node[0] = n7;
		mesh.elist[el]->Node[1] = n3;
		mesh.elist[el]->Node[2] = n5;
		mesh.elist[el]->Node[3] = n6;
		el++;
		mesh.elist[el]->Node[0] = n5;
		mesh.elist[el]->Node[1] = n3;
		mesh.elist[el]->Node[2] = n0;
		mesh.elist[el]->Node[3] = n6;
		el++;
	    }
	}
    }
    cout << "Generated element list (" << nels << " elements)" << endl;

    for (i = 0; i < mesh.nlist.Len(); i++)
        mesh.nlist[i].SetRegion (hreg);

    for (;;) {
        int nd, cmd;
	int xmin, xmax, ymin, ymax, zmin, zmax;

        cout << "Set region index on node subset (1/0) ? ";
        cin  >> cmd;
	if (cmd == 0) break;
	cout << "x-extension of perturbation:" << endl;
	cout << "xmin[0-" << dimx << "] xmax[0-" << dimx << "]: ";
	cin  >> xmin >> xmax;
	cout << "y-extension of perturbation:" << endl;
	cout << "ymin[0-" << dimy << "] ymax[0-" << dimy << "]: ";
	cin  >> ymin >> ymax;
	cout << "z-extension of perturbation:" << endl;
	cout << "zmin[0-" << dimz << "] zmax[0-" << dimz << "]: ";
	cin  >> zmin >> zmax;
	cout << "region index: ";
	cin  >> reg;

	for (k = zmin; k <= zmax; k++) {
	    if (k < 0 || k > dimz) continue;
	    for (j = ymin; j <= ymax; j++) {
	        if (j < 0 || j > dimy) continue;
		for (i = xmin; i <= xmax; i++) {
		    if (i < 0 || i > dimx) continue;
		    nd = k * (nnode_slice) + j * (dimy+1) + i;
		    mesh.nlist[nd].SetRegion (reg);
		}
	    }
	}
    }

    mesh.Setup();
    cout << "Mesh initialised" << endl;

    // write mesh
    ofstream ofs(fname);
    ofs << mesh;
    cout << "Mesh written to " << fname << endl;
    cout << "Bye." << endl;
    return 0;
}
