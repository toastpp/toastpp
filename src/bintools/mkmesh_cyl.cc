#include<fstream>
#include<mathlib.h>
#include<felib.h>

using namespace std;

void CreateNode (Node &node, double x, double y, double z, bool isBnd);
void CreateWedge (ElementList &elist, int &el, int n0, int n1, int n2, int n3,
		  int n4, int n5);
void CreateWedge2 (ElementList &elist, int &el, int n0, int n1, int n2, int n3,
		  int n4, int n5);

int main (void)
{
    double hmua, hmus, href;
    double rad, ringrad, zmin, zmax;
    double x, y, z, alpha;
    double *layerz, r, fac2, d0;
    int dimz, iz, ring, nring, sect, nsect, nnode_ring, nnode;
    int nwdg_sectring, nwdg_sect, nwdg_slice, nwdg;
    int nels_slice, nels;
    int nofs0, nofs1, ringofs00, ringofs01, ringofs10, ringofs11;
    int nnode_iring, nnode_oring, relnode0, relnode1;
    int n0, n1, n2, n3, n4, n5;
    int i, n, el, nnode_slice, hreg, cmd;
    char fname[256];
    Mesh mesh;

    cout << "mkmesh_cyl: Create a cylindrical mesh using tetrahedral elements"
	 << endl << endl;

    cout << "Mesh name: ";
    cin  >> fname;
    cout << "Physical dimensions:" << endl;
    cout << "radius: ";
    cin  >> rad;
    cout << "zmin zmax: ";
    cin  >> zmin >> zmax;

    cout << "Number of element layers along z: ";
    cin  >> dimz;
    layerz = new double[dimz+1];
    cout << "Layer thickness" << endl;
    cout << "(1) fixed" << endl;
    cout << "(2) linear from centre to top/bottom planes" << endl;
    cout << "(3) linear from top to bottom plane" << endl;
    cout << "(4) linear from top to bottom plane with fixed top layer" << endl;
    cin  >> cmd;
    switch (cmd) {
    case 1:
        for (i = 0; i <= dimz; i++)
	    layerz[i] = zmin + (zmax-zmin) * (double)i / (double)dimz;
	break;
    case 2:
        cout << "Ratio top/central layer thickness: ";
        cin  >> r;
	if (dimz % 2) { // odd number of layers
	    double sum = 1.0;
	    double fac = fac2 = pow (r, 2.0/(double)(dimz-1));
	    for (i = 1; i <= (dimz-1)/2; i++) {
	        sum += 2.0*fac2;
		fac2 *= fac;
	    }
	    d0 = (zmax-zmin)/sum; // central layer thickness
	    z = d0*0.5;
	    for (i = 0; i < (dimz+1)/2; i++) {
	        layerz[(dimz+1)/2+i] = z;
		layerz[(dimz+1)/2-i-1] = -z;
		d0 *= fac;
		z += d0;
	    }
	} else { // even number of layers
	    double sum = 1.0;
	    double fac = fac2 = pow (r, 1.0/(double)(dimz/2-1));
	    for (i = 1; i <= dimz/2-1; i++) {
	        sum += fac2;
		fac2 *= fac;
	    }
	    d0 = (zmax-zmin)/(2.0*sum); // thickness of the 2 central layers
	    z = 0.0;
	    layerz[dimz/2] = z;
	    for (i = 1; i <= dimz/2; i++) {
	        z += d0;
		layerz[dimz/2+i] = z;
		layerz[dimz/2-i] = -z;
		d0 *= fac;
	    }
	}
	//for (i = 0; i <= dimz; i++)
	//    layerz[i] += zmin;
	break;
    case 3: {
	cout << "Ratio bottom/top layer thickness: ";
	cin >> r;
	double sum = 1.0;
	double fac = fac2 = pow (r, 1.0/(double)(dimz-1));
	for (i = 1; i <= dimz-1; i++) {
	    sum += fac2;
	    fac2 *= fac;
	}
	d0 = (zmax-zmin)/sum;
	z = zmin;
	layerz[0] = z;
	for (i = 1; i <= dimz; i++) {
	    z += d0;
	    layerz[i] = z;
	    d0 *= fac;
	}
        } break;
    case 4: {
	cout << "Top layer thickness: ";
	double z0;
	cin >> z0;
	cout << "Ratio bottom/top layer thickness: ";
	cin >> r;
	double sum = 1.0;
	double fac = fac2 = pow (r, 1.0/(double)(dimz-2));
	for (i = 1; i <= dimz-2; i++) {
	    sum += fac2;
	    fac2 *= fac;
	}
	d0 = (zmax-zmin-z0)/sum;
	z = zmin;
	layerz[0] = z;
	z += z0;
	layerz[1] = z;
	for (i = 2; i <= dimz; i++) {
	    z += d0;
	    layerz[i] = z;
	    d0 *= fac;
	}
        } break;
    }
    cout << "Created node layers at z =" << endl;
    for (i = 0; i <= dimz; i++) cout << layerz[i] << ' ';
    cout << endl << endl;


    cout << "Number of rings per layer: ";
    cin  >> nring;
    cout << "Number of sectors per layer: ";
    cin  >> nsect;

    cout << "Background parameters:" << endl;
    cout << "mua mus refind region: " ;
    cin  >> hmua >> hmus >> href >> hreg;

    nnode_slice = 1;
    nnode_ring = nsect;
    for (i = 1; i <= nring; i++) {
        nnode_slice += nnode_ring;
	nnode_ring += nsect;
    }
    nnode = nnode_slice * (dimz+1);
    cout << nnode_slice << " nodes per slice" << endl;
    cout << nnode << " nodes total" << endl;

    nwdg_sectring = 1;
    nwdg_sect = 0;
    for (i = 1; i <= nring; i++) {
        nwdg_sect += nwdg_sectring;
	nwdg_sectring += 2;
    }
    nwdg_slice = nwdg_sect * nsect;
    nwdg = nwdg_slice * dimz;
    cout << nwdg_slice << " wedges per slice" << endl;
    cout << nwdg << " wedges total" << endl;

    nels_slice = nwdg_slice * 3;
    nels = nels_slice * dimz;
    cout << nels_slice << " tetrahedrals per slice" << endl;
    cout << nels << " tetrahedrals total" << endl;

    // create node list
    mesh.nlist.New (nnode);
    n = 0;
    for (iz = 0; iz <= dimz; iz++) { // loop over node slices
        bool isbnd = (iz == 0 || iz == dimz);
	z = layerz[iz];
        //z = zmin + (zmax-zmin)*(double)iz/(double)dimz;
	CreateNode (mesh.nlist[n++], 0.0, 0.0, z, isbnd);
	for (ring = 1; ring <= nring; ring++) { // loop over rings
	    isbnd = (isbnd || (ring == nring));
	    nnode_ring = nsect * ring;
	    ringrad = rad * (double)ring/(double)nring;
	    for (i = 0; i < nnode_ring; i++) { // loop over nodes in ring
	        alpha = 2.0*Pi * (double)i / (double)nnode_ring;
		x = ringrad * cos(alpha);
		y = ringrad * sin(alpha);
		CreateNode (mesh.nlist[n++], x, y, z, isbnd);
	    } // end loop nodes in ring
	} // end loop over rings
    } // end loop over slices
    cout << "Generated node list (" << n << " nodes)" << endl;

    // create element list
    mesh.elist.Clear();
    Element **list = new Element*[nels];
    for (el = 0; el < nels; el++)
        list[el] = new Tetrahedron4;
    mesh.elist.AppendList (nels, list);
    delete []list;

    el = 0;
    for (iz = 0; iz < dimz; iz++) { // loop over element slices
        nofs0 = nnode_slice * iz; // node no. offset for nodes at lower surface
	nofs1 = nofs0 + nnode_slice; // offset for nodes at upper surface

	ringofs00 = nofs0; // inner ring, lower surface
	ringofs10 = nofs1; // inner ring, upper surface
	ringofs01 = nofs0+1; // outer ring, lower surface
	ringofs11 = nofs1+1; // outer ring, upper surface

	for (ring = 0; ring < nring; ring++) { // loop over rings
	  
	    nnode_iring = nsect * ring;     // nodes on inner ring
	    nnode_oring = nsect * (ring+1); // nodes on outer ring

	    for (sect = 0; sect < nsect; sect++) { // loop over sectors

	        relnode0 = sect * ring;
		relnode1 = sect * (ring+1);

		// create first wedge in sector
		n0 = relnode0 + ringofs00;
		n1 = relnode1 + ringofs01;
		n2 = ((relnode1+1) % nnode_oring) + ringofs01;
		n3 = relnode0 + ringofs10;
		n4 = relnode1 + ringofs11;
		n5 = ((relnode1+1) % nnode_oring) + ringofs11;
		CreateWedge (mesh.elist, el, n0, n1, n2, n3, n4, n5);


		// create wedge pairs
		for (i = 0; i < ring; i++) {
		    n0 = i+relnode0 + ringofs00;
		    n1 = ((i+relnode0+1) % nnode_iring) + ringofs00;
		    n2 = (i+relnode1+1) + ringofs01;
		    n3 = i+relnode0 + ringofs10;
		    n4 = ((i+relnode0+1) % nnode_iring) + ringofs10;
		    n5 = (i+relnode1+1) + ringofs11;
		    CreateWedge2 (mesh.elist, el, n0, n1, n2, n3, n4, n5);

		    n0 = ((i+relnode0+1) % nnode_iring) + ringofs00;
		    n1 = (i+relnode1+1) + ringofs01;
		    n2 = ((i+relnode1+2) % nnode_oring) + ringofs01;
		    n3 = ((i+relnode0+1) % nnode_iring) + ringofs10;
		    n4 = (i+relnode1+1) + ringofs11;
		    n5 = ((i+relnode1+2) % nnode_oring) + ringofs11;
		    CreateWedge (mesh.elist, el, n0, n1, n2, n3, n4, n5);
		}

	    } // end loop over sectors

	    ringofs00 = ringofs01;
	    ringofs10 = ringofs11;
	    ringofs01 += nsect * (ring+1);
	    ringofs11 += nsect * (ring+1);

	} // end loop over rings
    } // end loop over element slices
    cout << "Generated element list (" << el << " elements)" << endl;

    for (i = 0; i < mesh.nlist.Len(); i++)
        mesh.nlist[i].SetRegion (hreg);

    for (;;) {
        int n, cmd, reg;
	double dx, dy, dz, dist, rad, xcnt, ycnt, zcnt, x1, y1, z1, x2, y2, z2;

	cout << endl << "Set region index on internal area?" << endl;
	cout << "(1) Rod" << endl;
	cout << "(2) Sphere" << endl;
	cout << "(3) Brick" << endl;
	cout << "(0) Done" << endl;
        cin  >> cmd;
	if (cmd == 0) break;
	switch (cmd) {
	case 1:
	    cout << "Centre (x y): ";
	    cin  >> xcnt >> ycnt;
	    cout << "Radius: ";
	    cin  >> rad;
	    cout << "region index: ";
	    cin  >> reg;
	    for (n = 0; n < mesh.nlen(); n++) {
	        dx = mesh.nlist[n][0] - xcnt;
		dy = mesh.nlist[n][1] - ycnt;
		dist = hypot (dx, dy);
		if (dist <= rad) {
		    mesh.nlist[n].SetRegion (reg);
		}
	    }
	    break;
	case 2:
	    cout << "Centre (x y z): ";
	    cin  >> xcnt >> ycnt >> zcnt;
	    cout << "Radius: ";
	    cin  >> rad;
	    cout << "region index: ";
	    cin  >> reg;
	    for (n = 0; n < mesh.nlen(); n++) {
	        dx = mesh.nlist[n][0] - xcnt;
		dy = mesh.nlist[n][1] - ycnt;
		dz = mesh.nlist[n][2] - zcnt;
		dist = sqrt (dx*dx + dy*dy + dz*dz);
		if (dist <= rad) {
		    mesh.nlist[n].SetRegion (reg);
		}
	    }
	    break;
	case 3:
	    cout << "Brick BB (xmin ymin zmin): ";
	    cin  >> x1 >> y1 >> z1;
	    cout << "Brick BB (xmax ymax zmax): ";
	    cin  >> x2 >> y2 >> z2;
	    cout << "region index: ";
	    cin  >> reg;
	    for (n = 0; n < mesh.nlen(); n++) {
	        if (x >= x1 && x <= x2 &&
		    y >= y1 && y <= y2 &&
		    z >= z1 && z <= z2) {
		  mesh.nlist[n].SetRegion (reg);
		}
	    }
	    break;
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

void CreateNode (Node &node, double x, double y, double z, bool isBnd)
{
  node.New (3);
  node[0] = x;
  node[1] = y;
  node[2] = z;
  node.SetBndTp (isBnd ? BND_DIRICHLET : BND_NONE);
}

void CreateWedge (ElementList &elist, int &el, int n0, int n1, int n2, int n3,
		  int n4, int n5)
{
    elist[el]->Node[0] = n2;
    elist[el]->Node[1] = n0;
    elist[el]->Node[2] = n1;
    elist[el]->Node[3] = n5;
    el++;
    elist[el]->Node[0] = n3;
    elist[el]->Node[1] = n5;
    elist[el]->Node[2] = n4;
    elist[el]->Node[3] = n0;
    el++;
    elist[el]->Node[0] = n1;
    elist[el]->Node[1] = n5;
    elist[el]->Node[2] = n0;
    elist[el]->Node[3] = n4;
    el++;
}

void CreateWedge2 (ElementList &elist, int &el, int n0, int n1, int n2, int n3,
		   int n4, int n5)
{
    elist[el]->Node[0] = n2;
    elist[el]->Node[1] = n1;
    elist[el]->Node[2] = n0;
    elist[el]->Node[3] = n5;
    el++;
    elist[el]->Node[0] = n3;
    elist[el]->Node[1] = n4;
    elist[el]->Node[2] = n5;
    elist[el]->Node[3] = n0;
    el++;
    elist[el]->Node[0] = n1;
    elist[el]->Node[1] = n0;
    elist[el]->Node[2] = n5;
    elist[el]->Node[3] = n4;
    el++;
}
