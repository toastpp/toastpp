#include <fstream.h>
#include <iostream.h>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

void AddObject (bool *bmesh, int nx, int ny, int nz);
void AddCylinder (bool *bmesh, int nx, int ny, int nz);
double dst (const Point &p, const Point &p0, const RVector &d);

int main (void)
{
    int nx, ny, nz, n, i;
    char cmd;
    bool *bmesh;

    cout << "Mesh bounding box (pixels):\n";
    cout << "[nx ny nz] >> ";
    cin >> nx >> ny >> nz;
    n = nx*ny*nz;
    bmesh = new bool[n];
    for (i = 0; i < n; i++) bmesh[i] = false;

    do {
	cout << "\n(1) Add object\n";
	cout << "(Q) Quit\n";
	cout << "[1|Q] >> ";
	cin >> cmd;
	switch (toupper (cmd)) {
	case '1':
	    AddObject (bmesh, nx, ny, nz);
	    break;
	}
    } while (toupper (cmd) != 'Q');

    // write pixel mesh to file
    FILE *f = fopen ("pixmesh.bin", "wb");
    fwrite (bmesh, sizeof(bool), n, f);
    fclose (f);

    Mesh mesh;
    CreateVoxelMesh (nx, ny, nz, bmesh, 1, 1, 1, mesh);
    ofstream ofs ("pixmesh.msh");
    ofs << mesh << endl;

    delete []bmesh;
    return 0;
}

void AddObject (bool *bmesh, int nx, int ny, int nz)
{
    char cmd;

    cout << "\nAdd object\n";
    cout << "(1) Add cylinder\n";
    cout << "(Q) Cancel\n";
    cout << "[1|Q] >> ";
    cin >> cmd;
    switch (toupper (cmd)) {
    case '1':
	AddCylinder (bmesh, nx, ny, nz);
	break;
    }
}

void AddCylinder (bool *bmesh, int nx, int ny, int nz)
{
    int i, j, k;
    double bx, by, bz, dx, dy, dz, len, rad;
    Point px(3), p0(3);
    RVector d(3);

    cout << "\nAdd cylinder\n";
    cout << "Enter base point:\n";
    cout << "[x y z] >> ";
    cin  >> bx >> by >> bz;
    cout << "Enter direction vector:\n";
    cout << "[dx dy dz] >> ";
    cin  >> dx >> dy >> dz;
    len = sqrt (dx*dx + dy*dy + dz*dz);
    dx /= len, dy /= len, dz /= len;
    cout << "Enter length:\n";
    cout << "[l] >> ";
    cin  >> len;
    cout << "Enter radius:\n";
    cout << "[r] >> ";
    cin  >> rad;
    
    p0[0] = bx; p0[1] = by; p0[2] = bz;
    d[0] = dx; d[1] = dy; d[2] = dz;

    for (k = 0; k < nz; k++) {
	px[2] = k+0.5;
	for (j = 0; j < ny; j++) {
	    px[1] = j+0.5;
	    for (i = 0; i < nx; i++) {
		px[0] = i+0.5;
		if (dst(px, p0, d) <= rad)
		    bmesh[i + j*nx + k*nx*ny] = true;
	    }
	}
    }
}

double dst (const Point &p, const Point &p0, const RVector &d)
{
    // helper function: distance of point p from straight line defined
    // by point p0 and direction d

    double f1 = (p[0]-p0[0])*d[1] - (p[1]-p0[1])*d[0];
    double f2 = (p[1]-p0[1])*d[2] - (p[2]-p0[2])*d[1];
    double f3 = (p[2]-p0[2])*d[0] - (p[0]-p0[0])*d[2];

    return sqrt ((f1*f1 + f2*f2 + f3*f3)/(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]));
}
