#include <mathlib.h>
#include <felib.h>
#include <iostream>
#include <fstream>
#include <ctype.h>

#define MAXQ 1000
#define MAXM 1000

using namespace std;

// a few ANSI terminal escape sequences
char clear[]  = "\033[2J\033[0;0f";
char normal[] = "\033[0m";
char bold[]   = "\033[1m";

Mesh mesh;
double meshsize;
Point meshcnt(3);

int nq = 0;
int nm = 0;

Point Q[MAXQ];
Point M[MAXM];

void PlacePlanar ();
void WriteQM ();

int main (void)
{
    char fmesh[256], cmd;

    cout << "makeqm3d: Generate source-detector description file for a 3D mesh"
	 << endl << endl;
    cout << "Mesh file: ";
    cin  >> fmesh;
    cout << "Reading and initialising mesh ..." << endl;
    ifstream ifs (fmesh);
    ifs >> mesh;
    mesh.Setup ();
    meshsize = mesh.Size (&meshcnt);

    do {
        cout << clear << bold << "makeqm3d" << normal << endl << endl;;
	cout << "(1) Generate a QM plane" << endl;
	cout << "(W) Write QM file" << endl;
	cout << "(Q) Quit" << endl;
	cin  >> cmd;
	cmd = toupper (cmd);
	switch (cmd) {
	case '1': PlacePlanar (); break;
	case 'W': WriteQM (); break;
	}
    } while (cmd != 'Q');
    return 0;
}

void PlacePlanar ()
{
    double a, b, c, d, rad, phi, dphi, phi_ofs, qofs;
    double scl, cosa, cosb, cosc, dst, theta, cost, sint, cosp, sinp;
    int i, nqm;
    Point p1(3), p2(3), pb(3);
    Point px(3), py(3), pz(3);

    cout << clear << bold << "Source-detector placement in a plane"
	 << normal << endl << endl;
    cout << "Define plane: ax + by + cz + d = 0" << endl;
    cout << "a b c d = "; 
    cin  >> a >> b >> c >> d;
    cout << "Number of optodes: ";
    cin  >> nqm;
    cout << "Angular offset of first optode [deg]: ";
    cin  >> phi_ofs;

    cout << "Sub-surface placement of source optodes (0=on surface) [mm]: ";
    cin  >> qofs;
    // this pushes source positions towards the origin!

    phi_ofs *= Pi/180.0;
    dphi = 2.0*Pi/nqm;
    rad = 2.0*meshsize;

    // direction cosines of plane normal
    scl = 1.0 / sqrt (a*a + b*b + c*c);
    cosa = a*scl, cosb = b*scl, cosc = c*scl;
    dst = -d*scl;

    // rotation angles for transform matrix
    theta = acos(cosc);
    phi   = atan2 (cosb, cosa);
    cerr << "theta=" << 180.0/Pi*theta << endl;
    cerr << "phi=" << 180.0/Pi*phi << endl;

    cost  = cos(theta), sint = sin(theta);
    cosp  = cos(phi),   sinp = sin(phi);

    // transform is as follows:
    // 1. move point by `dst' along z
    // 2. rotate by `theta' around y
    // 3. rotate by `phi' around z

    px[0] = px[1] = px[2] = 0.0;

    py[0] = px[0], py[1] = px[1], py[2] = px[2]+dst;

    pz[0] =  cost*py[0] + 0*py[1] + sint*py[2];
    pz[1] =     0*py[0] + 1*py[1] +    0*py[2];
    pz[2] = -sint*py[0] + 0*py[1] + cost*py[2];

    p1[0] =  cosp*pz[0] - sinp*pz[1] + 0*pz[2];
    p1[1] =  sinp*pz[0] + cosp*pz[1] + 0*pz[2];
    p1[2] =     0*pz[0] +    0*pz[1] + 1*pz[2];

    for (i = 0; i < nqm; i++) {
        cout << "\033[s" << i << endl << "\033[u" << flush;
        phi = dphi*i + phi_ofs;

	px[0] = rad * cos(phi), px[1] = rad*sin(phi), px[2] = 0.0;

	py[0] = px[0], py[1] = px[1], py[2] = px[2]+dst;

	pz[0] =  cost*py[0] + 0*py[1] + sint*py[2];
	pz[1] =     0*py[0] + 1*py[1] +    0*py[2];
	pz[2] = -sint*py[0] + 0*py[1] + cost*py[2];

	p2[0] =  cosp*pz[0] - sinp*pz[1] + 0*pz[2];
	p2[1] =  sinp*pz[0] + cosp*pz[1] + 0*pz[2];
	p2[2] =     0*pz[0] +    0*pz[1] + 1*pz[2];

	pb = mesh.BndIntersect (p1, p2);
	if (!(i & 1)) {
	    if (qofs) {
	        // move towards centre of optode plane
	        // this is not always guaranteed to work
	        RVector shift = pb-p1;
		shift *= qofs/length(shift);
		pb -= shift;
	    }
	    Q[nq].New(3); Q[nq++] = pb;
	} else {
	    M[nm].New(3); M[nm++] = pb;
	}
    }
	
}

void WriteQM ()
{
    char cbuf[256];
    int i, j;

    cout << clear << bold << "Write QM description to file" << normal
	 << endl << endl;
    cout << "File name: ";
    cin  >> cbuf;
    ofstream ofs (cbuf);
    ofs << "QM file 3D" << endl;
    ofs << "Dimension 3" << endl;
    ofs << endl;
    ofs.precision (8);
    ofs << "SourceList " << nq << endl;
    for (i = 0; i < nq; i++)
        for (j = 0; j < 3; j++)
	    ofs << Q[i][j] << (j == 2 ? '\n' : ' ');
    ofs << endl << "MeasurementList " << nm << endl;
    for (i = 0; i < nm; i++)
        for (j = 0; j < 3; j++)
	    ofs << M[i][j] << (j == 2 ? '\n' : ' ');
    ofs << endl << "LinkList" << endl;
    for (i = 0; i < nq; i++) {
        ofs << nm << ": ";
	for (j = 0; j < nm; j++)
	    ofs << j << (j == nm-1 ? '\n' : ' ');
    }
}
