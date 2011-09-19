#include "stoastlib.h"

using namespace std;
using namespace toast;

int main (void)
{
    QMMesh mesh;
    SourceMode qmode;
    int qprof;
    double qwidth;

    ifstream ifs ("cyl2.msh");
    ifs >> mesh;
    mesh.Setup();
    int n = mesh.nlen();
    int ne = mesh.elen();
    int i, j, k, ni;
    Point p1(3), p2(3);
    Point *lst;
    p1[0] = 0, p1[1] = 0, p1[2] = 0;
    p2[0] = 1, p2[1] = 1, p2[2] = 1;

    ofstream ofs ("intersection.dat");

    for (i = 0; i < ne; i++) {
	ni = mesh.elist[i]->GlobalIntersection (mesh.nlist, p1, p2, &lst);
	for (j = 0; j < ni; j++) {
	    ofs << i;
	    Point p = mesh.elist[i]->Global (mesh.nlist, lst[j]);
	    for (k = 0; k < 3; k++)
		ofs << ' ' << p[k];
	    ofs << endl;
	}
    }


    return 0;
}
