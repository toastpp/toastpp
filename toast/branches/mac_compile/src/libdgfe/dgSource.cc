#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
using namespace std;

RVector DGQVec_Gaussian (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode)
{
    int n = mesh.elen();
    RVector qvec(4*n); //for tetrahedral mesh;
    Element *pel;
    int i, j, is, js, el, nnode, *node;
    double d, q, w2 = w*w;
    double fac1 = 1.0 / (sqrt (2.0*Pi) * w);
    if (mesh.Dimension() == 3) fac1 /= (sqrt (2.0*Pi) * w);
    if (mode != SRCMODE_NEUMANN) fac1 /= (sqrt (2.0*Pi) * w);
    double fac2 = -0.5/w2;
    cout<<"Source coords: "<<cnt[0] <<"  "<< cnt[1]<<"  "<<cnt[2]<<endl;;

    // assemble source vector
    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    //if ((is = node[i]) >= n) continue;
	    is = node[i];
	    //if (mode == SRCMODE_NEUMANN && !mesh.nlist[is].isBnd()) continue;

	    // source strength at node is
	    d = cnt.Dist (mesh.nlist[is]);
	    q = exp (d*d*fac2) * fac1;

	    for (j = 0; j < nnode; j++) {
		//if ((js = node[j]) >= n) continue;
		js = node[j];
		/*if (mode == SRCMODE_NEUMANN) {
		    if (mesh.nlist[is].isBnd())
		        qvec[4*el+j] += q * pel->BndIntFF (i,j);
		} else*/
		    qvec[4*el+j] += q * pel->IntFF (i,j);
	    }
	}
    }
    return qvec;
}



