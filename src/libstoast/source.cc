#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

RVector QVec_Gaussian (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode)
{
    int n = mesh.nlen();
    RVector qvec(n);
    Element *pel;
    int i, j, is, js, el, nnode, *node;
    double d, q, w2 = w*w;
    double fac1 = 1.0 / (sqrt (2.0*Pi) * w);
    if (mesh.Dimension() == 3) fac1 /= (sqrt (2.0*Pi) * w);
    if (mode != SRCMODE_NEUMANN) fac1 /= (sqrt (2.0*Pi) * w);
    double fac2 = -0.5/w2;

    // assemble source vector
    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    if ((is = node[i]) >= n) continue;
	    if (mode == SRCMODE_NEUMANN && !mesh.nlist[is].isBnd()) continue;

	    // source strength at node is
	    d = cnt.Dist (mesh.nlist[is]);
	    q = exp (d*d*fac2) * fac1;

	    for (j = 0; j < nnode; j++) {
		if ((js = node[j]) >= n) continue;
		if (mode == SRCMODE_NEUMANN) {
		    if (mesh.nlist[js].isBnd())
		        qvec[js] += q * pel->BndIntFF (i,j);
		} else
		    qvec[js] += q * pel->IntFF (i,j);
	    }
	}
    }
    return qvec;
}


// ============================================================================


RVector QVec_Cosine (const Mesh &mesh, const Point &cnt, double w,
    SourceMode mode)
{
    int n = mesh.nlen();
    RVector qvec(n);
    Element *pel;
    int i, j, is, js, el, nnode, *node;
    double scale = 0.5*Pi/w;
    double d, q;
    double fac1 = 1.0 / (sqrt (2.0*Pi) * w);
    if (mesh.Dimension() == 3) fac1 /= (sqrt (2.0*Pi) * w);
    if (mode != SRCMODE_NEUMANN) fac1 /= (sqrt (2.0*Pi) * w);

    // assemble source vector
    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	for (i = 0; i < nnode; i++) {
	    if ((is = node[i]) >= n) continue;
	    if (mode == SRCMODE_NEUMANN && !mesh.nlist[is].isBnd()) continue;

	    // source strength at node is
	    d = cnt.Dist (mesh.nlist[is]);
	    q = (d < w ? cos (scale*d) : 0.0);
	    if (q) {
	        for (j = 0; j < nnode; j++) {
		    if ((js = node[j]) >= n) continue;
		    if (mode == SRCMODE_NEUMANN) {
		        if (mesh.nlist[is].isBnd())
			    qvec[js] += q * pel->BndIntFF (i,j);
		    } else
		        qvec[js] += q * pel->IntFF (i,j);
		}
	    }
	}
    }



    // The following doesn't appear to work (crashes with tetrahedral mesh
    // using NEUMANN, and freezes using ISOTROPIC.
    // Needs more thought.

#ifdef UNDEF
    double scale = 0.5*Pi/w;
    int n = mesh.nlen();
    int i, j, el, sd, nnode, *node, nabsc = 0;
    double *wght, val, tval = 0.0, d;
    Point *absc;
    RVector *F;
    RDenseMatrix *D;
    Element *pel;
    RVector qvec (n);
    double fac1 = 1.0 / sqrt (2.0*Pi * w*w);
    double fac2 = -0.5/(w*w);

    for (el = 0; el < mesh.elen(); el++) {
        pel = mesh.elist[el];
	nnode = pel->nNode();
	node  = pel->Node;
	if (mode == SRCMODE_NEUMANN) {
	    for (sd = 0; sd < pel->nSide(); sd++) {
		if (!pel->IsBoundarySide (sd)) continue;
//#define USE_SUBSAMPLING
#ifdef USE_SUBSAMPLING
		nabsc = pel->GetBndSubsampleFD (sd, nabsc, wght, absc, F, D,
						mesh.nlist);
		for (i = 0; i < nabsc; i++) {
		    d = cnt.Dist (absc[i]);
		    val = (d < w ? cos (scale*d) : 0.0);
		    if (val) {
			val *= wght[i];
			// scale with interval width
			val *= M_PI/(4.0*w);
			// normalise with integral over cosine
			for (j = 0; j < nnode; j++)
			    if (mesh.nlist[node[j]].isBnd())
				qvec[node[j]] += val * F[i][j];
		    }
		}
#else
		RVector ucos = pel->BndIntFCos (sd, cnt, w, mesh.nlist);
		ucos *= Pi/(4.0*w);
		// normalise with integral over cosine
		for (i = 0; i < pel->nSideNode (sd); i++)
		    qvec[node[pel->SideNode (sd, i)]] += ucos[i];
#endif
	    }
	} else {
	    nabsc = pel->GetSubsampleFD (nabsc, wght, absc, F, D, mesh.nlist);
	    for (i = 0; i < nabsc; i++) {
		d = cnt.Dist (absc[i]);
		val = (d < w ? cos (scale*d) : 0.0);
		if (val) {
		    val *= wght[i];
		    tval += val;
		    for (j = 0; j < nnode; j++)
			qvec[node[j]] += val * F[i][j];
		}
	    }
	}
    }
    if (mode == SRCMODE_ISOTROPIC) qvec /= tval;
#endif
    return qvec;
}

// ============================================================================

RVector QVec_Point (const Mesh &mesh, const Point &cnt, SourceMode mode)
{
    int i, nd, n = mesh.nlen();
    const NodeList &nlist = mesh.nlist;
    const ElementList &elist = mesh.elist;

    int el = mesh.ElFind (cnt);
    if (el < 0) { // point outside mesh
	double d, d2, dmin = 1e10, dmin2 = 1e10;
	const double eps = 1e-8;
	int j, k;
	for (i = 0; i < mesh.elen(); i++) {
	    for (j = 0; j < elist[i]->nNode(); j++) {
		nd = elist[i]->Node[j];
		d = nlist[nd].Dist(cnt);
		if (d < dmin+eps) {
		    for (k = 0; k < elist[i]->nNode(); k++) {
			if (k == j) continue;
			d2 = nlist[elist[i]->Node[k]].Dist(cnt);
			if (d2 < dmin2) {
			    dmin2 = d2;
			    dmin = d;
			    el = i;
			}
		    }
		}
	    }
	}
    }

    RVector qvec (n);
    Point loc = elist[el]->Local (mesh.nlist, cnt);
    RVector fun = elist[el]->LocalShapeF (loc);
    for (i = 0; i < elist[el]->nNode(); i++) {
	nd = elist[el]->Node[i];
	if (nlist[nd].isBnd() || mode != SRCMODE_NEUMANN)
	    qvec[nd] = fun[i];
    }
    return qvec;
}
