// ==========================================================================
// Module libfe
// File qmmesh.cc
// Definition of class QMMesh
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

QMMesh::QMMesh (): Mesh ()
{ 
    nQ = nM = 0;
    QMofs = 0;
    Mel = 0;
    Mcosingder = 0;
    source_profile = 0;
    meas_profile = 0;
    mptype = qptype = 0;
    mwidth = qwidth = 0;
    msup   = qsup = 0;
    fixed_q_pos = fixed_m_pos = false;
}

QMMesh::~QMMesh ()
{
    int i;

    if (nQ > 0) {
	for (i = 0; i < nQ; i++) {
	    if (nQMref[i]) delete []QMref[i];
	    if (QMofs && QMofs[i]) delete []QMofs[i];
	}
	delete []Q;
	delete []QMref;
	delete []nQMref;
	delete []Qofs;
	if (QMofs) delete []QMofs;
    }
    if (nM > 0) {
	delete []M;
	if (Mel) delete []Mel;
	if (Mcosingder) delete []Mcosingder;
    }
    if (meas_profile) delete []meas_profile;
    if (mptype) delete []mptype;
    if (mwidth) delete []mwidth;
    if (msup) delete []msup;
    if (source_profile) delete []source_profile;
    if (qptype) delete []qptype;
    if (qwidth) delete []qwidth;
    if (qsup) delete []qsup;
}

static struct {
    Surface *surf;
    RVector prm0;
    double sup;
    double fac;
} source_map_param;

void InitCosineWeight (Surface *surf, Point &p0, double width,
    double sup = -1.0)
{
    source_map_param.surf = surf;
    source_map_param.prm0.New (surf->ParamDim());
    surf->Point2Param (p0, source_map_param.prm0);
    source_map_param.fac  = Pi/(2.0*width);
    source_map_param.sup  = (sup > 0.0 ? sup : width);
}

void InitTophatWeight (Surface *surf, Point &p0, double width)
{
    source_map_param.surf = surf;
    source_map_param.prm0.New (surf->ParamDim());
    surf->Point2Param (p0, source_map_param.prm0);
    source_map_param.fac  = width;
}

double CosineWeight (const Point &glob)
{
    RVector prm(source_map_param.surf->ParamDim());
    source_map_param.surf->Point2Param (glob, prm);
    double d = source_map_param.surf->ChordDist (source_map_param.prm0, prm);
    if (d > source_map_param.sup) return 0.0; // outside support
    return cos (d * source_map_param.fac);
}

double TophatWeight (const Point &glob)
{
    RVector prm(source_map_param.surf->ParamDim());
    source_map_param.surf->Point2Param (glob, prm);
    double d = source_map_param.surf->ChordDist (source_map_param.prm0, prm);
    return (d > source_map_param.fac ? 0.0 : 1.0);
}

void QMMesh::InitM ()
{
    int i, j, k, is, in, el, sd, bndsd, nd, nnd, bnd;

    Mel = new int[nM];
    Mcosingder = new RDenseMatrix[nM];

    for (i = 0; i < nM; i++) {
        if (fixed_m_pos) {
	    el = ElFind (M[i]);
	    if (el < 0) { // hack wiggle positions a bit in case we are
		Point p(M[i]);
                const double eps = 1e-4;         // just outside the mesh
		for (int d = 0; d < p.Dim(); d++) {
		    for (int s = -1; s <= 1; s++) {
			M[i] = p;
			M[i][d] += s*eps;
			el = ElFind (M[i]);
			if (el >= 0) break;
		    }
		    if (el >= 0) break;
		}
	    }

	    xASSERT(el >= 0, No element found to place measurement.);
	    // the following searches a boundary side in el.
	    // This can go wrong if el contains more than 1 boundary side
	    bndsd = -1;
	    for (is = 0; is < elist[el]->nSide() && bndsd < 0; is++) {
	        bool isbnd = true;
		for (in = 0; in < elist[el]->nSideNode(is) && isbnd; in++) {
		    int nd = elist[el]->Node[elist[el]->SideNode(is,in)];
		    if (!nlist[nd].isBnd()) isbnd = false;
		}
		if (isbnd) bndsd = is;
	    }
	    xASSERT (bndsd >= 0, Element contains no boundary side);
	} else {
	    PullToBoundary (M[i], M[i], el, bndsd);
	    xASSERT(el >= 0, No element found to place measurement.);
	}
	Point loc = elist[el]->Local (nlist, M[i]);
	RDenseMatrix lder = elist[el]->LocalShapeD (loc);
	RDenseMatrix jacin = inverse (lder * ElGeom (el));
	RVector cosin = elist[el]->DirectionCosine (bndsd, jacin);
	RDenseMatrix gder = jacin * lder;
	Mel[i] = el;
	Mcosingder[i].New (gder.Dim(ROW), gder.Dim(COL));
	for (is = 0; is < gder.Dim(ROW); is++)
	    for (in = 0; in < gder.Dim(COL); in++)
		Mcosingder[i](is,in) = cosin[is] * gder(is,in);
    }

    for (nQM = 0, i = 0; i < nQ; i++) nQM += nQMref[i];
	// total number of measurements

    // calculate the source and measurement boundary profiles
    for (int which = 0; which < 2; which++) {
        RVector *profile = 0;
	Point *P = (which ? M : Q);
	int nP = (which ? nM : nQ);
	MProfileType *ptype = (which ? mptype : qptype);
	double *width = (which ? mwidth : qwidth);
	double *sup = (which ? msup : qsup);
	double (*weightfunc)(const Point&);

	if (Boundary()) {
	    profile = new RVector[nP];
	    int nb, *bndel, *bndsd;
	    nb = BoundaryList (&bndel, &bndsd);
	    RVector bndint;
	    for (i = 0; i < nP; i++) {
	        switch (ptype[i]) {
		case PROFILE_COSINE:
		    InitCosineWeight (Boundary(), P[i], width[i], sup[i]);
		    weightfunc = CosineWeight;
		    break;
		case PROFILE_TOPHAT:
		    InitTophatWeight (Boundary(), P[i], width[i]);
		    weightfunc = TophatWeight;
		    break;
		}
	        profile[i].New (nbnd());
		RVector pprm(Boundary()->ParamDim());
		Boundary()->Point2Param (P[i], pprm);
		for (j = 0; j < nb; j++) {
		    el = bndel[j];
		    sd = bndsd[j];
		    nnd = elist[el]->nSideNode(sd);
		    bndint.New(nnd);
		    switch (ptype[i]) {
		    case PROFILE_COSINE:
		    case PROFILE_TOPHAT:
		        bndint = elist[el]->BndIntFX (sd, weightfunc, nlist);
		        //bndint = elist[el]->BndIntFCos (sd, Boundary(), pprm,
			//     width[i], sup[i], nlist);
			break;
		    case PROFILE_POINT:
		    default: // assuming point as default profile
		        bndint = elist[el]->BndIntFDelta (sd, Boundary(), pprm,
			     nlist);
			break;
		    }
		    for (k = 0; k < nnd; k++) {
		        nd = elist[el]->Node[elist[el]->SideNode (sd, k)];
			bnd = IndexNode2Bnd[nd];
			profile[i][bnd] += bndint[k];
		    }
		}
	    }
	    delete []bndel;
	    delete []bndsd;
	}
	if (which) {
	    if (meas_profile) delete []meas_profile;
	    meas_profile = profile;
	} else {
	    if (source_profile) delete []source_profile;
	    source_profile = profile;
	}
    }
}

void QMMesh::SetupQM (Point *q, int nq, Point *m, int nm)
{
    // 1. assume all sources are linked to all measurements
    // 2. assume Point source and detector profiles

    int i, j;
    Q = q, nQ = nq;
    M = m, nM = nm;
    nQM = nq*nm;

    nQMref = new int[nQ];
    Qofs   = new int[nQ];
    QMref  = new int*[nQ];
    QMofs  = new int*[nQ];
    for (i = 0; i < nQ; i++) {
	nQMref[i] = nM;
	Qofs[i] = (i ? Qofs[i-1]+nQMref[i-1] : 0);
	QMref[i] = new int[nM];
	for (j = 0; j < nM; j++) QMref[i][j] = j;
	QMofs[i] = new int[nM];
	for (j = 0; j < nM; j++) QMofs[i][j] = Qofs[i]+j;
    }

    // reset source profile specs
    if (qptype) delete []qptype;     qptype = new MProfileType[nQ];
    if (qwidth) delete []qwidth;     qwidth = new double[nQ];
    if (qsup)   delete []qsup;       qsup = new double[nQ];
    for (i = 0; i < nQ; i++) {
        qptype[i] = PROFILE_POINT; // default
	qwidth[i] = 0.0;
	qsup[i] = 0.0;
    }

    // reset measurement profile specs
    if (mptype) delete []mptype;     mptype = new MProfileType[nM];
    if (mwidth) delete []mwidth;     mwidth = new double[nM];
    if (msup)   delete []msup;       msup = new double[nM];
    for (i = 0; i < nM; i++) {
        mptype[i] = PROFILE_POINT; // default
	mwidth[i] = 0.0;
	msup[i] = 0.0;
    }

    InitM ();
}

void QMMesh::SetupRegularQM (int nqm)
{
    int el, sd, dim = nlist[0].Dim();
    Point *newQ = new Point[nqm];
    Point *newM = new Point[nqm];
    Point cnt(dim);
    double radius = Size (&cnt);
    for (int i = 0; i < nqm; i++) {
	newQ[i].New (dim);
	newM[i].New (dim);
	newQ[i][0] = cnt[0] + radius * cos ((2.0 * Pi * i) / nqm);
	newQ[i][1] = cnt[1] + radius * sin ((2.0 * Pi * i) / nqm);
	if (dim == 3) newQ[i][2] = cnt[2];
	newM[i][0] = cnt[0] + radius * cos (2.0 * Pi * (i+0.5) / nqm);
	newM[i][1] = cnt[1] + radius * sin (2.0 * Pi * (i+0.5) / nqm);
	if (dim == 3) newM[i][2] = cnt[2];
	PullToBoundary (newQ[i], newQ[i], el, sd);
	PullToBoundary (newM[i], newM[i], el, sd);
    }
    SetupQM (newQ, nqm, newM, nqm);
}

void QMMesh::LoadQM (istream &is)
{
    char cbuf[256], flagstr[100], c;
    int i, j, dim, nitem;
    double param1, param2, crd[3];
    bool ok;

    // read header
    is.getline (cbuf, 256);
    if (strncmp (cbuf, "QM file", 7)) {
        xERROR(QM file not found or invalid format);
	return;
    }
    is.getline (cbuf, 256);
    xASSERT(!strncmp (cbuf, "Dimension", 9), Unknown QM file format.);
    sscanf (cbuf+9, "%d", &dim);
    xASSERT(dim == Dimension(), QM dimension does not match mesh);

    // read source list
    do {
	is.getline (cbuf, 256);
    } while ((ok = is.good ()) && strncmp (cbuf, "SourceList", 10));
    xASSERT(ok, SourceList not found.);
    nitem = sscanf (cbuf+10, "%d%s", &nQ, flagstr);
    xASSERT(nitem >= 1, Number of sources not found.);
    fixed_q_pos = (nitem > 1 && !strcasecmp (flagstr, "fixed"));
    Q = new Point[nQ];

    // reset source profile specs
    if (qptype) delete []qptype;     qptype = new MProfileType[nQ];
    if (qwidth) delete []qwidth;     qwidth = new double[nQ];
    if (qsup)   delete []qsup;       qsup = new double[nQ];
    for (i = 0; i < nQ; i++) {
        qptype[i] = PROFILE_POINT; // default
	qwidth[i] = 0.0;
	qsup[i] = 0.0;
    }

    for (i = 0; i < nQ; i++) {
	Q[i].New(dim);
	is.getline (cbuf, 256);
	if (dim == 2) {
	    nitem = sscanf (cbuf, "%lf%lf%s%lf%lf", crd+0, crd+1,
			    flagstr, &param1, &param2);
	} else if (dim == 3) {
	    nitem = sscanf (cbuf, "%lf%lf%lf%s%lf%lf", crd+0, crd+1, crd+2,
			    flagstr, &param1, &param2);
	}
	xASSERT(nitem >= dim, Parse error while reading source list);
	for (j = 0; j < dim; j++) Q[i][j] = crd[j];
	if (nitem == dim) continue; // use default profile
	nitem -= dim;
	if (!strcasecmp (flagstr, "POINT")) {
	    qptype[i] = PROFILE_POINT;
	} else if (!strcasecmp (flagstr, "COSINE")) {
	    xASSERT(nitem >= 2, Parse error while reading source list);
	    qptype[i] = PROFILE_COSINE;
	    qwidth[i] = param1;
	    if (nitem > 2) qsup[i] = param2;
	    else qsup[i] = Pi05/param1; // default support
	} else if (!strcasecmp (flagstr, "TOPHAT")) {
	    xASSERT(nitem == 2, Parse error while reading source list);
	    qptype[i] = PROFILE_TOPHAT;
	    qwidth[i] = param1;
	} else {
	    xERROR(Parse error while reading source list);
	}   
    }

    // read measurement list
    do {
	is.getline (cbuf, 256);
    } while ((ok = is.good ()) && strncmp (cbuf, "MeasurementList", 15));
    xASSERT(ok, MeasurementList not found.);
    nitem = sscanf (cbuf+15, "%d%s", &nM, flagstr);
    xASSERT(nitem >= 1, Number of measurements not found.);
    fixed_m_pos = (nitem > 1 && !strcasecmp (flagstr, "fixed"));
    M = new Point[nM];

    // reset measurement profile specs
    if (mptype) delete []mptype;     mptype = new MProfileType[nM];
    if (mwidth) delete []mwidth;     mwidth = new double[nM];
    if (msup)   delete []msup;       msup = new double[nM];
    for (i = 0; i < nM; i++) {
        mptype[i] = PROFILE_POINT; // default
	mwidth[i] = 0.0;
	msup[i] = 0.0;
    }

    for (i = 0; i < nM; i++) {
	M[i].New(dim);
	is.getline (cbuf, 256);
	if (dim == 2) {
	    nitem = sscanf (cbuf, "%lf%lf%s%lf%lf", crd+0, crd+1,
			    flagstr, &param1, &param2);
	} else if (dim == 3) {
	    nitem = sscanf (cbuf, "%lf%lf%lf%s%lf%lf", crd+0, crd+1, crd+2,
			    flagstr, &param1, &param2);
	}
	xASSERT(nitem >= dim, Parse error while reading measurement list);
	for (j = 0; j < dim; j++) M[i][j] = crd[j];
	if (nitem == dim) continue; // use default profile
	nitem -= dim;
	if (!strcasecmp (flagstr, "POINT")) {
	    mptype[i] = PROFILE_POINT;
	} else if (!strcasecmp (flagstr, "COSINE")) {
	    xASSERT(nitem >= 2, Parse error while reading measurement list);
	    mptype[i] = PROFILE_COSINE;
	    mwidth[i] = param1;
	    if (nitem > 2) msup[i] = param2;
	    else msup[i] = Pi05/param1; // default support
	} else if (!strcasecmp (flagstr, "TOPHAT")) {
	    xASSERT(nitem == 2, Parse error while reading measurement list);
	    mptype[i] = PROFILE_TOPHAT;
	    mwidth[i] = param1;
	} else {
	    xERROR(Parse error while reading measurement list);
	}   
    }

    // read link list
    do {
	is.getline (cbuf, 256);
    } while ((ok = is.good ()) && strncmp (cbuf, "LinkList", 8));
    xASSERT(ok, LinkList not found.);
    nQMref = new int[nQ];
    Qofs   = new int[nQ];
    QMref  = new int*[nQ];
    QMofs  = new int*[nQ];
    for (i = 0; i < nQ; i++) {
	is >> nQMref[i];
	Qofs[i] = (i ? Qofs[i-1]+nQMref[i-1] : 0);
	do { is.get(c); } while (c != ':');
	QMref[i] = new int[nQMref[i]];
	for (j = 0; j < nQMref[i]; j++) is >> QMref[i][j];
	QMofs[i] = new int[nM];
	for (j = 0; j < nM; j++) QMofs[i][j] = -1;
	for (j = 0; j < nQMref[i]; j++) QMofs[i][QMref[i][j]] = Qofs[i]+j;
    }

    InitM ();
}

void QMMesh::ScaleMesh (double scale)
{
    int i;
    Mesh::ScaleMesh (scale);
    for (i = 0; i < nQ; i++) Q[i] *= scale;
    for (i = 0; i < nM; i++) M[i] *= scale;
    InitM ();
}

void QMMesh::ScaleMesh (const RVector &scale)
{
    int i, j, dim = Dimension();
    Mesh::ScaleMesh (scale);
    for (i = 0; i < nQ; i++) 
        for (j = 0; j < dim; j++) Q[i][j] *= scale[j];
    for (i = 0; i < nM; i++)
        for (j = 0; j < dim; j++) M[i][j] *= scale[j];
    InitM ();
}

bool QMMesh::Connected (int q, int m) const
{
    for (int i = 0; i < nQMref[q]; i++)
	if (QMref[q][i] == m) return true;
    return false;
}

bool QMMesh::Connected (int qm) const
{
    return Connected (qm / nM, qm % nM);
}

int QMMesh::Meas (int q, int m) const
{
    for (int i = 0; i < nQMref[q]; i++)
	if (QMref[q][i] == m) return i;
    return -1;
}

void QMMesh::GetMesh (istream& i)
{
    i >> *(Mesh*)this;
}

void QMMesh::PutMesh (ostream& o)
{
    o << *(Mesh*)this;
}
