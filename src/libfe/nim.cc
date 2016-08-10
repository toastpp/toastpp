// -*-C++-*-
// Utilities for nodal images
#define FELIB_IMPLEMENTATION

#include "nim.h"

// ===========================================================================
// Local prototypes

int ProcessNeighbours (const Mesh *mesh, int **nbrs, int *nnbrs,
    RVector &phase, TVector<bool> &processed, TVector<int> &front, int &nfront);


// ===========================================================================

FELIB int ReadNim (const char *name, int idx, RVector &img,
    char *meshname)
{
    char cbuf[256];
    int i, j = 0, imgsize = 0;

    ifstream ifs (name);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM") && strcmp (cbuf, "RIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
	else if (!strncasecmp (cbuf, "Mesh", 4) && meshname)
	    sscanf (cbuf+6, "%s", meshname);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    for (;;) {
	do {
	    ifs.getline (cbuf, 256);
	} while (ifs.good() && strncasecmp (cbuf, "Image", 5));
	if (!ifs.good()) break;
	for (i = 0; i < imgsize; i++)
	    ifs >> img[i];
	if (++j == idx) break;
    }
    return true;
}

// ===========================================================================

FELIB bool ReadNimAll (char *name, RDenseMatrix &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (name);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM") && strcmp (cbuf, "RIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;

    img.New(imgsize,0);
    RDenseMatrix img_single(imgsize,1);
    for (;;) {
	do {
	    ifs.getline (cbuf, 256);
	} while (ifs.good() && strncasecmp (cbuf, "Image", 5));
	if (!ifs.good()) break;
	for (i = 0; i < imgsize; i++)
	    ifs >> img_single(i,0);
	img = cath(img,img_single);
    }
    return true;
}

// ===========================================================================

FELIB void WriteNim (const char *name, const char *meshname, const RVector &img)
{
    ofstream ofs(name);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = N/A" << endl;
    ofs << "ImageSize = " << img.Dim() << endl;
    ofs << "EndHeader" << endl;
    ofs << "Image 0" << endl;
    const double *val = img.data_buffer();
    ofs.precision(12);
    ofs.setf (ios::scientific);
    for (int i = 0; i < img.Dim(); i++)
	ofs << val[i] << ' ';
    ofs << endl;
}

// ===========================================================================
// Unwrap nodal phase vector 'phase', starting from point 'seed'
// Return the number of unwrapped nodes

FELIB int NimPhaseUnwrap (const Mesh *mesh, RVector &phase, Point seed)
{
    // Get the node neighbour list
    int *nnbrs, **nbrs;
    int i, nlen = mesh->nlen();
    mesh->NodeNeighbourList (&nnbrs, &nbrs);

    // Find node closest to seed
    double dstmin = seed.Dist (mesh->nlist[0]);
    int seednd = 0;
    for (i = 1; i < nlen; i++) {
	double dst = seed.Dist (mesh->nlist[i]);
	if (dst < dstmin) {
	    dstmin = dst;
	    seednd = i;
	}
    }

    TVector<bool> processed(nlen);
    TVector<int> frontline(nlen);
    frontline[0] = seednd;
    int nfront = 1;
    int nunwrap = 0;

    for (i = 0; i < nlen; i++)
	processed[i] = false;
    processed[seednd] = true;

    while (nfront) {
	nunwrap += ProcessNeighbours (mesh, nbrs, nnbrs, phase, processed,
				      frontline, nfront);
    }

    // cleanup
    for (i = 0; i < nlen; i++)
	delete []nbrs[i];
    delete []nbrs;
    delete []nnbrs;

    return nunwrap;
}

int ProcessNeighbours (const Mesh *mesh, int **nbrs, int *nnbrs,
    RVector &phase, TVector<bool> &processed, TVector<int> &front, int &nfront)
{
    int nunwrap = 0;
    int i, j, noldfront = nfront;
    TVector<int> oldfront(noldfront);
    for (i = 0; i < noldfront; i++)
	oldfront[i] = front[i];

    nfront = 0;
    for (i = 0; i < noldfront; i++) {
	int s = oldfront[i];
	for (j = 0; j < nnbrs[s]; j++) {
	    int nb = nbrs[s][j];
	    if (!processed[nb]) {
	        if (phase[nb] >= phase[s]+Pi) {
		    while (phase[nb] >= phase[s]+Pi)
		        phase[nb] -= Pi2;
		    nunwrap++;
		} else if (phase[nb] < phase[s]-Pi) {
		    while (phase[nb] < phase[s]-Pi)
		        phase[nb] += Pi2;
		    nunwrap++;
		}
		processed[nb] = true;
		front[nfront++] = nb;
	    }
	}
    }
    return nunwrap;
}
