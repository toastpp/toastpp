#include "stoastlib.h"
#include <iomanip>
#include <string.h>
#include "util.h"
#include "solver.h"
#include "source.h"
#include "timing.h"
#include <time.h>

using namespace std;
using namespace toast;

bool ReadNim (char *name, int idx, RVector &img)
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

int main (int argc, char *argv[])
{
    //char *meshname = "ellips_tri10.msh";
    //char *tgtname = "tgt_mua_ellips_tri10.nim";
    char *meshname = "circle25_32.msh";
    char *tgtname = "tmp.nim";

    Mesh mesh;
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();
    int dim = mesh.Dimension();

    IVector bdim(dim);
    bdim = 16;
    IVector gdim(bdim*4);

    Raster *raster;
    raster = new Raster_GaussBlob(bdim,bdim,&mesh,1,1);

    int n = mesh.nlen();
    RVector prm(n);
    n = 1;

    RVector bprm(raster->BLen());
    raster->Map_MeshToBasis(prm,bprm);

    bprm = 1;
    raster->Map_BasisToMesh(bprm,prm);

    return 0;
}
