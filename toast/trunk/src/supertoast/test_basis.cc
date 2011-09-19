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
    char *meshname = "cyl2.msh";
    char *tgtname = "mua_tgt_cyl2.nim";

    Mesh mesh;
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();
    int dim = mesh.Dimension();

    IVector bdim(dim);
    bdim = 16;
    IVector gdim(bdim*4);

    Raster *raster;
    raster = new Raster_CubicPixel (bdim, gdim, &mesh);

    int glen = raster->GLen();
    int blen = raster->BLen();
    int slen = raster->SLen();

    cout << "glen=" << glen << endl;
    cout << "blen=" << blen << endl;
    cout << "slen=" << slen << endl;

    int n = mesh.nlen();
    RVector mua(n);
    if (!ReadNim (tgtname, 0, mua))
	mua = 1;

    RVector gmua(glen);
    raster->Map_MeshToGrid (mua, gmua);
    ofstream ofs1 ("gmua.dat");
    ofs1 << gmua << endl;

    RVector bmua(blen);
    raster->Map_GridToBasis (gmua, bmua);
    ofstream ofs2 ("bmua.dat");
    ofs2 << bmua << endl;
    
    RVector gmua2(glen);
    raster->Map_BasisToGrid (bmua, gmua2);
    ofstream ofs3 ("gmua2.dat");
    ofs3 << gmua2 << endl;

    return 0;
}
