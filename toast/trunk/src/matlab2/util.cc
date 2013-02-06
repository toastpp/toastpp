#include "stoastlib.h"
#include "util.h"
#include <time.h>

using namespace std;

void WritePixmap (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname, bool colour)
{
    if (colour) WritePPM (img, gdim, scalemin, scalemax, fname);
    else        WritePGM (img, gdim, scalemin, scalemax, fname);
}

void WritePGM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname)
{
    int i, j, ii, dim = gdim[0]*gdim[1];
    double imgmin, imgmax, scale;
    unsigned char *pixmap = new unsigned char[dim];

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    // rescale image
    if (scalemin) {
        imgmin = *scalemin;
    } else {
	for (i = 0, imgmin = 1e100; i < dim; i++)
	    if (img[i+ii] < imgmin) imgmin = img[i+ii];
    }
    if (scalemax) {
        imgmax = *scalemax;
    } else {
      for (i = 0, imgmax = -1e100; i < dim; i++)
	  if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    scale = 256.0/(imgmax-imgmin);

    for (j = 0; j < gdim[1]; j++) {
	for (i = 0; i < gdim[0]; i++) {
	    int v = (int)((img[ii + (gdim[1]-j-1)*gdim[0] + i]-imgmin)*scale);
	    if      (v < 0  ) v = 0;
	    else if (v > 255) v = 255;
	    pixmap[j*gdim[0]+i] = (unsigned char)v;
	}
    }
    ofstream ofs(fname);
    ofs << "P5" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++) ofs << pixmap[i];
    ofs << endl;
    delete []pixmap;
}

void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname)
{
    typedef struct {
        unsigned char r,g,b;
    } RGB;

    int i, j, ii, dim = gdim[0]*gdim[1];
    unsigned char *pixmap = new unsigned char[dim];
    double imgmin, imgmax, scale;

    static RGB colmap[256];
    static bool have_colmap = false;

    if (!have_colmap) {
        int r, g, b;
        char colormap[256];
	strcpy (colormap, getenv ("TOASTDIR"));
	strcat (colormap, "/scales/fire2.pal");
        ifstream ifs (colormap);
	for (i = 0; i < 256; i++) {
	    ifs >> r >> g >> b;
	    colmap[i].r = (unsigned char)r;
	    colmap[i].g = (unsigned char)g;
	    colmap[i].b = (unsigned char)b;
	}
	have_colmap = true;
    }

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    // rescale image
    if (scalemin) {
        imgmin = *scalemin;
    } else {
	for (i = 0, imgmin = 1e100; i < dim; i++)
	    if (img[i+ii] < imgmin) imgmin = img[i+ii];
    }
    if (scalemax) {
        imgmax = *scalemax;
    } else {
      for (i = 0, imgmax = -1e100; i < dim; i++)
	  if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    scale = 256.0/(imgmax-imgmin);

    for (j = 0; j < gdim[1]; j++) {
	for (i = 0; i < gdim[0]; i++) {
	    int v = (int)((img[ii + (gdim[1]-j-1)*gdim[0] + i]-imgmin)*scale);
	    if      (v < 0  ) v = 0;
	    else if (v > 255) v = 255;
	    pixmap[j*gdim[0]+i] = (unsigned char)v;
	}
    }
    ofstream ofs(fname);
    ofs << "P6" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++)
        ofs << colmap[pixmap[i]].r
	    << colmap[pixmap[i]].g
	    << colmap[pixmap[i]].b;
    ofs << endl;

    delete []pixmap;
}

void WritePPMArray (const RVector *img, const IVector &gdim, int nimg,
    int scalemode, double *scalemin, double *scalemax, char *rootname)
{
    char fname[256];
    int i;
    double scmin, scmax, imin, imax;

    switch (scalemode) {
    case 0: // use provided image range
	break;
    case 1: // rescale each image individually
	scalemin = 0;
	scalemax = 0;
	break;
    case 2: // rescale to global range
	for (i = 0; i < nimg; i++) {
	    ImageScale (img[i], imin, imax, 0);
	    if (!i || imin < scmin) scmin = imin;
	    if (!i || imax > scmax) scmax = imax;
	}
	scalemin = &scmin;
	scalemax = &scmax;
	break;
    }
    for (int i = 0; i < nimg; i++) {
	sprintf (fname, "%s_%03d.ppm", rootname, i);
	WritePPM (img[i], gdim, scalemin, scalemax, fname);
    }
}

void WriteData (const RVector &data, char *fname)
{
    ofstream ofs (fname);
    ofs << data << endl;
}

void WriteData_pixmap (const RVector &data, const QMMesh &mesh, bool dolog,
    bool doinvert, char *fname)
{
    double vmn = (doinvert ? -vmax(data) : vmin(data));
    double vmx = (doinvert ? -vmin(data) : vmax(data));
    int nx = mesh.nM;
    int ny = mesh.nQ;
    int q, m, i;
    if (vmn <= 0) dolog = false;
    RVector fulldata (nx*ny);
    for (q = i = 0; q < ny; q++) {
	for (m = 0; m < nx; m++) {
	    if (mesh.Connected(q,m))
		fulldata[q*nx+m] = (doinvert ? -data[i++] : data[i++]);
	    else fulldata[q*nx+m] = vmn;
	    if (dolog) fulldata[q*nx+m] = log(fulldata[q*nx+m]);
	}
    }
    if (dolog) vmn = log(vmn), vmx = log(vmx);
    IVector gdim(2);
    gdim[0] = nx, gdim[1] = ny;
    WritePPM (fulldata, gdim, &vmn, &vmx, fname);
    cerr << "Data PPM written. Range: " << vmn << " to " << vmx << endl;
}

void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = " << type << endl;
    ofs << "ImageSize = " << imgsize << endl;
    ofs << "EndHeader" << endl;
}

void WriteRimHeader (const IVector &gdim, char *rimname)
{
    int i, isize = 1;
    ofstream ofs (rimname);
    ofs << "RIM" << endl;
    ofs << "Grid =";
    for (i = 0; i < gdim.Dim(); i++) {
        ofs << ' ' << gdim[i];
	isize *= gdim[i];
    }
    ofs << endl << "ImageSize = " << isize << endl;
    ofs << "EndHeader" << endl;
}

void WriteImage (const RVector &nim, int imgno, char *nimname)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << imgno << endl;
    for (int i = 0; i < nim.Dim(); i++)
        ofs << nim[i] << ' ';
    ofs << endl;
}

bool ReadNim (char *nimname, RVector &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (nimname);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    do {
        ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
        ifs >> img[i];
    return true;
}

void ImageScale (RVector img, double &imgmin, double &imgmax, int ofs)
{
    int i;
    imgmin = 1e100, imgmax = -1e100;
    for (i = 0; i < img.Dim(); i++) {
        if (img[i+ofs] < imgmin) imgmin = img[i+ofs];
	if (img[i+ofs] > imgmax) imgmax = img[i+ofs];
    }
}
