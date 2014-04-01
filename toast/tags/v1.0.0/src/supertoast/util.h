#ifndef __GBUTIL_H
#define __GBUTIL_H

#include "mathlib.h"
#include "felib.h"

#define IMGFMT_RAW 1
#define IMGFMT_NIM 2
#define IMGFMT_RAW_NIM 3
#define IMGFMT_PGM 4
#define IMGFMT_PPM 5

void WritePixmap (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, const char *fname, bool colour);

void WritePGM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, const char *fname);
bool ReadPGM (RVector &img, IVector &gdim, const char *fname);

void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, const char *fname);

void WritePPMArray (const RVector *img, const IVector &gdim, int nimg,
    int scalemode, double *scalemin, double *scalemax, const char *rootname);
// write array of images into sequence of PPM files
// scalemode: 0=use image range provided by scalemin,scalemax
//            1=rescale each image to its own range
//            2=rescale all images to global range

void WriteData (const RVector &data, const char *fname);
void WriteData_pixmap (const RVector &data, const QMMesh &mesh, bool dolog,
    bool doinvert, const char *fname);
void WriteNimHeader (const char *meshname, int imgsize, const char *nimname,
    const char *type);
void WriteEimHeader (const char *meshname, int imgsize, const char *eimname,
    const char *type);
void WriteRimHeader (const IVector &gdim, const char *rimname);
void WriteImage (const RVector &nim, int imgno, const char *nimname);
bool ReadNim (const char *nimname, RVector &img);
void ImageScale (RVector img, double &imgmin, double &imgmax, int ofs=0);

double tick();
double tock();
double tock (double tick);

#endif // !__GBUTIL_H
