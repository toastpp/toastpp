#ifndef __GBUTIL_H
#define __GBUTIL_H

#include "mathlib.h"
#include "felib.h"

#define IMGFMT_RAW 1
#define IMGFMT_NIM 2
#define IMGFMT_PGM 3
#define IMGFMT_PPM 4

void WritePixmap (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname, bool colour);
void WritePGM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname);
void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname);

void WritePPMArray (const RVector *img, const IVector &gdim, int nimg,
    int scalemode, double *scalemin, double *scalemax, char *rootname);
// write array of images into sequence of PPM files
// scalemode: 0=use image range provided by scalemin,scalemax
//            1=rescale each image to its own range
//            2=rescale all images to global range

void WriteData (const RVector &data, char *fname);
void WriteData_pixmap (const RVector &data, const QMMesh &mesh, bool dolog,
    bool doinvert, char *fname);
void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type);
void WriteRimHeader (const IVector &gdim, char *rimname);
void WriteImage (const RVector &nim, int imgno, char *nimname);
bool ReadNim (char *nimname, RVector &img);
void ImageScale (RVector img, double &imgmin, double &imgmax, int ofs=0);

#endif // !__GBUTIL_H
