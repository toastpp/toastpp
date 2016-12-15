/* 
---
Some utility functions for reading/writing data and parameters
from file or console input
---
*/

#ifndef UTIL_H
#define UTIL_H

#include "stoastlib.h"
#include "slu_zdefs.h"
#include <string>

#define MAXREGION 100
#define SOURCE_GAUSSIAN 0
#define SOURCE_COSINE 1

class Projector;
class Camera;
typedef enum { FSOLVER_MF, FSOLVER_MLEM, FSOLVER_NONE } FluoSolverType;

void OpenNIM (const char *nimname, const char *meshname, int size);
void WriteNIM (const char *nimname, const RVector &img, int size, int no);
bool ReadNim (char *nimname, RVector &img);
void WritePGM (const RVector &img, const IVector &gdim, char * fname, bool doScale=true);
void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname);
void WriteData (const RVector &data, char *fname);
void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname);
void SelectMesh (ParamParser &pp, char *meshname, QMMesh &mesh);
void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);
void SelectData (ParamParser &pp, double &freq);
void SelectInitialParams (ParamParser &pp, const Mesh &mesh, Solution &msol);
void SelectProjectors (ParamParser &pp, QMMesh & mesh, Projector ** projPtrs, Camera ** camPtrs, RVector * norms);
RVector ProjectToBoundaryNode(const QMMesh & mesh, const RVector & pos, const RVector & ray, int nBndElems, int * bndellist, int * bndsdlist);
void ReadData (char *fname, RVector &data);
void WriteBinaryData (const RVector &data, const char *fname);
bool SelectPriorImage(ParamParser &pp, char* pfname, IVector & gDim);
bool SelectLabelImage(ParamParser &pp, char * labelFname, IVector & gDim);
double SelectNoise(ParamParser &pp);
void SelectGrid (ParamParser &pp, QMMesh &mesh, 
	Raster ** rast);
void SelectRegularisation (ParamParser &pp, QMMesh &mesh, 
	Raster * rast, Regularisation ** reg, bool & linReg,
	const RVector * kref, double & tau);
bool SelectDataImage(ParamParser &pp, RVector & fluoData, RVector & excitData, int imageX, int imageY, int nImages);
void SelectSolvers (ParamParser &pp, bool & solveMua, bool & solveFluo, int & iters, bool & logSolve, double & fluoTol, IVector & dataWin);
void SelectOptical (ParamParser &pp, bool & setOptical);
bool SelectNormals(ParamParser &pp, RVector * srcNorms, RVector * detNorms, int nQ, int nM);
void ReadDoubleData (char *fname, RVector &data);
FluoSolverType SelectFluoSolver(ParamParser &pp);
#endif
