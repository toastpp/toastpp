#include "mathlib.h"
#include "felib.h"
#include "PN_incl.h"

void ReadMesh(const char* fname, QMMesh &qmmesh);

void ReadQM(const char* fname, QMMesh &qmmesh, int &nQ, int &nM,  RCompRowMatrix &qvec, RCompRowMatrix &mvec);

void ReadParams(const char* fname, RVector &muabs, RVector &muscat, RVector &ref);

void ReadSphOrder(const char* fname, IVector &sphorder);

void ReadDirections(const char* fname, const int nQ, RVector* &dirVec);

void WriteData (const RVector &data, char *fname);

void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname);

void OpenNIM (const char *nimname, const char *meshname, int size);

void WriteNIM (const char *nimname, const RVector &img, int size, int no);

bool ReadNim (char *nimname, RVector &img);




