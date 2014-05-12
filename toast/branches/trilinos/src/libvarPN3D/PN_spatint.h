#include "mathlib.h"
#include "felib.h"
/*Generating all the spatial matrices required by variable order PN method*/
void gen_spatint(const QMMesh& mesh, const RVector& muabs, const RVector& muscat, const RVector& ref, const RVector& delta, double w, double c, CCompRowMatrix& Sint, CCompRowMatrix& Sdx, CCompRowMatrix& Sdy, CCompRowMatrix& Sdz, RCompRowMatrix& Sx, RCompRowMatrix& Sy, RCompRowMatrix& Sz, RCompRowMatrix& Sdxx, RCompRowMatrix& Sdxy, RCompRowMatrix& Sdyx, RCompRowMatrix& Sdyy, RCompRowMatrix& Sdxz, RCompRowMatrix& Sdzx, RCompRowMatrix& Sdyz, RCompRowMatrix& Sdzy, RCompRowMatrix& Sdzz, RCompRowMatrix& spatA3_rte, RCompRowMatrix& spatA3_sdmx, RCompRowMatrix& spatA3_sdmy, RCompRowMatrix& spatA3_sdmz, RCompRowMatrix& SPS, RCompRowMatrix& SPSdx, RCompRowMatrix& SPSdy, RCompRowMatrix& SPSdz);

