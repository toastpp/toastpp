#include "mathdef.h"

//typedef int int;
//typedef bool logical;
typedef double doublereal;

idxtype perphn_(idxtype *n, idxtype *ja, idxtype *ia, idxtype *init,
	idxtype *iperm, idxtype *mask, idxtype *maskval, idxtype *nlev,
	idxtype *riord, idxtype *levels);
idxtype bfs_(idxtype *n, idxtype *ja, idxtype *ia, idxtype *
	nfirst, idxtype *iperm, idxtype *mask, idxtype *maskval, idxtype *
	riord, idxtype *levels, idxtype *nlev);
idxtype reversp_(idxtype *n, idxtype *riord);
