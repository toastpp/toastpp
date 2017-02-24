#include "mathdef.h"

//typedef int int;
// predclaration of functions
idxtype ordmmd_(idxtype *neqns, idxtype *xadj, idxtype *adjncy, idxtype *invp, 
	    idxtype *perm, idxtype *iwsiz, idxtype *iwork, idxtype *nofsub, 
	    idxtype *iflag);
idxtype genmmd_(idxtype *neqns, idxtype *xadj, idxtype *adjncy, idxtype *invp, 
	    idxtype *perm, idxtype *delta, idxtype *dhead, idxtype *qsize, 
	    idxtype *llist, idxtype *marker, idxtype *maxint, idxtype *nofsub);
idxtype mmdint_(idxtype *neqns, idxtype *xadj, idxtype *adjncy, idxtype *dhead, 
	    idxtype *dforw, idxtype *dbakw, idxtype *qsize, idxtype *llist, 
	    idxtype *marker);
idxtype mmdelm_(idxtype *mdnode, idxtype  *xadj, idxtype *adjncy, idxtype *dhead, 
	    idxtype *dforw, idxtype *dbakw, idxtype *qsize, idxtype *llist, 
	    idxtype *marker,idxtype  *maxint, idxtype *tag);
idxtype mmdnum_(idxtype *neqns, idxtype *perm, idxtype *invp, idxtype *qsize);
idxtype mmdupd_(idxtype *ehead, idxtype *neqns, idxtype *xadj, idxtype *adjncy, 
	    idxtype *delta, idxtype *mdeg, idxtype *dhead, idxtype *dforw, 
	    idxtype *dbakw, idxtype *qsize, idxtype *llist, idxtype *marker, 
	    idxtype *maxint, idxtype *tag);
