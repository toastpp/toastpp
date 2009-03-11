typedef int integer;
// predclaration of functions
int ordmmd_(integer *neqns, integer *xadj, integer *adjncy, integer *invp, 
	    integer *perm, integer *iwsiz, integer *iwork, integer *nofsub, 
	    integer *iflag);
int genmmd_(integer *neqns, integer *xadj, integer *adjncy, integer *invp, 
	    integer *perm, integer *delta, integer *dhead, integer *qsize, 
	    integer *llist, integer *marker, integer *maxint, integer *nofsub);
int mmdint_(integer *neqns, integer *xadj, integer *adjncy, integer *dhead, 
	    integer *dforw, integer *dbakw, integer *qsize, integer *llist, 
	    integer *marker);
int mmdelm_(integer *mdnode, integer  *xadj, integer *adjncy, integer *dhead, 
	    integer *dforw, integer *dbakw, integer *qsize, integer *llist, 
	    integer *marker,integer  *maxint, integer *tag);
int mmdnum_(integer *neqns, integer *perm, integer *invp, integer *qsize);
int mmdupd_(integer *ehead, integer *neqns, integer *xadj, integer *adjncy, 
	    integer *delta, integer *mdeg, integer *dhead, integer *dforw, 
	    integer *dbakw, integer *qsize, integer *llist, integer *marker, 
	    integer *maxint, integer *tag);
