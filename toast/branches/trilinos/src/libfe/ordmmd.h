//typedef int int;
// predclaration of functions
int ordmmd_(int *neqns, int *xadj, int *adjncy, int *invp, 
	    int *perm, int *iwsiz, int *iwork, int *nofsub, 
	    int *iflag);
int genmmd_(int *neqns, int *xadj, int *adjncy, int *invp, 
	    int *perm, int *delta, int *dhead, int *qsize, 
	    int *llist, int *marker, int *maxint, int *nofsub);
int mmdint_(int *neqns, int *xadj, int *adjncy, int *dhead, 
	    int *dforw, int *dbakw, int *qsize, int *llist, 
	    int *marker);
int mmdelm_(int *mdnode, int  *xadj, int *adjncy, int *dhead, 
	    int *dforw, int *dbakw, int *qsize, int *llist, 
	    int *marker,int  *maxint, int *tag);
int mmdnum_(int *neqns, int *perm, int *invp, int *qsize);
int mmdupd_(int *ehead, int *neqns, int *xadj, int *adjncy, 
	    int *delta, int *mdeg, int *dhead, int *dforw, 
	    int *dbakw, int *qsize, int *llist, int *marker, 
	    int *maxint, int *tag);
