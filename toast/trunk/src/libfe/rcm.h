typedef int integer;
typedef bool logical;
typedef double doublereal;

int perphn_(integer *n, integer *ja, integer *ia, integer *init,
	integer *iperm, integer *mask, integer *maskval, integer *nlev,
	integer *riord, integer *levels);
int bfs_(integer *n, integer *ja, integer *ia, integer *
	nfirst, integer *iperm, integer *mask, integer *maskval, integer *
	riord, integer *levels, integer *nlev);
int reversp_(integer *n, integer *riord);
