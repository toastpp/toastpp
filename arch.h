// =======================================================================
// Architecture-dependent features
// =======================================================================

#ifndef __ARCH_H
#define __ARCH_H

#if defined(WIN32)||defined(WIN64)

#include <process.h>
#include <direct.h>

// work around language incompatibilities
#define strcasecmp(str1,str2) _stricmp((str1),(str2))
#define strncasecmp(str1,str2,n) _strnicmp((str1),(str2),(n))

// avoid annoying warnings
#define hypot(x,y) _hypot((x),(y))
#define getpid() _getpid()
#define getcwd(buffer,maxlen) _getcwd(buffer,maxlen)
#define unlink(fname) _unlink(fname)

inline double drand48(void) { return (double)rand()/(double)RAND_MAX; }
// quick fix of missing function. this should be improved upon!

#endif

#endif // !__ARCH_H