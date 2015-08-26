// =======================================================================
// Architecture-dependent features
// =======================================================================

#ifndef __ARCH_H
#define __ARCH_H

#if defined(_WIN32) || defined(_WIN64)

#include <process.h>
#include <direct.h>
#include <float.h>
#include <time.h>

#define __func__ __FUNCTION__

// work around language incompatibilities
#define strcasecmp(str1,str2) _stricmp((str1),(str2))
#define strncasecmp(str1,str2,n) _strnicmp((str1),(str2),(n))
#define isnan(arg) _isnan((arg))

// avoid annoying warnings
#define hypot _hypot
#define isnan(x) _isnan((x))
#define getpid() _getpid()
#define getcwd(buffer,maxlen) _getcwd(buffer,maxlen)
#define unlink(fname) _unlink(fname)

inline double drand48(void) { return (double)rand()/(double)RAND_MAX; }
// quick fix of missing function. this should be improved upon!

#endif //_WIN32

#endif // !__ARCH_H