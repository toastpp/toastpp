// =======================================================================
// Architecture-dependent features
// =======================================================================

#ifndef __TOASTARCH_H
#define __TOASTARCH_H

#if defined(_WIN32) || defined(_WIN64)

#include <process.h>
#include <direct.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>

#define __func__ __FUNCTION__

// work around language incompatibilities
#define strcasecmp(str1,str2) _stricmp((str1),(str2))
#define strncasecmp(str1,str2,n) _strnicmp((str1),(str2),(n))

#if _MSC_VER <= 1500 // version limit may be higher
namespace std {
	inline bool isnan(double arg) { return _isnan(arg); }
}
#endif
//#define std::isnan(arg) _isnan((arg))

// avoid annoying warnings
#define hypot _hypot
#define getpid() _getpid()
#define getcwd(buffer,maxlen) _getcwd(buffer,maxlen)
#define unlink(fname) _unlink(fname)

inline double drand48(void) { return (double)rand()/(double)RAND_MAX; }
// quick fix of missing function. this should be improved upon!

#endif // WINxx

#endif // !__TOASTARCH_H
