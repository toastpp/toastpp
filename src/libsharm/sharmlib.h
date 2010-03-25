#ifndef __SHARMLIB_H
#define __SHARMLIB_H

#include "mathlib.h"
#include "felib.h"

// Symbol import/export direction
#ifdef SHARMLIB_IMPLEMENTATION
#define SHARMLIB DLLEXPORT
#else
#define SHARMLIB DLLIMPORT
#endif


 #include "diffusphere.h"
//#include "optimizzz.h"
//#include "usefulfan.h"
//#include "legendre.h"
//#include "pplot3d.h"


const double M_PI=3.141592653589793;


#endif // !__SHARMLIB_H
