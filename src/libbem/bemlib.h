#ifndef __BEMLIB_H
#define __BEMLIB_H

#include "felib.h"

#ifdef BEMLIB_IMPLEMENTATION
#define BEMLIB DLLEXPORT
#else
#define BEMLIB DLLIMPORT
#endif

#include "bem_element.h"
#include "bem_region.h"
#include "bem_surface.h"
#include "bem_kernel.h"

const double eps=1.E-10;
const double pi=3.141592653589793;
const double pi1=1./(2.*pi);
const double pi4=1./(4.*pi);

#endif
