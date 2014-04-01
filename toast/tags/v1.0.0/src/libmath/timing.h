#ifndef __TIMING_H
#define __TIMING_H

MATHLIB double tic();
MATHLIB double toc();
MATHLIB double toc(double tic);

MATHLIB double walltic();
MATHLIB double walltoc();
MATHLIB double walltoc(double walltic);

#ifdef DBG_TIMING
#define TIC (tic())
#define TOC(t) (t=toc())
#define TOCADD(t) (t+=toc())
#else
#define TIC
#define TOC(t)
#define TOCADD(t)
#endif

#endif // !__TIMING_H
