#define MATHLIB_IMPLEMENTATION

#include "mathlib.h"
#include "timing.h"

#if defined(WIN32) || defined(WIN64)
#include <windows.h>
#else
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>
#endif
#include <time.h>

// ==========================================================================
// prototypes

static double tick_per_sec();

// ==========================================================================
// globals

static double tps = tick_per_sec();
static double g_tic = 0;
static struct timeval g_tv = {0,0};

// ==========================================================================
// Returns clock ticks per second (system dependent)

static double tick_per_sec ()
{
    clock_t ticks_per_second;

#if defined (WIN32) || defined(WIN64)
    ticks_per_second = 100;
#else
#ifndef CLK_TCK
    ticks_per_second = sysconf(_SC_CLK_TCK);
#else
    ticks_per_second = CLK_TCK;
#endif
#endif
    return ticks_per_second;
}

// ==========================================================================
// Converts clock ticks to seconds

inline double cpu2sec (clock_t cputime)
{
    return (double)cputime / tps;
}

inline double getclock ()
{
#if defined (WIN32) || defined (WIN64)
    return cpu2sec (clock());
#else
    struct tms tm;
    times (&tm);
    return cpu2sec (tm.tms_utime);
#endif
}

// ==========================================================================
// 

double tic ()
{
    return g_tic = getclock();
}

double toc ()
{
    return getclock() - g_tic;
}

double toc (double tic)
{
    return getclock() - tic;
}

// ==========================================================================
// 

double walltic ()
{
#if defined (WIN32) || defined (WIN64)
	return tic();
#else
	gettimeofday(&g_tv, NULL);
    return g_tv.tv_sec + g_tv.tv_usec*1e-6;
#endif
}

double walltoc ()
{
#if defined (WIN32) || defined (WIN64)
	return toc();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)(tv.tv_sec-g_tv.tv_sec) +
	(double)(tv.tv_usec-g_tv.tv_usec)*1e-6;
#endif
}

double walltoc (double walltic)
{
#if defined (WIN32) || defined (WIN64)
	return toc(walltic);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec*1e-6 - walltic;
#endif
}
