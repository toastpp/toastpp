// -*-C++-*-
// ==========================================================================
// Module mathlib
// File error.h
// Error handling routines
// ==========================================================================

#ifndef __ERROR_H
#define __ERROR_H

#include "toastdef.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdarg.h>

#ifndef VERBOSE_LEVEL
#ifdef FEM_DEBUG
#define VERBOSE_LEVEL 0
#else
#define VERBOSE_LEVEL 0
#endif // FEM_DEBUG
#endif // !VERBOSE_LEVEL

#ifndef __FUNCTION__
#define __FUNCTION__ ""
#endif

#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ ""
#endif

#define ERROR_UNDEF Error_Undef(__PRETTY_FUNCTION__,__FILE__,__LINE__)

// unconditional error checking
#define xERROR(msg, ...) Error(__PRETTY_FUNCTION__,__FILE__,__LINE__, msg, ##__VA_ARGS__)
#define xASSERT(cond, msg, ...) if (!(cond)) xERROR(msg, ##__VA_ARGS__ )
#define xRANGE_CHECK(cond) if (!(cond)) xERROR("Index out of range.");

// error checking in debug environment
#ifdef FEM_DEBUG
#define dERROR(msg, ...) Error(__PRETTY_FUNCTION__,__FILE__,__LINE__, msg, ##__VA_ARGS__)
#define dASSERT(cond, msg, ...) if (!(cond)) dERROR(msg, ##__VA_ARGS__)
#define RANGE_CHECK(cond) if (!(cond)) dERROR("Index out of range.");
#else
#define dERROR(msg,...)
#define dASSERT(cond,msg,...)
#define RANGE_CHECK(cond)
#endif

// output to log file for different verbose levels
#define LOGOUT_ON LogoutOn()
#define LOGOUT_OFF LogoutOff()
#define LOGOUT(msg, ...) LogOut(msg, ##__VA_ARGS__)
#define LOGOUT_ENTER LogOut_Enter(__FUNCTION__)
#define LOGOUT_OPEN(msg) LogOut_Enter(msg)
#define LOGOUT_EXIT LogOut_Exit()
#define LOGOUT_INIT_PROGRESSBAR(name,len,maxcount) \
    LogOut_InitProgressbar(name,len,maxcount)
#define LOGOUT_PROGRESS(count) LogOut_Progress(count)

#if (VERBOSE_LEVEL >= 1)
#define LOGOUT1(msg,...) LogOut(msg, ##__VA_ARGS__)
#define LOGOUT1_ENTER LogOut_Enter(__FUNCTION__)
#define LOGOUT1_OPEN(msg) LogOut_Enter(msg)
#define LOGOUT1_EXIT LogOut_Exit()
#define LOGOUT1_INIT_PROGRESSBAR(name,len,maxcount) \
    LogOut_InitProgressbar(name,len,maxcount)
#define LOGOUT1_PROGRESS(count) LogOut_Progress(count)
#else
#define LOGOUT1(msg,...)
#define LOGOUT1_ENTER
#define LOGOUT1_OPEN(msg)
#define LOGOUT1_EXIT
#define LOGOUT1_INIT_PROGRESSBAR(name,len,maxcount)
#define LOGOUT1_PROGRESS(count)
#endif

#if (VERBOSE_LEVEL >= 2)
#define LOGOUT2(msg,...) LogOut(msg, ##__VA_ARGS__)
#define LOGOUT2_ENTER LogOut_Enter(__FUNCTION__)
#define LOGOUT2_OPEN(msg) LogOut_Enter(msg)
#define LOGOUT2_EXIT LogOut_Exit()
#define LOGOUT2_INIT_PROGRESSBAR(name,len,maxcount) \
    LogOut_InitProgressbar(name,len,maxcount)
#define LOGOUT2_PROGRESS(count) LogOut_Progress(count)
#else
#define LOGOUT2(msg,...)
#define LOGOUT2_ENTER
#define LOGOUT2_OPEN(msg)
#define LOGOUT2_EXIT
#define LOGOUT2_INIT_PROGRESSBAR(name,len,maxcount)
#define LOGOUT2_PROGRESS(count)
#endif

// version string macros and functions
void SetVersion (const char *vstr);
MATHLIB const char *Version (const char *type, const char *date);
void OutputProgramInfo (std::ostream &os);

#ifdef FEM_DEBUG
#define VERSION_STRING Version("DEBUG", __DATE__)
#else
#define VERSION_STRING Version("RELEASE", __DATE__)
#endif

// handling of template instantiation
#ifdef __BORLANDC__
#define INSTTEMPLATE
#else
#define INSTTEMPLATE template
#endif

// error routines
MATHLIB void SetErrorhandler (void (*ErrorFunc)(char*));
MATHLIB void Error (const char *name, const char *file, int line);
MATHLIB void Error (const char *name, const char *file, int line, const char *msg, ...);
MATHLIB void Error_Undef (const char *name, const char *file, int line);

// log file output routines
void LogoutOn();                   // turn on output to log file
void LogoutOff();                  // turn off output to log file
MATHLIB void LogfileOpen (const char *fname, bool rewind = false);
void LogfileClose ();              // close an open log file
MATHLIB void LogOut (const char *msg, ...);  // write message `msg' to log file
void LogOut_Enter (const char *routine); // write `enter' note to log file
void LogOut_Exit ();               // write `leave routine' note to log file
MATHLIB void LogOut_InitProgressbar (const char *name, int len, int maxcount);
    // initialise output of progress bar of `len' characters, to monitor
    // progress up to `maxcount'
MATHLIB void LogOut_Progress (int count);
    // set progress bar to show state `count' of `maxcount'
#ifndef __ERROR_CC
MATHLIB extern int toastVerbosity;  // global verbosity flag
MATHLIB extern char logbuf[256];
MATHLIB extern std::ofstream logfile;
#endif

// expiry-checking routines
void SetExpiryhandler (void (*ExpiryFunc)());
MATHLIB void Check_Expired (int month, int year);
#ifdef EXPIRY_YEAR
#ifndef EXPIRY_MONTH
#define EXPIRY_MONTH 1
#endif // !EXPIRY_MONTH
#define CHECK_EXPIRED() Check_Expired(EXPIRY_MONTH,EXPIRY_YEAR)
#else
#define CHECK_EXPIRED()
#endif // EXPIRY_YEAR

// flop counting routines
void ResetFlops();
unsigned int FlopsAdd();
unsigned int FlopsMul();
void Inc_FlopsAdd (unsigned int n);
void Inc_FlopsMul (unsigned int n);
#ifdef COMPUTE_FLOPS
#define INC_FLOPS_ADD(n) Inc_FlopsAdd(n)
#define INC_FLOPS_MUL(n) Inc_FlopsMul(n)
#else
#define INC_FLOPS_ADD(n)
#define INC_FLOPS_MUL(n)
#endif

#endif // !__ERROR_H
