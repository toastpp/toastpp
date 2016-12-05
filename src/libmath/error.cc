// ==========================================================================
// Module mathlib
// File error.cpp
// Error handling routines
// ==========================================================================

#define MATHLIB_IMPLEMENTATION

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__GNUC__)
#include <sys/time.h>
#endif
#include <time.h>

#include "mathlib.h"

#define __ERROR_CC

using namespace std;

const char *TOAST_VERSION = "15";

#define ECHOLOG // echo log file entries on screen

#define EXPIRED_MSG "Sorry, the license for this program has expired.\n\
Please contact S.Arridge@cs.ucl.ac.uk or martins@medphys.ucl.ac.uk\n"

MATHLIB ofstream logfile;   // global log file handler
MATHLIB char logbuf[256];

MATHLIB int toastVerbosity = 0;

void DefaultErrorhandler (char *msg)
{
    if (msg)
        cerr << msg << endl;
    exit (1);
}
static void (*Errorhandler)(char*) = DefaultErrorhandler;

void SetErrorhandler (void (*ErrorFunc)(char*))
{
    Errorhandler = ErrorFunc;
}

MATHLIB void Error (const char *name, const char *file, int line)
{
    char cbuf[500];
    sprintf (cbuf, "**** Error in module mathlib ****\nIn %s\nFile %s line %d\n"
	, name, file, line);
    Errorhandler (cbuf);
}

MATHLIB void Error (const char *name, const char *file, int line, const char *msg, ...)
{
    cerr << "#########################" << endl;
    cerr << "ERROR IN LIBFE3:" << endl;
    va_list ap;
    va_start (ap,msg);
    vsnprintf (logbuf, 255, msg, ap);
    va_end(ap);
    cerr << logbuf << endl;
    cerr << name << endl << file << ", " << line << endl;
    cerr << "#########################" << endl;
    Errorhandler (NULL);
}

MATHLIB void Error_Undef (const char *name, const char *file, int line)
{
    char cbuf[500];
    sprintf (cbuf, "#########################\nERROR IN LIBFE3:\nFunction not implemented:\n%s\n%s, %d\n#########################\n",
	     name, file, line);
    Errorhandler (cbuf);
}
// ==========================================================================
// log file output routines

const int MAXLOGLEVELDEPTH = 10;
static char levelstr[MAXLOGLEVELDEPTH][50];
static int logleveldepth = 0;
static int progress_len;
static int progress_maxlen;
static int progress_maxcount;
static int logout_on = 1;

void LogoutOn()
{
    logout_on = 1;
}

void LogoutOff()
{
    logout_on = 0;
}

void LogOut (const char *msg, ...)
{
    if (logout_on) {
        va_list ap;
	va_start(ap,msg);
	vsnprintf (logbuf, 255, msg, ap);
	va_end(ap);
	for (int i = 0; i < logleveldepth; i++) logfile << "| ";
	logfile << logbuf << endl;
#ifdef ECHOLOG
	cout << logbuf << endl;
#endif
    }
}

void LogOut_Enter (const char *msg)
{
    if (logout_on) {
	for (int i = 0; i < logleveldepth; i++) logfile << "| ";
	logfile << "+-Enter: " << msg << endl;
	if (logleveldepth < MAXLOGLEVELDEPTH)
	    strcpy (levelstr[logleveldepth], msg);
    }
    logleveldepth++;
}

void LogOut_Exit()
{
    if (logleveldepth <= 0) return;  // something wrong
    logleveldepth--;
    if (logout_on) {
	for (int i = 0; i < logleveldepth; i++) logfile << "| ";
	logfile << "+-Leave: " << levelstr[logleveldepth] << endl;
    }
}

void LogOut_InitProgressbar (const char *name, int len, int maxcount)
{
    int i;
    size_t j;
    if (logout_on) {
	for (i = 0; i < logleveldepth; i++) logfile << "| ";
	logfile << name;
	for (j = strlen(name) + 2*logleveldepth; j < 25; j++) logfile << ' ';
	logfile << '[' << flush;
    }
    progress_len = 0;
    progress_maxlen = len;
    progress_maxcount = maxcount;
}

void LogOut_Progress (int count)
{
    int len = (progress_maxlen * (count+1)) / progress_maxcount;
    if (len > progress_len) {
	if (logout_on) {
	    for (int i = 0; i < len-progress_len; i++) logfile << '=';
	    logfile << flush;
	}
	progress_len = len;
    }
    if (logout_on && count == progress_maxcount-1) logfile << ']' << endl;
}

// ==========================================================================
// expiry-handling routines

void DefaultExpiryhandler ()
{
    cerr << EXPIRED_MSG;
    exit (1);
}
static void (*Expiryhandler)() = DefaultExpiryhandler;

void SetExpiryhandler (void (*ExpiryFunc)())
{
    Expiryhandler = ExpiryFunc;
}

int Expired (int month, int year)
{
    if (month < 0) return 0;
    struct tm expiry_date = {0};
    expiry_date.tm_mday = 1;
    expiry_date.tm_mon = month-1;
    expiry_date.tm_year = year-1900;
    return (time(0) > mktime (&expiry_date));
}

MATHLIB void Check_Expired (int month, int year)
{
    if (Expired(month, year)) Expiryhandler();
}

static char versionstr[256] = "";
static char fullversionstr[256] = "";

void SetVersion (const char *vstr)
{
    strcpy (versionstr, vstr);
}

const char *Version (const char *type, const char *date)
{
    sprintf (fullversionstr, "TOAST%s distribution [%s] - Build %s\n(c) Martin Schweiger and Simon Arridge\n", TOAST_VERSION, type, date);
#ifdef TOAST_PARALLEL
    strcat (fullversionstr, "Running parallel (pthreads) version\n");
#endif
    //if (versionstr[0]) sprintf (fullversionstr, "%s version %s, build %s",
    //	type, versionstr, date);
    //else sprintf (fullversionstr, "%s version, build %s", type, date);
    return fullversionstr;
}

void OutputProgramInfo (ostream &os)
{
    time_t tm = time(0);
    os << "Executed " << ctime(&tm);
    char *host = getenv ("HOST");
    if (host) os << "on host " << host << endl;
}

// ==========================================================================
// log file routines

static int logfile_isopen = 0;
static int logfile_isfile = 0;

void LogfileClose ()
{
    if (logfile_isopen) {
        if (logfile_isfile) logfile.close();
	logfile_isopen = 0;
    }
}

void LogfileOpen (const char *fname, bool rewind)
{
    if (logfile_isopen) {
        if (logfile_isfile && rewind)
	    LogfileClose();
	else
	    return;  // append to existing log file
    }
    if (fname[0]) {
        logfile.open (fname);
	logfile_isfile = 1;
    } else {
        xERROR("Function not supported");
        //logfile.attach (1);  // log to stdout
	logfile_isfile = 0;
    }
    logfile_isopen = 1;
}

// ==========================================================================
// flop counting routines

static unsigned int flop_add = 0;
static unsigned int flop_mul = 0;

void ResetFlops ()
{
    flop_add = 0;
    flop_mul = 0;
}

unsigned int FlopsAdd() { return flop_add; }
unsigned int FlopsMul() { return flop_mul; }
void Inc_FlopsAdd (unsigned int n) { flop_add += n; }
void Inc_FlopsMul (unsigned int n) { flop_mul += n; }
