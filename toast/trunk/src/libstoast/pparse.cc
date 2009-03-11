// ============================================================================
// ParamParser: file parser for parameter files
// ============================================================================

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include <fstream>

using namespace std;

ParamParser::ParamParser ()
{
    bfile = false;
}

bool ParamParser::Open (char *fname)
{
    prmf.open (fname);
    return bfile = prmf.good();
}

void ParamParser::LogOpen (char *fname)
{
    LogfileOpen (fname);
    LogOut ("================== TOAST LOG ==================");
    LogOut (fname);
}

char *ParamParser::trim_string (char *cbuf)
{
    char *c;

    // strip comments starting with ';'
    for (c = cbuf; *c; c++) {
        if (*c == ';') {
	    *c = '\0';
	    break;
	}
    }
    // strip trailing white space
    for (--c; c >= cbuf; c--) {
        if (*c == ' ' || *c == '\t') *c = '\0';
	else break;
    }
    // skip leading white space
    for (c = cbuf; *c; c++)
        if (*c != ' ' && *c != '\t') return c;

    // should never get here
    return c;
}

void ParamParser::Lineout (char *str)
{
    LogOut (str);
}

bool ParamParser::GetString (const char *cat, char *str)
{
    char cbuf[256];
    int i, ok;
    
    if (!bfile) return false;
    prmf.clear();
    prmf.seekg (0, ios::beg);
    while ((ok = (prmf.getline (cbuf, 256) != 0)) &&
	   strncasecmp (cbuf, cat, strlen(cat)));
    if (!ok) {
        prmf.clear();
	return false;
    }
    for (i = 0; cbuf[i] && cbuf[i] != '='; i++);
    if (!cbuf[i]) return false;
    strcpy (str, trim_string (cbuf+i+1));
    return true;
}

void ParamParser::PutString (const char *cat, const char *str)
{
    char cbuf[256];
    sprintf (cbuf, "%s = %s", cat, str);
    Lineout (cbuf);
}

bool ParamParser::GetReal (const char *cat, double &val)
{
    char cbuf[256];
    if (!GetString (cat, cbuf)) return false;
    return (sscanf (cbuf, "%lf", &val) == 1);
}

void ParamParser::PutReal (const char *cat, double val)
{
    char cbuf[256];
    sprintf (cbuf, "%s = %g", cat, val);
    Lineout (cbuf);
}

bool ParamParser::GetInt (const char *cat, int &val)
{
    char cbuf[256];
    if (!GetString (cat, cbuf)) return false;
    return (sscanf (cbuf, "%d", &val) == 1);
}

void ParamParser::PutInt (const char *cat, int val)
{
    char cbuf[256];
    sprintf (cbuf, "%s = %d", cat, val);
    Lineout (cbuf);
}

bool ParamParser::GetBool(const char *cat, bool &val)
{
    char cbuf[256];
    if (!GetString (cat, cbuf)) return false;
    if (!strcasecmp (cbuf, "true")) { val = true; return true; }
    else if (!strcasecmp (cbuf, "false")) { val = false; return true; }
    return false;
}

void ParamParser::PutBool (const char *cat, bool val)
{
    char cbuf[256];
    sprintf (cbuf, "%s = %s", cat, val ? "TRUE":"FALSE");
    Lineout (cbuf);
}
