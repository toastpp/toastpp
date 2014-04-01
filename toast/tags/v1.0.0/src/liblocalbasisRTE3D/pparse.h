// =========================================================================
// ParamParser: file parser for parameter files
// ============================================================================

#ifndef __PPARSE_H
#define __PPARSE_H

#include <fstream>
#include <iostream>

using namespace std;
class Raster;

class ParamParser {
public:
    ParamParser ();
    bool Open (char *fname);
    void LogOpen (char *fname);
    void Lineout (char *str);
    bool GetString (const char *cat, char *str);
    void PutString (const char *cat, const char *str);
    bool GetReal (const char *cat, double &val);
    void PutReal (const char *cat, double val);
    bool GetInt  (const char *cat, int &val);
    void PutInt  (const char *cat, int val);
    bool GetBool (const char *cat, bool &val);
    void PutBool (const char *cat, bool val);

private:
    char *trim_string (char *cbuf);
    bool bfile;
    ifstream prmf;
};

#endif // !__PPARSE_H
