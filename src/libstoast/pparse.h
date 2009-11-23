// =========================================================================
// ParamParser: file parser for parameter files
// ============================================================================

#ifndef __PPARSE_H
#define __PPARSE_H

#include <fstream>
#include <iostream>

class STOASTLIB ParamParser {
public:
    ParamParser ();
    bool Open (const char *fname);
    void LogOpen (const char *fname);
    void Lineout (const char *str);
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
    std::ifstream prmf;
};

#endif // !__PPARSE_H
