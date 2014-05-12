// Generate file containing initial boundary displacements
// defined by user input

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "felib.h"

using namespace std;
const char *WS = " \t";
typedef enum {
  OP_GT,
  OP_LT,
  OP_GE,
  OP_LE,
  OP_AND,
  OP_OR,
  OP_NIL
} Operator;

void ParseError();
void SelectNodes (const Mesh &mesh, bool *sel, int var, Operator op,
		  double val, Operator bop, bool bnd_only = true);

int main (void)
{
    char fname[256], cbuf[256], *s;
    int i, j, n, nn, cmd, expect, var, res;
    double val;
    bool *sel;
    Operator op, bop;
    Mesh mesh;

    cout << "initdisp: Generate initial boundary displacements\n";
    cout << "          or boundary forces\n\n";
    cout << "Mesh file: ";
    cin >> fname;
    ifstream ifs (fname);
    ifs >> mesh;
    mesh.Setup();
    n = mesh.nlen();
    sel = new bool[n];

    cout << "Output file: ";
    cin >> fname;
    ofstream ofs (fname);

    for (;;) {
        cout << "Add boundary displacements/forces (1/0): ";
	cin >> cmd;
	if (!cmd) break;
	cin.getline (cbuf, 256);
	cout << "Define surface region (tokens separated by w.s.):\n";
	cout << "(e.g. x >= 0 and y < 1)\n";
	cin.getline (cbuf, 256);
	expect = 0; // expect a variable
	s = strtok (cbuf, WS);
	bop = OP_NIL;
	while (s) {
	    switch (expect) {
	    case 0: // read variable
	        var = toupper (s[0]) - 'X';
		if (var < 0 || var > 2) ParseError();
		expect = 1; // expect an operator
		cout << (char)(var + 'X') << ' ';
		break;
	    case 1: // read operator
	        if      (!strcmp (s, ">"))  op = OP_GT;
		else if (!strcmp (s, ">=")) op = OP_GE;
		else if (!strcmp (s, "<"))  op = OP_LT;
		else if (!strcmp (s, "<=")) op = OP_LE;
		else    ParseError();
		expect = 2; // expect a value;
		cout << s << ' ';
		break;
	    case 2: // read value
	        res = sscanf (s, "%lf", &val);
		if (res < 1) ParseError();
		SelectNodes (mesh, sel, var, op, val, bop);
		expect = 3; // expect a boolean operator
		cout << val << ' ';
		break;
	    case 3: // read boolean operator
	        if      (!strcasecmp (s, "AND")) bop = OP_AND;
		else if (!strcasecmp (s, "OR")) bop = OP_OR;
		else    ParseError();
		expect = 0;
		cout << s << ' ';
		break;
	    }
	    s = strtok (NULL, WS);
	}
	cout << endl;
	cout << "Define displacement [A] [x y z]\n";
	cout << "Use '-' to indicate no displacement for an axis\n";
	cout << "Use leading 'A' to indicate absolute position\n";
	cin.getline (cbuf, 256);
	bool isabsolute;
	char *str, dsp[3][64];
	int res;
	if (toupper (cbuf[0]) == 'A') {
	    isabsolute = true;
	    str = cbuf+2;
	} else {
	    isabsolute = false;
	    str = cbuf;
	}
	res = sscanf (str, "%s%s%s", dsp[0], dsp[1], dsp[2]);
	nn = 0;
	for (i = 0; i < n; i++)
	    if (sel[i]) {
	        ofs << i << ':';
		for (j = 0; j < res; j++) {
		    if (!strcmp(dsp[j], "-")) ofs << " -";
		    else if (isabsolute)
		        ofs << ' ' << atof(dsp[j])-mesh.nlist[i][j];
		    else
		        ofs << ' ' << atof(dsp[j]);
		}
		ofs << endl;
		      
	        //ofs << i << ": " << cbuf << endl;
		nn++;
	    }
	cout << "Displaced " << nn << " nodes\n";
    }
    delete []sel;
    return 0;
}

void ParseError ()
{
    cerr << "Parse error" << endl;
    exit (1);
}

void SelectNodes (const Mesh &mesh, bool *sel, int var, Operator op,
		  double val, Operator bop, bool bnd_only)
{
    int i, n = mesh.nlen();
    bool match;
    for (i = 0; i < n; i++) {
        if (bnd_only && !mesh.nlist[i].isBnd()) continue;
        switch (op) {
	case OP_GT: match = (mesh.nlist[i][var] >  val); break;
	case OP_GE: match = (mesh.nlist[i][var] >= val); break;
	case OP_LT: match = (mesh.nlist[i][var] <  val); break;
	case OP_LE: match = (mesh.nlist[i][var] <= val); break;
	}
	switch (bop) {
	case OP_NIL: sel[i] = match;  break;
	case OP_AND: sel[i] &= match; break;
	case OP_OR:  sel[i] |= match; break;
	}
    }
}
