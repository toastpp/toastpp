// ==========================================================================
// TOAST package v14 - FEM library
// Martin Schweiger and Simon R. Arridge
// param.cc
//
// Description:
// Definition of class Parameter: holds a single set of optical parameters
// (mua, kappa, n) and defines related methods
// ==========================================================================

#define FELIB_IMPLEMENTATION
#define __PARAM_CC

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

// ==========================================================================

double A_Keijzer (double n)
{
    // Source: M. Keijzer et al. AO 27, 1820-1824 (1988)

    double th    = asin (1.0/n);
    double costh = fabs (cos (th));
    double R0    = ((n-1.0)*(n-1.0)) / ((n+1.0)*(n+1.0));
    return (2.0/(1.0-R0) - 1.0 + costh*costh*costh) / (1.0 - costh*costh);
}

double A_Contini (double n)
{
    // Source: Contini et al. AO 36, 4587-4599 (1997)
    // Polynomial fit as described in appendix

    if (n > 1.0) {
        return 504.332889 + n * (-2641.00214 + n * (5923.699064 + n * (
	       -7376.355814 + n * (5507.53041 + n * (-2463.357945 + n * (
	       610.956547 - n * 64.8047))))));
    } else {
        return 3.084635 + n * (-6.531194 + n * (8.357854 + n * (
               -5.082751 + n * 1.171382)));
    }
}

double Afixed = 1.0;
double A_Const (double /*n*/)
{
    // Preset constant reflection term
    return Afixed;
}

double (*ReflectionTerm)(double) = A_Keijzer; // default

void SetReflectionFunction (ReflectionType rt, double A)
{
    switch (rt) {
    case REFLECTION_KEIJZER:
        ReflectionTerm = A_Keijzer;
	break;
    case REFLECTION_CONTINI:
        ReflectionTerm = A_Contini;
	break;
    case REFLECTION_CONST:
        Afixed = A;
        ReflectionTerm = A_Const;
	break;
    default:
        xERROR("Panic");
    }
}

double Parameter::C2A (ReflectionType rt, double n)
{
    typedef double(*A_FUNC)(double);
    A_FUNC A_func[3] = {A_Keijzer, A_Contini, A_Const};
    return c0/(2.0*n*A_func[rt](n));
}

Parameter &Parameter::operator= (const Parameter &prm)
{
    mua = prm.mua;
    kappa = prm.kappa;
    n = prm.n;
    a = prm.a;
    return *this;
}
    
int Parameter::operator== (const Parameter &prm) const
{
    return (mua == prm.mua && kappa == prm.kappa && n == prm.n);
}

int Parameter::operator!= (const Parameter &prm) const
{
    return (mua != prm.mua || kappa != prm.kappa || n != prm.n);
}

double Parameter::A() const { return (*ReflectionTerm)(n); }
double Parameter::bndA(double n2) const { return (*ReflectionTerm)(n/n2); }

double Parameter::Param (ParameterType prmtp) const
{
    switch (prmtp) {
    case PRM_MUA:    return Mua();
    case PRM_KAPPA:  return Kappa();
    case PRM_N:      return N();
    case PRM_CMUA:   return CMua();
    case PRM_CKAPPA: return CKappa();
    case PRM_MUS:    return Mus();
    case PRM_CMUS:   return CMus();
    case PRM_C:      return C();
    case PRM_ETA:    return CMua()/CKappa();
    case PRM_XI:     return 1.0/CKappa();
    case PRM_A:      return A();
    case PRM_C2A:    return C2A();
    default:         xERROR("Function not supported for this parameter type");
                     return 0.0; // dummy
    }
}

void Parameter::SetParam (ParameterType prmtp, const double val)
{
    switch (prmtp) {
    case PRM_MUA:    SetMua (val); return;
    case PRM_KAPPA:  SetKappa (val); return;
    case PRM_N:      SetN (val); return;
    case PRM_CMUA:   SetCMua (val); return;
    case PRM_CKAPPA: SetCKappa (val); return;
    case PRM_MUS:    SetMus (val); return;
    case PRM_CMUS:   xERROR("Set c_mus not implemented."); return;
    case PRM_C:      SetC (val); return;
    case PRM_A:      SetA (val); return;
    default:         xERROR("Function not supported for this parameter type");
    }
}

istream& operator>> (istream& is, Parameter& prm)
{
    // need to read A
    is >> prm.mua >> prm.kappa >> prm.n;
    return is;
}

ostream& operator<< (ostream& os, Parameter& prm)
{
    // need to write A
    os << prm.mua << ' ' << prm.kappa << ' ' << prm.n;
    return os;
}

void Parameter::get (istream &is, ParameterType p1, ParameterType p2,
		     ParameterType p3)
{
    dASSERT(p1 == PRM_MUA || p1 == PRM_CMUA, "Invalid parameter type 1.");
    dASSERT(p2 == PRM_KAPPA || p2 == PRM_CKAPPA || p2 == PRM_MUS ||
	   p2 == PRM_CMUS, "Invalid parameter type 2.");
    dASSERT(p3 == PRM_N || p3 == PRM_C, "Invalid parameter type 3.");
    double v1, v2, v3;
    is >> v1 >> v2 >> v3;
    SetParam (p1, v1);
    SetParam (p2, v2);
    SetParam (p3, v3);
}

void Parameter::put (ostream &os, ParameterType p1, ParameterType p2,
		     ParameterType p3)
{
    dASSERT(p1 == PRM_MUA || p1 == PRM_CMUA, "Invalid parameter type 1.");
    dASSERT(p2 == PRM_KAPPA || p2 == PRM_CKAPPA || p2 == PRM_MUS ||
	   p2 == PRM_CMUS, "Invalid parameter type 2.");
    dASSERT(p3 == PRM_N || p3 == PRM_C, "Invalid parameter type 3.");
    os << Param(p1) << ' ' << Param(p2) << ' ' << Param(p3);
}

void Swap (Parameter& prm1, Parameter& prm2)
{
    double tmp;
    tmp = prm1.mua, prm1.mua = prm2.mua, prm2.mua = tmp;
    tmp = prm1.kappa, prm1.kappa = prm2.kappa, prm2.kappa = tmp;
    tmp = prm1.n, prm1.n = prm2.n, prm2.n = tmp;
}

// ==========================================================================

ParameterList::ParameterList ()
{
    size = 0;
    list = 0;
    output_prmtp_p1 = PRM_MUA;
    output_prmtp_p2 = PRM_KAPPA;
    output_prmtp_p3 = PRM_N;
}

ParameterList::ParameterList (const int _size)
{
    list = new Parameter[size = _size];
    output_prmtp_p1 = PRM_MUA;
    output_prmtp_p2 = PRM_KAPPA;
    output_prmtp_p3 = PRM_N;
}

void ParameterList::New (const int _size)
{
    if (size) delete []list;
    if ((size = _size) != 0) {
	list = new Parameter[size];
	dASSERT(list, "Memory allocation failed.");
    } else list = 0;
}

void ParameterList::Clear ()
{
    New (0);
}

void ParameterList::SetList (int no, Parameter *prms)
{
    Clear();
    size = no;
    list = prms;
}

void ParameterList::Swap (const int rec1, const int rec2)
{
    dASSERT(rec1 >= 0 && rec1 < size && rec2 >= 0 && rec2 < size,
	   "Index out of range.");
    ::Swap (list[rec1], list[rec2]);
}

void ParameterList::Append (int number)
{
    Parameter *tmplist = new Parameter[size+number];
    dASSERT(tmplist, "Memory allocation failed");
    for (int i = 0; i < size; i++) tmplist[i] = list[i];
    if (size) delete []list;
    list = tmplist;
    size += number;
}

void ParameterList::Remove (const int rec)
{
    int i;

    dASSERT (rec >= 0 && rec < size, "Index out of range");
    Parameter *tmplist = new Parameter[size-1];
    for (i = 0; i < rec; i++) tmplist[i] = list[i];
    for (i = rec+1; i < size; i++) tmplist[i-1] = list[i];
    delete []list;
    list = tmplist;
    size--;
}

ParameterList &ParameterList::operator= (const ParameterList &plist)
{
    New (plist.Len());
    for (int i = 0; i < size; i++) {
	list[i] = plist.list[i];
    }
    output_prmtp_p1 = plist.output_prmtp_p1;
    output_prmtp_p2 = plist.output_prmtp_p2;
    output_prmtp_p3 = plist.output_prmtp_p3;
    plist_type = plist.plist_type;
    return *this;
}

RVector ParameterList::Param (ParameterType prmtp) const
{
    RVector prmvec(size);
    for (int i = 0; i < size; i++) prmvec[i] = list[i].Param (prmtp);
    return prmvec;
}

void ParameterList::Param (ParameterType prmtp, RVector &prmvec) const
{
    dASSERT(prmvec.Dim() == size, "Wrong size vector.");
    for (int i = 0; i < size; i++) prmvec[i] = list[i].Param (prmtp);
}

void ParameterList::SetParam (ParameterType prmtp, const RVector &prmvec)
{
    dASSERT(prmvec.Dim() == size, "Vector has wrong size.");
    for (int i = 0; i < size; i++) list[i].SetParam (prmtp, prmvec[i]);
}

void ParameterList::SetParam (ParameterType prmtp, double prm)
{
    for (int i = 0; i < size; i++) list[i].SetParam (prmtp, prm);
}

#ifdef FEM_DEBUG
Parameter& ParameterList::operator[] (const int rec) const
{
    dASSERT(rec >= 0 && rec < size, "Index out of range.");
    return list[rec];
}
#endif

istream& operator>> (istream& is, ParameterList& plist)
{
    char cbuf[256], valstr[30];
    int i, size = 0;
    ParameterType p1 = PRM_MUA, p2 = PRM_KAPPA, p3 = PRM_N;
    PlistType plist_type;
    do {   // find parameter list
	is.getline (cbuf, 256);
    } while (strncmp (cbuf, "[ParameterList]", 15));
    if(!strncmp(cbuf+15,"E",1))
      plist_type = PLIST_BY_ELEMENT;
    else
      plist_type = PLIST_BY_NODE;
    do {   // read parameter list header
	is.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "Size", 4)) {
	    xASSERT(sscanf (cbuf+4, "%d", &size) == 1, "Parse error.");
	} else if (!strncasecmp (cbuf, "Param1", 6)) {
	    xASSERT(sscanf (cbuf+6, "%s", valstr) == 1, "Parse error.");
	    if (!strcasecmp (valstr, "MUA")) p1 = PRM_MUA;
	    else if (!strcasecmp (valstr, "CMUA")) p1 = PRM_CMUA;
	    else xERROR("Parse error.");
	} else if (!strncasecmp (cbuf, "Param2", 6)) {
	    xASSERT(sscanf (cbuf+6, "%s", valstr) == 1, "Parse error.");
	    if (!strcasecmp (valstr, "KAPPA")) p2 = PRM_KAPPA;
	    else if (!strcasecmp (valstr, "CKAPPA")) p2 = PRM_CKAPPA;
	    else if (!strcasecmp (valstr, "MUS")) p2 = PRM_MUS;
	    else if (!strcasecmp (valstr, "CMUS")) p2 = PRM_CMUS;
	    else xERROR("Parse error.");
	} else if (!strncasecmp (cbuf, "Param3", 6)) {
	    xASSERT(sscanf (cbuf+6, "%s", valstr) == 1, "Parse error.");
	    if (!strcasecmp (valstr, "N")) p3 = PRM_N;
	    else if (!strcasecmp (valstr, "C")) p3 = PRM_C;
	    else xERROR("Parse error.");
	}
    } while (strncasecmp (cbuf, "Data", 4));

    xASSERT(size > 0, "No parameter list found or zero list size");
    plist.New (size);
    plist.plist_type = plist_type;
    for (i = 0; i < size; i++) plist[i].get (is, p1, p2, p3);  // read data
    return is;
}

void ParameterList::SetOutputParameterTypes (const ParameterType p1,
    const ParameterType p2, const ParameterType p3)
{
    if (p1 == PRM_MUA || p1 == PRM_CMUA)
	output_prmtp_p1 = p1;
    else xERROR("Invalid parameter type for P1");

    if (p2 == PRM_KAPPA || p2 == PRM_CKAPPA || p2 == PRM_MUS || p2 == PRM_CMUS)
	output_prmtp_p2 = p2;
    else xERROR("Invalid parameter type for P2");

    if (p3 == PRM_N || p3 == PRM_C)
	output_prmtp_p3 = p3;
    else xERROR("Invalid parameter type for P3");
}

ostream& operator<< (ostream& os, const ParameterList& plist)
{
    os << "[ParameterList]" << (plist.plist_type == PLIST_BY_ELEMENT ? 'E' : 'N') << endl;
    os << "Size " << plist.Len() << endl;
    os << "Param1 ";
    switch (plist.output_prmtp_p1) {
        case PRM_MUA: os << "MUA" << endl; break;
        case PRM_CMUA: os << "CMUA" << endl; break;
        default: xERROR("Invalid parameter 1 setting");
    }
    os << "Param2 ";
    switch (plist.output_prmtp_p2) {
        case PRM_KAPPA: os << "KAPPA" << endl; break;
        case PRM_CKAPPA: os << "CKAPPA" << endl; break;
        case PRM_MUS: os << "MUS" << endl; break;
        case PRM_CMUS: os << "CMUS" << endl; break;
        default: xERROR("Invalid parameter 2 setting");
    }
    os << "Param3 ";
    switch (plist.output_prmtp_p3) {
        case PRM_N: os << "N" << endl; break;
        case PRM_C: os << "C" << endl; break;
        default: xERROR("Invalid parameter 3 setting");
    }
    os << "Data" << endl;
    for (int i = 0; i < plist.Len(); i++) {
	plist.list[i].put (os, plist.output_prmtp_p1, plist.output_prmtp_p2,
	    plist.output_prmtp_p3);
	os << endl;
    }
    return os;
}
