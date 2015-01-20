// -*-C++-*-
// ==========================================================================
// TOAST package v14 - FEM library
// Martin Schweiger and Simon R. Arridge
// param.h
//
// Description:
// * Declaration of class Parameter: holds a single set of optical parameters
//   (mua, kappa, n) and defines related methods
// * Declaration of class ParameterList: holds a list of parameter sets, e.g.
//   for all nodes (or elements) of a mesh. Contains methods to read and write
//   a list from/to a stream.
// ==========================================================================

#ifndef __PARAM_H
#define __PARAM_H

#ifdef UNDEF
// some typedefs
typedef enum {
    PRM_MUA,
    PRM_KAPPA,
    PRM_N,
    PRM_CMUA,
    PRM_CKAPPA,
    PRM_MUS,
    PRM_CMUS,
    PRM_C,
    PRM_ETA,   // canonical parameter: eta = cmua/ckappa
    PRM_XI,    // canonical parameter: xi  = 1/ckappa
    PRM_A,     // boundary reflection term A
    PRM_C2A,   // c/2A
    PRM_ZERO   // dummy
} ParameterType;
#endif
// some typedefs
typedef enum {
    PRM_CMUA,   //    
    PRM_CKAPPA, // primary parameters
    PRM_N,      //

    PRM_MUA,    //
    PRM_KAPPA,  //
    PRM_MUS,    // derived parameters
    PRM_CMUS,   //
    PRM_C,      //

    PRM_ETA,    // canonical parameter: eta = cmua/ckappa
    PRM_XI,     // canonical parameter: xi  = 1/ckappa
    PRM_A,      // boundary reflection term A
    PRM_C2A,    // c/2A
    PRM_ZERO    // dummy
} ParameterType;


typedef enum {
    REFLECTION_KEIJZER,
    REFLECTION_CONTINI,
    REFLECTION_CONST
} ReflectionType;


typedef enum {
    PLIST_BY_NODE,
    PLIST_BY_ELEMENT
} PlistType;


// general conversion mus <-> kappa
inline double MuaMusToKappa (const double mua, const double mus)
{ return 1.0 / (3.0*(mua+mus)); }

inline double MuaKappaToMus (const double mua, const double kappa)
{ return 1.0/(3.0*kappa) - mua; }

// default parameter settings
const double c0            = 0.3;   // speed of light in vacuum [mm/ps]
const double default_n     = 1.4;   // default refractive index of medium
const double default_mua   = 0.025; // default absorption coefficient
const double default_mus   = 2.0;   // default (reduced) scatter coefficient
const double default_kappa = MuaMusToKappa (default_mua, default_mus);
const double default_a     = -1.0;  // "undefined"

// some useful conversions
inline double lspeed (double n) { return c0/n; }
void SetReflectionFunction (ReflectionType rt, double A = 1.0);
// select reflection calculation method
// A is only used for reflection type REFLECTION_CONST

// some algorithms for reflection term calculations
FELIB double A_Keijzer (double n);
FELIB double A_Contini (double n);
FELIB double A_Const (double /*n*/);

// ==========================================================================

class FELIB Parameter {
public:
    Parameter ()
    { mua = default_mua, kappa = default_kappa, n = default_n, a = default_a; }

    Parameter (const double _mua, const double _kappa, const double _n,
	       const double _a = default_a)
    { mua = _mua, kappa = _kappa, n = _n, a = _a; }

    Parameter (const Parameter &prm)
    { mua = prm.mua, kappa = prm.kappa, n = prm.n, a = prm.a; }

    // basic conversions mus <-> kappa
    friend double MuaMusToKappa (const double mua, const double mus);
    friend double MuaKappaToMus (const double mua, const double kappa);

    // some operators
    Parameter &operator= (const Parameter &prm);
    int operator== (const Parameter &prm) const;
    int operator!= (const Parameter &prm) const;

    // read parameters
    double Mua()    const { return mua; }
    double Kappa()  const { return kappa; }
    double N()      const { return n; }
    double C() const { return lspeed(n); }
    double CMua()   const { return C()*Mua(); }
    double CKappa() const { return C()*Kappa(); }
    double Mus()    const { return MuaKappaToMus (mua, kappa); }
    double CMus()   const { return C()*Mus(); }
    double A() const;
    double bndA (double n2 = 1.0) const;
    double C2A() const { return C()/(2.0*A()); }
    double Param (ParameterType prmtp) const;

    static double C2A (ReflectionType rt, double n);

    // set parameters
    void SetMua    (const double _mua) { mua = _mua; }
    void SetKappa  (const double _kappa) { kappa = _kappa; }
    void SetN      (const double _n) { n = _n; }
    void SetC      (const double _c) { n = c0/_c; }
    void SetMus    (const double _mus) { kappa = MuaMusToKappa (mua, _mus); }
    void SetCMua   (const double _cmua) { mua = _cmua/C(); }
    void SetCKappa (const double _ckappa) { kappa = _ckappa/C(); }
    void SetA      (const double _a) { a = _a; }
    void SetMuaKappa (const double _mua, const double _kappa)
        { mua = _mua, kappa = _kappa; }
    void SetMuaMus (const double _mua, const double _mus)
        { mua = _mua, kappa = MuaMusToKappa (_mua, _mus); }
    void SetParam (ParameterType prmtp, const double val);

    // I/O
    friend std::istream& operator>> (std::istream& is, Parameter& prm);
    friend std::ostream& operator<< (std::ostream& os, Parameter& prm);
    void get (std::istream &is, ParameterType p1 = PRM_MUA,
	      ParameterType p2 = PRM_KAPPA, ParameterType p3 = PRM_N);
    void put (std::ostream &os, ParameterType p1 = PRM_MUA,
	      ParameterType p2 = PRM_KAPPA, ParameterType p3 = PRM_N);

    // miscellaneous
    friend void Swap (Parameter& prm1, Parameter& prm2);

private:
    double mua;    // mua
    double kappa;  // kappa
    double n;      // refractive index
    double a;      // reflection term
};

// ==========================================================================

class FELIB ParameterList {
public:
    ParameterList ();
    ParameterList (const int _size);
    ~ParameterList() { if (size) delete []list; }

    void New (const int _size);
    // create new list with '_size' entries and set to default values

    void Clear();
    // clear the list

    int Len() const { return size; }

    void SetList (int no, Parameter *prms);
    // replace parameter list with prms of length no

    void Swap (const int rec1, const int rec2);
    // swap entries rec1 and rec2 in the list

    void Append (int number);
    // Append 'number' empty to the list. The parameters of the new entries
    // are set to default values.

    void Remove (const int rec);
    // remove entry `rec' from the list

    ParameterList &operator= (const ParameterList &plist);
    // make *this a copy of plist

    // list type 
    PlistType pltype() const { return plist_type;}
    void SetListType(PlistType p) { plist_type = p;}
    // read parameter vectors
    RVector Mua()    const { return Param (PRM_MUA); }
    RVector Kappa()  const { return Param (PRM_KAPPA); }
    RVector N()      const { return Param (PRM_N); }
    RVector C()      const { return Param (PRM_C); }
    RVector CMua()   const { return Param (PRM_CMUA); }
    RVector CKappa() const { return Param (PRM_CKAPPA); }
    RVector Mus()    const { return Param (PRM_MUS); }
    RVector CMus()   const { return Param (PRM_CMUS); }
    RVector A()      const { return Param (PRM_A); }
    RVector C2A()    const { return Param (PRM_C2A); }
    RVector Param (ParameterType prmtp) const;

    // as above but use vector supplied in parameter list. This avoids
    // creation of a temporary vector
    void Mua (RVector &prmvec) const { Param (PRM_MUA, prmvec); }
    void Kappa (RVector &prmvec) const { Param (PRM_KAPPA, prmvec); }
    void N (RVector &prmvec) const { Param (PRM_N, prmvec); }
    void C (RVector &prmvec) const { Param (PRM_C, prmvec); }
    void CMua (RVector &prmvec) const { Param (PRM_CMUA, prmvec); }
    void CKappa (RVector &prmvec) const { Param (PRM_CKAPPA, prmvec); }
    void Mus (RVector &prmvec) const { Param (PRM_MUS, prmvec); }
    void CMus (RVector &prmvec) const { Param (PRM_CMUS, prmvec); }
    void A (RVector &prmvec) const { Param (PRM_A, prmvec); }
    void C2A (RVector &prmvec) const { Param (PRM_C2A, prmvec); }
    void Param (ParameterType prmtp, RVector &prmvec) const;

    // set parameter vectors. Vector must be of the same size as parameter list
    void SetMua    (const RVector& prmvec) { SetParam (PRM_MUA, prmvec); }
    void SetKappa  (const RVector& prmvec) { SetParam (PRM_KAPPA, prmvec); }
    void SetN      (const RVector& prmvec) { SetParam (PRM_N, prmvec); }
    void SetC      (const RVector& prmvec) { SetParam (PRM_C, prmvec); }
    void SetCMua   (const RVector& prmvec) { SetParam (PRM_CMUA, prmvec); }
    void SetCKappa (const RVector& prmvec) { SetParam (PRM_CKAPPA, prmvec); }
    void SetMus    (const RVector& prmvec) { SetParam (PRM_MUS, prmvec); }
    void SetA      (const RVector& prmvec) { SetParam (PRM_A, prmvec); }
    void SetParam  (ParameterType prmtp, const RVector& prmvec);

    // set parameter vector to homogeneous value
    void SetMua    (double prm) { SetParam (PRM_MUA, prm); }
    void SetKappa  (double prm) { SetParam (PRM_KAPPA, prm); }
    void SetN      (double prm) { SetParam (PRM_N, prm); }
    void SetC      (double prm) { SetParam (PRM_C, prm); }
    void SetCMua   (double prm) { SetParam (PRM_CMUA, prm); }
    void SetCKappa (double prm) { SetParam (PRM_CKAPPA, prm); }
    void SetMus    (double prm) { SetParam (PRM_MUS, prm); }
    void SetA      (double prm) { SetParam (PRM_A, prm); }
    void SetParam  (ParameterType prmtp, double prm);

#ifdef FEM_DEBUG
    Parameter& operator[] (const int rec) const;
#else
    Parameter& operator[] (const int rec) const { return list[rec]; }
#endif

    // I/O
    void SetOutputParameterTypes (const ParameterType p1,
	const ParameterType p2, const ParameterType p3);

    friend std::istream& operator>> (std::istream& is, ParameterList& plist);
    friend std::ostream& operator<< (std::ostream& os,
        const ParameterList& plist);

private:
    int size;         // list length
    Parameter *list;  // array of parameters
    ParameterType output_prmtp_p1, output_prmtp_p2, output_prmtp_p3;
    PlistType plist_type;  // type of list (by node (default), or by element)
};

#endif // !__PARAM_H
