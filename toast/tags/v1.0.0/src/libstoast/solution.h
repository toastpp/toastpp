// -*-C++-*-
// ==========================================================================
// Solution class
// ==========================================================================

#ifndef __SOLUTION_H
#define __SOLUTION_H

#include "mathlib.h"

// definition of parameters for standard OT
#define OT_NPARAM     4
#define OT_CMUA       0
#define OT_CKAPPA     1
#define OT_C2A        2
#define OT_N          3

// definition of parameters for fluorescence imaging
#define FLU_NPARAM    7
#define FLU_CMUA_EX   0
#define FLU_CKAPPA_EX 1
#define FLU_C2A_EX    2
#define FLU_CMUA_FL   3
#define FLU_CKAPPA_FL 4
#define FLU_C2A_FL    5
#define FLU_YIELD     6

// =========================================================================
// Prototypes

class STOASTLIB Solution;

STOASTLIB double l2norm (const Solution &sol);
STOASTLIB Solution exp (const Solution &sol);
STOASTLIB Solution log (const Solution &sol);

// =========================================================================

class STOASTLIB Solution {
    friend class Raster;

public:
    //enum PARAM { CMUA, CKAPPA, C2A };

    Solution (int nparam, int n);
    Solution (const Solution &sol);
    ~Solution();

    int nParam() const { return nprm; }
    void SetParam (int which, const RVector &prm);
    void SetParam (int which, double prm) { param[which] = prm; }
    void Set (const RVector &prm);
    const RVector GetParam (int which) const;
  
    const RVector GetActiveParams () const;
    // returns concatenated vector of all active parameters

    void SetActiveParams (const RVector &prm);
    // set active parameters from concatenated vector

    void SetActive (int which, bool act);
    inline bool IsActive (int which) const { return active[which]; }
    inline int nActive() const { return nactive; }

    int ActiveDim() const;
    // returns the summed dimensions of all active components

    void Extents (int which, double &vmin, double &vmax) const;
    bool Bounded (int which, double vmin, double vmax);
    bool Valid ();

    //void SetParameterScale (PARAM which, double scale);
    //void ScaleParams ();
    //void UnscaleParams ();

    void Scale (int which, double factor);
    void Add (const Solution &sol);
    void vmul (const Solution &sol);
    void vmul (const RVector &v) const;
    Solution &operator= (const Solution &sol);

    Solution operator* (double s) const;
    Solution &operator*= (double s);
    // multiplication of _active_ parameters with s

    Solution operator+ (const Solution &sol) const;
    Solution operator+ (const RVector &v) const;
    Solution &operator*= (const RVector &v);
    Solution &operator+= (const Solution &sol);
    Solution &operator+= (const RVector &v);

    void WriteImgGeneric (int imgno, const char *filename, 
        int prmind, bool append = true);
    // write param[prmind] to image file

    static void WriteImgGeneric (int imgno, const char *filename,
        const RVector &img, bool append = true);
    // write img to image file

    void WriteImg_mua (int imgno, const char *nimname, bool append = true);
    void WriteImg_mus (int imgno, const char *nimname, bool append = true);

    friend STOASTLIB double l2norm (const Solution &sol);
    friend STOASTLIB Solution exp (const Solution &sol);
    friend STOASTLIB Solution log (const Solution &sol);

protected:
    int nprm;          // number of components
    RVector *param;    // parameter components
    bool *active;      // activation flag
    int nactive;       // number of active components
};

inline Solution log (const Solution &sol)
{
    Solution res(sol);
    for (int i = 0; i < sol.nprm; i++)
        if (res.active[i]) res.param[i] = log(sol.param[i]);
    return res;
}

typedef Solution(*ParamTransform)(const Solution &sol);

#endif // !__SOLUTION_H
