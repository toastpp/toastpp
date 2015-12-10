// -*-C++-*-
// ==========================================================================
// Objective function for parameter optimisation
// ==========================================================================

#ifndef __OF_H
#define __OF_H

#include "stoastlib.h"

class Raster;
class Solution;
template<class T>class TFwdSolver;

typedef enum {
    PRIOR_NONE_OLD,
    PRIOR_LAPLACIAN_OLD,
    PRIOR_DIAG_OLD
} PRIOR_OLD;

// =========================================================================

class STOASTLIB ObjectiveFunction {
public:
    ObjectiveFunction (const Raster *_raster);
    void SetData (const RVector *_data);
    void SetSD (const RVector *_sd);
    void SetScalingMatrix (const RMatrix *mat);
    void SetGMRFParam (double _gmrf_exp, double _hypermua, double _hyperkappa,
		       const int *_nbg_rowptr, const int *_nbg_colidx);
    void SetTVParam (double _tv_beta2, double _hypermua, double _hyperkappa);
    const RVector *Data () { return data; }
    const RVector *SD () { return sd; }
    double get_posterior (const RVector *proj) const;
    //double get_prior (const Solution *sol) const;
    static double get_prior (PRIOR_OLD prior, const RVector &x, const RVector &xp,
        const RVector &x0, double tau, const RCompRowMatrix &LTL);
    static double get_value (const RVector &data, const RVector &proj,
        const RVector &sd, const RMatrix *scale = 0);
    double get_value (const RVector *proj, const Solution *sol = 0) const;
    void add_gradient_data (RVector &grad, const Raster &raster,
	const CFwdSolver &FWS, const RVector &proj,
        CVector *dphi, const CCompRowMatrix &mvec) const;
    void add_gradient_prior (const Solution &so, RVector &grad) const;
    void get_gradient (const Raster &raster,
        const CFwdSolver &FWS, const RVector &proj,
        CVector *dphi, const CCompRowMatrix &mvec, const Solution *sol,
        RVector &grad) const;
    RVector get_gradient ();

private:
    const RVector *data;
    const RVector *sd;
    const RMatrix *scale;
    bool apply_prior;
    double gmrf_exp;
    double tv_beta2;
    double hyperp[2];
    const int *nbg_rowptr, *nbg_colidx;
    const Raster *raster;
};

inline void ObjectiveFunction::SetData (const RVector *_data)
{ data = _data; }

inline void ObjectiveFunction::SetSD (const RVector *_sd)
{ sd = _sd; }

inline void ObjectiveFunction::SetScalingMatrix (const RMatrix *mat)
{ scale = mat; }

inline void ObjectiveFunction::SetGMRFParam (double _gmrf_exp,
    double _hypermua, double _hyperkappa,
    const int *_nbg_rowptr, const int *_nbg_colidx)
{
    gmrf_exp = _gmrf_exp;
    hyperp[0] = _hypermua;
    hyperp[1] = _hyperkappa;
    nbg_rowptr = _nbg_rowptr;
    nbg_colidx = _nbg_colidx;
    apply_prior = true;
}

// =========================================================================
// this should probably go into ObjectiveFunction, or into a new class
// 'Prior'

STOASTLIB double penalty_tikh0 (const RVector &x, const RVector &x0,
    const RVector &xs);
STOASTLIB void penalty_gradient_tikh0_add (const RVector &x, const RVector &x0,
    const RVector &xs, RVector &grad, const double tau);
STOASTLIB void penalty_gradient_tikh0_rescale_add (const RVector &x, const RVector &x0,
    RVector &grad, const double tau);
STOASTLIB double penalty_tikh1 (const RVector &x, const RVector &x0,
    const double tau1, const double tau2,const RCompRowMatrix &LTL);
STOASTLIB void penalty_gradient_tikh1_add (const RVector &x, const RVector &x0,
    RVector &grad, const double tau1, const double tau2, 
    const RCompRowMatrix &LTL);

#endif // !__OF_H
