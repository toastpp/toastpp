#ifndef __SUPERTOAST_CW_MW_H
#define __SUPERTOAST_CW_MW_H

// ==========================================================================
// Context data for line search callback function

struct OF_CLBK_DATA {
    RFwdSolverMW *fws;
    const Raster *raster;
    const Scaler *pscaler;
    MWsolution *meshsol;
    const Regularisation *reg;
    const RCompRowMatrix *qvec;
    const RCompRowMatrix *mvec;
    const RVector *data;
    const RVector *sd;
};

// ==========================================================================
// The callback function for obtaining the objective function during
// line search

double of_clbk (const RVector &x, double *of_sub, void *context);

#endif // !__SUPERTOAST_CW_MW_H
