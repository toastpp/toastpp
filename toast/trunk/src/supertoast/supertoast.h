#ifndef __SUPERTOAST_H
#define __SUPERTOAST_H

// ==========================================================================
// Context data for line search callback function

struct OF_CLBK_DATA {
    CFwdSolver *fws;
    const Raster *raster;
    const Scaler *pscaler;
    Solution *meshsol;
    const Regularisation *reg;
    const CCompRowMatrix *qvec;
    const CCompRowMatrix *mvec;
    double omega;
    const RVector *data;
    const RVector *sd;
};

// ==========================================================================
// The callback function for obtaining the objective function during
// line search

double of_clbk (const RVector &x, double *of_sub, void *context);

#endif // !__SUPERTOAST_H
