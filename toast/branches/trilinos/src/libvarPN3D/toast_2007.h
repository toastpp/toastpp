// -*-C++-*-
// ============================================================================
// M. Schweiger and S. R. Arridge					28.5.96
// general TOAST header file
// ============================================================================

#ifndef __TOAST_H
#define __TOAST_H

#include <fstream.h>
#include <mathlib.h>
#include <felib.h>

// some global constants

#define MAXREG 50	// maximum admissible number of regions
#define NREGS 3		// number of regularisation entries
#define MAXMTP 20	// maximum number of measurement types
#define START_TAU (1e-3)	/* regularisation start parameter */

// enums and typedefs

typedef enum {
    ENORM_CHISQUARED,
    ENORM_UNNORMALISED,
    ENORM_NORMALISED,
    ENORM_CALCULATE,
    ENORM_CONSTANTRECEIVED,
    ENORM_READ
} EnormType;

typedef enum {
    RECON_MUA,
    RECON_KAPPA,
    RECON_BOTH
} ReconType;

typedef enum {
    RESET_MESH,
    RESET_HOMOG,
    RESET_NIM,
    RESET_PROBE
} ResetType;

typedef enum {
    KERNEL_SVD,
    KERNEL_GJ,
    KERNEL_ART,
    KERNEL_BLOCKART,
    KERNEL_AUGMENTART,
    KERNEL_AUGMENTBLOCKART,
    KERNEL_CG_FR,
    KERNEL_CG_PR,
    KERNEL_BLOCK_KACZMARZ,
    KERNEL_TRUNCATED_NEWTON,
    KERNEL_FILTBACKPROP,
    KERNEL_QUASI_NEWTON
} KernelType;

typedef enum {
    GRADIENT_EXPLICIT,
    GRADIENT_PARTSYS,
    GRADIENT_PMDF,
    GRADIENT_DIRECT
} GradientEvalType;

//typedef enum {
//    CG_PRECON_NONE,
//    CG_PRECON_DIAG,
//    CG_PRECON_EXP_R,
//    CG_PRECON_CH
//} PreconType;

typedef enum {
    EACH_BACKPROJECTION,
    EACH_ART_ITERATION,
    EACH_UPDATE
} ImageFilterAction;

typedef enum {
    CONSTRAINT_NONE,
    CONSTRAINT_EACH_BACKPROJECTION,
    CONSTRAINT_EACH_ART_ITERATION,
    CONSTRAINT_EACH_UPDATE
} ConstraintAction;

typedef enum {			/*** SRA 21-4-97 ***/
    BLOCK_SOLVE_MINNORM,
    BLOCK_SOLVE_PARALLEL_ART
} BlockSolveType;

typedef enum {
    BC_DIRICHLET,
    BC_ROBIN
} BCType;

typedef enum {
    ISOTROPIC,
    NEUMANN,
    EXPLICIT
} SourceType;

typedef enum {
    SRC_POINT,
    SRC_GAUSSIAN
} SourceProfType;

typedef enum {
    MEAS_POINT,
    MEAS_GAUSSIAN
} MeasurementProfType;

typedef enum {
    SAS,
    RAS
} SourceOrderType;

typedef enum {
    ALTERNATING,
    SIMULTANEOUS
} DualSolveType;

// ============================================================================
// class forward declarations

class DataMap;
class Environment;
class ForwardSolver;

class DataType;
    class Data;
    class Projection;
    class SD;

class DataElement;
class Dvec;
class Mvec;
class QMvec;
class MDvec;

class Basis;
    class NodeBasis;
    class PixelBasis;
    class CubicPixelBasis2D;
    class FourierBasis;
    class WaveletBasis;
    class BlobBasis2D;
        class BesselBlobBasis;
        class GaussBlobBasis;
        class SplineBlobBasis;
        class HanningBlobBasis;
        class RampBlobBasis;

class Image;
    class PixelImage;
    class NodalImage;
    class ElementImage;

class pmdf;
class Dpmdf;
class Mpmdf;
class QMpmdf;
class MDpmdf;

class Solution;
    class TOAST_Solution;
    class State;
    class Update;

class MType;
    class Intensity;
    class Moment;
    class CentralMoment;
    class LTransform;
    class MellinLaplace;

class Kernel;
    class ExplicitKernel;
        class SVDKernel;
        class GJKernel;
    class ARTKernel;
    class BlockARTKernel;
    class AugmentedARTKernel;

class Solver;
    class JacobianSolver;
    class GradientSolver;
        class CGSolver;
        class TNSolver;
        class QNSolver;

class ObjectiveFunction;

class Regularisation;
    class Regularisation_Cluster;

// local includes

#include "mtype.h"
#include "basis.h"
#include "cpbasis.h"
#include "fbasis.h"
#include "bbasis.h"
#include "gbbasis.h"
#include "sbbasis.h"
#include "bbbasis.h"
#include "hbbasis.h"
#include "rbbasis.h"
#include "wvbasis.h"
#include "image.h"
#include "environ.h"
#include "fwdsolve.h"
#include "solution.h"
#include "TOASTsolution.h"
#include "update.h"
#include "state.h"
#include "tpsf.h"
#include "datamap.h"
#include "datatype.h"
#include "project.h"
#include "sd.h"
#include "data.h"
#include "globlopt.h"
#include "pmdf.h"
#include "kernel.h"
#include "regular.h"
#include "reg_cluster.h"
#include "gradient.h"
#include "solver.h"
#include "jcsolver.h"
#include "gsolver.h"
#include "cgsolver.h"
#include "bksolver.h"
#include "tnsolver.h"
#include "qnsolver.h"
//#include "seq.h"
#include "stack.h"


#endif // !__TOAST_H
