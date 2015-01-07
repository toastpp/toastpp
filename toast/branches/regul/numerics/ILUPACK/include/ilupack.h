/*! \file ilupack.h
   \brief main header for ILUPACK

   This header contains all definitions of functions as well as those of the
   constants
*/
#ifndef _ILU_PACK_H
#define _ILU_PACK_H
#define USE_LAPACK_DRIVER


#include <stdlib.h>

#include "long_integer.h"
#include "namesilupack.h"



#define LONG_INT integer
#define MEDIUM_INT int

/*! switch to indicate inverse-based dropping. It is used in AMGINIT
    AMGGETPARAMS and AMGSETPARAMS. The parameter "flag" is bitwise
    modified by flag|=DROP_INVERSE to set and flag&=~DROP_INVERSE to
    turn off inverse-based dropping.
    In AMGINIT, DROP_INVERSE is set by default
 */
#define DROP_INVERSE                     1

/*! switch for not shifting away zero pivots. This switch is used in ILUC
    which does not have pivoting prevent small diagonal entries.
    The parameter "param" is bitwise modified by param|=NO_SHIFT to 
    suppress shifts and param&=~NO_SHIFT to allow shifts.
 */
#define NO_SHIFT                         2

/*! switch for using Tismenetsky update. 
 */
#define TISMENETSKY_SC                   4
/* switch for repeated ILU */
#define REPEAT_FACT                      8
/* switch for enhanced estimate for the norm of the inverses */
#define IMPROVED_ESTIMATE               16
/* switch for using diagonal compensation */
#define DIAGONAL_COMPENSATION           32
/* switch for reducing the partial factorization to the non-coarse part */
#define COARSE_REDUCE                   64

/* switch for using a different pivoting strategy, if the regular reordering
   fails and before we switch to ILUTP
*/
#define FINAL_PIVOTING                 128
/* enforce the positve definite property */
#define ENSURE_SPD                     256

/* switch for the most simple Schur complement update */
#define SIMPLE_SC                      512


#define PREPROCESS_INITIAL_SYSTEM     1024
#define PREPROCESS_SUBSYSTEMS         2048
#define MULTI_PILUC                   4096


#define RE_FACTOR                     8192
#define AGGRESSIVE_DROPPING          16384
#define DISCARD_MATRIX               32768
#define SYMMETRIC_STRUCTURE          65536

#define STATIC_TESTVECTOR           131072
#define DYNAMIC_TESTVECTOR          262144
#define ADAPT_CONDEST               524288

#define SADDLE_POINT               1048576
#define BLOCK_STRUCTURE            2097152
/*
                                   4194304
                                   8388608
                                  16777216
                                  33554432
                                  67108864
                                 134217728 
                                 268435456 
                                 536870912
                                1073741824
                                2147483648
*/



#define _D_REAL_MAX_        1.7e+308
#define _S_REAL_MAX_        1.7e+38

/* ***************************************************** */
/* ******      Definitions for preconditioners     ***** */
typedef struct {
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   void *a;
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
} SPARSEmat;

typedef struct {
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   doubleprecision *a;
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
} Dmat;

typedef struct {
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   real *a;
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
} Smat;

typedef struct {
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   ilu_doublecomplex *a;
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
} Zmat;

typedef struct {
   integer nr;
   integer nc;
   integer nnz;
   integer *ia;
   integer *ja;
   ilu_complex *a;
   integer issymmetric;
   integer isdefinite;
   integer ishermitian;
   integer isskew;
   integer isreal;
   integer issingle;
} Cmat;



#define ILUPACK_NIPAR   50
#define ILUPACK_NFPAR   50




typedef struct  AMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer nlev;                  
  integer n;                  
  integer nB;
  SPARSEmat A; 
  SPARSEmat LU;
  integer *LUperm;
  SPARSEmat E;
  SPARSEmat F;
  integer *p;
  integer *invq;
  void *rowscal;
  void *colscal;
  void *absdiag;
  struct AMGLM *prev;
  struct AMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  doubleprecision errorL;
  doubleprecision errorU;
  doubleprecision errorS;
  /*
     size_t *pardiso_pt[64];
     int     pardiso_iparm[64];
  */
} AMGlevelmat; 

typedef struct  DAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer nlev;                  
  integer n;                  
  integer nB;
  Dmat A; 
  Dmat LU;
  integer *LUperm;
  Dmat E;
  Dmat F;
  integer *p;
  integer *invq;
  doubleprecision *rowscal;
  doubleprecision *colscal;
  doubleprecision *absdiag;
  struct DAMGLM *prev;
  struct DAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  doubleprecision errorL;
  doubleprecision errorU;
  doubleprecision errorS;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} DAMGlevelmat; 

typedef struct  SAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer nlev;                  
  integer n;                  
  integer nB; 
  Smat A; 
  Smat LU;
  integer *LUperm;
  Smat E;
  Smat F;
  integer *p;
  integer *invq;
  real *rowscal;
  real *colscal;
  real *absdiag;
  struct SAMGLM *prev;
  struct SAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  real errorL;
  real errorU;
  real errorS;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} SAMGlevelmat; 

typedef struct  ZAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer nlev;                  
  integer n;                  
  integer nB; 
  Zmat A; 
  Zmat LU;
  integer *LUperm;
  Zmat E;
  Zmat F;
  integer *p;
  integer *invq;
  ilu_doublecomplex *rowscal;
  ilu_doublecomplex *colscal;
  ilu_doublecomplex *absdiag;
  struct ZAMGLM *prev;
  struct ZAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  doubleprecision errorL;
  doubleprecision errorU;
  doubleprecision errorS;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} ZAMGlevelmat; 

typedef struct CAMGLM {
  integer issymmetric;
  integer isdefinite;
  integer ishermitian;
  integer isskew;
  integer isreal;
  integer issingle;
  integer isblock;
  integer nlev;                  
  integer n;                  
  integer nB; 
  Cmat A; 
  Cmat LU;
  integer *LUperm;
  Cmat E;
  Cmat F;
  integer *p;
  integer *invq;
  ilu_complex *rowscal;
  ilu_complex *colscal;
  ilu_complex *absdiag;
  struct CAMGLM *prev;
  struct CAMGLM *next;
  integer *nextblock;
  integer *blocksize;
  integer maxblocksize;
  real errorL;
  real errorU;
  real errorS;
  /*
    size_t *pardiso_pt[64];
    int     pardiso_iparm[64];
  */
} CAMGlevelmat; 


typedef struct {
   integer       ipar[ILUPACK_NIPAR];
   real          fpar[ILUPACK_NFPAR];
   integer       type;
   integer       *ibuff;
   integer       *iaux;
   real          *dbuff;
   real          *daux;
   integer       *ju;
   integer       *jlu;
   real          *alu;
   real          *testvector;
   size_t        nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer       rcomflag, returnlabel;
   real          *tv;
   integer       *ind;
   integer       nindicator;
   integer       *indicator;
   Smat          A;
   integer       istack[30], *pistack[20];
   real          rstack[30], *prstack[10];
   real          fstack[30], *pfstack[10];
   size_t        ststack[5], *pststack[5];
   Smat          mstack[1];
   SAMGlevelmat  *amglmstack[1];
   integer       (*intfctstack[3])();
   integer       matching;
   char          *ordering;
   real          droptol;
   real          droptolS;
   real          condest;
   real          restol;
   integer       maxit;
   real          elbow;
   integer       lfil;
   integer       lfilS;
   char          *typetv;
   char          *amg;
   integer       npresmoothing;
   integer       npostsmoothing;
   integer       ncoarse;
   char          *presmoother;
   char          *postsmoother;
   char          *FCpart;
   char          *typecoarse;
   integer       nrestart;
   integer       flags;
   char          *solver;
   real          damping;
   integer       mixedprecision;
   integer       (*perm0)();
   integer       (*perm)();
   integer       (*permf)();
   real          shift0;
   real          shiftmax;
   integer       nshifts;
   real          *shifts;
   Smat          *shiftmatrix;
} SILUPACKparam;

typedef struct {
   integer          ipar[ILUPACK_NIPAR];
   doubleprecision  fpar[ILUPACK_NFPAR];
   integer          type;
   integer          *ibuff;
   integer          *iaux;
   doubleprecision  *dbuff;
   doubleprecision  *daux;
   integer          *ju;
   integer          *jlu;
   doubleprecision  *alu;
   doubleprecision  *testvector;
   size_t           nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer          rcomflag, returnlabel;
   doubleprecision  *tv;
   integer          *ind;
   integer          nindicator;
   integer          *indicator;
   Dmat             A;
   integer          istack[30], *pistack[20];
   doubleprecision  rstack[30], *prstack[10];
   doubleprecision  fstack[30], *pfstack[10];
   size_t           ststack[5], *pststack[5];
   Dmat             mstack[1];
   DAMGlevelmat     *amglmstack[1];
   integer          (*intfctstack[3])();
   integer          matching;
   char             *ordering;
   doubleprecision  droptol;
   doubleprecision  droptolS;
   doubleprecision  condest;
   doubleprecision  restol;
   integer          maxit;
   doubleprecision  elbow;
   integer          lfil;
   integer          lfilS;
   char             *typetv;
   char             *amg;
   integer          npresmoothing;
   integer          npostsmoothing;
   integer          ncoarse;
   char             *presmoother;
   char             *postsmoother;
   char             *FCpart;
   char             *typecoarse;
   integer          nrestart;
   integer          flags;
   char             *solver;
   doubleprecision  damping;
   integer          mixedprecision;
   integer          (*perm0)();
   integer          (*perm)();
   integer          (*permf)();
   doubleprecision  shift0;
   doubleprecision  shiftmax;
   integer          nshifts;
   doubleprecision  *shifts;
   Dmat             *shiftmatrix;
} DILUPACKparam;

typedef struct {
   integer     ipar[ILUPACK_NIPAR];
   real    fpar[ILUPACK_NFPAR];
   integer     type;
   integer     *ibuff;
   integer     *iaux;
   ilu_complex *dbuff;
   ilu_complex *daux;
   integer     *ju;
   integer     *jlu;
   ilu_complex *alu;
   ilu_complex *testvector;
   size_t  nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer          rcomflag, returnlabel;
   ilu_complex          *tv;
   integer          *ind;
   integer              nindicator;
   integer          *indicator;
   Cmat             A;
   integer          istack[30], *pistack[20];
   real             rstack[30], *prstack[10];
   ilu_complex          fstack[30], *pfstack[10];
   size_t           ststack[5], *pststack[5];
   Cmat             mstack[1];
   CAMGlevelmat     *amglmstack[1];
   integer          (*intfctstack[3])();
   integer      matching;
   char         *ordering;
   real         droptol;
   real         droptolS;
   real         condest;
   real         restol;
   integer      maxit;
   real         elbow;
   integer      lfil;
   integer      lfilS;
   char         *typetv;
   char         *amg;
   integer      npresmoothing;
   integer      npostsmoothing;
   integer      ncoarse;
   char         *presmoother;
   char         *postsmoother;
   char         *FCpart;
   char         *typecoarse;
   integer       nrestart;
   integer       flags;
   char         *solver;
   ilu_complex       damping;
   integer       mixedprecision;
   integer       (*perm0)();
   integer       (*perm)();
   integer       (*permf)();
   ilu_complex       shift0;
   ilu_complex       shiftmax;
   integer       nshifts;
   ilu_complex       *shifts;
   Cmat          *shiftmatrix;
} CILUPACKparam;

typedef struct {
   integer              ipar[ILUPACK_NIPAR];
   doubleprecision  fpar[ILUPACK_NFPAR];
   integer              type;
   integer              *ibuff;
   integer              *iaux;
   ilu_doublecomplex    *dbuff;
   ilu_doublecomplex    *daux;
   integer              *ju;
   integer              *jlu;
   ilu_doublecomplex    *alu;
   ilu_doublecomplex    *testvector;
   size_t           nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer          rcomflag, returnlabel;
   ilu_doublecomplex    *tv;
   integer          *ind;
   integer              nindicator;
   integer          *indicator;
   Zmat             A;
   integer          istack[30], *pistack[20];
   doubleprecision  rstack[30], *prstack[10];
   ilu_doublecomplex    fstack[30], *pfstack[10];
   size_t           ststack[5], *pststack[5];
   Zmat             mstack[1];
   ZAMGlevelmat     *amglmstack[1];
   integer          (*intfctstack[3])();
   integer      matching;
   char         *ordering;
   doubleprecision       droptol;
   doubleprecision       droptolS;
   doubleprecision       condest;
   doubleprecision       restol;
   integer      maxit;
   doubleprecision       elbow;
   integer      lfil;
   integer      lfilS;
   char         *typetv;
   char         *amg;
   integer      npresmoothing;
   integer      npostsmoothing;
   integer      ncoarse;
   char         *presmoother;
   char         *postsmoother;
   char         *FCpart;
   char         *typecoarse;
   integer       nrestart;
   integer       flags;
   char         *solver;
   ilu_doublecomplex damping;
   integer       mixedprecision;
   integer       (*perm0)();
   integer       (*perm)();
   integer       (*permf)();
   ilu_doublecomplex shift0;
   ilu_doublecomplex shiftmax;
   integer       nshifts;
   ilu_doublecomplex *shifts;
   Zmat          *shiftmatrix;
} ZILUPACKparam;


typedef struct {
   integer              ipar[ILUPACK_NIPAR];
   doubleprecision      fpar[ILUPACK_NFPAR];
   integer              type;
   integer              *ibuff;
   integer              *iaux;
   void                 *dbuff;
   void                 *daux;
   integer              *ju;
   integer              *jlu;
   void                 *alu;
   void                 *testvector;
   size_t               nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
   integer              rcomflag, returnlabel;
   void                 *tv;
   integer              *ind;
   integer              nindicator;
   integer              *indicator;
   SPARSEmat            A;
   integer              istack[30], *pistack[20];
   doubleprecision      rstack[30], *prstack[10];
   ilu_doublecomplex        fstack[30], *pfstack[10];
   size_t               ststack[5], *pststack[5];
   SPARSEmat            mstack[1];
   AMGlevelmat          *amglmstack[1];
   integer              (*intfctstack[3])();
   integer              matching;
   char                 *ordering;
   doubleprecision      droptol;
   doubleprecision      droptolS;
   doubleprecision      condest;
   doubleprecision      restol;
   integer              maxit;
   doubleprecision      elbow;
   integer              lfil;
   integer              lfilS;
   char                 *typetv;
   char                 *amg;
   integer              npresmoothing;
   integer              npostsmoothing;
   integer              ncoarse;
   char                 *presmoother;
   char                 *postsmoother;
   char                 *FCpart;
   char                 *typecoarse;
   integer              nrestart;
   integer              flags;
   char                 *solver;
   ilu_doublecomplex        damping;
   integer              mixedprecision;
   integer              (*perm0)();
   integer              (*perm)();
   integer              (*permf)();
} ILUPACKparam;



void DGNLAMGsol1(DAMGlevelmat *, integer, integer, integer, integer, integer, 
		 DILUPACKparam *, doubleprecision *, doubleprecision *);
void SGNLAMGsol1(SAMGlevelmat *, integer, integer, integer, integer, integer, 
		 SILUPACKparam *, real *, real *);
void ZGNLAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CGNLAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);

void DGNLAMGsol2(DAMGlevelmat *, integer, integer, integer, integer, integer, 
		 DILUPACKparam *, doubleprecision *, doubleprecision *);
void SGNLAMGsol2(SAMGlevelmat *, integer, integer, integer, integer, integer, 
		 SILUPACKparam *, real *, real *);
void ZGNLAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CGNLAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);


void DSPDAMGsol1(DAMGlevelmat *, integer, integer, integer, integer, integer, 
		 DILUPACKparam *, doubleprecision *, doubleprecision *);
void SSPDAMGsol1(SAMGlevelmat *, integer, integer, integer, integer, integer, 
		 SILUPACKparam *, real *, real *);
void ZHPDAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CHPDAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);

void DSPDAMGsol2(DAMGlevelmat *, integer, integer, integer, integer, integer, 
		 DILUPACKparam *, doubleprecision *, doubleprecision *);
void SSPDAMGsol2(SAMGlevelmat *, integer, integer, integer, integer, integer, 
		 SILUPACKparam *, real *, real *);
void ZHPDAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CHPDAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);


void DSYMAMGsol1(DAMGlevelmat *, integer, integer, integer, integer, integer, 
		 DILUPACKparam *, doubleprecision *, doubleprecision *);
void SSYMAMGsol1(SAMGlevelmat *, integer, integer, integer, integer, integer, 
		 SILUPACKparam *, real *, real *);
void ZHERAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CHERAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);
void ZSYMAMGsol1(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CSYMAMGsol1(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);

void DSYMAMGsol2(DAMGlevelmat *, integer, integer, integer, integer, integer, 
		 DILUPACKparam *, doubleprecision *, doubleprecision *);
void SSYMAMGsol2(SAMGlevelmat *, integer, integer, integer, integer, integer, 
		 SILUPACKparam *, real *, real *);
void ZHERAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CHERAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);
void ZSYMAMGsol2(ZAMGlevelmat *, integer, integer, integer, integer, integer, 
		 ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void CSYMAMGsol2(CAMGlevelmat *, integer, integer, integer, integer, integer, 
		 CILUPACKparam *, ilu_complex *, ilu_complex *);





void   AMGsol(AMGlevelmat *,ILUPACKparam *, void *, void *, void *);
void   AMGtsol(AMGlevelmat *,ILUPACKparam *, void *, void *, void *);
void   AMGhsol(AMGlevelmat *,ILUPACKparam *, void *, void *, void *);


void    DGNLAMGsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void    DGNLAMGsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
void  DGNLAMGdlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void  DGNLAMGdusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void   DGNLAMGlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void   DGNLAMGusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLAMGtdlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLAMGtdusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void  DGNLAMGtlsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void  DGNLAMGtusol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void   DGNLAMGtsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void   DGNLAMGtsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
void    DSPDAMGsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void    DSPDAMGsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
void    DSYMAMGsol_internal(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void    DSYMAMGsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *);
void    DSYMAMGbsol(DAMGlevelmat *,DILUPACKparam *, doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLAMGextract(Dmat *,Dmat *, Dmat, integer *,integer *,  integer);
void DSYMAMGextract(Dmat *, Dmat, integer *,integer *,  integer);
void DSSMAMGextract(Dmat *, Dmat, integer *,integer *,  integer);

void    SGNLAMGsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void    SGNLAMGsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
void  SGNLAMGdlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void  SGNLAMGdusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void   SGNLAMGusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void   SGNLAMGlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void SGNLAMGtdlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void SGNLAMGtdusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void  SGNLAMGtusol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void  SGNLAMGtlsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void   SGNLAMGtsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void   SGNLAMGtsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
void    SSPDAMGsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void    SSPDAMGsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
void    SSYMAMGsol_internal(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void    SSYMAMGsol(SAMGlevelmat *,SILUPACKparam *, real *, real *);
void    SSYMAMGbsol(SAMGlevelmat *,SILUPACKparam *, real *, real *, real *);
void SGNLAMGextract(Smat *,Smat *, Smat, integer *,integer *,  integer);
void SSYMAMGextract(Smat *, Smat, integer *,integer *,  integer);
void SSSMAMGextract(Smat *, Smat, integer *,integer *,  integer);

void    ZGNLAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
#ifdef __cplusplus
extern "C" {
#endif
void    ZGNLAMGsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
#ifdef __cplusplus
}
#endif

void  ZGNLAMGdlsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void  ZGNLAMGdusol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void   ZGNLAMGlsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void   ZGNLAMGusol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void ZGNLAMGtdlsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void ZGNLAMGtdusol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void  ZGNLAMGtlsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void  ZGNLAMGtusol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void   ZGNLAMGtsol_internal(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void   ZGNLAMGtsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZHPDAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZHPDAMGsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZHERAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZHERAMGsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZSYMAMGsol_internal(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZSYMAMGsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZHERAMGbsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void    ZSYMAMGbsol(ZAMGlevelmat *,ZILUPACKparam *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *);
void ZGNLAMGextract(Zmat *,Zmat *, Zmat, integer *,integer *,  integer);
void ZHERAMGextract(Zmat *, Zmat, integer *,integer *,  integer);
void ZSHRAMGextract(Zmat *, Zmat, integer *,integer *,  integer);
void ZSYMAMGextract(Zmat *, Zmat, integer *,integer *,  integer);
void ZSSMAMGextract(Zmat *, Zmat, integer *,integer *,  integer);

void    CGNLAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void    CGNLAMGsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *);
void  CGNLAMGdlsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void  CGNLAMGdusol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void   CGNLAMGlsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void   CGNLAMGusol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void CGNLAMGtdlsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void CGNLAMGtdusol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void  CGNLAMGtusol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void  CGNLAMGtlsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void   CGNLAMGtsol_internal(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void   CGNLAMGtsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *);
void    CHPDAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void    CHPDAMGsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *);
void    CHERAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void    CHERAMGsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *);
void    CSYMAMGsol_internal(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void    CSYMAMGsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *);
void    CHERAMGbsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void    CSYMAMGbsol(CAMGlevelmat *,CILUPACKparam *, ilu_complex *, ilu_complex *, ilu_complex *);
void CGNLAMGextract(Cmat *,Cmat *, Cmat, integer *,integer *,  integer);
void CHERAMGextract(Cmat *, Cmat, integer *,integer *,  integer);
void CSHRAMGextract(Cmat *, Cmat, integer *,integer *,  integer);
void CSYMAMGextract(Cmat *, Cmat, integer *,integer *,  integer);
void CSSMAMGextract(Cmat *, Cmat, integer *,integer *,  integer);


void DGNLlupq(integer *, doubleprecision *, integer *, integer *, integer *);
void DGNLlupqsol (integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqtsol(integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqlsol (integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqtlsol(integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqusol (integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqtusol(integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqdlsol (integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqtdlsol(integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqdusol (integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);
void DGNLlupqtdusol(integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, doubleprecision *, doubleprecision *);

void SGNLlupq(integer *, real *, integer *, integer *, integer *);
void SGNLlupqsol (integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqtsol(integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqlsol (integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqtlsol(integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqusol (integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqtusol(integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqdlsol (integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqtdlsol(integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqdusol (integer *, real *, integer *, integer *, real *, real *, 
		  real *);
void SGNLlupqtdusol(integer *, real *, integer *, integer *, real *, real *, 
		  real *);

void ZGNLlupq(integer *, ilu_doublecomplex *, integer *, integer *, integer *);
void ZGNLlupqsol (integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqtsol(integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqlsol (integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqtlsol(integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqusol (integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqtusol(integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqdlsol (integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqtdlsol(integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqdusol (integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);
void ZGNLlupqtdusol(integer *, ilu_doublecomplex *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, 
		  ilu_doublecomplex *);

void CGNLlupq(integer *, ilu_complex *, integer *, integer *, integer *);
void CGNLlupqsol (integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqtsol(integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqlsol (integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqtlsol(integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqusol (integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqtusol(integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqdlsol (integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqtdlsol(integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqdusol (integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);
void CGNLlupqtdusol(integer *, ilu_complex *, integer *, integer *, ilu_complex *, ilu_complex *,
		  ilu_complex *);


void DSPDldlp(integer *, doubleprecision *, integer *, integer *, integer *);
void DSPDldlpsol (integer *, doubleprecision *, integer *,
		  doubleprecision *, doubleprecision *, integer *);

void SSPDldlp(integer *, real *, integer *, integer *, integer *);
void SSPDldlpsol (integer *, real *, integer *,
		  real *, real *, integer *);

void ZHPDldlp(integer *, ilu_doublecomplex *, integer *, integer *, integer *);
void ZHPDldlpsol (integer *, ilu_doublecomplex *, integer *,
		  ilu_doublecomplex *, ilu_doublecomplex *, integer *);

void CHPDldlp(integer *, ilu_complex *, integer *, integer *, integer *);
void CHPDldlpsol (integer *, ilu_complex *, integer *,
		  ilu_complex *, ilu_complex *, integer *);

integer AMGfactor(SPARSEmat *, AMGlevelmat *, ILUPACKparam *);

#ifdef __cplusplus
extern "C" {
#endif
integer DGNLAMGfactor(Dmat *, DAMGlevelmat *, DILUPACKparam *);
#ifdef __cplusplus
}
#endif

integer DSPDAMGfactor(Dmat *, DAMGlevelmat *, DILUPACKparam *);
integer DSYMAMGfactor(Dmat *, DAMGlevelmat *, DILUPACKparam *);

#ifdef __cplusplus
extern "C" {
#endif
integer SGNLAMGfactor(Smat *, SAMGlevelmat *, SILUPACKparam *);
#ifdef __cplusplus
}
#endif

integer SSPDAMGfactor(Smat *, SAMGlevelmat *, SILUPACKparam *);
integer SSYMAMGfactor(Smat *, SAMGlevelmat *, SILUPACKparam *);

#ifdef __cplusplus
extern "C" {
#endif
integer ZGNLAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
#ifdef __cplusplus
}
#endif

integer ZHPDAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
integer ZHERAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
#ifdef __cplusplus
extern "C" {
#endif
integer ZSYMAMGfactor(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
integer CGNLAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);
#ifdef __cplusplus
}
#endif

integer CHPDAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);
integer CHERAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);
integer CSYMAMGfactor(Cmat *, CAMGlevelmat *, CILUPACKparam *);


void       DGNLilut(integer *,doubleprecision *,integer *,integer *,integer *,
		    doubleprecision *,
		    doubleprecision *,integer *,integer *,integer *, 
		    doubleprecision *,integer *, integer *);
void       DGNLilutp(integer *,doubleprecision *,integer *,integer *,integer *,
		     doubleprecision *,doubleprecision *,integer *,
		     doubleprecision *,integer *,integer *,integer *, 
		     doubleprecision *,integer *,integer *,integer *);
void       DGNLlusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLlutsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLludlsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLlutdlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLludusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLlutdusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLlulsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLlutlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLluusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);
void       DGNLlutusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,integer *);

void       SGNLilut(integer *,real *,integer *,integer *,integer *,
		    real *,
		    real *,integer *,integer *,integer *, 
		    real *,integer *, integer *);
void       SGNLilutp(integer *,real *,integer *,integer *,integer *,
		     real *,real *,integer *,
		     real *,integer *,integer *,integer *, 
		     real *,integer *,integer *,integer *);
void       SGNLlusol (integer *,real *,real *,real *,integer *,integer *);
void       SGNLlutsol(integer *,real *,real *,real *,integer *,integer *);
void       SGNLlulsol (integer *,real *,real *,real *,integer *,integer *);
void       SGNLlutlsol(integer *,real *,real *,real *,integer *,integer *);
void       SGNLluusol (integer *,real *,real *,real *,integer *,integer *);
void       SGNLlutusol(integer *,real *,real *,real *,integer *,integer *);
void       SGNLludlsol (integer *,real *,real *,real *,integer *,integer *);
void       SGNLlutdlsol(integer *,real *,real *,real *,integer *,integer *);
void       SGNLludusol (integer *,real *,real *,real *,integer *,integer *);
void       SGNLlutdusol(integer *,real *,real *,real *,integer *,integer *);

void       ZGNLilut (integer *,ilu_doublecomplex *,integer *,integer *,integer *,
		     doubleprecision *,
		     ilu_doublecomplex *,integer *,integer *,integer *, 
		     ilu_doublecomplex *,integer *,integer *);
void       ZGNLilutp(integer *,ilu_doublecomplex *,integer *,integer *,integer *,
		     doubleprecision *,doubleprecision *,integer *,
		     ilu_doublecomplex *,integer *,integer *,integer *, 
		     ilu_doublecomplex *,integer *,integer *,integer *);
void       ZGNLlusol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLlutsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLlulsol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLlutlsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLluusol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLlutusol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLludlsol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLlutdlsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLludusol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);
void       ZGNLlutdusol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *);

void       CGNLilut (integer *,ilu_complex *,integer *,integer *,integer *,
		     real *,
		     ilu_complex *,integer *,integer *,integer *,
		     ilu_complex *,integer *,integer *);
void       CGNLilutp(integer *,ilu_complex *,integer *,integer *,integer *,
		     real *,real *,integer *,
		     ilu_complex *,integer *,integer *,integer *,
		     ilu_complex *,integer *,integer *,integer *);
void       CGNLlusol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLlutsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLlulsol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLlutlsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLluusol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLlutusol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLludlsol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLlutdlsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLludusol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);
void       CGNLlutdusol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,integer *);


void DGNLiluc(integer *,doubleprecision *,integer *,integer *,integer *,
	      doubleprecision *,integer *,doubleprecision *,integer *,integer *,
	      integer *,doubleprecision *,integer *,integer *);
void DGNLilucsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLiluctsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLilucdlsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLiluctdlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLilucdusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLiluctdusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLiluclsol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLiluctlsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLilucusol (integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);
void DGNLiluctusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,integer *,
		  integer *);

void DGNLpiluclsol  (integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpilucdlsol (integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpilucusol  (integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpilucdusol (integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpiluctlsol (integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpiluctdlsol(integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpiluctusol (integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);
void DGNLpiluctdusol(integer *,integer *, doubleprecision *,doubleprecision *,
		     doubleprecision *,integer *,integer *);

void DSYMildlc(integer *,doubleprecision *,integer *,integer *,integer *,
	       doubleprecision *,integer *,doubleprecision *,integer *,integer *,
	       doubleprecision *,integer *,integer *,integer *);
void DSYMildlcsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		  integer *);
void DSYMpildlcdlsol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSYMpildlcdusol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSYMpildlclsol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSYMpildlcusol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSYMpilucsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		 integer *);
void DSYMbpilucsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		   integer *,integer *,doubleprecision *);
void DSSMildlc(integer *,doubleprecision *,integer *,integer *,integer *,
	       doubleprecision *,integer *,doubleprecision *,integer *,integer *,
	       doubleprecision *,integer *,integer *);
void DSSMildlcsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		 integer *);
void DSSMpildlcdlsol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSSMpildlcdusol(integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSSMpildlclsol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DSSMpildlcusol (integer *,integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		     integer *);
void DGNLpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       doubleprecision *,integer *,integer *,integer *,doubleprecision *,integer *,
	       integer *, doubleprecision *, doubleprecision *, integer*, integer*);
void DGNLspiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,integer *,integer *,integer *,integer *,
		doubleprecision *,integer *,integer *,integer *,doubleprecision *,integer *,
		integer *, doubleprecision *, doubleprecision *, integer*, integer*);
void DGNLmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		doubleprecision *,integer *,integer *,integer *,doubleprecision *,integer *,
		integer *, doubleprecision *, doubleprecision *, integer *, integer*);
void DSPDpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       doubleprecision *,integer *,integer *,doubleprecision *,integer *,
	       integer *, doubleprecision *, doubleprecision *, integer *, integer *);
void DSPDmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		doubleprecision *,integer *,integer *,doubleprecision *,integer *,
		integer *, doubleprecision *, doubleprecision *, integer *, integer *);
void DSYMpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       doubleprecision *,integer *,integer *,doubleprecision *,integer *,
	       integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void DSYMbpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       doubleprecision *,integer *,integer *,doubleprecision *,integer *,
	       integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *, integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void DSYMmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       doubleprecision *,integer *,integer *,doubleprecision *,integer *,
	       integer *, doubleprecision *,  doubleprecision *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
/* void DSYMmpiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		doubleprecision *,integer *,integer *,doubleprecision *,integer *,
		integer *, doubleprecision *,  doubleprecision *, integer *, integer *);
*/

void DSYMiluc(integer *,doubleprecision *,integer *,integer *,integer *,doubleprecision *,
	       integer *,integer *,integer *,
	       doubleprecision *,integer *,integer *,doubleprecision *,integer *,
	       integer *,integer *);
void SGNLiluc(integer *,real *,integer *,integer *,integer *,
	      real *,integer *,real *,integer *,integer *,
	      integer *,real *,integer *,integer *);
void SGNLilucsol (integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLiluctsol(integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLilucdlsol (integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLiluctdlsol(integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLilucdusol (integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLiluctdusol(integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLiluclsol (integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLiluctlsol(integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLilucusol (integer *,real *,real *,real *,integer *,
		  integer *);
void SGNLiluctusol(integer *,real *,real *,real *,integer *,
		  integer *);

void SGNLpiluclsol  (integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpilucdlsol (integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpilucusol  (integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpilucdusol (integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpiluctlsol (integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpiluctdlsol(integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpiluctusol (integer *,integer *, real *,real *,
		     real *,integer *,integer *);
void SGNLpiluctdusol(integer *,integer *, real *,real *,
		     real *,integer *,integer *);

void SSYMildlc(integer *,real *,integer *,integer *,integer *,
	       real *,integer *,real *,integer *,integer *,
	       real *,integer *,integer *,integer *);
void SSYMildlcsol(integer *,real *,real *,real *,
		  integer *);
void SSYMpildlcdlsol(integer *,integer *,real *,real *,real *,
		     integer *);
void SSYMpildlcdusol(integer *,integer *,real *,real *,real *,
		     integer *);
void SSYMpildlclsol (integer *,integer *,real *,real *,real *,
		     integer *);
void SSYMpildlcusol (integer *,integer *,real *,real *,real *,
		     integer *);
void SSYMpilucsol(integer *,real *,real *,real *,
		 integer *);
void SSYMbpilucsol(integer *,real *,real *,real *,
		   integer *,integer *,real *);
void SSSMildlc(integer *,real *,integer *,integer *,integer *,
	       real *,integer *,real *,integer *,integer *,
	       real *,integer *,integer *);
void SSSMildlcsol(integer *,real *,real *,real *,
		 integer *);
void SSSMpildlcdlsol(integer *,integer *,real *,real *,real *,
		     integer *);
void SSSMpildlcdusol(integer *,integer *,real *,real *,real *,
		     integer *);
void SSSMpildlclsol (integer *,integer *,real *,real *,real *,
		     integer *);
void SSSMpildlcusol (integer *,integer *,real *,real *,real *,
		     integer *);
void SGNLpiluc(integer *,real *,integer *,integer *,integer *,real *,
	       real *,integer *,integer *,integer *,integer *,
	       real *,integer *,integer *,integer *,real *,integer *,
	       integer *, real *, real *, integer *, integer*);
void SGNLspiluc(integer *,real *,integer *,integer *,integer *,real *,
		real *,integer *,integer *,integer *,integer *,
		real *,integer *,integer *,integer *,real *,integer *,
		integer *, real *, real *, integer *, integer*);
void SGNLmpiluc(integer *,real *,integer *,integer *,integer *,real *,
		real *,real *,integer *,integer *,integer *,integer *,
		real *,integer *,integer *,integer *,real *,integer *,
		integer *, real *, real *, integer *, integer*);
void SSPDpiluc(integer *,real *,integer *,integer *,integer *,real *,
	       real *,integer *,integer *,integer *,integer *,
	       real *,integer *,integer *,real *,integer *,
	       integer *, real *, real *, integer *, integer *);
void SSPDmpiluc(integer *,real *,integer *,integer *,integer *,real *,
		real *,real *,integer *,integer *,integer *,integer *,
		real *,integer *,integer *,real *,integer *,
		integer *, real *, real *, integer *, integer *);
void SSYMpiluc(integer *,real *,integer *,integer *,integer *,real *,
	       real *,integer *,integer *,integer *,integer *,
	       real *,integer *,integer *,real *,integer *,
	       integer *, real *, real *, integer *, integer *, integer *,
	       real *, real *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void SSYMbpiluc(integer *,real *,integer *,integer *,integer *,real *,
	       real *,integer *,integer *,integer *,integer *,
	       real *,integer *,integer *,real *,integer *,
	       integer *, real *, real *, integer *, integer *, integer *,
	       real *, real *, integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
/* void SSYMmpiluc(integer *,real *,integer *,integer *,integer *,real *,
		real *,real *,integer *,integer *,integer *,integer *,
		real *,integer *,integer *,real *,integer *,
		integer *, real *, real *, integer *, integer *);
*/
void SSYMiluc(integer *,real *,integer *,integer *,integer *,real *,
	       integer *,integer *,integer *,
	       real *,integer *,integer *,real *,integer *,
	       integer *, integer *);

void ZGNLiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,
	      doubleprecision *,integer *,ilu_doublecomplex *,integer *,integer *,
	      integer *,ilu_doublecomplex *,integer *,integer *);
void ZGNLilucsol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLiluctsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLilucdlsol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLiluctdlsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLilucdusol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLiluctdusol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLiluclsol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLiluctlsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLilucusol (integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);
void ZGNLiluctusol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,
		  integer *);

void ZGNLpiluclsol  (integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpilucdlsol (integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpilucusol  (integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpilucdusol (integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpiluctlsol (integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpiluctdlsol(integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpiluctusol (integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);
void ZGNLpiluctdusol(integer *,integer *, ilu_doublecomplex *,ilu_doublecomplex *,
		     ilu_doublecomplex *,integer *,integer *);

void ZHERildlc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,
	       doubleprecision *,integer *,ilu_doublecomplex *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,integer *);
void ZHERildlcsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		  integer *);
void ZHERpildlcdlsol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZHERpildlcdusol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZHERpildlclsol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZHERpildlcusol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZHERpilucsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		 integer *);
void ZHERbpilucsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		   integer *,integer *,ilu_doublecomplex *);
void ZSYMpilucsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		 integer *);
void ZSYMbpilucsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		   integer *, integer *,ilu_doublecomplex *);
void ZSYMildlc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,
	       doubleprecision *,integer *,ilu_doublecomplex *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,integer *);
void ZSYMildlcsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		  integer *);
void ZSYMpildlcdlsol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSYMpildlcdusol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSYMpildlclsol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSYMpildlcusol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSHRildlc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,
	       doubleprecision *,integer *,ilu_doublecomplex *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *);
void ZSHRildlcsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		 integer *);
void ZSHRpildlcdlsol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSHRpildlcdusol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSHRpildlclsol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSHRpildlcusol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSSMildlc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,
	       doubleprecision *,integer *,ilu_doublecomplex *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *);
void ZSSMildlcsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		 integer *);
void ZSSMpildlcdlsol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSSMpildlcdusol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSSMpildlclsol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZSSMpildlcusol (integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *,
		     integer *);
void ZGNLpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
void ZGNLspiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,integer *,integer *,integer *,integer *,
		ilu_doublecomplex *,integer *,integer *,integer *,ilu_doublecomplex *,integer *,
		integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
void ZGNLmpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		ilu_doublecomplex *,integer *,integer *,integer *,ilu_doublecomplex *,integer *,
		integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
void ZHPDpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
void ZHPDmpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
		integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
void ZHERpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void ZHERbpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *, integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
/* void ZHERmpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
		integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
*/
void ZHERiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *,integer *);

void ZSYMpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void ZSYMbpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       doubleprecision *,integer *,integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *, integer *,
	       doubleprecision *, doubleprecision *, integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
/*void ZSYMmpiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
		doubleprecision *,doubleprecision *,integer *,integer *,integer *,integer *,
		ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
		integer *, doubleprecision *, ilu_doublecomplex *, integer *, integer *);
*/
void ZSYMiluc(integer *,ilu_doublecomplex *,integer *,integer *,integer *,doubleprecision *,
	       integer *,integer *,integer *,
	       ilu_doublecomplex *,integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *, integer *);

void CGNLiluc(integer *,ilu_complex *,integer *,integer *,integer *,
	      real *,integer *,ilu_complex *,integer *,integer *,
	      integer *,ilu_complex *,integer *,integer *);
void CGNLilucsol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLiluctsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLilucdlsol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLiluctdlsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLilucdusol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLiluctdusol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLilucusol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLiluctusol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLiluclsol (integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);
void CGNLiluctlsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,integer *,
		  integer *);

void CGNLpiluclsol  (integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpilucdlsol (integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpilucusol  (integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpilucdusol (integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpiluctlsol (integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpiluctdlsol(integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpiluctusol (integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);
void CGNLpiluctdusol(integer *,integer *, ilu_complex *,ilu_complex *,
		ilu_complex *,integer *,integer *);

void CHERildlc(integer *,ilu_complex *,integer *,integer *,integer *,
	       real *,integer *,ilu_complex *,integer *,integer *,
	       ilu_complex *,integer *,integer *,integer *);
void CHERildlcsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		  integer *);
void CHERpildlcdlsol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CHERpildlcdusol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CHERpildlclsol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CHERpildlcusol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CHERpilucsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		  integer *);
void CSYMpilucsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		  integer *);
void CHERbpilucsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		   integer *,integer *,ilu_complex *);
void CSYMbpilucsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		   integer *, integer *,ilu_complex *);
void CSYMpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
               real *,integer *,integer *,integer *,integer *,
               ilu_complex *,integer *,integer *,ilu_complex *,integer *,
               integer *, real *, ilu_complex *, integer *, integer *, integer *,
	       real *, real *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void CSYMbpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
               real *,integer *,integer *,integer *,integer *,
               ilu_complex *,integer *,integer *,ilu_complex *,integer *,
               integer *, real *, ilu_complex *, integer *, integer *, integer *,
	       real *, real *, integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void CHERpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
               real *,integer *,integer *,integer *,integer *,
               ilu_complex *,integer *,integer *,ilu_complex *,integer *,
               integer *, real *, ilu_complex *, integer *, integer *, integer *,
	       real *, real *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void CHERbpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
               real *,integer *,integer *,integer *,integer *,
               ilu_complex *,integer *,integer *,ilu_complex *,integer *,
               integer *, real *, ilu_complex *, integer *, integer *, integer *,
	       real *, real *, integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,
	       integer *,integer *,integer *,integer *,integer *,integer *,integer *,integer *);
void CSYMiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
	       integer *,integer *,integer *,
	       ilu_complex *,integer *,integer *,ilu_complex *,integer *,
	       integer *, integer *);
void CHERiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
	       integer *,integer *,integer *,
	       ilu_complex *,integer *,integer *,ilu_complex *,integer *,
	       integer *, integer *);

void CSYMildlc(integer *,ilu_complex *,integer *,integer *,integer *,
	       real *,integer *,ilu_complex *,integer *,integer *,
	       ilu_complex *,integer *,integer *,integer *);
void CSYMildlcsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		 integer *);
void CSYMpildlcdlsol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSYMpildlcdusol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSYMpildlclsol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSYMpildlcusol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSHRildlc(integer *,ilu_complex *,integer *,integer *,integer *,
	       real *,integer *,ilu_complex *,integer *,integer *,
	       ilu_complex *,integer *,integer *);
void CSHRildlcsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		  integer *);
void CSHRpildlcdlsol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSHRpildlcdusol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSHRpildlclsol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSHRpildlcusol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSSMildlc(integer *,ilu_complex *,integer *,integer *,integer *,
	       real *,integer *,ilu_complex *,integer *,integer *,
	       ilu_complex *,integer *,integer *);
void CSSMildlcsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		 integer *);
void CSSMpildlcdlsol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSSMpildlcdusol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSSMpildlclsol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CSSMpildlclsol (integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *,
		     integer *);
void CGNLpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
	       real *,integer *,integer *,integer *,integer *,
	       ilu_complex *,integer *,integer *,integer *,ilu_complex *,integer *,
	       integer *, real *, ilu_complex *, integer *, integer *);
void CGNLspiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
		real *,integer *,integer *,integer *,integer *,
		ilu_complex *,integer *,integer *,integer *,ilu_complex *,integer *,
		integer *, real *, ilu_complex *, integer *, integer *);
void CGNLmpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
		real *,real *,integer *,integer *,integer *,integer *,
		ilu_complex *,integer *,integer *,integer *,ilu_complex *,integer *,
		integer *, real *, ilu_complex *, integer *, integer *);
void CHPDpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
	       real *,integer *,integer *,integer *,integer *,
	       ilu_complex *,integer *,integer *,ilu_complex *,integer *,
	       integer *, real *, ilu_complex *, integer *, integer *);
void CHPDmpiluc(integer *,ilu_complex *,integer *,integer *,integer *,real *,
		real *,real *,integer *,integer *,integer *,integer *,
		ilu_complex *,integer *,integer *,ilu_complex *,integer *,
		integer *, real *, ilu_complex *, integer *, integer *);


/* *********************************************** */
/* ******      Definitions for orderings     ***** */

integer    DGNLperm_null        (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_nd          (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_nd_fc       (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_nd_fcv      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_nd_fc       (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_nd_fcv      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_rcm         (Dmat, doubleprecision *,doubleprecision *,
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_rcm_fc      (Dmat, doubleprecision *,doubleprecision *,
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_rcm_fcv     (Dmat, doubleprecision *,doubleprecision *,
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mmd         (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mmd_fc      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mmd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mmd_fc      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mmd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_amf         (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_amf_fc      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_amf_fcv     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_amf_fc      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_amf_fcv     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_amd         (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_amd_fc      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_amd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_amd_fc      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_amd_fcv     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_metis_e     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_metis_e_fc  (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_metis_e_fcv (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_metis_e_fc  (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_metis_e_fcv (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_metis_n     (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_metis_n_fc  (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_metis_n_fcv (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_metis_n_fc  (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_metis_n_fcv (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_pq          (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_fc          (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_fc          (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_fcv         (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_fcv         (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_p           (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_indset      (Dmat, doubleprecision *,doubleprecision *, 
				 integer *,integer *, integer *, DILUPACKparam *);


integer    DGNLperm_mwm_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_amf        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_amd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwm_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);

integer    DSYMperm_mwm_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_rcm_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_mmd_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amf        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amf_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amd_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_n_sp (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mwm_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);


integer    DGNLperm_matching_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_amf        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_amd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_matching_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);

integer    DSYMperm_matching_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_rcm_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_mmd_sp        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amf        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amf_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amd        (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amd_sp     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_e_sp (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_n_sp (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_matching_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				    integer *,integer *, integer *, DILUPACKparam *);


integer    DSYMperm_mc64_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_rcm_sp        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_mmd_sp        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amf        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amf_sp        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amd        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amd_sp        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_e_sp    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_n_sp    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DSYMperm_mc64_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);



integer    DGNLperm_mc64_null       (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_amf        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_amd        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mc64_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);

integer    DGNLperm_mwa_rcm        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_rcm_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_rcm_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_mmd        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_mmd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_mmd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_amf        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_amf_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_amf_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_amd        (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_amd_fc     (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_amd_fcv    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_metis_e    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_metis_e_fc (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_metis_e_fcv(Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_metis_n    (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_metis_n_fc (Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);
integer    DGNLperm_mwa_metis_n_fcv(Dmat, doubleprecision *,doubleprecision *, 
				     integer *,integer *, integer *, DILUPACKparam *);


integer    SGNLperm_null        (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_nd          (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_nd_fc       (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_nd_fcv      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_nd_fc       (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_nd_fcv      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_rcm         (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_rcm_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_rcm_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_mmd         (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_mmd_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_mmd_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_mmd_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_mmd_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_amf         (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_amf_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_amf_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_amf_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_amf_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_amd         (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_amd_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_amd_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_amd_fc      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_amd_fcv     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_metis_e     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_metis_e_fc  (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_metis_e_fcv (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_metis_e_fc  (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_metis_e_fcv (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_metis_n     (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_metis_n_fc  (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_metis_n_fcv (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_metis_n_fc  (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_metis_n_fcv (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_pq          (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_fc          (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_fc          (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_fcv         (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SSYMperm_fcv         (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_indset      (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);
integer    SGNLperm_p           (Smat, real *,real *, integer *,integer *,
				 integer *, SILUPACKparam *);

integer    SGNLperm_mwm_rcm        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_rcm_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_rcm_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_mmd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_mmd_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_mmd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_mmd_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amf        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amf_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amf_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amf_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amd_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_amd_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_e    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_e_sp (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_e_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_n    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_n_sp (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_n_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_mwm_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);

integer    SSYMperm_mwm_rcm        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_rcm_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_rcm_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_mmd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_mmd_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_mmd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_mmd_fcv    (Smat, real *,real *, integer *,integer *,
			            integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amf        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amf_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amf_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amf_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amd_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_amd_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_e    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_e_sp (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_e_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_n    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_n_sp (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_n_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_mwm_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);


integer    SGNLperm_matching_rcm        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_rcm_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_mmd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_mmd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_mmd_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_amf        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_amf_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_amf_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_amd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_amd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_amd_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_metis_e    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_metis_e_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_metis_n    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_metis_n_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SGNLperm_matching_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);

integer    SSYMperm_matching_rcm        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_rcm_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_rcm_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_mmd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_mmd_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_mmd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_mmd_fcv    (Smat, real *,real *, integer *,integer *,
			            integer *, SILUPACKparam *);
integer    SSYMperm_matching_amf        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amf_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amf_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amf_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amd        (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amd_sp     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amd_fc     (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_amd_fcv    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_e    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_e_sp (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_e_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_n    (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_n_sp (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_n_fc (Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);
integer    SSYMperm_matching_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				    integer *, SILUPACKparam *);


integer    SSYMperm_mc64_rcm        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_rcm_sp        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_rcm_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_mmd        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_mmd_sp        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_mmd_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_mmd_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amf        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amf_sp        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amf_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amf_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amd        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amd_sp        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amd_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_amd_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_e    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_e_sp    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_e_fc (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_n    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_n_sp (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_n_fc (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SSYMperm_mc64_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);



integer    SGNLperm_mc64_null       (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_rcm        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_rcm_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_mmd        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_mmd_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_mmd_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_amf        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_amf_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_amf_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_amd        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_amd_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_amd_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_metis_e    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_metis_e_fc (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_metis_n    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_metis_n_fc (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mc64_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);

integer    SGNLperm_mwa_rcm        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_rcm_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_rcm_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_mmd        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_mmd_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_mmd_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_amf        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_amf_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_amf_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_amd        (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_amd_fc     (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_amd_fcv    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_metis_e    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_metis_e_fc (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_metis_e_fcv(Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_metis_n    (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_metis_n_fc (Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);
integer    SGNLperm_mwa_metis_n_fcv(Smat, real *,real *, integer *,integer *,
				     integer *, SILUPACKparam *);


integer    ZGNLperm_null        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_nd          (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_nd_fc       (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_nd_fcv      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_nd_fc       (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_nd_fcv      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_rcm         (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_rcm_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_rcm_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mmd         (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mmd_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mmd_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mmd_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mmd_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_amf         (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_amf_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_amf_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_amf_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_amf_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_amd         (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_amd_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_amd_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_amd_fc      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_amd_fcv     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_metis_e     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_metis_e_fc  (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_metis_e_fcv (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_metis_e_fc  (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_metis_e_fcv (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_metis_n     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_metis_n_fc  (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_metis_n_fcv (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_metis_n_fc  (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_metis_n_fcv (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_pq          (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_fc          (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_fc          (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_fcv         (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_fcv         (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_indset      (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				 integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_p           (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				 integer *,integer *, integer *, ZILUPACKparam *);

integer    ZGNLperm_mwm_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwm_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);

integer    ZHERperm_mwm_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_rcm_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_mmd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amf_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_e_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_n_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mwm_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);

integer    ZSYMperm_mwm_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_rcm_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_mmd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amf_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_e_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_n_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mwm_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);


integer    ZGNLperm_matching_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_matching_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);

integer    ZHERperm_matching_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_rcm_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_mmd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amf_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_e_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_n_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_matching_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);

integer    ZSYMperm_matching_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_rcm_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_mmd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amf_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_e_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_n_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_matching_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				    integer *,integer *, integer *, ZILUPACKparam *);


integer    ZHERperm_mc64_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_rcm_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_mmd_sp        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amf_sp        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amd_sp        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_e_sp    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_n_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZHERperm_mc64_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);

integer    ZSYMperm_mc64_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_rcm_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_mmd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *, 
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amf_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amd_sp     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_e_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_n_sp (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZSYMperm_mc64_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);


integer    ZGNLperm_mc64_null       (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mc64_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);


integer    ZGNLperm_mwa_rcm        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_rcm_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_rcm_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_mmd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_mmd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_mmd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_amf        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_amf_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_amf_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_amd        (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_amd_fc     (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_amd_fcv    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_metis_e    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_metis_e_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_metis_e_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_metis_n    (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_metis_n_fc (Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);
integer    ZGNLperm_mwa_metis_n_fcv(Zmat, ilu_doublecomplex *,ilu_doublecomplex *,
				     integer *,integer *, integer *, ZILUPACKparam *);


integer    CGNLperm_null        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_nd          (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_nd_fc       (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_nd_fcv      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_nd_fc       (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_nd_fcv      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_rcm         (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_rcm_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_rcm_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_mmd         (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_mmd_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_mmd_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_mmd_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_mmd_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_amf         (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_amf_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_amf_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_amf_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_amf_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_amd         (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_amd_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_amd_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_amd_fc      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_amd_fcv     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_metis_e     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_metis_e_fc  (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_metis_e_fcv (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_metis_e_fc  (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_metis_e_fcv (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_metis_n     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_metis_n_fc  (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_metis_n_fcv (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_metis_n_fc  (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_metis_n_fcv (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_pq          (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_fc          (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_fc          (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_fcv         (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CSYMperm_fcv         (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_indset      (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);
integer    CGNLperm_p           (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				 integer *, CILUPACKparam *);

integer    CGNLperm_mwm_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_mwm_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);

integer    CHERperm_mwm_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_rcm_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_mmd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amf_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_e_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_n_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_mwm_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);

integer    CSYMperm_mwm_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_rcm_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_mmd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amf_sp    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_e_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_n_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_mwm_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);

integer    CGNLperm_matching_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CGNLperm_matching_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);

integer    CHERperm_matching_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_rcm_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_mmd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amf_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_e_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_n_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CHERperm_matching_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);

integer    CSYMperm_matching_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_rcm_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_mmd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amf_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_e_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_n_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);
integer    CSYMperm_matching_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				    integer *, CILUPACKparam *);


integer    CHERperm_mc64_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_rcm_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_mmd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amf_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_e_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_n_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CHERperm_mc64_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);

integer    CSYMperm_mc64_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_rcm_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_mmd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amf_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amd_sp     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_e_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_n_sp (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CSYMperm_mc64_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);


integer    CGNLperm_mc64_null       (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mc64_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);


integer    CGNLperm_mwa_rcm        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_rcm_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_rcm_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_mmd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_mmd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_mmd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_amf        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_amf_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_amf_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_amd        (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_amd_fc     (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_amd_fcv    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_metis_e    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_metis_e_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_metis_e_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_metis_n    (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_metis_n_fc (Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);
integer    CGNLperm_mwa_metis_n_fcv(Cmat, ilu_complex *,ilu_complex *, integer *,integer *,
				     integer *, CILUPACKparam *);




#define DSPDperm_null  DGNLperm_null
#define SSPDperm_null  SGNLperm_null
#define ZHPDperm_null  ZGNLperm_null
#define CHPDperm_null  CGNLperm_null

#define DSPDpermnull  DGNLperm_null
#define SSPDpermnull  SGNLperm_null
#define ZHPDpermnull  ZGNLperm_null
#define CHPDpermnull  CGNLperm_null

#define DSYMperm_null  DGNLperm_null
#define SSYMperm_null  SGNLperm_null
#define ZSYMperm_null  ZGNLperm_null
#define CSYMperm_null  CGNLperm_null
#define ZHERperm_null  ZGNLperm_null
#define CHERperm_null  CGNLperm_null


#define DSPDperm_nd    DGNLperm_nd
#define SSPDperm_nd    SGNLperm_nd
#define ZHPDperm_nd    ZGNLperm_nd
#define CHPDperm_nd    CGNLperm_nd

#define DSPDpermnd    DGNLperm_nd
#define SSPDpermnd    SGNLperm_nd
#define ZHPDpermnd    ZGNLperm_nd
#define CHPDpermnd    CGNLperm_nd

#define DSYMperm_nd    DGNLperm_nd
#define SSYMperm_nd    SGNLperm_nd
#define ZSYMperm_nd    ZGNLperm_nd
#define CSYMperm_nd    CGNLperm_nd
#define ZHERperm_nd    ZGNLperm_nd
#define CHERperm_nd    CGNLperm_nd


#define DSPDperm_amf   DGNLperm_amf
#define SSPDperm_amf   SGNLperm_amf
#define ZHPDperm_amf   ZGNLperm_amf
#define CHPDperm_amf   CGNLperm_amf

#define DSPDpermamf   DGNLperm_amf
#define SSPDpermamf   SGNLperm_amf
#define ZHPDpermamf   ZGNLperm_amf
#define CHPDpermamf   CGNLperm_amf

#define DSYMperm_amf   DGNLperm_amf
#define SSYMperm_amf   SGNLperm_amf
#define ZSYMperm_amf   ZGNLperm_amf
#define CSYMperm_amf   CGNLperm_amf
#define ZHERperm_amf   ZGNLperm_amf
#define CHERperm_amf   CGNLperm_amf


#define DSPDperm_amd   DGNLperm_amd
#define SSPDperm_amd   SGNLperm_amd
#define ZHPDperm_amd   ZGNLperm_amd
#define CHPDperm_amd   CGNLperm_amd

#define DSPDpermamd   DGNLperm_amd
#define SSPDpermamd   SGNLperm_amd
#define ZHPDpermamd   ZGNLperm_amd
#define CHPDpermamd   CGNLperm_amd

#define DSYMperm_amd   DGNLperm_amd
#define SSYMperm_amd   SGNLperm_amd
#define ZSYMperm_amd   ZGNLperm_amd
#define CSYMperm_amd   CGNLperm_amd
#define ZHERperm_amd   ZGNLperm_amd
#define CHERperm_amd   CGNLperm_amd


#define DSPDperm_metis_e   DGNLperm_metis_e
#define SSPDperm_metis_e   SGNLperm_metis_e
#define ZHPDperm_metis_e   ZGNLperm_metis_e
#define CHPDperm_metis_e   CGNLperm_metis_e

#define DSPDpermmetis_e   DGNLperm_metis_e
#define SSPDpermmetis_e   SGNLperm_metis_e
#define ZHPDpermmetis_e   ZGNLperm_metis_e
#define CHPDpermmetis_e   CGNLperm_metis_e

#define DSYMperm_metis_e   DGNLperm_metis_e
#define SSYMperm_metis_e   SGNLperm_metis_e
#define ZSYMperm_metis_e   ZGNLperm_metis_e
#define CSYMperm_metis_e   CGNLperm_metis_e
#define ZHERperm_metis_e   ZGNLperm_metis_e
#define CHERperm_metis_e   CGNLperm_metis_e


#define DSPDperm_metis_n   DGNLperm_metis_n
#define SSPDperm_metis_n   SGNLperm_metis_n
#define ZHPDperm_metis_n   ZGNLperm_metis_n
#define CHPDperm_metis_n   CGNLperm_metis_n

#define DSPDpermmetis_n   DGNLperm_metis_n
#define SSPDpermmetis_n   SGNLperm_metis_n
#define ZHPDpermmetis_n   ZGNLperm_metis_n
#define CHPDpermmetis_n   CGNLperm_metis_n

#define DSYMperm_metis_n   DGNLperm_metis_n
#define SSYMperm_metis_n   SGNLperm_metis_n
#define ZSYMperm_metis_n   ZGNLperm_metis_n
#define CSYMperm_metis_n   CGNLperm_metis_n
#define ZHERperm_metis_n   ZGNLperm_metis_n
#define CHERperm_metis_n   CGNLperm_metis_n


#define DSPDperm_fc    DSYMperm_fc
#define SSPDperm_fc    SSYMperm_fc
#define ZHPDperm_fc    ZSYMperm_fc
#define CHPDperm_fc    CSYMperm_fc

#define DSPDperm_fcv    DSYMperm_fcv
#define SSPDperm_fcv    SSYMperm_fcv
#define ZHPDperm_fcv    ZSYMperm_fcv
#define CHPDperm_fcv    CSYMperm_fcv

#define DSPDpermfc    DSYMperm_fc
#define SSPDpermfc    SSYMperm_fc
#define ZHPDpermfc    ZSYMperm_fc
#define CHPDpermfc    CSYMperm_fc

#define DSSMperm_fc    DSYMperm_fc
#define SSSMperm_fc    SSYMperm_fc
#define ZSSMperm_fc    ZSYMperm_fc
#define CSSMperm_fc    CSYMperm_fc

#define DSSMpermfc    DSYMperm_fc
#define SSSMpermfc    SSYMperm_fc
#define ZSSMpermfc    ZSYMperm_fc
#define CSSMpermfc    CSYMperm_fc

#define ZHERperm_fc    ZSYMperm_fc
#define CHERperm_fc    CSYMperm_fc

#define ZHERpermfc    ZSYMperm_fc
#define CHERpermfc    CSYMperm_fc

#define ZSHRperm_fc    ZSYMperm_fc
#define CSHRperm_fc    CSYMperm_fc

#define ZSHRpermfc    ZSYMperm_fc
#define CSHRpermfc    CSYMperm_fc

#define DSPDperm_amd_fc    DSYMperm_amd_fc
#define SSPDperm_amd_fc    SSYMperm_amd_fc
#define ZHPDperm_amd_fc    ZSYMperm_amd_fc
#define CHPDperm_amd_fc    CSYMperm_amd_fc

#define DSPDperm_amd_fcv    DSYMperm_amd_fcv
#define SSPDperm_amd_fcv    SSYMperm_amd_fcv
#define ZHPDperm_amd_fcv    ZSYMperm_amd_fcv
#define CHPDperm_amd_fcv    CSYMperm_amd_fcv


#define DSPDperm_amf_fc    DSYMperm_amf_fc
#define SSPDperm_amf_fc    SSYMperm_amf_fc
#define ZHPDperm_amf_fc    ZSYMperm_amf_fc
#define CHPDperm_amf_fc    CSYMperm_amf_fc

#define DSPDperm_amf_fcv    DSYMperm_amf_fcv
#define SSPDperm_amf_fcv    SSYMperm_amf_fcv
#define ZHPDperm_amf_fcv    ZSYMperm_amf_fcv
#define CHPDperm_amf_fcv    CSYMperm_amf_fcv


#define DSPDperm_mmd_fc    DSYMperm_mmd_fc
#define SSPDperm_mmd_fc    SSYMperm_mmd_fc
#define ZHPDperm_mmd_fc    ZSYMperm_mmd_fc
#define CHPDperm_mmd_fc    CSYMperm_mmd_fc

#define DSPDperm_mmd_fcv    DSYMperm_mmd_fcv
#define SSPDperm_mmd_fcv    SSYMperm_mmd_fcv
#define ZHPDperm_mmd_fcv    ZSYMperm_mmd_fcv
#define CHPDperm_mmd_fcv    CSYMperm_mmd_fcv


#define DSPDperm_rcm_fc    DSYMperm_rcm_fc
#define SSPDperm_rcm_fc    SSYMperm_rcm_fc
#define ZHPDperm_rcm_fc    ZSYMperm_rcm_fc
#define CHPDperm_rcm_fc    CSYMperm_rcm_fc

#define DSPDperm_rcm_fcv    DSYMperm_rcm_fcv
#define SSPDperm_rcm_fcv    SSYMperm_rcm_fcv
#define ZHPDperm_rcm_fcv    ZSYMperm_rcm_fcv
#define CHPDperm_rcm_fcv    CSYMperm_rcm_fcv


#define DSPDperm_metis_n_fc    DSYMperm_metis_n_fc
#define SSPDperm_metis_n_fc    SSYMperm_metis_n_fc
#define ZHPDperm_metis_n_fc    ZSYMperm_metis_n_fc
#define CHPDperm_metis_n_fc    CSYMperm_metis_n_fc

#define DSPDperm_metis_n_fcv    DSYMperm_metis_n_fcv
#define SSPDperm_metis_n_fcv    SSYMperm_metis_n_fcv
#define ZHPDperm_metis_n_fcv    ZSYMperm_metis_n_fcv
#define CHPDperm_metis_n_fcv    CSYMperm_metis_n_fcv


#define DSPDperm_metis_e_fc    DSYMperm_metis_e_fc
#define SSPDperm_metis_e_fc    SSYMperm_metis_e_fc
#define ZHPDperm_metis_e_fc    ZSYMperm_metis_e_fc
#define CHPDperm_metis_e_fc    CSYMperm_metis_e_fc

#define DSPDperm_metis_e_fcv    DSYMperm_metis_e_fcv
#define SSPDperm_metis_e_fcv    SSYMperm_metis_e_fcv
#define ZHPDperm_metis_e_fcv    ZSYMperm_metis_e_fcv
#define CHPDperm_metis_e_fcv    CSYMperm_metis_e_fcv


#define DSPDperm_nd_fc    DSYMperm_nd_fc
#define SSPDperm_nd_fc    SSYMperm_nd_fc
#define ZHPDperm_nd_fc    ZSYMperm_nd_fc
#define CHPDperm_nd_fc    CSYMperm_nd_fc

#define DSPDperm_nd_fcv    DSYMperm_nd_fcv
#define SSPDperm_nd_fcv    SSYMperm_nd_fcv
#define ZHPDperm_nd_fcv    ZSYMperm_nd_fcv
#define CHPDperm_nd_fcv    CSYMperm_nd_fcv




integer    DSPDperm_rcm   (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
			   integer *, DILUPACKparam *);

integer    SSPDperm_rcm   (Smat, real *, real *, integer *,integer *,
			   integer *, SILUPACKparam *);

integer    ZHPDperm_rcm   (Zmat, ilu_doublecomplex *,   ilu_doublecomplex *,   integer *,integer *, 
			   integer *, ZILUPACKparam *);

integer    CHPDperm_rcm   (Cmat, ilu_complex *,   ilu_complex *,   integer *,integer *,
			   integer *, CILUPACKparam *);


integer    DSYMperm_rcm_fc (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
			   integer *, DILUPACKparam *);

integer    SSYMperm_rcm_fc (Smat, real *, real *, integer *,integer *,
			   integer *, SILUPACKparam *);

integer    ZSYMperm_rcm_fc (Zmat, ilu_doublecomplex *,   ilu_doublecomplex *,   integer *,integer *, 
			   integer *, ZILUPACKparam *);

integer    CSYMperm_rcm_fc (Cmat, ilu_complex *,   ilu_complex *,   integer *,integer *,
			   integer *, CILUPACKparam *);


integer    DSYMperm_rcm_fcv(Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
			   integer *, DILUPACKparam *);

integer    SSYMperm_rcm_fcv(Smat, real *, real *, integer *,integer *,
			   integer *, SILUPACKparam *);

integer    ZSYMperm_rcm_fcv(Zmat, ilu_doublecomplex *,   ilu_doublecomplex *,   integer *,integer *, 
			   integer *, ZILUPACKparam *);

integer    CSYMperm_rcm_fcv(Cmat, ilu_complex *,   ilu_complex *,   integer *,integer *,
			   integer *, CILUPACKparam *);


#define DSPDpermrcm    DSPDperm_rcm   
#define SSPDpermrcm    SSPDperm_rcm   
#define ZHPDpermrcm    ZHPDperm_rcm   
#define CHPDpermrcm    CHPDperm_rcm   

#define DSYMperm_rcm   DSPDperm_rcm
#define SSYMperm_rcm   SSPDperm_rcm
#define CSYMperm_rcm   CHPDperm_rcm
#define ZSYMperm_rcm   ZHPDperm_rcm
#define CHERperm_rcm   CHPDperm_rcm
#define ZHERperm_rcm   ZHPDperm_rcm



#define DSPDperm_mmd   DGNLperm_mmd
#define SSPDperm_mmd   SGNLperm_mmd
#define ZHPDperm_mmd   ZGNLperm_mmd
#define CHPDperm_mmd   CGNLperm_mmd

#define DSPDpermmmd   DGNLperm_mmd
#define SSPDpermmmd   SGNLperm_mmd
#define ZHPDpermmmd   ZGNLperm_mmd
#define CHPDpermmmd   CGNLperm_mmd

#define DSYMperm_mmd   DGNLperm_mmd
#define SSYMperm_mmd   SGNLperm_mmd
#define ZSYMperm_mmd   ZGNLperm_mmd
#define CSYMperm_mmd   CGNLperm_mmd
#define ZHERperm_mmd   ZGNLperm_mmd
#define CHERperm_mmd   CGNLperm_mmd


#define DSPDperm_indset DGNLperm_indset
#define SSPDperm_indset SGNLperm_indset
#define ZHPDperm_indset ZGNLperm_indset
#define CHPDperm_indset CGNLperm_indset

#define DSPDpermindset DGNLperm_indset
#define SSPDpermindset SGNLperm_indset
#define ZHPDpermindset ZGNLperm_indset
#define CHPDpermindset CGNLperm_indset

#define DSYMperm_indset DGNLperm_indset
#define SSYMperm_indset SGNLperm_indset
#define ZSYMperm_indset ZGNLperm_indset
#define CSYMperm_indset CGNLperm_indset
#define ZHERperm_indset ZGNLperm_indset
#define CHERperm_indset CGNLperm_indset


integer    DSPDperm_pp    (Dmat, doubleprecision *, doubleprecision *, integer *,integer *,
			   integer *, DILUPACKparam *);

integer    SSPDperm_pp    (Smat, real *, real *, integer *,integer *,
			   integer *, SILUPACKparam *);

integer    ZHPDperm_pp    (Zmat, ilu_doublecomplex *,   ilu_doublecomplex *, integer *,integer *,
			   integer *, ZILUPACKparam *);

integer    CHPDperm_pp    (Cmat, ilu_complex *,   ilu_complex *, integer *,integer *,
			   integer *, CILUPACKparam *);

#define DSPDpermpp     DSPDperm_pp    
#define SSPDpermpp     SSPDperm_pp    
#define ZHPDpermpp     ZHPDperm_pp    
#define CHPDpermpp     CHPDperm_pp    




#define DGNLpermnull      DGNLperm_null   
#define DGNLpermnd        DGNLperm_nd     
#define DGNLpermrcm       DGNLperm_rcm    
#define DGNLpermamf       DGNLperm_amf    
#define DGNLpermamd       DGNLperm_amd 
#define DGNLpermmmd       DGNLperm_mmd    
#define DGNLpermpq        DGNLperm_pq     
#define DGNLpermfc        DGNLperm_fc     
#define DSYMpermfc        DSYMperm_fc     
#define DGNLpermp         DGNLperm_p      
#define DGNLpermindset    DGNLperm_indset         


#define DGNLpermmwm_rcm          DGNLperm_mwm_rcm       
#define DGNLpermmwm_mmd          DGNLperm_mwm_mmd       
#define DGNLpermmwm_amf          DGNLperm_mwm_amf       
#define DGNLpermmwm_amd          DGNLperm_mwm_amd
#define DGNLpermmwm_metis_e      DGNLperm_mwm_metis_e    
#define DGNLpermmwm_metis_n      DGNLperm_mwm_metis_n    

#define DGNLpermmatching_rcm          DGNLperm_matching_rcm       
#define DGNLpermmatching_mmd          DGNLperm_matching_mmd       
#define DGNLpermmatching_amf          DGNLperm_matching_amf       
#define DGNLpermmatching_amd          DGNLperm_matching_amd
#define DGNLpermmatching_metis_e      DGNLperm_matching_metis_e    
#define DGNLpermmatching_metis_n      DGNLperm_matching_metis_n    
                                                       
#define DGNLpermmc64_rcm         DGNLperm_mc64_rcm      
#define DGNLpermmc64_mmd         DGNLperm_mc64_mmd      
#define DGNLpermmc64_amf         DGNLperm_mc64_amf      
#define DGNLpermmc64_amd         DGNLperm_mc64_amd
#define DGNLpermmc64_metis_e     DGNLperm_mc64_metis_e  
#define DGNLpermmc64_metis_n     DGNLperm_mc64_metis_n  


#define SGNLpermnull    SGNLperm_null    
#define SGNLpermnd      SGNLperm_nd      
#define SGNLpermrcm     SGNLperm_rcm     
#define SGNLpermamf     SGNLperm_amf     
#define SGNLpermamd     SGNLperm_amd
#define SGNLpermmmd     SGNLperm_mmd     
#define SGNLpermpq      SGNLperm_pq      
#define SGNLpermfc      SGNLperm_fc      
#define SSYMpermfc      SSYMperm_fc      
#define SGNLpermp       SGNLperm_p       
#define SGNLpermindset  SGNLperm_indset  

#define SSPDpermrcm	    SSPDperm_rcm	
#define SSPDpermpp	    SSPDperm_pp	

#define SGNLpermmwm_rcm       SGNLperm_mwm_rcm    
#define SGNLpermmwm_mmd       SGNLperm_mwm_mmd    
#define SGNLpermmwm_amf       SGNLperm_mwm_amf    
#define SGNLpermmwm_amd       SGNLperm_mwm_amd
#define SGNLpermmwm_metis_e   SGNLperm_mwm_metis_e
#define SGNLpermmwm_metis_n   SGNLperm_mwm_metis_n

#define SGNLpermmatching_rcm       SGNLperm_matching_rcm    
#define SGNLpermmatching_mmd       SGNLperm_matching_mmd    
#define SGNLpermmatching_amf       SGNLperm_matching_amf    
#define SGNLpermmatching_amd       SGNLperm_matching_amd
#define SGNLpermmatching_metis_e   SGNLperm_matching_metis_e
#define SGNLpermmatching_metis_n   SGNLperm_matching_metis_n


#define SGNLpermmc64_rcm          SGNLperm_mc64_rcm    
#define SGNLpermmc64_mmd	  SGNLperm_mc64_mmd    
#define SGNLpermmc64_amf	  SGNLperm_mc64_amf    
#define SGNLpermmc64_amd	  SGNLperm_mc64_amd
#define SGNLpermmc64_metis_e	  SGNLperm_mc64_metis_e
#define SGNLpermmc64_metis_n	  SGNLperm_mc64_metis_n


#define CGNLpermnull       CGNLperm_null  
#define CGNLpermnd	   CGNLperm_nd    
#define CGNLpermrcm	   CGNLperm_rcm   
#define CGNLpermamf	   CGNLperm_amf   
#define CGNLpermamd	   CGNLperm_amd
#define CGNLpermmmd	   CGNLperm_mmd   
#define CGNLpermpq	   CGNLperm_pq    
#define CGNLpermfc	   CGNLperm_fc    
#define CSYMpermfc	   CSYMperm_fc    
#define CGNLpermp	   CGNLperm_p     
#define CGNLpermindset	   CGNLperm_indset

#define CGNLpermmwm_rcm      CGNLperm_mwm_rcm	   
#define CGNLpermmwm_mmd	     CGNLperm_mwm_mmd	   
#define CGNLpermmwm_amf	     CGNLperm_mwm_amf	   
#define CGNLpermmwm_amd	     CGNLperm_mwm_amd
#define CGNLpermmwm_metis_e  CGNLperm_mwm_metis_e
#define CGNLpermmwm_metis_n  CGNLperm_mwm_metis_n

#define CGNLpermmatching_rcm      CGNLperm_matching_rcm	   
#define CGNLpermmatching_mmd	     CGNLperm_matching_mmd	   
#define CGNLpermmatching_amf	     CGNLperm_matching_amf	   
#define CGNLpermmatching_amd	     CGNLperm_matching_amd
#define CGNLpermmatching_metis_e  CGNLperm_matching_metis_e
#define CGNLpermmatching_metis_n  CGNLperm_matching_metis_n


#define CGNLpermmc64_rcm      CGNLperm_mc64_rcm    
#define CGNLpermmc64_mmd      CGNLperm_mc64_mmd    
#define CGNLpermmc64_amf      CGNLperm_mc64_amf    
#define CGNLpermmc64_amd      CGNLperm_mc64_amd 
#define CGNLpermmc64_metis_e  CGNLperm_mc64_metis_e
#define CGNLpermmc64_metis_n  CGNLperm_mc64_metis_n

#define CHPDpermrcm	 CHPDperm_rcm   
#define CHPDpermpp	 CHPDperm_pp    


#define ZGNLpermnull      ZGNLperm_null	 
#define ZGNLpermnd	  ZGNLperm_nd	 
#define ZGNLpermrcm	  ZGNLperm_rcm	 
#define ZGNLpermamf	  ZGNLperm_amf	 
#define ZGNLpermamd	  ZGNLperm_amd
#define ZGNLpermmmd	  ZGNLperm_mmd	 
#define ZGNLpermpq	  ZGNLperm_pq	 
#define ZGNLpermfc	  ZGNLperm_fc	 
#define ZSYMpermfc	  ZSYMperm_fc	 
#define ZGNLpermp	  ZGNLperm_p	 
#define ZGNLpermindset    ZGNLperm_indset   

#define ZGNLpermmwm_rcm         ZGNLperm_mwm_rcm	       
#define ZGNLpermmwm_mmd		ZGNLperm_mwm_mmd	       
#define ZGNLpermmwm_amf		ZGNLperm_mwm_amf	       
#define ZGNLpermmwm_amd		ZGNLperm_mwm_amd
#define ZGNLpermmwm_metis_e	ZGNLperm_mwm_metis_e    
#define ZGNLpermmwm_metis_n	ZGNLperm_mwm_metis_n    
				                       
#define ZGNLpermmatching_rcm         ZGNLperm_matching_rcm	       
#define ZGNLpermmatching_mmd		ZGNLperm_matching_mmd	       
#define ZGNLpermmatching_amf		ZGNLperm_matching_amf	       
#define ZGNLpermmatching_amd		ZGNLperm_matching_amd
#define ZGNLpermmatching_metis_e	ZGNLperm_matching_metis_e    
#define ZGNLpermmatching_metis_n	ZGNLperm_matching_metis_n    
				                       
#define ZGNLpermmc64_rcm	ZGNLperm_mc64_rcm       
#define ZGNLpermmc64_mmd	ZGNLperm_mc64_mmd       
#define ZGNLpermmc64_amf	ZGNLperm_mc64_amf       
#define ZGNLpermmc64_amd	ZGNLperm_mc64_amd
#define ZGNLpermmc64_metis_e	ZGNLperm_mc64_metis_e   
#define ZGNLpermmc64_metis_n	ZGNLperm_mc64_metis_n   
				                       
#define ZHPDpermrcm		ZHPDperm_rcm	       
#define ZHPDpermpp		ZHPDperm_pp	       


#define ZHERindfc  ZSYMindfc               
#define CHERindfc  CSYMindfc               



void swapj(integer *, integer, integer);
void dswapm(doubleprecision *, integer, integer);
void sswapm(real *, integer, integer);
integer DPQpermF(Dmat, integer, integer *, integer *, integer *, doubleprecision,doubleprecision *,    integer *);
integer ZPQpermF(Zmat, integer, integer *, integer *, integer *, doubleprecision,ilu_doublecomplex *, integer *);
integer indAMF(Dmat, integer, integer *, integer *, doubleprecision);

integer Dindset(Dmat, integer, integer *, integer *, doubleprecision);
integer Sindset(Smat, integer, integer *, integer *, real);
integer Zindset(Zmat, integer, integer *, integer *, doubleprecision);
integer Cindset(Cmat, integer, integer *, integer *, real);

void Dindfc(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void Sindfc(Smat, integer *, integer *, integer *, real, real *, integer *);
void Zindfc(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void Cindfc(Cmat, integer *, integer *, integer *, real, real *, integer *);

void Dindfc_rs(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void Sindfc_rs(Smat, integer *, integer *, integer *, real, real *, integer *);
void Zindfc_rs(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void Cindfc_rs(Cmat, integer *, integer *, integer *, real, real *, integer *);

void DSYMindfc(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void SSYMindfc(Smat, integer *, integer *, integer *, real, real *, integer *);
void ZSYMindfc(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void CSYMindfc(Cmat, integer *, integer *, integer *, real, real *, integer *);

void DSYMindfc_rs(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void SSYMindfc_rs(Smat, integer *, integer *, integer *, real, real *, integer *);
void ZSYMindfc_rs(Zmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
void CSYMindfc_rs(Cmat, integer *, integer *, integer *, real, real *, integer *);


void SSYMpindfc(Smat, integer *,integer *, integer *,integer *,
		integer *, real, real *,integer *);
void DSYMpindfc(Dmat, integer *,integer *, integer *,integer *,
		integer *, doubleprecision, doubleprecision *,integer *);
void CSYMpindfc(Cmat, integer *,integer *, integer *,integer *,
		integer *, real, real *,integer *);
void ZSYMpindfc(Zmat, integer *,integer *, integer *,integer *,
		integer *, doubleprecision, doubleprecision *,integer *);
void CHERpindfc(Cmat, integer *,integer *, integer *,integer *,
		integer *, real, real *,integer *);
void ZHERpindfc(Zmat, integer *,integer *, integer *,integer *,
		integer *, doubleprecision, doubleprecision *,integer *);

void SSYMpindfc_rs(Smat, integer *,integer *, integer *,integer *,
		integer *, real, real *,integer *);
void DSYMpindfc_rs(Dmat, integer *,integer *, integer *,integer *,
		integer *, doubleprecision, doubleprecision *,integer *);
void CSYMpindfc_rs(Cmat, integer *,integer *, integer *,integer *,
		integer *, real, real *,integer *);
void ZSYMpindfc_rs(Zmat, integer *,integer *, integer *,integer *,
		integer *, doubleprecision, doubleprecision *,integer *);
void CHERpindfc_rs(Cmat, integer *,integer *, integer *,integer *,
		integer *, real, real *,integer *);
void ZHERpindfc_rs(Zmat, integer *,integer *, integer *,integer *,
		integer *, doubleprecision, doubleprecision *,integer *);


void SSYMpindfcv(Smat, integer *,integer *, integer *,integer *,
		 integer *, real, real *, integer, 
		 real *,integer *,real *);
void DSYMpindfcv(Dmat, integer *,integer *, integer *,integer *,
		 integer *, doubleprecision, doubleprecision *, integer, 
		 doubleprecision *,integer *,doubleprecision *);
void CSYMpindfcv(Cmat, integer *,integer *, integer *,integer *,
		 integer *, real, ilu_complex *, integer,
		 real *,integer *,ilu_complex *);
void ZSYMpindfcv(Zmat, integer *,integer *, integer *,integer *,
		 integer *, doubleprecision, ilu_doublecomplex *, integer, 
		 doubleprecision *,integer *,ilu_doublecomplex *);
void CHERpindfcv(Cmat, integer *,integer *, integer *,integer *,
		 integer *, real, ilu_complex *, integer,
		 real *,integer *,ilu_complex *);
void ZHERpindfcv(Zmat, integer *,integer *, integer *,integer *,
		 integer *, doubleprecision, ilu_doublecomplex *, integer, 
		 doubleprecision *,integer *,ilu_doublecomplex *);

void SSYMpindfcv_rs(Smat, integer *,integer *, integer *,integer *,
		 integer *, real, real *, integer, 
		 real *,integer *,real *);
void DSYMpindfcv_rs(Dmat, integer *,integer *, integer *,integer *,
		 integer *, doubleprecision, doubleprecision *, integer, 
		 doubleprecision *,integer *,doubleprecision *);
void CSYMpindfcv_rs(Cmat, integer *,integer *, integer *,integer *,
		 integer *, real, ilu_complex *, integer,
		 real *,integer *,ilu_complex *);
void ZSYMpindfcv_rs(Zmat, integer *,integer *, integer *,integer *,
		 integer *, doubleprecision, ilu_doublecomplex *, integer, 
		 doubleprecision *,integer *,ilu_doublecomplex *);
void CHERpindfcv_rs(Cmat, integer *,integer *, integer *,integer *,
		 integer *, real, ilu_complex *, integer,
		 real *,integer *,ilu_complex *);
void ZHERpindfcv_rs(Zmat, integer *,integer *, integer *,integer *,
		 integer *, doubleprecision, ilu_doublecomplex *, integer, 
		 doubleprecision *,integer *,ilu_doublecomplex *);


void SSYMbuildblock(Smat, integer *,integer *, integer *, real *);
void DSYMbuildblock(Dmat, integer *,integer *, integer *, doubleprecision *);
void CSYMbuildblock(Cmat, integer *,integer *, integer *, real *);
void ZSYMbuildblock(Zmat, integer *,integer *, integer *, doubleprecision *);
void CHERbuildblock(Cmat, integer *,integer *, integer *, real *);
void ZHERbuildblock(Zmat, integer *,integer *, integer *, doubleprecision *);







void Dindfcv(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer, doubleprecision *, integer *);
void Sindfcv(Smat, integer *, integer *, integer *, real,            real *,            integer, real *,            integer *);
void Zindfcv(Zmat, integer *, integer *, integer *, doubleprecision, ilu_doublecomplex *,   integer, doubleprecision *, integer *);
void Cindfcv(Cmat, integer *, integer *, integer *, real,            ilu_complex *,         integer, real *,            integer *);

void DSYMindfcv(Dmat, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer, doubleprecision *, integer *);
void SSYMindfcv(Smat, integer *, integer *, integer *, real,            real *,            integer, real *,            integer *);
void ZSYMindfcv(Zmat, integer *, integer *, integer *, doubleprecision, ilu_doublecomplex *,   integer, doubleprecision *, integer *);
void CSYMindfcv(Cmat, integer *, integer *, integer *, real,            ilu_complex *,         integer, real *,            integer *);
void ZHERindfcv(Zmat, integer *, integer *, integer *, doubleprecision, ilu_doublecomplex *,   integer, doubleprecision *, integer *);
void CHERindfcv(Cmat, integer *, integer *, integer *, real,            ilu_complex *,         integer, real *,            integer *);


void dqsortr2i(doubleprecision *, integer *, integer *, integer, integer);
void sqsortr2i(real *, integer *, integer *, integer, integer);

void Dclear(integer, doubleprecision *, integer);
void Sclear(integer, real *,            integer);
void Zclear(integer, ilu_doublecomplex *,   integer);
void Cclear(integer, ilu_complex *,         integer);

void IP_etree(integer *, integer *, integer, integer *, 
	      integer *, integer *, integer *);
void IP_post_order(integer *, integer, integer *, integer *);
integer IP_tdfs(integer, integer, integer *, integer *, 
		integer *, integer *);


/* ********************************************* */
/* ******      Definitions for solvers     ***** */
void      Dpcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
void      Dfpcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
void      Dbcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
void      DSYMbcg(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
void      DSYMqmr(integer *,doubleprecision *,   doubleprecision *,   integer *,doubleprecision *,doubleprecision *);
void      Dgmres(integer *,doubleprecision *,doubleprecision *,integer *,doubleprecision *,doubleprecision *);
void      Dfgmres(integer *,doubleprecision *,doubleprecision *,integer *,doubleprecision *,doubleprecision *);

void      Spcg(integer *,real *,   real *,   integer *,real *,real *);
void      Sfpcg(integer *,real *,   real *,   integer *,real *,real *);
void      Sbcg(integer *,real *,   real *,   integer *,real *,real *);
void      SSYMbcg(integer *,real *,   real *,   integer *,real *,real *);
void      SSYMqmr(integer *,real *,   real *,   integer *,real *,real *);
void      Sgmres(integer *,real *,real *,integer *,real *,real *);
void      Sfgmres(integer *,real *,real *,integer *,real *,real *);

void      Zpcg(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      Zfpcg(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      Zbcg(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      ZSYMbcg(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      ZHERbcg(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      ZSYMqmr(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      ZHERqmr(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      Zgmres(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);
void      Zfgmres(integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,doubleprecision *,ilu_doublecomplex *);

void      Cpcg(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      Cfpcg(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      Cbcg(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      CSYMbcg(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      CHERbcg(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      CSYMqmr(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      CHERqmr(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      Cgmres(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);
void      Cfgmres(integer *,ilu_complex *,ilu_complex *,integer *,real *,ilu_complex *);


doubleprecision Ddistdot(integer *,doubleprecision *,integer *,doubleprecision *,integer *);

real            Sdistdot(integer *,real *,integer *,real *,integer *);

ilu_doublecomplex   Zdistdotc(integer *,ilu_doublecomplex *,integer *,ilu_doublecomplex *,integer *);
ilu_doublecomplex   Zdistdotu(integer *,ilu_doublecomplex *,integer *,ilu_doublecomplex *,integer *);

ilu_complex         Cdistdotc(integer *,ilu_complex *,integer *,ilu_complex *,integer *);
ilu_complex         Cdistdotu(integer *,ilu_complex *,integer *,ilu_complex *,integer *);

integer AMGsolver(SPARSEmat *, AMGlevelmat *, ILUPACKparam *, 
		  void *, void *);
#ifdef __cplusplus
extern "C" {
#endif
integer DGNLAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
		      doubleprecision *, doubleprecision *);
#ifdef __cplusplus
}
#endif

integer DGNLSYMAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
		      doubleprecision *, doubleprecision *);
integer DGNLSPDAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
		      doubleprecision *, doubleprecision *);
integer DSPDAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
		      doubleprecision *, doubleprecision *);
integer DSYMAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
		      doubleprecision *, doubleprecision *);
integer DSYMSPDAMGsolver(Dmat *, DAMGlevelmat *, DILUPACKparam *, 
			 doubleprecision *, doubleprecision *);


integer SGNLAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
		      real *, real *);
integer SGNLSPDAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
			 real *, real *);
integer SGNLSYMAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
			 real *, real *);
integer SSPDAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
		      real *, real *);
integer SSYMAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
		      real *, real *);
integer SSYMSPDAMGsolver(Smat *, SAMGlevelmat *, SILUPACKparam *, 
			 real *, real *);
#ifdef __cplusplus
extern "C" {
#endif
integer ZGNLAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
		      ilu_doublecomplex *, ilu_doublecomplex *);
#ifdef __cplusplus
}
#endif

integer ZGNLDGNLAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZGNLHPDAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
			 ilu_doublecomplex *, ilu_doublecomplex *);
integer ZGNLHERAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
			 ilu_doublecomplex *, ilu_doublecomplex *);
integer ZGNLSYMAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
			 ilu_doublecomplex *, ilu_doublecomplex *);
integer ZGNLDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZGNLDSYMAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZHPDAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
		      ilu_doublecomplex *, ilu_doublecomplex *);
integer ZHPDDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZHERHPDAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
			 ilu_doublecomplex *, ilu_doublecomplex *);
integer ZHERDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZHERDSYMAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZHERAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
		      ilu_doublecomplex *, ilu_doublecomplex *);
#ifdef __cplusplus
extern "C" {
#endif
integer ZSYMAMGsolver(Zmat *, ZAMGlevelmat *, ZILUPACKparam *, 
		      ilu_doublecomplex *, ilu_doublecomplex *);
#ifdef __cplusplus
}
#endif
integer ZSYMDSPDAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);
integer ZSYMDSYMAMGsolver(Zmat *, DAMGlevelmat *, DILUPACKparam *, 
			  ilu_doublecomplex *, ilu_doublecomplex *);

integer CGNLAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
		      ilu_complex *, ilu_complex *);
integer CGNLSGNLAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
			  ilu_complex *, ilu_complex *);
integer CGNLSYMAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
			 ilu_complex *, ilu_complex *);
integer CGNLHERAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
			 ilu_complex *, ilu_complex *);
integer CGNLHPDAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
			 ilu_complex *, ilu_complex *);
integer CGNLDSPDAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
			  ilu_complex *, ilu_complex *);
integer CGNLDSYMAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
			  ilu_complex *, ilu_complex *);
integer CHPDAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
		      ilu_complex *, ilu_complex *);
integer CHERAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
		      ilu_complex *, ilu_complex *);
integer CSYMAMGsolver(Cmat *, CAMGlevelmat *, CILUPACKparam *, 
		      ilu_complex *, ilu_complex *);
integer CSYMSSPDAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
			  ilu_complex *, ilu_complex *);
integer CSYMSSYMAMGsolver(Cmat *, SAMGlevelmat *, SILUPACKparam *, 
			  ilu_complex *, ilu_complex *);

void AMGinit(SPARSEmat *, ILUPACKparam *);

#ifdef __cplusplus
extern "C" {
#endif
void DGNLAMGinit(Dmat *, DILUPACKparam *);
void DGNLAMGgetparams(DILUPACKparam *,
		      integer *, integer *, doubleprecision *, doubleprecision *,
		      doubleprecision *, integer *, integer *);
void DGNLAMGsetparams(Dmat *, DILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer , integer );
void SGNLAMGinit(Smat *, SILUPACKparam *);
#ifdef __cplusplus
}
#endif

void SGNLAMGgetparams(SILUPACKparam *,
		     integer *, integer *, real *, real *,
		     real *, integer *, integer *);
void SGNLAMGsetparams(Smat *, SILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer , integer );

#ifdef __cplusplus
extern "C" {
#endif
void ZGNLAMGinit(Zmat *, ZILUPACKparam *);

void ZGNLAMGgetparams(ZILUPACKparam *,
		     integer *, integer *, doubleprecision *, doubleprecision *,
		     doubleprecision *, integer *, integer *);
void ZGNLAMGsetparams(Zmat *, ZILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer , integer );

void CGNLAMGinit(Cmat *, CILUPACKparam *);
#ifdef __cplusplus
}
#endif


void CGNLAMGgetparams(CILUPACKparam *,
		     integer *, integer *, real *, real *,
		     real *, integer *, integer *);
void CGNLAMGsetparams(Cmat *, CILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer , integer );


void DSPDAMGinit(Dmat *, DILUPACKparam *);
void DSPDAMGgetparams(DILUPACKparam *,
		     integer *, integer *, doubleprecision *, doubleprecision *,
		     doubleprecision *, integer *);
void DSPDAMGsetparams(Dmat *, DILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer );

void DSYMAMGinit(Dmat *, DILUPACKparam *);
void DSYMAMGgetparams(DILUPACKparam *,
		     integer *, integer *, doubleprecision *, doubleprecision *,
		     doubleprecision *, integer *);
void DSYMAMGsetparams(Dmat *, DILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer );

void SSPDAMGinit(Smat *, SILUPACKparam *);
void SSPDAMGgetparams(SILUPACKparam *,
		     integer *, integer *, real *, real *,
		     real *, integer *);
void SSPDAMGsetparams(Smat *, SILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer );

void SSYMAMGinit(Smat *, SILUPACKparam *);
void SSYMAMGgetparams(SILUPACKparam *,
		     integer *, integer *, real *, real *,
		     real *, integer *);
void SSYMAMGsetparams(Smat *, SILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer );

void ZHPDAMGinit(Zmat *, ZILUPACKparam *);
void ZHERAMGinit(Zmat *, ZILUPACKparam *);

#ifdef __cplusplus
extern "C" {
#endif
void ZSYMAMGinit(Zmat *, ZILUPACKparam *);
#ifdef __cplusplus
}
#endif
void ZHPDAMGgetparams(ZILUPACKparam *,
		      integer *, integer *, doubleprecision *, doubleprecision *,
		      doubleprecision *, integer *);
void ZHPDAMGsetparams(Zmat *, ZILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer );
void ZHERAMGgetparams(ZILUPACKparam *,
		      integer *, integer *, doubleprecision *, doubleprecision *,
		      doubleprecision *, integer *);
void ZHERAMGsetparams(Zmat *, ZILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer );
#ifdef __cplusplus
extern "C" {
#endif
void ZSYMAMGgetparams(ZILUPACKparam *,
		      integer *, integer *, doubleprecision *, doubleprecision *,
		      doubleprecision *, integer *);
void ZSYMAMGsetparams(Zmat *, ZILUPACKparam *,
		      integer , integer , doubleprecision *, doubleprecision ,
		      doubleprecision , integer );
#ifdef __cplusplus
}
#endif

void CHPDAMGinit(Cmat *, CILUPACKparam *);
void CHERAMGinit(Cmat *, CILUPACKparam *);
void CSYMAMGinit(Cmat *, CILUPACKparam *);
void CHPDAMGgetparams(CILUPACKparam *,
		      integer *, integer *, real *, real *,
		      real *, integer *);
void CHPDAMGsetparams(Cmat *, CILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer );
void CHERAMGgetparams(CILUPACKparam *,
		      integer *, integer *, real *, real *,
		      real *, integer *);
void CHERAMGsetparams(Cmat *, CILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer );
void CSYMAMGgetparams(CILUPACKparam *,
		      integer *, integer *, real *, real *,
		      real *, integer *);
void CSYMAMGsetparams(Cmat *, CILUPACKparam *,
		      integer , integer , real *, real ,
		      real , integer );



/* ********************************************* */
/* ******      Definitions for sparskit    ***** */
void      Dcsrcsc(integer *, integer *, integer *, doubleprecision *, integer *, integer *, 
		  doubleprecision *, integer *, integer *);
void      Dcsrcsc(integer *,integer *,integer *,doubleprecision *,
		  integer *,integer *,doubleprecision *,integer *,integer *);

void      Scsrcsc(integer *, integer *, integer *, real *, integer *, integer *, 
		  real *, integer *, integer *);
void      Scsrcsc(integer *,integer *,integer *,real *,
		  integer *,integer *,real *,integer *,integer *);

void      Zcsrcsc(integer *, integer *, integer *, ilu_doublecomplex *, integer *, integer *, 
		  ilu_doublecomplex *, integer *, integer *);
void      Zcsrcsc (integer *,integer *,integer *,ilu_doublecomplex *,
		   integer *,integer *,ilu_doublecomplex *,integer *,integer *);

void      Ccsrcsc(integer *, integer *, integer *, ilu_complex *, integer *, integer *,
		  ilu_complex *, integer *, integer *);
void      Ccsrcsc(integer *,integer *,integer *,ilu_complex *,
		  integer *,integer *,ilu_complex *,integer *,integer *);




/* ******************************************* */
/* ******      Definitions for tools     ***** */
void DGNLmatvec(Dmat, doubleprecision *, doubleprecision *);
void DGNLmattvec(Dmat, doubleprecision *, doubleprecision *);
void DGNLmathvec(Dmat, doubleprecision *, doubleprecision *);

void SGNLmatvec(Smat, real *, real *);
void SGNLmattvec(Smat, real *, real *);
void SGNLmathvec(Smat, real *, real *);

void ZGNLmatvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);
void ZGNLmattvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);
void ZGNLmathvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);

void CGNLmatvec(Cmat, ilu_complex *, ilu_complex *);
void CGNLmattvec(Cmat, ilu_complex *, ilu_complex *);
void CGNLmathvec(Cmat, ilu_complex *, ilu_complex *);


void DSYMmatvec(Dmat, doubleprecision *, doubleprecision *);
void DSSMmatvec(Dmat, doubleprecision *, doubleprecision *);

void SSYMmatvec(Smat, real *, real *);
void SSSMmatvec(Smat, real *, real *);

void ZHERmatvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);
void ZSHRmatvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);
void ZSYMmatvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);
void ZSSMmatvec(Zmat, ilu_doublecomplex *, ilu_doublecomplex *);

void CHERmatvec(Cmat, ilu_complex *, ilu_complex *);
void CSHRmatvec(Cmat, ilu_complex *, ilu_complex *);
void CSYMmatvec(Cmat, ilu_complex *, ilu_complex *);
void CSSMmatvec(Cmat, ilu_complex *, ilu_complex *);


void Sqsort(real *,            integer *, integer *, integer *); 
void Dqsort(doubleprecision *, integer *, integer *, integer *); 
void Cqsort(ilu_complex *,         integer *, integer *, integer *);
void Zqsort(ilu_doublecomplex *,   integer *, integer *, integer *); 

void Sbqsort(real *,           integer *,integer *,integer *,integer *); 
void Dbqsort(doubleprecision *,integer *,integer *,integer *,integer *); 
void Cbqsort(ilu_complex *,        integer *,integer *,integer *,integer *);
void Zbqsort(ilu_doublecomplex *,  integer *,integer *,integer *,integer *); 


void qqsorti(integer *, integer *, integer *); 

integer  Dspartran(Dmat, Dmat *, integer, integer);
void Dsetupgraph(Dmat, Dmat *, integer *, integer *, size_t);
void Dsetupgraph_epsilon(Dmat, Dmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t);
void Dsetupgraph_epsilon_sp(Dmat, Dmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t, integer *);
void Dqqsort(doubleprecision *, integer *, integer *, integer *, integer *); 
void Dqqsort2(doubleprecision *, integer *, integer *, integer *, integer *); 
void Dqqsorts(doubleprecision *, integer *, integer *, integer *); 
void Dqqsorts2(doubleprecision *, integer *, integer *, integer *); 
void Dcperm(Dmat *, integer *);
void Drperm(Dmat *, integer *);

integer  Sspartran(Smat, Smat *, integer, integer);
void Ssetupgraph(Smat, Smat *, integer *, integer *, size_t);
void Ssetupgraph_epsilon(Smat, Smat *, real, real *, integer *, integer *, size_t);
void Ssetupgraph_epsilon_sp(Smat, Smat *, real, real *, integer *, integer *, size_t, integer *);
void Sqqsort(real *, integer *, integer *, integer *, integer *); 
void Scperm(Smat *, integer *);
void Srperm(Smat *, integer *);

integer  Zspartran(Zmat, Zmat *, integer, integer);
void Zsetupgraph(Zmat, Zmat *, integer *, integer *, size_t);
void Zsetupgraph_epsilon(Zmat, Zmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t);
void Zsetupgraph_epsilon_sp(Zmat, Zmat *, doubleprecision, doubleprecision *, integer *, integer *, size_t, integer *);
void Zqqsort(ilu_doublecomplex *, integer *, integer *, integer *, integer *); 
void Zcperm(Zmat *, integer *);
void Zrperm(Zmat *, integer *);

integer  Cspartran(Cmat, Cmat *, integer, integer);
void Csetupgraph(Cmat, Cmat *, integer *, integer *, size_t);
void Csetupgraph_epsilon(Cmat, Cmat *, real, real *, integer *, integer *, size_t);
void Csetupgraph_epsilon_sp(Cmat, Cmat *, real, real *, integer *, integer *, size_t, integer *);
void Cqqsort(ilu_complex *, integer *, integer *, integer *, integer *);
void Ccperm(Cmat *, integer *);
void Crperm(Cmat *, integer *);


void *MAlloc(size_t, char *);
void *ReAlloc(void *, size_t, char *);
void FRee(void *);
doubleprecision dgeteps();
real            sgeteps();
float evaluate_time(float *, float *);


void Droscal(integer *,integer *,integer *,doubleprecision *,integer *,
	    integer *,doubleprecision *,doubleprecision *,integer *,integer *,
	    integer *);
void Dcoscal(integer *,integer *,integer *,doubleprecision *,integer *,
	    integer *,doubleprecision *,doubleprecision *,integer *,integer *,
	    integer *);
void Drowscale(integer *,integer *,doubleprecision *,integer *,
	       integer *,doubleprecision *,integer *);
void Dcolscale(integer *,integer *,doubleprecision *,integer *,
	       integer *,doubleprecision *,integer *);
void DSPDscale(integer *, doubleprecision *,integer *,
	       integer *,doubleprecision *,integer *);
void DSYMscale(integer *, doubleprecision *,integer *,
	       integer *,doubleprecision *,doubleprecision *,integer *);

void Sroscal(integer *,integer *,integer *,real *,integer *,
	    integer *,real *,real *,integer *,integer *,
	    integer *);
void Scoscal(integer *,integer *,integer *,real *,integer *,
	    integer *,real *,real *,integer *,integer *,
	    integer *);
void Srowscale(integer *,integer *,real *,integer *,
	       integer *,real *,integer *);
void Scolscale(integer *,integer *,real *,integer *,
	       integer *,real *,integer *);
void SSPDscale(integer *, real *,integer *,
	       integer *,real *,integer *);
void SSYMscale(integer *, real *,integer *,
	       integer *,real *,real *,integer *);

void Zroscal(integer *,integer *,integer *,ilu_doublecomplex *,integer *,
	    integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *,
	    integer *);
void Zcoscal(integer *,integer *,integer *,ilu_doublecomplex *,integer *,
	    integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *,integer *,
	    integer *);
void Zrowscale(integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *,ilu_doublecomplex *,integer *);
void Zcolscale(integer *,integer *,ilu_doublecomplex *,integer *,
	       integer *,ilu_doublecomplex *,integer *);
void ZHPDscale(integer *, ilu_doublecomplex *,integer *,
	       integer *,ilu_doublecomplex *,integer *);
void ZSYMscale(integer *, ilu_doublecomplex *,integer *,
	       integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *);
void ZHERscale(integer *, ilu_doublecomplex *,integer *,
	       integer *,ilu_doublecomplex *,ilu_doublecomplex *,integer *);

void Croscal(integer *,integer *,integer *,ilu_complex *,integer *,
	    integer *,ilu_complex *,ilu_complex *,integer *,integer *,
	    integer *);
void Ccoscal(integer *,integer *,integer *,ilu_complex *,integer *,
	    integer *,ilu_complex *,ilu_complex *,integer *,integer *,
	    integer *);
void Crowscale(integer *,integer *,ilu_complex *,integer *,
	       integer *,ilu_complex *,integer *);
void Ccolscale(integer *,integer *,ilu_complex *,integer *,
	       integer *,ilu_complex *,integer *);
void CHPDscale(integer *, ilu_complex *,integer *,
	       integer *,ilu_complex *,integer *);
void CSYMscale(integer *, ilu_complex *,integer *,
	       integer *, ilu_complex *,ilu_complex *,integer *);
void CHERscale(integer *, ilu_complex *,integer *,
	       integer *,ilu_complex *,ilu_complex *,integer *);


integer DPQpermF(Dmat, integer, integer *, integer *, integer *, doubleprecision, doubleprecision *, integer *);
integer SPQpermF(Smat, integer, integer *, integer *, integer *, real,            real *,            integer *);
integer ZPQpermF(Zmat, integer, integer *, integer *, integer *, doubleprecision, ilu_doublecomplex *,   integer *);
integer CPQpermF(Cmat, integer, integer *, integer *, integer *, real,            ilu_complex *,         integer *);

integer DPPpermF(Dmat, integer, integer *, integer *, doubleprecision, doubleprecision *, integer *);
integer SPPpermF(Smat, integer, integer *, integer *, real,            real *,            integer *);
integer ZPPpermF(Zmat, integer, integer *, integer *, doubleprecision, ilu_doublecomplex *,   integer *);
integer CPPpermF(Cmat, integer, integer *, integer *, real,            ilu_complex *,         integer *);





void      Dreadmtc(integer *,integer *,integer *,character *,doubleprecision *,integer *,
		   integer *,doubleprecision *,integer *,character *,integer *,integer *,
		   integer *,character *,character *,character *,
		   integer *, integer *, integer *, doubleprecision *,
		   integer *,ftnlen,ftnlen, ftnlen,ftnlen,ftnlen);
void      Sreadmtc(integer *,integer *,integer *,character *,real *,integer *,
		   integer *,real *,integer *,character *,integer *,integer *,
		   integer *,character *,character *,character *,
		   integer *, integer *, integer *, real *,
		   integer *,ftnlen,ftnlen, ftnlen,ftnlen,ftnlen);
void      Zreadmtc(integer *,integer *,integer *,character *,ilu_doublecomplex *,integer *,
		   integer *,ilu_doublecomplex *,integer *,character *,integer *,integer *,
		   integer *,character *,character *,character *,
		   integer *, integer *, integer *, ilu_doublecomplex *,
		   integer *,ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
void      Creadmtc(integer *,integer *,integer *,character *,ilu_complex *,integer *,
		   integer *,ilu_complex *,integer *,character *,integer *,integer *,
		   integer *,character *,character *,character *,
		   integer *, integer *, integer *, ilu_complex *,
		   integer *,ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
void      Dwritemtc(character *, doubleprecision *,integer *,integer *,
		    doubleprecision *,integer *, character *,
		    integer *,integer *,integer *,
		    character *,character *,character *,
		    ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
void      Swritemtc(character *, real *,integer *,integer *,
		    real *,integer *, character *,
		    integer *,integer *,integer *,
		    character *,character *,character *,
		    ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
void      Zwritemtc(character *, ilu_doublecomplex *,integer *,integer *,
		    ilu_doublecomplex *,integer *, character *,
		    integer *,integer *,integer *,
		    character *,character *,character *,
		    ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);
void      Cwritemtc(character *, ilu_complex *,integer *,integer *,
		    ilu_complex *,integer *, character *,
		    integer *,integer *,integer *,
		    character *,character *,character *,
		    ftnlen,ftnlen,ftnlen,ftnlen,ftnlen);

void      Dreadvectors(character *,doubleprecision *,integer *,integer *,
		       character *,character *, ftnlen,ftnlen,ftnlen);
void      Sreadvectors(character *,real *,integer *,integer *,
		       character *,character *, ftnlen,ftnlen,ftnlen);
void      Zreadvectors(character *,ilu_doublecomplex *,integer *,integer *,
		       character *,character *, ftnlen,ftnlen,ftnlen);
void      Creadvectors(character *,ilu_complex *,integer *,integer *,
		       character *,character *, ftnlen,ftnlen,ftnlen);

void      Dwritevectors(character *,doubleprecision *,integer *,integer *,
			character *,character *, ftnlen,ftnlen,ftnlen);
void      Swritevectors(character *,real *,integer *,integer *,
			character *,character *, ftnlen,ftnlen,ftnlen);
void      Zwritevectors(character *,ilu_doublecomplex *,integer *,integer *,
			character *,character *, ftnlen,ftnlen,ftnlen);
void      Cwritevectors(character *,ilu_complex *,integer *,integer *,
			character *,character *, ftnlen,ftnlen,ftnlen);


integer  DSSMsmwm(Dmat, integer *, integer *, doubleprecision *),
     SSSMsmwm(Smat, integer *, integer *, real *),
     DSYMsmwm(Dmat, integer *, integer *, doubleprecision *),
     SSYMsmwm(Smat, integer *, integer *, real *),
     CSSMsmwm(Cmat, integer *, integer *, ilu_complex *),
     ZSSMsmwm(Zmat, integer *, integer *, ilu_doublecomplex *),
     CSYMsmwm(Cmat, integer *, integer *, ilu_complex *),
     ZSYMsmwm(Zmat, integer *, integer *, ilu_doublecomplex *),
     CSHRsmwm(Cmat, integer *, integer *, ilu_complex *),
     ZSHRsmwm(Zmat, integer *, integer *, ilu_doublecomplex *),
     CHERsmwm(Cmat, integer *, integer *, ilu_complex *),
     ZHERsmwm(Zmat, integer *, integer *, ilu_doublecomplex *),
     DGNLsmwm(Dmat, integer *, integer *, doubleprecision *),
     SGNLsmwm(Smat, integer *, integer *, real *),
     CGNLsmwm(Cmat, integer *, integer *, ilu_complex *),
     ZGNLsmwm(Zmat, integer *, integer *, ilu_doublecomplex *);




void SSYMAMGdelete(Smat *, SAMGlevelmat *, SILUPACKparam *);
void DSYMAMGdelete(Dmat *, DAMGlevelmat *, DILUPACKparam *);
void CSYMAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
#ifdef __cplusplus
extern "C" {
#endif
void ZSYMAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
#ifdef __cplusplus
}
#endif

void CHERAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
void ZHERAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);


void SSPDAMGdelete(Smat *, SAMGlevelmat *, SILUPACKparam *);
void DSPDAMGdelete(Dmat *, DAMGlevelmat *, DILUPACKparam *);
void CHPDAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
void ZHPDAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
		        
#ifdef __cplusplus
extern "C" {
#endif
void SGNLAMGdelete(Smat *, SAMGlevelmat *, SILUPACKparam *);
void DGNLAMGdelete(Dmat *, DAMGlevelmat *, DILUPACKparam *);
void CGNLAMGdelete(Cmat *, CAMGlevelmat *, CILUPACKparam *);
void ZGNLAMGdelete(Zmat *, ZAMGlevelmat *, ZILUPACKparam *);
#ifdef __cplusplus
}
#endif


integer dsymilupack   (integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
		   doubleprecision *, doubleprecision *, integer *, integer *, 
		   doubleprecision *, integer *);
integer dsymilupackfac(long *, long *, 
		   integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);
integer dsymilupacksol(long *, long *, 
		   integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);
integer dsymilupackdel(long *, long *, 
		   integer *, integer *, integer *, doubleprecision *, doubleprecision *, doubleprecision *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);


integer ssymilupack   (integer *, integer *, integer *, real *, real *, real *,
		   real *, real *, integer *, integer *, 
		   real *, integer *);
integer ssymilupackfac(long *, long *, 
		   integer *, integer *, integer *, real *, real *, real *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer ssymilupacksol(long *, long *, 
		   integer *, integer *, integer *, real *, real *, real *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer ssymilupackdel(long *, long *, 
		   integer *, integer *, integer *, real *, real *, real *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);


integer csymilupack   (integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, 
		   real *, integer *);
integer csymilupackfac(long *, long *, 
		   integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer csymilupacksol(long *, long *, 
		   integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer csymilupackdel(long *, long *, 
		   integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer chermilupack   (integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, 
		   real *, integer *);
integer cherilupackfac(long *, long *, 
		   integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer cherilupacksol(long *, long *, 
		   integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);
integer cherilupackdel(long *, long *, 
		   integer *, integer *, integer *, ilu_complex *, ilu_complex *, ilu_complex *,
		   real *, real *, integer *, integer *, integer *, 
		   real *, integer *);


integer zsymilupack   (integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, 
		   doubleprecision *, integer *);
integer zsymilupackfac(long *, long *, 
		   integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);
integer zsymilupacksol(long *, long *, 
		   integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);
integer zsymilupackdel(long *, long *, 
		   integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);

integer zherilupack   (integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, 
		   doubleprecision *, integer *);
integer zherilupackfac(long *, long *, 
		   integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);
integer zherilupacksol(long *, long *, 
		   integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);
integer zherilupackdel(long *, long *, 
		   integer *, integer *, integer *, ilu_doublecomplex *, ilu_doublecomplex *, ilu_doublecomplex *,
		   doubleprecision *, doubleprecision *, integer *, integer *, integer *, 
		   doubleprecision *, integer *);




void sgnlamginit(integer *, integer *, integer *,
		 real *, integer *, character *, real *,real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void dgnlamginit(integer *, integer *, integer *,
		 double *, integer *, character *, double *,double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);
void cgnlamginit(integer *, integer *, integer *,
		 ilu_complex *, integer *, character *, real *, real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void zgnlamginit(integer *, integer *, integer *,
		 ilu_doublecomplex *, integer *, character *, double *, double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);

void sspdamginit(integer *, integer *, integer *,
		 real *, integer *, character *, real *, real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void dspdamginit(integer *, integer *, integer *,
		 double *, integer *, character *, double *, double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);
void chpdamginit(integer *, integer *, integer *,
		 ilu_complex *, integer *, character *, real *, real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void zhpdamginit(integer *, integer *, integer *,
		 ilu_doublecomplex *, integer *, character *, double *, double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);

void ssymamginit(integer *, integer *, integer *,
		 real *, integer *, character *, real *, real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void dsymamginit(integer *, integer *, integer *,
		 double *, integer *, character *, double *, double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);
void cheramginit(integer *, integer *, integer *,
		 ilu_complex *, integer *, character *, real *, real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void zheramginit(integer *, integer *, integer *,
		 ilu_doublecomplex *, integer *, character *, double *, double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);
void csymamginit(integer *, integer *, integer *,
		 ilu_complex *, integer *, character *, real *, real *,
		 real *, real *, integer *, real *,
		 integer *, integer *, integer *, integer *, integer *);
void zsymamginit(integer *, integer *, integer *,
		 ilu_doublecomplex *, integer *, character *, double *, double *,
		 double *, double *, integer *, double *,
		 integer *, integer *, integer *, integer *, integer *);


int sgnlamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  real *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int dgnlamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  double *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int cgnlamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zgnlamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);

int sspdamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  real *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int dspdamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  double *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int chpdamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zhpdamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);

int ssymamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  real *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int dsymamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  double *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int cheramgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zheramgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int csymamgfactor(size_t *, size_t *, 
		  integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zsymamgfactor(size_t *, size_t *,
		  integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);


int sgnlamgsolver(size_t *, size_t *, 
		  real *, real *, integer *, integer *, integer *,
		  real *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int dgnlamgsolver(size_t *, size_t *, 
		  double *, double *, integer *, integer *, integer *,
		  double *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int cgnlamgsolver(size_t *, size_t *, 
		  ilu_complex *, ilu_complex *, integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);

int zgnlamgsolver(size_t *, size_t *, 
		  ilu_doublecomplex *, ilu_doublecomplex *, integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);

int sspdamgsolver(size_t *, size_t *, 
		  real *, real *, integer *, integer *, integer *,
		  real *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int dspdamgsolver(size_t *, size_t *, 
		  double *, double *, integer *, integer *, integer *,
		  double *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int chpdamgsolver(size_t *, size_t *, 
		  ilu_complex *, ilu_complex *, integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zhpdamgsolver(size_t *, size_t *, 
		  ilu_doublecomplex *, ilu_doublecomplex *, integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);

int ssymamgsolver(size_t *, size_t *, 
		  real *, real *, integer *, integer *, integer *,
		  real *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int dsymamgsolver(size_t *, size_t *, 
		  double *, double *, integer *, integer *, integer *,
		  double *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int cheramgsolver(size_t *, size_t *, 
		  ilu_complex *, ilu_complex *, integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zheramgsolver(size_t *, size_t *, 
		  ilu_doublecomplex *, ilu_doublecomplex *, integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);
int csymamgsolver(size_t *, size_t *, 
		  ilu_complex *, ilu_complex *, integer *, integer *, integer *,
		  ilu_complex *, integer *, character *, real *, real *,
		  real *, real *, integer *, real *,
		  integer *, integer *, integer *, integer *, integer *);
int zsymamgsolver(size_t *, size_t *, 
		  ilu_doublecomplex *, ilu_doublecomplex *, integer *, integer *, integer *,
		  ilu_doublecomplex *, integer *, character *, double *, double *,
		  double *, double *, integer *, double *,
		  integer *, integer *, integer *, integer *, integer *);


void sgnlamgsol(size_t *, size_t *, 
		real *, real *, integer *);
void dgnlamgsol(size_t *, size_t *, 
		double *, double *, integer *);
void cgnlamgsol(size_t *, size_t *, 
		ilu_complex *, ilu_complex *, integer *);
void zgnlamgsol(size_t *, size_t *, 
		ilu_doublecomplex *, ilu_doublecomplex *, integer *);

void sspdamgsol(size_t *, size_t *, 
		real *, real *, integer *);
void dspdamgsol(size_t *, size_t *, 
		double *, double *, integer *);
void chpdamgsol(size_t *, size_t *, 
		ilu_complex *, ilu_complex *, integer *);
void zhpdamgsol(size_t *, size_t *, 
		ilu_doublecomplex *, ilu_doublecomplex *, integer *);

void ssymamgsol(size_t *, size_t *, 
		real *, real *, integer *);
void dsymamgsol(size_t *, size_t *, 
		double *, double *, integer *);
void csymamgsol(size_t *, size_t *, 
		ilu_complex *, ilu_complex *, integer *);
void zsymamgsol(size_t *, size_t *, 
		ilu_doublecomplex *, ilu_doublecomplex *, integer *);
void cheramgsol(size_t *, size_t *, 
		ilu_complex *, ilu_complex *, integer *);
void zheramgsol(size_t *, size_t *, 
		ilu_doublecomplex *, ilu_doublecomplex *, integer *);
void ssymamgsol(size_t *, size_t *, 
		real *, real *, integer *);
void dsymamgsol(size_t *, size_t *, 
		double *, double *, integer *);
void csymamgsol(size_t *, size_t *, 
		ilu_complex *, ilu_complex *, integer *);
void zsymamgsol(size_t *, size_t *, 
		ilu_doublecomplex *, ilu_doublecomplex *, integer *);
void cheramgsol(size_t *, size_t *, 
		ilu_complex *, ilu_complex *, integer *);
void zheramgsol(size_t *, size_t *, 
		ilu_doublecomplex *, ilu_doublecomplex *, integer *);


void sgnlamgdelete(size_t *, size_t *);
void dgnlamgdelete(size_t *, size_t *);
void cgnlamgdelete(size_t *, size_t *);
void zgnlamgdelete(size_t *, size_t *);

void sspdamgdelete(size_t *, size_t *);
void dspdamgdelete(size_t *, size_t *);
void chpdamgdelete(size_t *, size_t *);
void zhpdamgdelete(size_t *, size_t *);

void ssymamgdelete(size_t *, size_t *);
void dsymamgdelete(size_t *, size_t *);
void cheramgdelete(size_t *, size_t *);
void zheramgdelete(size_t *, size_t *);
void csymamgdelete(size_t *, size_t *);
void zsymamgdelete(size_t *, size_t *);



void sspdamginfo(size_t *, size_t *, integer *, integer *, 
	         integer *, real *);
void dspdamginfo(size_t *, size_t *, 
		 integer *, integer *, integer *, double *);
void chpdamginfo(size_t *, size_t *, 
		 integer *, integer *, integer *, ilu_complex *);
void zhpdamginfo(size_t *, size_t *, 
		 integer *, integer *, integer *, ilu_doublecomplex *);

void ssymamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, real *);
void dsymamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, double *);
void cheramginfo(size_t *, size_t *,
		 integer *, integer *, integer *, ilu_complex *);
void zheramginfo(size_t *, size_t *,
		 integer *, integer *, integer *, ilu_doublecomplex *);
void csymamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, ilu_complex *);
void zsymamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, ilu_doublecomplex *);

void sgnlamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, real *);
void dgnlamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, double *);
void cgnlamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, ilu_complex *);
void zgnlamginfo(size_t *, size_t *,
		 integer *, integer *, integer *, ilu_doublecomplex *);


size_t sspdamgnnz(size_t *, size_t *);
size_t dspdamgnnz(size_t *, size_t *);
size_t chpdamgnnz(size_t *, size_t *);
size_t zhpdamgnnz(size_t *, size_t *);

size_t ssymamgnnz(size_t *, size_t *);
size_t dsymamgnnz(size_t *, size_t *);
size_t cheramgnnz(size_t *, size_t *);
size_t zheramgnnz(size_t *, size_t *);
size_t csymamgnnz(size_t *, size_t *);
size_t zsymamgnnz(size_t *, size_t *);

size_t sgnlamgnnz(size_t *, size_t *);
size_t dgnlamgnnz(size_t *, size_t *);
size_t cgnlamgnnz(size_t *, size_t *);
size_t zgnlamgnnz(size_t *, size_t *);

void ssymspdamgconvert(size_t *, size_t *);
void dsymspdamgconvert(size_t *, size_t *);
void cherhpdamgconvert(size_t *, size_t *);
void zherhpdamgconvert(size_t *, size_t *);

void samgundoscaling(size_t *, size_t *, integer *, integer *, integer *, real *);
void damgundoscaling(size_t *, size_t *, integer *, integer *, integer *, doubleprecision *);
void camgundoscaling(size_t *, size_t *, integer *, integer *, integer *, ilu_complex *);
void zamgundoscaling(size_t *, size_t *, integer *, integer *, integer *, ilu_doublecomplex *);


/*
void DGNLSYM(Dmat , Dmat *, integer *);
void SGNLSYM(Smat , Smat *, integer *);
void CGNLSYM(Cmat , Cmat *, integer *);
void ZGNLSYM(Zmat , Zmat *, integer *);
void CGNLHER(Cmat , Cmat *, integer *);
void ZGNLHER(Zmat , Zmat *, integer *);
*/
void SSYMSPDAMGconvert(SAMGlevelmat *);
void DSYMSPDAMGconvert(DAMGlevelmat *);
void CHERHPDAMGconvert(CAMGlevelmat *);
void ZHERHPDAMGconvert(ZAMGlevelmat *);



void SSYMpiluclsol(integer *,real *,real *,real *, integer *);
void SSYMpilucusol(integer *,real *,real *,real *, integer *);
void DSYMpiluclsol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		   integer *);
void DSYMpilucusol(integer *,doubleprecision *,doubleprecision *,doubleprecision *,
		   integer *);
void CSYMpiluclsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void CSYMpilucusol(integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void CHERpiluclsol(integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void CHERpilucusol(integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void ZSYMpiluclsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *, integer *);
void ZSYMpilucusol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *, integer *);
void ZHERpiluclsol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *, integer *);
void ZHERpilucusol(integer *,ilu_doublecomplex *,ilu_doublecomplex *,ilu_doublecomplex *, integer *);


void SSYMppiluclsol(integer *,integer *,real *,real *,real *, integer *);
void SSYMppilucusol(integer *,integer *,real *,real *,real *, integer *);
void DSYMppiluclsol(integer *,integer *,doubleprecision *,doubleprecision *,
		    doubleprecision *, integer *);
void DSYMppilucusol(integer *,integer *,doubleprecision *,doubleprecision *,
		    doubleprecision *, integer *);
void CSYMppiluclsol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void CSYMppilucusol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void CHERppiluclsol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void CHERppilucusol(integer *,integer *,ilu_complex *,ilu_complex *,ilu_complex *, integer *);
void ZSYMppiluclsol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,
		    ilu_doublecomplex *, integer *);
void ZSYMppilucusol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,
		    ilu_doublecomplex *, integer *);
void ZHERppiluclsol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,
		    ilu_doublecomplex *, integer *);
void ZHERppilucusol(integer *,integer *,ilu_doublecomplex *,ilu_doublecomplex *,
		    ilu_doublecomplex *, integer *);


integer Smps_arms(Smat, integer *, real *, real *,  
		  integer *, real *);
integer Dmps_arms(Dmat, integer *, doubleprecision *, doubleprecision *, 
		  integer *, doubleprecision *);
integer Cmps_arms(Cmat, integer *, real *, real *,  
		  integer *, real *);
integer Zmps_arms(Zmat, integer *, doubleprecision *, doubleprecision *, 
		  integer *, doubleprecision *);

void sprivatesptrs(character *uplo, integer *n, integer *nrhs, real *ap, integer *ipiv, real *b, integer *ldb, integer *info, ftnlen uplolen);
void dprivatesptrs(character *uplo, integer *n, integer *nrhs, doubleprecision *ap, integer *ipiv, doubleprecision *b, integer *ldb, integer *info, ftnlen uplolen);
void cprivatehptrs(character *uplo, integer *n, integer *nrhs, ilu_complex *ap, integer *ipiv, ilu_complex *b, integer *ldb, integer *info, ftnlen uplolen);
void zprivatehptrs(character *uplo, integer *n, integer *nrhs, ilu_doublecomplex *ap, integer *ipiv, ilu_doublecomplex *b, integer *ldb, integer *info, ftnlen uplolen);


void SGNLAMGsetupparameters(Smat *, SILUPACKparam *, integer);
void DGNLAMGsetupparameters(Dmat *, DILUPACKparam *, integer);
void CGNLAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
void ZGNLAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);

void SSPDAMGsetupparameters(Smat *, SILUPACKparam *, integer);
void DSPDAMGsetupparameters(Dmat *, DILUPACKparam *, integer);
void CHPDAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
void ZHPDAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);

void SSYMAMGsetupparameters(Smat *, SILUPACKparam *, integer);
void DSYMAMGsetupparameters(Dmat *, DILUPACKparam *, integer);
void CHERAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
void ZHERAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);
void CSYMAMGsetupparameters(Cmat *, CILUPACKparam *, integer);
void ZSYMAMGsetupparameters(Zmat *, ZILUPACKparam *, integer);



void SILUPACKparamdelete(SILUPACKparam *);
void DILUPACKparamdelete(DILUPACKparam *);
void CILUPACKparamdelete(CILUPACKparam *);
void ZILUPACKparamdelete(ZILUPACKparam *);


void Smergematrices(Smat *C, Smat *A, Smat *B, 
		    integer *p, integer *invq, integer *lenB, 
		    real shift, integer *buff);
void Dmergematrices(Dmat *C, Dmat *A, Dmat *B, 
		    integer *p, integer *invq, integer *lenB, 
		    doubleprecision shift, integer *buff);
void Cmergematrices(Cmat *C, Cmat *A, Cmat *B, 
		    integer *p, integer *invq, integer *lenB, 
		    ilu_complex shift, integer *buff);
void Zmergematrices(Zmat *C, Zmat *A, Zmat *B, 
		    integer *p, integer *invq, integer *lenB, 
		    ilu_doublecomplex shift, integer *buff);

void SSYMmergematrices(Smat *C, Smat *A, Smat *B, 
		       real *rowscale,
		       integer *p, integer *invq, integer *lenB, 
		       real shift, integer *buff);
void DSYMmergematrices(Dmat *C, Dmat *A, Dmat *B, 
		       doubleprecision *rowscale,
		       integer *p, integer *invq, integer *lenB, 
		       doubleprecision shift, integer *buff);
void CSYMmergematrices(Cmat *C, Cmat *A, Cmat *B, 
		       ilu_complex *rowscale,
		       integer *p, integer *invq, integer *lenB, 
		       ilu_complex shift, integer *buff);
void ZSYMmergematrices(Zmat *C, Zmat *A, Zmat *B, 
		       ilu_doublecomplex *rowscale,
		       integer *p, integer *invq, integer *lenB, 
		       ilu_doublecomplex shift, integer *buff);
void CHERmergematrices(Cmat *C, Cmat *A, Cmat *B, 
		       ilu_complex *rowscale,
		       integer *p, integer *invq, integer *lenB, 
		       ilu_complex shift, integer *buff);
void ZHERmergematrices(Zmat *C, Zmat *A, Zmat *B, 
		       ilu_doublecomplex *rowscale,
		       integer *p, integer *invq, integer *lenB, 
		       ilu_doublecomplex shift, integer *buff);



/* global variable to measure the timings of the components within ILUPACK */
void ilupack_dummy_000();
#define ILUPACK_secnds_length   10
#define ILUPACK_mem_length      12

#ifdef  _DECLARE_ILUPACK_GLOBALS_
    double ILUPACK_secnds[ILUPACK_secnds_length];
    size_t ILUPACK_mem[ILUPACK_mem_length];
#else
    extern double ILUPACK_secnds[ILUPACK_secnds_length];
    extern size_t ILUPACK_mem[ILUPACK_mem_length];
#endif

#endif /* _ILU_PACK_H */





