#ifndef _ILU_PACK_MACROS_H
#define _ILU_PACK_MACROS_H

#ifdef _DOUBLE_REAL_
#define weightsC          DweightsC
#define add2is            Dadd2is
#define add2com           Dadd2com
#define CSRMERGEMATRICES  Dmergematrices

#define ILUPACKPARAMDELETE DILUPACKparamdelete
#define ILUCUPDATE        DILUCUPDAT
#define PILUCUPDATE       DPILUCUPDAT
#define SILUCUPDATE       DSILUCUPDAT
#define ILUCLIST          DILUCLIS
#define PILUCLIST         DPILUCLIS
#define BILUCLIST         DBILUCLIST
#define ILDLCUPDATE       DILDLCUPDAT
#define PILDLCUPDATE      DPILDLCUPDAT
#define SILDLCUPDATE      DSILDLCUPDAT

#define ILUPACKINIT       dgnlamginit
#define ILUPACKFACTOR     dgnlamgfactor
#define ILUPACKSOLVER     dgnlamgsolver
#define ILUPACKSOL        dgnlamgsol
#define ILUPACKDELETE     dgnlamgdelete
#define ILUPACKINFO       dgnlamginfo
#define ILUPACKNNZ        dgnlamgnnz
#define ILUPACKUNDOSCALING damgundoscaling

#define SPDILUPACKINIT    dspdamginit
#define SPDILUPACKFACTOR  dspdamgfactor
#define SPDILUPACKSOLVER  dspdamgsolver
#define SPDILUPACKSOL     dspdamgsol
#define SPDILUPACKDELETE  dspdamgdelete
#define SPDILUPACKINFO    dspdamginfo
#define SPDILUPACKNNZ     dspdamgnnz

#define SYMSPDILUPACKCONVERT      dsymspdamgconvert
#define SYMSPDCONVERT             DSYMSPDAMGconvert

#define GETEPS       dgeteps
#define _REAL_MAX_   _D_REAL_MAX_
#define REALS        doubleprecision
#define FLOAT        doubleprecision
#define CSRMAT       Dmat
#define RCSRMAT      Dmat
#define PRO2SUM      Dpro2sum
#define INITARMS     Dinitarms
#define MPSARMS      Dmps_arms
#define AMGLEVELMAT  DAMGlevelmat
#define ILUPACKPARAM DILUPACKparam
#define GNLAGGRESSIVEDROPPING  DGNLaggressivedropping 
#define GNLFAGGRESSIVEDROPPING DGNLfaggressivedroppin
#define SPDAGGRESSIVEDROPPING  DSPDaggressivedropping 
#define SPDFAGGRESSIVEDROPPING DSPDfaggressivedroppin
#define SPTRF        dsptrf
#define SPTRS        dsptrs
#define PPTRF        dpptrf
#define PPTRS        dpptrs
#define GETRF        dgetrf
#define GETRS        dgetrs
#define GNLSYMMWM    DGNLsmwm

#define ROWSCALE     Drowscale
#define COLSCALE     Dcolscale
#define SPDSCALE     DSPDscale
#define SYMSCALE     DSYMscale
#define CNRMS        Dcnrms
#define RNRMS        Drnrms
#define DIAMUA       DDiagMultA
#define AMUDIA       DAMultDiag

#define CSRMATVEC    DGNLmatvec
#define CSRMATTVEC   DGNLmattvec
#define CSRMATHVEC   DGNLmathvec
#define SYMMATVEC    DSYMmatvec
#define SSMMATVEC    DSSMmatvec
#define HERMATVEC    DSYMmatvec
#define SETUPGRAPH   Dsetupgraph
#define SETUPGRAPH_EPSILON   Dsetupgraph_epsilon
#define SETUPGRAPH_EPSILON_SP   Dsetupgraph_epsilon_sp

#define AMGDELETE       DGNLAMGdelete
#define SPDAMGDELETE    DSPDAMGdelete

#define AMGSETUP        DGNLAMGfactor
#define SPDAMGSETUP     DSPDAMGfactor
#define AMGFACTOR       DGNLAMGfactor
#define SPDAMGFACTOR    DSPDAMGfactor
#define AMGEXTRACT      DGNLAMGextract
#define SPDAMGEXTRACT   DSYMAMGextract
#define SYMAMGEXTRACT   DSYMAMGextract

#define AMGSOL1         DGNLAMGsol1
#define AMGTSOL1        DGNLAMGtsol1
#define AMGSOL2         DGNLAMGsol2
#define AMGTSOL2        DGNLAMGtsol2
#define SPDAMGSOL1      DSPDAMGsol1
#define SPDAMGSOL2      DSPDAMGsol2
#define SYMAMGSOL1      DSYMAMGsol1
#define SYMAMGSOL2      DSYMAMGsol2

#define AMGSOL          DGNLAMGsol_internal
#define AMGSOL_EXT      DGNLAMGsol
#define AMGTSOL         DGNLAMGtsol_internal
#define AMGTSOL_EXT     DGNLAMGtsol
#define AMGDLSOL        DGNLAMGdlsol
#define AMGDUSOL        DGNLAMGdusol
#define AMGLSOL         DGNLAMGlsol
#define AMGUSOL         DGNLAMGusol
#define AMGTDLSOL       DGNLAMGtdlsol
#define AMGTDUSOL       DGNLAMGtdusol
#define AMGTLSOL        DGNLAMGtlsol
#define AMGTUSOL        DGNLAMGtusol
#define AMGDRIVER       DGNLAMGdriver
#define AMGSOLVER       DGNLAMGsolver
#define SPDAMGSOLVER    DSPDAMGsolver
#define SYMAMGSOLVER    DSYMAMGsolver
#define AMGINIT         DGNLAMGinit
#define AMGGETPARAMS    DGNLAMGgetparams
#define AMGSETPARAMS    DGNLAMGsetparams
#define AMGSETUPPARAMETERS    DGNLAMGsetupparameters
#define SPDAMGSETUPPARAMETERS    DSPDAMGsetupparameters
#define SPDAMGSOL       DSPDAMGsol_internal
#define SPDAMGSOL_EXT   DSPDAMGsol
#define SPDAMGINIT      DSPDAMGinit
#define SPDAMGGETPARAMS DSPDAMGgetparams
#define SPDAMGSETPARAMS DSPDAMGsetparams
#define SYMAMGINIT      DSYMAMGinit
#define SYMAMGGETPARAMS DSYMAMGgetparams
#define SYMAMGSETPARAMS DSYMAMGsetparams

#define LDLP         DSPDldlp
#define LDLPSOL      DSPDldlpsol
#define LUPQ         DGNLlupq
#define LUPQSOL      DGNLlupqsol
#define LUPQLSOL     DGNLlupqlsol
#define LUPQUSOL     DGNLlupqusol
#define LUPQDLSOL    DGNLlupqdlsol
#define LUPQDUSOL    DGNLlupqdusol
#define LUPQTSOL     DGNLlupqtsol
#define LUPQTLSOL    DGNLlupqtlsol
#define LUPQTUSOL    DGNLlupqtusol
#define LUPQTDLSOL   DGNLlupqtdlsol
#define LUPQTDUSOL   DGNLlupqtdusol


#define ILDLC        DSYMildlc
#define ILDLCSOL     DSYMildlcsol

#define PILDLCDLSOL  DSYMpildlcdlsol
#define PILDLCDUSOL  DSYMpildlcdusol
#define PILDLCLSOL   DSYMpildlclsol
#define PILDLCUSOL   DSYMpildlcusol

#define ILUC         DGNLiluc
#define PILUC        DGNLpiluc
#define SPILUC       DGNLspiluc
#define MPILUC       DGNLmpiluc
#define PILDLC       DSPDpiluc
#define MPILDLC      DSPDmpiluc
#define ILUCSOL      DGNLilucsol
#define ILUCTSOL     DGNLiluctsol
#define ILUCDLSOL    DGNLilucdlsol
#define ILUCTDLSOL   DGNLiluctdlsol
#define ILUCDUSOL    DGNLilucdusol
#define ILUCTDUSOL   DGNLiluctdusol
#define ILUCLSOL     DGNLiluclsol
#define ILUCTLSOL    DGNLiluctlsol
#define ILUCUSOL     DGNLilucusol
#define ILUCTUSOL    DGNLiluctusol

#define PILUCDLSOL    DGNLpilucdlsol
#define PILUCTDLSOL   DGNLpiluctdlsol
#define PILUCDUSOL    DGNLpilucdusol
#define PILUCTDUSOL   DGNLpiluctdusol
#define PILUCLSOL     DGNLpiluclsol
#define PILUCTLSOL    DGNLpiluctlsol
#define PILUCUSOL     DGNLpilucusol
#define PILUCTUSOL    DGNLpiluctusol

#define ILUTP        DGNLilutp
#define ILUT         DGNLilut
#define LUSOL        DGNLlusol
#define LUTSOL       DGNLlutsol
#define LULSOL       DGNLlulsol
#define LUTLSOL      DGNLlutlsol
#define LUUSOL       DGNLluusol
#define LUTUSOL      DGNLlutusol
#define LUDLSOL      DGNLludlsol
#define LUTDLSOL     DGNLlutdlsol
#define LUDUSOL      DGNLludusol
#define LUTDUSOL     DGNLlutdusol

#define CPERM        Dcperm
#define RPERM        Drperm

#define PERMNULL     DGNLperm_null
#define PERMND       DGNLperm_nd
#define PERMRCM      DGNLperm_rcm
#define PERMAMF      DGNLperm_amf
#define PERMAMD      DGNLperm_amd
#define PERMMMD      DGNLperm_mmd
#define PERMMETISE   DGNLperm_metis_e
#define PERMMETISN   DGNLperm_metis_n

#define PERMNDFC       DGNLperm_nd_fc
#define PERMRCMFC      DGNLperm_rcm_fc
#define PERMAMFFC      DGNLperm_amf_fc
#define PERMAMDFC      DGNLperm_amd_fc
#define PERMMMDFC      DGNLperm_mmd_fc
#define PERMMETISEFC   DGNLperm_metis_e_fc
#define PERMMETISNFC   DGNLperm_metis_n_fc
#define PERMNDFCV      DGNLperm_nd_fcv
#define PERMRCMFCV     DGNLperm_rcm_fcv
#define PERMAMFFCV     DGNLperm_amf_fcv
#define PERMAMDFCV     DGNLperm_amd_fcv
#define PERMMMDFCV     DGNLperm_mmd_fcv
#define PERMMETISEFCV  DGNLperm_metis_e_fcv
#define PERMMETISNFCV  DGNLperm_metis_n_fcv
#define SYMPERMNDFC       DSYMperm_nd_fc
#define SYMPERMRCMFC      DSYMperm_rcm_fc
#define SYMPERMAMFFC      DSYMperm_amf_fc
#define SYMPERMAMDFC      DSYMperm_amd_fc
#define SYMPERMMMDFC      DSYMperm_mmd_fc
#define SYMPERMMETISEFC   DSYMperm_metis_e_fc
#define SYMPERMMETISNFC   DSYMperm_metis_n_fc
#define SYMPERMNDFCV      DSYMperm_nd_fcv
#define SYMPERMRCMFCV     DSYMperm_rcm_fcv
#define SYMPERMAMFFCV     DSYMperm_amf_fcv
#define SYMPERMAMDFCV     DSYMperm_amd_fcv
#define SYMPERMMMDFCV     DSYMperm_mmd_fcv
#define SYMPERMMETISEFCV  DSYMperm_metis_e_fcv
#define SYMPERMMETISNFCV  DSYMperm_metis_n_fcv

#define PERMPQ       DGNLperm_pq
#define PERMFC       DGNLperm_fc
#define INDFC        Dindfc
#define INDFC_RS     Dindfc_rs
#define SYMPERMFC    DSYMperm_fc
#define SYMINDFC     DSYMindfc
#define SYMINDFC_RS  DSYMindfc_rs
#define SPDPERMFC    DSPDperm_fc
#define SPDINDFC     DSPDindfc
#define SPDINDFC_RS  DSPDindfc_rs

#define PERMFCV      DGNLperm_fcv
#define INDFCV       Dindfcv
#define INDFCV_RS    Dindfcv_rs
#define SYMPERMFCV   DSYMperm_fcv
#define SYMINDFCV    DSYMindfcv
#define SYMINDFCV_RS DSYMindfcv_rs
#define SPDPERMFCV   DSPDperm_fcv
#define SPDINDFCV    DSPDindfcv
#define SPDINDFCV_RS DSPDindfcv_rs

#define PERMP        DGNLperm_p
#define PERMINDSET   DGNLperm_indset

#define SPDPERMNULL     DSPDperm_null
#define SPDPERMND       DSPDperm_nd
#define SPDPERMRCM      DSPDperm_rcm
#define SPDPERMAMF      DSPDperm_amf
#define SPDPERMAMD      DSPDperm_amd
#define SPDPERMMMD      DSPDperm_mmd
#define SPDPERMMETISE   DSPDperm_metis_e
#define SPDPERMMETISN   DSPDperm_metis_n
#define SPDPERMPP       DSPDperm_pp
#define SPDPERMINDSET   DSPDperm_indset


#define SPDPERMNDFC       DSPDperm_nd_fc
#define SPDPERMRCMFC      DSPDperm_rcm_fc
#define SPDPERMAMFFC      DSPDperm_amf_fc
#define SPDPERMAMDFC      DSPDperm_amd_fc
#define SPDPERMMMDFC      DSPDperm_mmd_fc
#define SPDPERMMETISEFC   DSPDperm_metis_e_fc
#define SPDPERMMETISNFC   DSPDperm_metis_n_fc

#define SPDPERMNDFCV      DSPDperm_nd_fcv
#define SPDPERMRCMFCV     DSPDperm_rcm_fcv
#define SPDPERMAMFFCV     DSPDperm_amf_fcv
#define SPDPERMAMDFCV     DSPDperm_amd_fcv
#define SPDPERMMMDFCV     DSPDperm_mmd_fcv
#define SPDPERMMETISEFCV  DSPDperm_metis_e_fcv
#define SPDPERMMETISNFCV  DSPDperm_metis_n_fcv


#define SYMPERMNULL     DSYMperm_null
#define SYMPERMND       DSYMperm_nd
#define SYMPERMRCM      DSYMperm_rcm
#define SYMPERMAMF      DSYMperm_amf
#define SYMPERMAMD      DSYMperm_amd
#define SYMPERMMMD      DSYMperm_mmd
#define SYMPERMMETISE   DSYMperm_metis_e
#define SYMPERMMETISN   DSYMperm_metis_n
#define SYMPERMINDSET   DSYMperm_indset

#define PERMMWMRCM    DGNLperm_mwm_rcm
#define PERMMWMMMD    DGNLperm_mwm_mmd
#define PERMMWMAMF    DGNLperm_mwm_amf
#define PERMMWMAMD    DGNLperm_mwm_amd
#define PERMMWMMETISE DGNLperm_mwm_metis_e
#define PERMMWMMETISN DGNLperm_mwm_metis_n

#define SYMPERMMWMRCM    DSYMperm_mwm_rcm
#define SYMPERMMWMMMD    DSYMperm_mwm_mmd
#define SYMPERMMWMAMF    DSYMperm_mwm_amf
#define SYMPERMMWMAMD    DSYMperm_mwm_amd
#define SYMPERMMWMMETISE DSYMperm_mwm_metis_e
#define SYMPERMMWMMETISN DSYMperm_mwm_metis_n

#define SYMPERMMATCHINGRCM    DSYMperm_matching_rcm
#define SYMPERMMATCHINGMMD    DSYMperm_matching_mmd
#define SYMPERMMATCHINGAMF    DSYMperm_matching_amf
#define SYMPERMMATCHINGAMD    DSYMperm_matching_amd
#define SYMPERMMATCHINGMETISE DSYMperm_matching_metis_e
#define SYMPERMMATCHINGMETISN DSYMperm_matching_metis_n

#define PERMMC64NULL   DGNLperm_mc64_null
#define PERMMC64RCM    DGNLperm_mc64_rcm
#define PERMMC64MMD    DGNLperm_mc64_mmd
#define PERMMC64AMF    DGNLperm_mc64_amf
#define PERMMC64AMD    DGNLperm_mc64_amd
#define PERMMC64METISE DGNLperm_mc64_metis_e
#define PERMMC64METISN DGNLperm_mc64_metis_n

#define PERMMATCHINGNULL   DGNLperm_matching_null
#define PERMMATCHINGRCM    DGNLperm_matching_rcm
#define PERMMATCHINGMMD    DGNLperm_matching_mmd
#define PERMMATCHINGAMF    DGNLperm_matching_amf
#define PERMMATCHINGAMD    DGNLperm_matching_amd
#define PERMMATCHINGMETISE DGNLperm_matching_metis_e
#define PERMMATCHINGMETISN DGNLperm_matching_metis_n

#define PERMMWARCM    DGNLperm_mwa_rcm
#define PERMMWAMMD    DGNLperm_mwa_mmd
#define PERMMWAAMF    DGNLperm_mwa_amf
#define PERMMWAAMD    DGNLperm_mwa_amd
#define PERMMWAMETISE DGNLperm_mwa_metis_e
#define PERMMWAMETISN DGNLperm_mwa_metis_n

#define SPARTRAN     Dspartran
#define ETREE        Detree
#define QSORT        Dqsort
#define QSORT2       Dqsort2
#define QQSORT       Dqqsort
#define QQSORT2      Dqqsort2
#define QQSORTS      Dqqsorts
#define QQSORTS2     Dqqsorts2
#define BQSORT       Dbqsort
#define QSPLIT       Dqsplit
#define QSPLIT2      Dqsplit2
#define QSPLIT3      Dqsplit3
#define BQSPLIT      Dbqsplit
#define QSEPARATE    Dqseparate
#define CLEAR        Dclear
#define INDSET       Dindset
#define PQPERM       DPQpermF
#define PPERM        DPpermF
#define PPPERM       DPPpermF
#define WDIAG        DwDiagRF
#define WPDIAG       DwPDiagRF
#define QSORTR2I     dqsortr2i
#define SWAPM        dswapm
#define FSWAP        dswap
#define SPDWDIAG     DSPDwDiagRF
#define BISINIT      Dbisinit
#define STOPBIS      Dstopbis
#define TIDYCG       Dtidycg
#define MGSRO        Dmgsro
#define GIVENS       Dgivens
#define BRKDN        Dbrkdn

#define CSRCSC       Dcsrcsc
#define CSRCSC2      Dcsrcsc2
#define READMTC      Dreadmtc
#define WRITEMTC     Dwritemtc
#define READVECTORS  Dreadvectors
#define WRITEVECTORS Dwritevectors
#define PCG          Dpcg
#define FPCG         Dfpcg
#define BCG          Dbcg
#define SBCG         DSYMbcg
#define SQMR         DSYMqmr
#define GMRES        Dgmres
#define FGMRES       Dfgmres

#define DISTDOT      ddot
#define DISTDOTU     ddot
#define DISTDOT2     ddot2
#define GEMM         dgemm
#define COPY         dcopy
#define AXPY         daxpy
#define ROTG         drotg
#define ROT          drot
#define NRM          dnrm2
#define I_AMAX       idamax
#define ASUM         dasum
#define SCAL         dscal

#define MC64I        mc64id
#define MC64A        mc64ad
#define MC64IR       mc64id
#define MC64AR       mc64ad

#define MUMPS_MATCH   dmumps_match
#define MUMPS_MATCHR  dmumps_match

#define FABS(A)      (((A)>=0)?(A):(-(A)))
#define ABS(A)       dabs(A)
#define FLOATABS(A)  dabs(A)
#define FLOATABSNAME dabs
#define SQRT(A)      dsqrt(A)
#define LOG(A)       dlog(A)
#define LOG10(A)     dlog10(A)
#define IREAL(A)     dble(A)
#define ABSFNAME     dabs
#define SQRTFNAME    dsqrt
#define LOGFNAME     dlog
#define LOGTENFNAME  dlog10
#define CTOD(A)      (A)
#define CONJ(A)      (A)
#define FNULL        0.0d0
#define FONE         1.0d0
#define RONE         1.0d0
#define RTWO         2.0d0
#define RFOUR        4.0d0
#define RZERO        0.0d0
#define REAL_MAX     1.0d300
#define SIGNUM       1.0d0
#define SIGNUM2      1.0d0
#define TWICE        1



#elif defined _SINGLE_REAL_
#define weightsC          SweightsC
#define add2is            Sadd2is
#define add2com           Sadd2com
#define CSRMERGEMATRICES  Smergematrices

#define ILUPACKPARAMDELETE SILUPACKparamdelete
#define ILUCUPDATE        SILUCUPDAT
#define PILUCUPDATE       SPILUCUPDAT
#define SILUCUPDATE       SSILUCUPDAT
#define ILUCLIST          SILUCLIS
#define PILUCLIST         SPILUCLIS
#define BILUCLIST         SBILUCLIST
#define ILDLCUPDATE       SILDLCUPDAT
#define PILDLCUPDATE      SPILDLCUPDAT
#define SILDLCUPDATE      SSILDLCUPDAT

#define ILUPACKINIT       sgnlamginit
#define ILUPACKFACTOR     sgnlamgfactor
#define ILUPACKSOLVER     sgnlamgsolver
#define ILUPACKSOL        sgnlamgsol
#define ILUPACKDELETE     sgnlamgdelete
#define ILUPACKINFO       sgnlamginfo
#define ILUPACKNNZ        sgnlamgnnz
#define ILUPACKUNDOSCALING samgundoscaling

#define SPDILUPACKINIT    sspdamginit
#define SPDILUPACKFACTOR  sspdamgfactor
#define SPDILUPACKSOLVER  sspdamgsolver
#define SPDILUPACKSOL     sspdamgsol
#define SPDILUPACKDELETE  sspdamgdelete
#define SPDILUPACKINFO    sspdamginfo
#define SPDILUPACKNNZ     sspdamgnnz

#define SYMSPDILUPACKCONVERT      ssymspdamgconvert
#define SYMSPDCONVERT             SSYMSPDAMGconvert

#define GETEPS       sgeteps
#define _REAL_MAX_   _S_REAL_MAX_
#define REALS        real
#define FLOAT        real
#define CSRMAT       Smat
#define RCSRMAT      Smat
#define PRO2SUM      Spro2sum
#define INITARMS     Sinitarms
#define MPSARMS      Smps_arms
#define AMGLEVELMAT  SAMGlevelmat
#define ILUPACKPARAM SILUPACKparam
#define GNLAGGRESSIVEDROPPING  SGNLaggressivedropping 
#define GNLFAGGRESSIVEDROPPING SGNLfaggressivedroppin
#define SPDAGGRESSIVEDROPPING  SSPDaggressivedropping 
#define SPDFAGGRESSIVEDROPPING SSPDfaggressivedroppin
#define SPTRF        ssptrf
#define SPTRS        ssptrs
#define PPTRF        spptrf
#define PPTRS        spptrs
#define GETRF        sgetrf
#define GETRS        sgetrs
#define GNLSYMMWM    SGNLsmwm

#define ROWSCALE     Srowscale
#define COLSCALE     Scolscale
#define SPDSCALE     SSPDscale
#define SYMSCALE     SSYMscale
#define CNRMS        Scnrms
#define RNRMS        Srnrms
#define DIAMUA       SDiagMultA
#define AMUDIA       SAMultDiag

#define CSRMATVEC    SGNLmatvec
#define CSRMATTVEC   SGNLmattvec
#define CSRMATHVEC   SGNLmathvec
#define SYMMATVEC    SSYMmatvec
#define SSMMATVEC    SSSMmatvec
#define HERMATVEC    SSYMmatvec
#define SETUPGRAPH   Ssetupgraph
#define SETUPGRAPH_EPSILON   Ssetupgraph_epsilon
#define SETUPGRAPH_EPSILON_SP   Ssetupgraph_epsilon_sp

#define AMGDELETE       SGNLAMGdelete
#define SPDAMGDELETE    SSPDAMGdelete

#define AMGSETUP        SGNLAMGfactor
#define SPDAMGSETUP     SSPDAMGfactor
#define AMGFACTOR       SGNLAMGfactor
#define SPDAMGFACTOR    SSPDAMGfactor
#define AMGEXTRACT      SGNLAMGextract
#define SPDAMGEXTRACT   SSYMAMGextract
#define SYMAMGEXTRACT   SSYMAMGextract

#define AMGSOL1         SGNLAMGsol1
#define AMGTSOL1        SGNLAMGtsol1
#define AMGSOL2         SGNLAMGsol2
#define AMGTSOL2        SGNLAMGtsol2
#define SPDAMGSOL1      SSPDAMGsol1
#define SPDAMGSOL2      SSPDAMGsol2
#define SYMAMGSOL1      SSYMAMGsol1
#define SYMAMGSOL2      SSYMAMGsol2

#define AMGSOL          SGNLAMGsol_internal
#define AMGSOL_EXT      SGNLAMGsol
#define AMGDLSOL        SGNLAMGdlsol
#define AMGDUSOL        SGNLAMGdusol
#define AMGLSOL         SGNLAMGlsol
#define AMGUSOL         SGNLAMGusol
#define AMGTDLSOL       SGNLAMGtdlsol
#define AMGTDUSOL       SGNLAMGtdusol
#define AMGTLSOL        SGNLAMGtlsol
#define AMGTUSOL        SGNLAMGtusol
#define AMGTSOL         SGNLAMGtsol_internal
#define AMGTSOL_EXT     SGNLAMGtsol
#define SPDAMGSOL       SSPDAMGsol_internal
#define SPDAMGSOL_EXT   SSPDAMGsol
#define AMGDRIVER       SGNLAMGdriver
#define AMGSOLVER       SGNLAMGsolver
#define SPDAMGSOLVER    SSPDAMGsolver
#define SYMAMGSOLVER    SSYMAMGsolver
#define AMGINIT         SGNLAMGinit
#define AMGGETPARAMS    SGNLAMGgetparams
#define AMGSETPARAMS    SGNLAMGsetparams
#define AMGSETUPPARAMETERS    SGNLAMGsetupparameters
#define SPDAMGSETUPPARAMETERS    SSPDAMGsetupparameters
#define SPDAMGINIT      SSPDAMGinit
#define SPDAMGGETPARAMS SSPDAMGgetparams
#define SPDAMGSETPARAMS SSPDAMGsetparams
#define SYMAMGINIT      SSYMAMGinit
#define SYMAMGGETPARAMS SSYMAMGgetparams
#define SYMAMGSETPARAMS SSYMAMGsetparams

#define LDLP         SSPDldlp
#define LDLPSOL      SSPDldlpsol
#define LUPQ         SGNLlupq
#define LUPQSOL      SGNLlupqsol
#define LUPQLSOL     SGNLlupqlsol
#define LUPQUSOL     SGNLlupqusol
#define LUPQDLSOL    SGNLlupqdlsol
#define LUPQDUSOL    SGNLlupqdusol
#define LUPQTSOL     SGNLlupqtsol
#define LUPQTLSOL    SGNLlupqtlsol
#define LUPQTUSOL    SGNLlupqtusol
#define LUPQTDLSOL   SGNLlupqtdlsol
#define LUPQTDUSOL   SGNLlupqtdusol



#define ILDLC        SSYMildlc
#define ILDLCSOL     SSYMildlcsol

#define PILDLCDLSOL  SSYMpildlcdlsol
#define PILDLCDUSOL  SSYMpildlcdusol
#define PILDLCLSOL   SSYMpildlclsol
#define PILDLCUSOL   SSYMpildlcusol

#define ILUC         SGNLiluc
#define PILUC        SGNLpiluc
#define SPILUC       SGNLspiluc
#define MPILUC       SGNLmpiluc
#define PILDLC       SSPDpiluc
#define MPILDLC      SSPDmpiluc
#define ILUCSOL      SGNLilucsol
#define ILUCTSOL     SGNLiluctsol
#define ILUCDLSOL    SGNLilucdlsol
#define ILUCTDLSOL   SGNLiluctdlsol
#define ILUCDUSOL    SGNLilucdusol
#define ILUCTDUSOL   SGNLiluctdusol
#define ILUCLSOL     SGNLiluclsol
#define ILUCTLSOL    SGNLiluctlsol
#define ILUCUSOL     SGNLilucusol
#define ILUCTUSOL    SGNLiluctusol

#define PILUCDLSOL    SGNLpilucdlsol
#define PILUCTDLSOL   SGNLpiluctdlsol
#define PILUCDUSOL    SGNLpilucdusol
#define PILUCTDUSOL   SGNLpiluctdusol
#define PILUCLSOL     SGNLpiluclsol
#define PILUCTLSOL    SGNLpiluctlsol
#define PILUCUSOL     SGNLpilucusol
#define PILUCTUSOL    SGNLpiluctusol

#define ILUTP        SGNLilutp
#define ILUT         SGNLilut
#define LUSOL        SGNLlusol
#define LUTSOL       SGNLlutsol
#define LUUSOL       SGNLluusol
#define LUTUSOL      SGNLlutusol
#define LULSOL       SGNLlulsol
#define LUTLSOL      SGNLlutlsol
#define LUDUSOL      SGNLludusol
#define LUTDUSOL     SGNLlutdusol
#define LUDLSOL      SGNLludlsol
#define LUTDLSOL     SGNLlutdlsol

#define CPERM        Scperm
#define RPERM        Srperm

#define PERMNULL     SGNLperm_null
#define PERMND       SGNLperm_nd
#define PERMRCM      SGNLperm_rcm
#define PERMAMF      SGNLperm_amf
#define PERMAMD      SGNLperm_amd
#define PERMMMD      SGNLperm_mmd
#define PERMMETISE   SGNLperm_metis_e
#define PERMMETISN   SGNLperm_metis_n
#define PERMPQ       SGNLperm_pq

#define PERMNDFC       SGNLperm_nd_fc
#define PERMRCMFC      SGNLperm_rcm_fc
#define PERMAMFFC      SGNLperm_amf_fc
#define PERMAMDFC      SGNLperm_amd_fc
#define PERMMMDFC      SGNLperm_mmd_fc
#define PERMMETISEFC   SGNLperm_metis_e_fc
#define PERMMETISNFC   SGNLperm_metis_n_fc
#define PERMNDFCV      SGNLperm_nd_fcv
#define PERMRCMFCV     SGNLperm_rcm_fcv
#define PERMAMFFCV     SGNLperm_amf_fcv
#define PERMAMDFCV     SGNLperm_amd_fcv
#define PERMMMDFCV     SGNLperm_mmd_fcv
#define PERMMETISEFCV  SGNLperm_metis_e_fcv
#define PERMMETISNFCV  SGNLperm_metis_n_fcv
#define SYMPERMNDFC       SSYMperm_nd_fc
#define SYMPERMRCMFC      SSYMperm_rcm_fc
#define SYMPERMAMFFC      SSYMperm_amf_fc
#define SYMPERMAMDFC      SSYMperm_amd_fc
#define SYMPERMMMDFC      SSYMperm_mmd_fc
#define SYMPERMMETISEFC   SSYMperm_metis_e_fc
#define SYMPERMMETISNFC   SSYMperm_metis_n_fc
#define SYMPERMNDFCV      SSYMperm_nd_fcv
#define SYMPERMRCMFCV     SSYMperm_rcm_fcv
#define SYMPERMAMFFCV     SSYMperm_amf_fcv
#define SYMPERMAMDFCV     SSYMperm_amd_fcv
#define SYMPERMMMDFCV     SSYMperm_mmd_fcv
#define SYMPERMMETISEFCV  SSYMperm_metis_e_fcv
#define SYMPERMMETISNFCV  SSYMperm_metis_n_fcv


#define PERMFC       SGNLperm_fc
#define INDFC        Sindfc
#define INDFC_RS     Sindfc_rs
#define SYMPERMFC    SSYMperm_fc
#define SYMINDFC     SSYMindfc
#define SYMINDFC_RS  SSYMindfc_rs
#define SPDPERMFC    SSPDperm_fc
#define SPDINDFC     SSPDindfc
#define SPDINDFC_RS  SSPDindfc_rs

#define PERMFCV      SGNLperm_fcv
#define INDFCV       Sindfcv
#define INDFCV_RS    Sindfcv_rs
#define SYMPERMFCV   SSYMperm_fcv
#define SYMINDFCV    SSYMindfcv
#define SYMINDFCV_RS SSYMindfcv_rs
#define SPDPERMFCV   SSPDperm_fcv
#define SPDINDFCV    SSPDindfcv
#define SPDINDFCV_RS SSPDindfcv_rs

#define PERMP        SGNLperm_p
#define PERMINDSET   SGNLperm_indset

#define SPDPERMNULL     SSPDperm_null
#define SPDPERMND       SSPDperm_nd
#define SPDPERMRCM      SSPDperm_rcm
#define SPDPERMAMF      SSPDperm_amf
#define SPDPERMAMD      SSPDperm_amd
#define SPDPERMMMD      SSPDperm_mmd
#define SPDPERMMETISE   SSPDperm_metis_e
#define SPDPERMMETISN   SSPDperm_metis_n
#define SPDPERMPP       SSPDperm_pp
#define SPDPERMINDSET   SSPDperm_indset


#define SPDPERMNDFC       SSPDperm_nd_fc
#define SPDPERMRCMFC      SSPDperm_rcm_fc
#define SPDPERMAMFFC      SSPDperm_amf_fc
#define SPDPERMAMDFC      SSPDperm_amd_fc
#define SPDPERMMMDFC      SSPDperm_mmd_fc
#define SPDPERMMETISEFC   SSPDperm_metis_e_fc
#define SPDPERMMETISNFC   SSPDperm_metis_n_fc

#define SPDPERMNDFCV      SSPDperm_nd_fcv
#define SPDPERMRCMFCV     SSPDperm_rcm_fcv
#define SPDPERMAMFFCV     SSPDperm_amf_fcv
#define SPDPERMAMDFCV     SSPDperm_amd_fcv
#define SPDPERMMMDFCV     SSPDperm_mmd_fcv
#define SPDPERMMETISEFCV  SSPDperm_metis_e_fcv
#define SPDPERMMETISNFCV  SSPDperm_metis_n_fcv


#define SYMPERMNULL     SSYMperm_null
#define SYMPERMND       SSYMperm_nd
#define SYMPERMRCM      SSYMperm_rcm
#define SYMPERMAMF      SSYMperm_amf
#define SYMPERMAMD      SSYMperm_amd
#define SYMPERMMMD      SSYMperm_mmd
#define SYMPERMMETISE   SSYMperm_metis_e
#define SYMPERMMETISN   SSYMperm_metis_n
#define SYMPERMINDSET   SSYMperm_indset

#define PERMMWMRCM    SGNLperm_mwm_rcm
#define PERMMWMMMD    SGNLperm_mwm_mmd
#define PERMMWMAMF    SGNLperm_mwm_amf
#define PERMMWMAMD    SGNLperm_mwm_amd
#define PERMMWMMETISE SGNLperm_mwm_metis_e
#define PERMMWMMETISN SGNLperm_mwm_metis_n

#define SYMPERMMWMRCM    SSYMperm_mwm_rcm
#define SYMPERMMWMMMD    SSYMperm_mwm_mmd
#define SYMPERMMWMAMF    SSYMperm_mwm_amf
#define SYMPERMMWMAMD    SSYMperm_mwm_amd
#define SYMPERMMWMMETISE SSYMperm_mwm_metis_e
#define SYMPERMMWMMETISN SSYMperm_mwm_metis_n

#define SYMPERMMATCHINGRCM    SSYMperm_matching_rcm
#define SYMPERMMATCHINGMMD    SSYMperm_matching_mmd
#define SYMPERMMATCHINGAMF    SSYMperm_matching_amf
#define SYMPERMMATCHINGAMD    SSYMperm_matching_amd
#define SYMPERMMATCHINGMETISE SSYMperm_matching_metis_e
#define SYMPERMMATCHINGMETISN SSYMperm_matching_metis_n

#define PERMMC64NULL   SGNLperm_mc64_null
#define PERMMC64RCM    SGNLperm_mc64_rcm
#define PERMMC64MMD    SGNLperm_mc64_mmd
#define PERMMC64AMF    SGNLperm_mc64_amf
#define PERMMC64AMD    SGNLperm_mc64_amd
#define PERMMC64METISE SGNLperm_mc64_metis_e
#define PERMMC64METISN SGNLperm_mc64_metis_n

#define PERMMATCHINGNULL   SGNLperm_matching_null
#define PERMMATCHINGRCM    SGNLperm_matching_rcm
#define PERMMATCHINGMMD    SGNLperm_matching_mmd
#define PERMMATCHINGAMF    SGNLperm_matching_amf
#define PERMMATCHINGAMD    SGNLperm_matching_amd
#define PERMMATCHINGMETISE SGNLperm_matching_metis_e
#define PERMMATCHINGMETISN SGNLperm_matching_metis_n

#define PERMMWARCM    SGNLperm_mwa_rcm
#define PERMMWAMMD    SGNLperm_mwa_mmd
#define PERMMWAAMF    SGNLperm_mwa_amf
#define PERMMWAAMD    SGNLperm_mwa_amd
#define PERMMWAMETISE SGNLperm_mwa_metis_e
#define PERMMWAMETISN SGNLperm_mwa_metis_n


#define SPARTRAN     Sspartran
#define ETREE        Setree
#define QQSORT       Sqqsort
#define QSORT        Sqsort
#define QQSORT2      Sqqsort2
#define QQSORTS      Sqqsorts
#define QQSORTS2     Sqqsorts2
#define QSORT2       Sqsort2
#define BQSORT       Sbqsort
#define QSPLIT       Sqsplit
#define QSPLIT2      Sqsplit2
#define QSPLIT3      Sqsplit3
#define BQSPLIT      Sbqsplit
#define QSEPARATE    Sqseparate
#define CLEAR        Sclear
#define INDSET       Sindset
#define PQPERM       SPQpermF
#define PPERM        SPpermF
#define PPPERM       SPPpermF
#define WDIAG        SwDiagRF
#define QSORTR2I     sqsortr2i
#define SWAPM        sswapm
#define FSWAP        sswap
#define WPDIAG       SwPDiagRF
#define SPDWDIAG     SSPDwDiagRF
#define BISINIT      Sbisinit
#define STOPBIS      Sstopbis
#define TIDYCG       Stidycg
#define MGSRO        Smgsro
#define GIVENS       Sgivens
#define BRKDN        Sbrkdn

#define CSRCSC       Scsrcsc
#define CSRCSC2      Scsrcsc2
#define READMTC      Sreadmtc
#define WRITEMTC     Swritemtc
#define READVECTORS  Sreadvectors
#define WRITEVECTORS Swritevectors
#define PCG          Spcg
#define FPCG         Sfpcg
#define BCG          Sbcg
#define SBCG         SSYMbcg
#define SQMR         SSYMqmr
#define GMRES        Sgmres
#define FGMRES       Sfgmres

#define DISTDOT      sdot
#define DISTDOTU     sdot
#define DISTDOT2     sdot2
#define GEMM         sgemm
#define COPY         scopy
#define AXPY         saxpy
#define ROTG         srotg
#define ROT          srot
#define NRM          snrm2
#define I_AMAX       isamax
#define ASUM         sasum
#define SCAL         sscal

#define MC64I        mc64is
#define MC64A        mc64as
#define MC64IR       mc64is
#define MC64AR       mc64as

#define MUMPS_MATCH   smumps_match
#define MUMPS_MATCHR  smumps_match

#define FABS(A)      (((A)>=0)?(A):(-(A)))
#define ABS(A)       abs(A)
#define FLOATABS(A)  abs(A)
#define FLOATABSNAME abs
#define SQRT(A)      sqrt(A)
#define LOG(A)       log(A)
#define LOG10(A)     log10(A)
#define IREAL(A)     real(A)
#define ABSFNAME     abs
#define SQRTFNAME    sqrt
#define LOGFNAME     log
#define LOGTENFNAME  log10
#define CTOD(A)      (A)
#define CONJ(A)      (A)
#define FNULL        0.0e0
#define FONE         1.0e0
#define RONE         1.0e0
#define RTWO         2.0e0
#define RFOUR        4.0e0
#define RZERO        0.0e0
#define REAL_MAX     1.0e38
#define SIGNUM       1.0e0
#define SIGNUM2      1.0e0
#define TWICE        1




#elif defined _SINGLE_COMPLEX_
#define weightsC          CweightsC
#define add2is            Cadd2is
#define add2com           Cadd2com
#define CSRMERGEMATRICES  Cmergematrices

#define ILUPACKPARAMDELETE CILUPACKparamdelete
#define ILUCUPDATE        CILUCUPDAT
#define PILUCUPDATE       CPILUCUPDAT
#define SILUCUPDATE       CSILUCUPDAT
#define ILUCLIST          CILUCLIS
#define PILUCLIST         CPILUCLIS
#define BILUCLIST         CBILUCLIST
#define ILDLCUPDATE       CILDLCUPDAT
#define PILDLCUPDATE      CPILDLCUPDAT
#define SILDLCUPDATE      CSILDLCUPDAT

#define ILUPACKINIT       cgnlamginit
#define ILUPACKFACTOR     cgnlamgfactor
#define ILUPACKSOLVER     cgnlamgsolver
#define ILUPACKSOL        cgnlamgsol
#define ILUPACKDELETE     cgnlamgdelete
#define ILUPACKINFO       cgnlamginfo
#define ILUPACKNNZ        cgnlamgnnz
#define ILUPACKUNDOSCALING camgundoscaling

#define SPDILUPACKINIT    chpdamginit
#define SPDILUPACKFACTOR  chpdamgfactor
#define SPDILUPACKSOLVER  chpdamgsolver
#define SPDILUPACKSOL     chpdamgsol
#define SPDILUPACKDELETE  chpdamgdelete
#define SPDILUPACKINFO    chpdamginfo
#define SPDILUPACKNNZ     chpdamgnnz

#define SYMSPDILUPACKCONVERT      cherhpdamgconvert
#define SYMSPDCONVERT             CHERHPDAMGconvert

#define GETEPS       sgeteps
#define _REAL_MAX_   _S_REAL_MAX_
#define REALS        real
#define FLOAT        complex
#define CSRMAT       Cmat
#define RCSRMAT      Smat
#define PRO2SUM      Cpro2sum
#define INITARMS     Cinitarms
#define MPSARMS      Cmps_arms
#define AMGLEVELMAT  CAMGlevelmat
#define ILUPACKPARAM CILUPACKparam
#define GNLAGGRESSIVEDROPPING  CGNLaggressivedropping 
#define GNLFAGGRESSIVEDROPPING CGNLfaggressivedroppin
#define SPDAGGRESSIVEDROPPING  CHPDaggressivedropping 
#define SPDFAGGRESSIVEDROPPING CHPDfaggressivedroppin

#define SPTRF        csptrf
#define SPTRS        csptrs
#define PPTRF        cpptrf
#define PPTRS        cpptrs
#define GETRF        cgetrf
#define GETRS        cgetrs
#define GNLSYMMWM    CGNLsmwm

#define ROWSCALE     Crowscale
#define COLSCALE     Ccolscale
#define SPDSCALE     CHPDscale
#define SYMSCALE     CSYMscale
#define CNRMS        Ccnrms
#define RNRMS        Crnrms
#define DIAMUA       CDiagMultA
#define AMUDIA       CAMultDiag

#define CSRMATVEC    CGNLmatvec
#define CSRMATTVEC   CGNLmattvec
#define CSRMATHVEC   CGNLmathvec
#define SYMMATVEC    CSYMmatvec
#define SSMMATVEC    CSSMmatvec
#define HERMATVEC    CHERmatvec
#define SHRMATVEC    CSHRmatvec
#define SETUPGRAPH   Csetupgraph
#define SETUPGRAPH_EPSILON   Csetupgraph_epsilon
#define SETUPGRAPH_EPSILON_SP   Csetupgraph_epsilon_sp

#define AMGDELETE       CGNLAMGdelete
#define SPDAMGDELETE    CHPDAMGdelete

#define AMGSETUP        CGNLAMGfactor
#define SPDAMGSETUP     CHPDAMGfactor
#define AMGFACTOR       CGNLAMGfactor
#define SPDAMGFACTOR    CHPDAMGfactor
#define AMGEXTRACT      CGNLAMGextract
#define SPDAMGEXTRACT   CHERAMGextract
#define SYMAMGEXTRACT   CHERAMGextract

#define AMGSOL1         CGNLAMGsol1
#define AMGTSOL1        CGNLAMGtsol1
#define AMGSOL2         CGNLAMGsol2
#define AMGTSOL2        CGNLAMGtsol2
#define SPDAMGSOL1      CHPDAMGsol1
#define SPDAMGSOL2      CHPDAMGsol2
#define SYMAMGSOL1      CHERAMGsol1
#define SYMAMGSOL2      CHERAMGsol2

#define AMGSOL          CGNLAMGsol_internal
#define AMGSOL_EXT      CGNLAMGsol
#define AMGDLSOL        CGNLAMGdlsol
#define AMGDUSOL        CGNLAMGdusol
#define AMGLSOL         CGNLAMGlsol
#define AMGUSOL         CGNLAMGusol
#define AMGTDLSOL       CGNLAMGtdlsol
#define AMGTDUSOL       CGNLAMGtdusol
#define AMGTLSOL        CGNLAMGtlsol
#define AMGTUSOL        CGNLAMGtusol
#define AMGTSOL         CGNLAMGtsol_internal
#define AMGTSOL_EXT     CGNLAMGtsol
#define SPDAMGSOL       CHPDAMGsol_internal
#define SPDAMGSOL_EXT   CHPDAMGsol
#define AMGDRIVER       CGNLAMGdriver
#define AMGSOLVER       CGNLAMGsolver
#define SPDAMGSOLVER    CHPDAMGsolver
#define SYMAMGSOLVERS   CSYMAMGsolver
#define SYMAMGSOLVER    CHERAMGsolver
#define AMGINIT         CGNLAMGinit
#define AMGGETPARAMS    CGNLAMGgetparams
#define AMGSETPARAMS    CGNLAMGsetparams
#define AMGSETUPPARAMETERS    CGNLAMGsetupparameters
#define SPDAMGSETUPPARAMETERS    CHPDAMGsetupparameters
#define SPDAMGINIT      CHPDAMGinit
#define SPDAMGGETPARAMS CHPDAMGgetparams
#define SPDAMGSETPARAMS CHPDAMGsetparams
#define SYMAMGINITS      CSYMAMGinit
#define SYMAMGGETPARAMSS CSYMAMGgetparams
#define SYMAMGSETPARAMSS CSYMAMGsetparams
#define SYMAMGINIT      CHERAMGinit
#define SYMAMGGETPARAMS CHERAMGgetparams
#define SYMAMGSETPARAMS CHERAMGsetparams

#define LDLP         CHPDldlp
#define LDLPSOL      CHPDldlpsol
#define LUPQ         CGNLlupq
#define LUPQSOL      CGNLlupqsol
#define LUPQLSOL     CGNLlupqlsol
#define LUPQUSOL     CGNLlupqusol
#define LUPQDLSOL    CGNLlupqdlsol
#define LUPQDUSOL    CGNLlupqdusol
#define LUPQTSOL     CGNLlupqtsol
#define LUPQTLSOL    CGNLlupqtlsol
#define LUPQTUSOL    CGNLlupqtusol
#define LUPQTDLSOL   CGNLlupqtdlsol
#define LUPQTDUSOL   CGNLlupqtdusol



#define ILDLC        CHERildlc
#define ILDLCSOL     CHERildlcsol

#define PILDLCDLSOL  CHERpildlcdlsol
#define PILDLCDUSOL  CHERpildlcdusol
#define PILDLCLSOL   CHERpildlclsol
#define PILDLCUSOL   CHERpildlcusol

#define ILUC         CGNLiluc
#define PILUC        CGNLpiluc
#define SPILUC       CGNLspiluc
#define MPILUC       CGNLmpiluc
#define PILDLC       CHPDpiluc
#define MPILDLC      CHPDmpiluc
#define ILUCSOL      CGNLilucsol
#define ILUCTSOL     CGNLiluctsol
#define ILUCDLSOL    CGNLilucdlsol
#define ILUCTDLSOL   CGNLiluctdlsol
#define ILUCDUSOL    CGNLilucdusol
#define ILUCTDUSOL   CGNLiluctdusol
#define ILUCLSOL     CGNLiluclsol
#define ILUCTLSOL    CGNLiluctlsol
#define ILUCUSOL     CGNLilucusol
#define ILUCTUSOL    CGNLiluctusol

#define PILUCDLSOL    CGNLpilucdlsol
#define PILUCTDLSOL   CGNLpiluctdlsol
#define PILUCDUSOL    CGNLpilucdusol
#define PILUCTDUSOL   CGNLpiluctdusol
#define PILUCLSOL     CGNLpiluclsol
#define PILUCTLSOL    CGNLpiluctlsol
#define PILUCUSOL     CGNLpilucusol
#define PILUCTUSOL    CGNLpiluctusol

#define ILUTP        CGNLilutp
#define ILUT         CGNLilut
#define LUSOL        CGNLlusol
#define LUTSOL       CGNLlutsol
#define LULSOL       CGNLlulsol
#define LUTLSOL      CGNLlutlsol
#define LUUSOL       CGNLluusol
#define LUTUSOL      CGNLlutusol
#define LUDUSOL      CGNLludusol
#define LUTDUSOL     CGNLlutdusol
#define LUDLSOL      CGNLludlsol
#define LUTDLSOL     CGNLlutdlsol

#define CPERM        Ccperm
#define RPERM        Crperm

#define PERMNULL     CGNLperm_null
#define PERMND       CGNLperm_nd
#define PERMRCM      CGNLperm_rcm
#define PERMAMF      CGNLperm_amf
#define PERMAMD      CGNLperm_amd
#define PERMMMD      CGNLperm_mmd
#define PERMMETISE   CGNLperm_metis_e
#define PERMMETISN   CGNLperm_metis_n
#define PERMPQ       CGNLperm_pq

#define PERMNDFC       CGNLperm_nd_fc
#define PERMRCMFC      CGNLperm_rcm_fc
#define PERMAMFFC      CGNLperm_amf_fc
#define PERMAMDFC      CGNLperm_amd_fc
#define PERMMMDFC      CGNLperm_mmd_fc
#define PERMMETISEFC   CGNLperm_metis_e_fc
#define PERMMETISNFC   CGNLperm_metis_n_fc
#define PERMNDFCV      CGNLperm_nd_fcv
#define PERMRCMFCV     CGNLperm_rcm_fcv
#define PERMAMFFCV     CGNLperm_amf_fcv
#define PERMAMDFCV     CGNLperm_amd_fcv
#define PERMMMDFCV     CGNLperm_mmd_fcv
#define PERMMETISEFCV  CGNLperm_metis_e_fcv
#define PERMMETISNFCV  CGNLperm_metis_n_fcv
#define SYMPERMNDFC       CSYMperm_nd_fc
#define SYMPERMRCMFC      CSYMperm_rcm_fc
#define SYMPERMAMFFC      CSYMperm_amf_fc
#define SYMPERMAMDFC      CSYMperm_amd_fc
#define SYMPERMMMDFC      CSYMperm_mmd_fc
#define SYMPERMMETISEFC   CSYMperm_metis_e_fc
#define SYMPERMMETISNFC   CSYMperm_metis_n_fc
#define SYMPERMNDFCV      CSYMperm_nd_fcv
#define SYMPERMRCMFCV     CSYMperm_rcm_fcv
#define SYMPERMAMFFCV     CSYMperm_amf_fcv
#define SYMPERMAMDFCV     CSYMperm_amd_fcv
#define SYMPERMMMDFCV     CSYMperm_mmd_fcv
#define SYMPERMMETISEFCV  CSYMperm_metis_e_fcv
#define SYMPERMMETISNFCV  CSYMperm_metis_n_fcv


#define PERMFC       CGNLperm_fc
#define INDFC        Cindfc
#define INDFC_RS     Cindfc_rs
#define SYMPERMFC    CSYMperm_fc
#define SYMINDFC     CSYMindfc
#define SYMINDFC_RS  CSYMindfc_rs
#define SPDPERMFC    CHPDperm_fc
#define SPDINDFC     CHPDindfc
#define SPDINDFC_RS  CHPDindfc_rs

#define PERMFCV      CGNLperm_fcv
#define INDFCV       Cindfcv
#define INDFCV_RS    Cindfcv_rs
#define SYMPERMFCV   CSYMperm_fcv
#define SYMINDFCV    CSYMindfcv
#define SYMINDFCV_RS CSYMindfcv_rs
#define SPDPERMFCV   CHPDperm_fcv
#define SPDINDFCV    CHPDindfcv
#define SPDINDFCV_RS CHPDindfcv_rs

#define PERMP        CGNLperm_p
#define PERMINDSET   CGNLperm_indset

#define PERMMWMRCM    CGNLperm_mwm_rcm
#define PERMMWMMMD    CGNLperm_mwm_mmd
#define PERMMWMAMF    CGNLperm_mwm_amf
#define PERMMWMAMD    CGNLperm_mwm_amd
#define PERMMWMMETISE CGNLperm_mwm_metis_e
#define PERMMWMMETISN CGNLperm_mwm_metis_n

#define SYMPERMMWMRCM    CHERperm_mwm_rcm
#define SYMPERMMWMMMD    CHERperm_mwm_mmd
#define SYMPERMMWMAMF    CHERperm_mwm_amf
#define SYMPERMMWMAMD    CHERperm_mwm_amd
#define SYMPERMMWMMETISE CHERperm_mwm_metis_e
#define SYMPERMMWMMETISN CHERperm_mwm_metis_n

#define SYMPERMMATCHINGRCM    CHERperm_matching_rcm
#define SYMPERMMATCHINGMMD    CHERperm_matching_mmd
#define SYMPERMMATCHINGAMF    CHERperm_matching_amf
#define SYMPERMMATCHINGAMD    CHERperm_matching_amd
#define SYMPERMMATCHINGMETISE CHERperm_matching_metis_e
#define SYMPERMMATCHINGMETISN CHERperm_matching_metis_n

#define PERMMC64NULL   CGNLperm_mc64_null
#define PERMMC64RCM    CGNLperm_mc64_rcm
#define PERMMC64MMD    CGNLperm_mc64_mmd
#define PERMMC64AMF    CGNLperm_mc64_amf
#define PERMMC64AMD    CGNLperm_mc64_amd
#define PERMMC64METISE CGNLperm_mc64_metis_e
#define PERMMC64METISN CGNLperm_mc64_metis_n

#define PERMMATCHINGNULL   CGNLperm_matching_null
#define PERMMATCHINGRCM    CGNLperm_matching_rcm
#define PERMMATCHINGMMD    CGNLperm_matching_mmd
#define PERMMATCHINGAMF    CGNLperm_matching_amf
#define PERMMATCHINGAMD    CGNLperm_matching_amd
#define PERMMATCHINGMETISE CGNLperm_matching_metis_e
#define PERMMATCHINGMETISN CGNLperm_matching_metis_n

#define PERMMWARCM    CGNLperm_mwa_rcm
#define PERMMWAMMD    CGNLperm_mwa_mmd
#define PERMMWAAMF    CGNLperm_mwa_amf
#define PERMMWAAMD    CGNLperm_mwa_amd
#define PERMMWAMETISE CGNLperm_mwa_metis_e
#define PERMMWAMETISN CGNLperm_mwa_metis_n


#define SPDPERMNULL     CHPDperm_null
#define SPDPERMND       CHPDperm_nd
#define SPDPERMRCM      CHPDperm_rcm
#define SPDPERMAMF      CHPDperm_amf
#define SPDPERMAMD      CHPDperm_amd
#define SPDPERMMMD      CHPDperm_mmd
#define SPDPERMMETISE   CHPDperm_metis_e
#define SPDPERMMETISN   CHPDperm_metis_n
#define SPDPERMPP       CHPDperm_pp
#define SPDPERMINDSET   CHPDperm_indset


#define SPDPERMNDFC       CHPDperm_nd_fc
#define SPDPERMRCMFC      CHPDperm_rcm_fc
#define SPDPERMAMFFC      CHPDperm_amf_fc
#define SPDPERMAMDFC      CHPDperm_amd_fc
#define SPDPERMMMDFC      CHPDperm_mmd_fc
#define SPDPERMMETISEFC   CHPDperm_metis_e_fc
#define SPDPERMMETISNFC   CHPDperm_metis_n_fc

#define SPDPERMNDFCV      CHPDperm_nd_fcv
#define SPDPERMRCMFCV     CHPDperm_rcm_fcv
#define SPDPERMAMFFCV     CHPDperm_amf_fcv
#define SPDPERMAMDFCV     CHPDperm_amd_fcv
#define SPDPERMMMDFCV     CHPDperm_mmd_fcv
#define SPDPERMMETISEFCV  CHPDperm_metis_e_fcv
#define SPDPERMMETISNFCV  CHPDperm_metis_n_fcv


#define SYMPERMNULL     CHERperm_null
#define SYMPERMND       CHERperm_nd
#define SYMPERMRCM      CHERperm_rcm
#define SYMPERMAMF      CHERperm_amf
#define SYMPERMAMD      CHERperm_amd
#define SYMPERMMMD      CHERperm_mmd
#define SYMPERMMETISE   CHERperm_metis_e
#define SYMPERMMETISN   CHERperm_metis_n
#define SYMPERMINDSET   CHERperm_indset

#define SPARTRAN     Cspartran
#define ETREE        Cetree
#define QQSORT       Cqqsort
#define QSORT        Cqsort
#define QQSORT2      Cqqsort2
#define QSORT2       Cqsort2
#define QQSORTS      Cqqsorts
#define QQSORTS2     Cqqsorts2
#define BQSORT       Cbqsort
#define CLEAR        Cclear
#define QSPLIT       Cqsplit
#define QSPLIT2      Cqsplit2
#define QSPLIT3      Cqsplit3
#define BQSPLIT      Cbqsplit
#define QSEPARATE    Cqseparate
#define INDSET       Cindset
#define PQPERM       CPQpermF
#define PPERM        CPpermF
#define PPPERM       CPPpermF
#define WDIAG        CwDiagRF
#define WPDIAG       CwPDiagRF
#define SPDWDIAG     CHPDwDiagRF
#define QSORTR2I     sqsortr2i
#define SWAPM        sswapm
#define FSWAP        cswap
#define BISINIT      Cbisinit
#define STOPBIS      Cstopbis
#define TIDYCG       Ctidycg
#define MGSRO        Cmgsro
#define GIVENS       Cgivens
#define BRKDN        Cbrkdn

#define CSRCSC       Ccsrcsc
#define CSRCSC2      Ccsrcsc2
#define READMTC      Creadmtc
#define WRITEMTC     Cwritemtc
#define READVECTORS  Creadvectors
#define WRITEVECTORS Cwritevectors
#define PCG          Cpcg
#define FPCG         Cfpcg
#define BCG          Cbcg
#define SBCGS        CSYMbcg
#define SBCG         CHERbcg
#define SQMR         CHERqmr
#define SQMRS        CSYMqmr
#define GMRES        Cgmres
#define FGMRES       Cfgmres

#define DISTDOT      cdotc
#define DISTDOTU     cdotu
#define DISTDOT2     cdotc2
#define GEMM         cgemm
#define COPY         ccopy
#define AXPY         caxpy
#define ROTG         crotg
#define ROT          crot
#define NRM          scnrm2
#define I_AMAX       icamax
#define ASUM         scasum
#define SCAL         cscal

#define MC64I        mc64ic
#define MC64A        mc64ac
#define MC64IR       mc64is
#define MC64AR       mc64as

#define MUMPS_MATCH   cmumps_match
#define MUMPS_MATCHR  smumps_match

#define FABS(A)      sqrt((double)((A).r*(A).r+(A).i*(A).i))
#define ABS(A)       cabs(A)
#define FLOATABS(A)  abs(A)
#define FLOATABSNAME abs
#define SQRT(A)      sqrt(A)
#define LOG(A)       log(A)
#define LOG10(A)     log10(A)
#define IREAL(A)     real(A)
#define ABSFNAME     cabs
#define SQRTFNAME    sqrt
#define LOGFNAME     log
#define LOGTENFNAME  log10
#define CTOD(A)      real(A)
#define CONJ(A)      conjg(A)
#define FNULL        cmplx(0.0e0,0.0e0)
#define FONE         cmplx(1.0e0,0.0e0)
#define RONE         1.0e0
#define RTWO         2.0e0
#define RFOUR        4.0e0
#define RZERO        0.0e0
#define REAL_MAX     1.0e38
#define SIGNUM       signum
#define SIGNUM2      signum2
#define TWICE        2



#else
#define weightsC          ZweightsC
#define add2is            Zadd2is
#define add2com           Zadd2com
#define CSRMERGEMATRICES  Zmergematrices

#define ILUPACKPARAMDELETE ZILUPACKparamdelete
#define ILUCUPDATE        ZILUCUPDAT
#define PILUCUPDATE       ZPILUCUPDAT
#define SILUCUPDATE       ZSILUCUPDAT
#define ILUCLIST          ZILUCLIS
#define PILUCLIST         ZPILUCLIS
#define BILUCLIST         ZBILUCLIST
#define ILDLCUPDATE       ZILDLCUPDAT
#define PILDLCUPDATE      ZPILDLCUPDAT
#define SILDLCUPDATE      ZSILDLCUPDAT

#define ILUPACKINIT       zgnlamginit
#define ILUPACKFACTOR     zgnlamgfactor
#define ILUPACKSOLVER     zgnlamgsolver
#define ILUPACKSOL        zgnlamgsol
#define ILUPACKDELETE     zgnlamgdelete
#define ILUPACKINFO       zgnlamginfo
#define ILUPACKNNZ        zgnlamgnnz
#define ILUPACKUNDOSCALING zamgundoscaling

#define SPDILUPACKINIT    zhpdamginit
#define SPDILUPACKFACTOR  zhpdamgfactor
#define SPDILUPACKSOLVER  zhpdamgsolver
#define SPDILUPACKSOL     zhpdamgsol
#define SPDILUPACKDELETE  zhpdamgdelete
#define SPDILUPACKINFO    zhpdamginfo
#define SPDILUPACKNNZ     zhpdamgnnz

#define SYMSPDILUPACKCONVERT      zherhpdamgconvert
#define SYMSPDCONVERT             ZHERHPDAMGconvert


#define GETEPS       dgeteps
#define _REAL_MAX_   _D_REAL_MAX_
#define REALS        doubleprecision
#define FLOAT        ilu_doublecomplex
#define CSRMAT       Zmat
#define RCSRMAT      Dmat
#define PRO2SUM      Zpro2sum
#define INITARMS     Zinitarms
#define MPSARMS      Zmps_arms
#define AMGLEVELMAT  ZAMGlevelmat
#define ILUPACKPARAM ZILUPACKparam
#define GNLAGGRESSIVEDROPPING  ZGNLaggressivedropping 
#define GNLFAGGRESSIVEDROPPING ZGNLfaggressivedroppin
#define SPDAGGRESSIVEDROPPING  ZHPDaggressivedropping 
#define SPDFAGGRESSIVEDROPPING ZHPDfaggressivedroppin

#define SPTRF        zsptrf
#define SPTRS        zsptrs
#define PPTRF        zpptrf
#define PPTRS        zpptrs
#define GETRF        zgetrf
#define GETRS        zgetrs
#define GNLSYMMWM    ZGNLsmwm

#define ROWSCALE     Zrowscale
#define COLSCALE     Zcolscale
#define SPDSCALE     ZHPDscale
#define SYMSCALE     ZSYMscale
#define CNRMS        Zcnrms
#define RNRMS        Zrnrms
#define DIAMUA       ZDiagMultA
#define AMUDIA       ZAMultDiag

#define CSRMATVEC    ZGNLmatvec
#define CSRMATTVEC   ZGNLmattvec
#define CSRMATHVEC   ZGNLmathvec
#define SYMMATVEC    ZSYMmatvec
#define SSMMATVEC    ZSSMmatvec
#define HERMATVEC    ZHERmatvec
#define SHRMATVEC    ZSHRmatvec
#define SETUPGRAPH   Zsetupgraph
#define SETUPGRAPH_EPSILON   Zsetupgraph_epsilon
#define SETUPGRAPH_EPSILON_SP   Zsetupgraph_epsilon_sp

#define AMGDELETE       ZGNLAMGdelete
#define SPDAMGDELETE    ZHPDAMGdelete

#define AMGSETUP        ZGNLAMGfactor
#define SPDAMGSETUP     ZHPDAMGfactor
#define AMGFACTOR       ZGNLAMGfactor
#define SPDAMGFACTOR    ZHPDAMGfactor
#define AMGEXTRACT      ZGNLAMGextract
#define SPDAMGEXTRACT   ZHERAMGextract
#define SYMAMGEXTRACT   ZHERAMGextract

#define AMGSOL1         ZGNLAMGsol1
#define AMGTSOL1        ZGNLAMGtsol1
#define AMGSOL2         ZGNLAMGsol2
#define AMGTSOL2        ZGNLAMGtsol2
#define SPDAMGSOL1      ZHPDAMGsol1
#define SPDAMGSOL2      ZHPDAMGsol2
#define SYMAMGSOL1      ZHERAMGsol1
#define SYMAMGSOL2      ZHERAMGsol2

#define AMGSOL          ZGNLAMGsol_internal
#define AMGSOL_EXT      ZGNLAMGsol
#define AMGDLSOL        ZGNLAMGdlsol
#define AMGDUSOL        ZGNLAMGdusol
#define AMGLSOL         ZGNLAMGlsol
#define AMGUSOL         ZGNLAMGusol
#define AMGTDLSOL       ZGNLAMGtdlsol
#define AMGTDUSOL       ZGNLAMGtdusol
#define AMGTLSOL        ZGNLAMGtlsol
#define AMGTUSOL        ZGNLAMGtusol
#define AMGTSOL         ZGNLAMGtsol_internal
#define AMGTSOL_EXT     ZGNLAMGtsol
#define SPDAMGSOL       ZHPDAMGsol_internal
#define SPDAMGSOL_EXT   ZHPDAMGsol
#define AMGDRIVER       ZGNLAMGdriver
#define AMGSOLVER       ZGNLAMGsolver
#define SPDAMGSOLVER    ZHPDAMGsolver
#define SYMAMGSOLVERS   ZSYMAMGsolver
#define SYMAMGSOLVER    ZHERAMGsolver
#define AMGINIT         ZGNLAMGinit
#define AMGGETPARAMS    ZGNLAMGgetparams
#define AMGSETPARAMS    ZGNLAMGsetparams
#define AMGSETUPPARAMETERS    ZGNLAMGsetupparameters
#define SPDAMGSETUPPARAMETERS    ZHPDAMGsetupparameters
#define SPDAMGINIT      ZHPDAMGinit
#define SPDAMGGETPARAMS ZHPDAMGgetparams
#define SPDAMGSETPARAMS ZHPDAMGsetparams
#define SYMAMGINITS     ZSYMAMGinit
#define SYMAMGGETPARAMSS ZSYMAMGgetparams
#define SYMAMGSETPARAMSS ZSYMAMGsetparams
#define SYMAMGINIT      ZHERAMGinit
#define SYMAMGGETPARAMS ZHERAMGgetparams
#define SYMAMGSETPARAMS ZHERAMGsetparams

#define LDLP         ZHPDldlp
#define LDLPSOL      ZHPDldlpsol
#define LUPQ         ZGNLlupq
#define LUPQSOL      ZGNLlupqsol
#define LUPQLSOL     ZGNLlupqlsol
#define LUPQUSOL     ZGNLlupqusol
#define LUPQDLSOL    ZGNLlupqdlsol
#define LUPQDUSOL    ZGNLlupqdusol
#define LUPQTSOL     ZGNLlupqtsol
#define LUPQTLSOL    ZGNLlupqtlsol
#define LUPQTUSOL    ZGNLlupqtusol
#define LUPQTDLSOL   ZGNLlupqtdlsol
#define LUPQTDUSOL   ZGNLlupqtdusol


#define ILDLC        ZHERildlc
#define ILDLCSOL     ZHERildlcsol

#define PILDLCDLSOL  ZHERpildlcdlsol
#define PILDLCDUSOL  ZHERpildlcdusol
#define PILDLCLSOL   ZHERpildlclsol
#define PILDLCUSOL   ZHERpildlcusol

#define ILUC         ZGNLiluc
#define PILUC        ZGNLpiluc
#define SPILUC       ZGNLspiluc
#define MPILUC       ZGNLmpiluc
#define PILDLC       ZHPDpiluc
#define MPILDLC      ZHPDmpiluc
#define ILUCSOL      ZGNLilucsol
#define ILUCTSOL     ZGNLiluctsol
#define ILUCDLSOL    ZGNLilucdlsol
#define ILUCTDLSOL   ZGNLiluctdlsol
#define ILUCDUSOL    ZGNLilucdusol
#define ILUCTDUSOL   ZGNLiluctdusol
#define ILUCLSOL     ZGNLiluclsol
#define ILUCTLSOL    ZGNLiluctlsol
#define ILUCUSOL     ZGNLilucusol
#define ILUCTUSOL    ZGNLiluctusol

#define PILUCDLSOL    ZGNLpilucdlsol
#define PILUCTDLSOL   ZGNLpiluctdlsol
#define PILUCDUSOL    ZGNLpilucdusol
#define PILUCTDUSOL   ZGNLpiluctdusol
#define PILUCLSOL     ZGNLpiluclsol
#define PILUCTLSOL    ZGNLpiluctlsol
#define PILUCUSOL     ZGNLpilucusol
#define PILUCTUSOL    ZGNLpiluctusol

#define ILUTP        ZGNLilutp
#define ILUT         ZGNLilut
#define LUSOL        ZGNLlusol
#define LUTSOL       ZGNLlutsol
#define LUUSOL       ZGNLluusol
#define LUTUSOL      ZGNLlutusol
#define LULSOL       ZGNLlulsol
#define LUTLSOL      ZGNLlutlsol
#define LUDUSOL      ZGNLludusol
#define LUTDUSOL     ZGNLlutdusol
#define LUDLSOL      ZGNLludlsol
#define LUTDLSOL     ZGNLlutdlsol

#define CPERM        Zcperm
#define RPERM        Zrperm

#define PERMNULL     ZGNLperm_null
#define PERMND       ZGNLperm_nd
#define PERMRCM      ZGNLperm_rcm
#define PERMAMF      ZGNLperm_amf
#define PERMAMD      ZGNLperm_amd
#define PERMMMD      ZGNLperm_mmd
#define PERMMETISE   ZGNLperm_metis_e
#define PERMMETISN   ZGNLperm_metis_n
#define PERMPQ       ZGNLperm_pq

#define PERMNDFC       ZGNLperm_nd_fc
#define PERMRCMFC      ZGNLperm_rcm_fc
#define PERMAMFFC      ZGNLperm_amf_fc
#define PERMAMDFC      ZGNLperm_amd_fc
#define PERMMMDFC      ZGNLperm_mmd_fc
#define PERMMETISEFC   ZGNLperm_metis_e_fc
#define PERMMETISNFC   ZGNLperm_metis_n_fc
#define PERMNDFCV      ZGNLperm_nd_fcv
#define PERMRCMFCV     ZGNLperm_rcm_fcv
#define PERMAMFFCV     ZGNLperm_amf_fcv
#define PERMAMDFCV     ZGNLperm_amd_fcv
#define PERMMMDFCV     ZGNLperm_mmd_fcv
#define PERMMETISEFCV  ZGNLperm_metis_e_fcv
#define PERMMETISNFCV  ZGNLperm_metis_n_fcv
#define SYMPERMNDFC       ZSYMperm_nd_fc
#define SYMPERMRCMFC      ZSYMperm_rcm_fc
#define SYMPERMAMFFC      ZSYMperm_amf_fc
#define SYMPERMAMDFC      ZSYMperm_amd_fc
#define SYMPERMMMDFC      ZSYMperm_mmd_fc
#define SYMPERMMETISEFC   ZSYMperm_metis_e_fc
#define SYMPERMMETISNFC   ZSYMperm_metis_n_fc
#define SYMPERMNDFCV      ZSYMperm_nd_fcv
#define SYMPERMRCMFCV     ZSYMperm_rcm_fcv
#define SYMPERMAMFFCV     ZSYMperm_amf_fcv
#define SYMPERMAMDFCV     ZSYMperm_amd_fcv
#define SYMPERMMMDFCV     ZSYMperm_mmd_fcv
#define SYMPERMMETISEFCV  ZSYMperm_metis_e_fcv
#define SYMPERMMETISNFCV  ZSYMperm_metis_n_fcv

#define PERMFC       ZGNLperm_fc
#define INDFC        Zindfc
#define INDFC_RS     Zindfc_rs
#define SYMPERMFC    ZSYMperm_fc
#define SYMINDFC     ZSYMindfc
#define SYMINDFC_RS  ZSYMindfc_rs
#define SPDPERMFC    ZHPDperm_fc
#define SPDINDFC     ZHPDindfc
#define SPDINDFC_RS  ZHPDindfc_rs

#define PERMFCV      ZGNLperm_fcv
#define INDFCV       Zindfcv
#define INDFCV_RS    Zindfcv_rs
#define SYMPERMFCV   ZSYMperm_fcv
#define SYMINDFCV    ZSYMindfcv
#define SYMINDFCV_RS ZSYMindfcv_rs
#define SPDPERMFCV   ZHPDperm_fcv
#define SPDINDFCV_RS ZHPDindfcv_rs

#define PERMP        ZGNLperm_p
#define PERMINDSET   ZGNLperm_indset

#define PERMMWMRCM    ZGNLperm_mwm_rcm
#define PERMMWMMMD    ZGNLperm_mwm_mmd
#define PERMMWMAMF    ZGNLperm_mwm_amf
#define PERMMWMAMD    ZGNLperm_mwm_amd
#define PERMMWMMETISE ZGNLperm_mwm_metis_e
#define PERMMWMMETISN ZGNLperm_mwm_metis_n

#define SYMPERMMWMRCM    ZHERperm_mwm_rcm
#define SYMPERMMWMMMD    ZHERperm_mwm_mmd
#define SYMPERMMWMAMF    ZHERperm_mwm_amf
#define SYMPERMMWMAMD    ZHERperm_mwm_amd
#define SYMPERMMWMMETISE ZHERperm_mwm_metis_e
#define SYMPERMMWMMETISN ZHERperm_mwm_metis_n

#define SYMPERMMATCHINGRCM    ZHERperm_matching_rcm
#define SYMPERMMATCHINGMMD    ZHERperm_matching_mmd
#define SYMPERMMATCHINGAMF    ZHERperm_matching_amf
#define SYMPERMMATCHINGAMD    ZHERperm_matching_amd
#define SYMPERMMATCHINGMETISE ZHERperm_matching_metis_e
#define SYMPERMMATCHINGMETISN ZHERperm_matching_metis_n

#define PERMMC64NULL   ZGNLperm_mc64_null
#define PERMMC64RCM    ZGNLperm_mc64_rcm
#define PERMMC64MMD    ZGNLperm_mc64_mmd
#define PERMMC64AMF    ZGNLperm_mc64_amf
#define PERMMC64AMD    ZGNLperm_mc64_amd
#define PERMMC64METISE ZGNLperm_mc64_metis_e
#define PERMMC64METISN ZGNLperm_mc64_metis_n

#define PERMMATCHINGNULL   ZGNLperm_matching_null
#define PERMMATCHINGRCM    ZGNLperm_matching_rcm
#define PERMMATCHINGMMD    ZGNLperm_matching_mmd
#define PERMMATCHINGAMF    ZGNLperm_matching_amf
#define PERMMATCHINGAMD    ZGNLperm_matching_amd
#define PERMMATCHINGMETISE ZGNLperm_matching_metis_e
#define PERMMATCHINGMETISN ZGNLperm_matching_metis_n

#define PERMMWARCM    ZGNLperm_mwa_rcm
#define PERMMWAMMD    ZGNLperm_mwa_mmd
#define PERMMWAAMF    ZGNLperm_mwa_amf
#define PERMMWAAMD    ZGNLperm_mwa_amd
#define PERMMWAMETISE ZGNLperm_mwa_metis_e
#define PERMMWAMETISN ZGNLperm_mwa_metis_n


#define SPDPERMNULL     ZHPDperm_null
#define SPDPERMND       ZHPDperm_nd
#define SPDPERMRCM      ZHPDperm_rcm
#define SPDPERMAMF      ZHPDperm_amf
#define SPDPERMAMD      ZHPDperm_amd
#define SPDPERMMMD      ZHPDperm_mmd
#define SPDPERMMETISE   ZHPDperm_metis_e
#define SPDPERMMETISN   ZHPDperm_metis_n
#define SPDPERMPP       ZHPDperm_pp
#define SPDPERMINDSET   ZHPDperm_indset


#define SPDPERMNDFC       ZHPDperm_nd_fc
#define SPDPERMRCMFC      ZHPDperm_rcm_fc
#define SPDPERMAMFFC      ZHPDperm_amf_fc
#define SPDPERMAMDFC      ZHPDperm_amd_fc
#define SPDPERMMMDFC      ZHPDperm_mmd_fc
#define SPDPERMMETISEFC   ZHPDperm_metis_e_fc
#define SPDPERMMETISNFC   ZHPDperm_metis_n_fc

#define SPDPERMNDFCV      ZHPDperm_nd_fcv
#define SPDPERMRCMFCV     ZHPDperm_rcm_fcv
#define SPDPERMAMFFCV     ZHPDperm_amf_fcv
#define SPDPERMAMDFCV     ZHPDperm_amd_fcv
#define SPDPERMMMDFCV     ZHPDperm_mmd_fcv
#define SPDPERMMETISEFCV  ZHPDperm_metis_e_fcv
#define SPDPERMMETISNFCV  ZHPDperm_metis_n_fcv


#define SYMPERMNULL     ZHERperm_null
#define SYMPERMND       ZHERperm_nd
#define SYMPERMRCM      ZHERperm_rcm
#define SYMPERMAMF      ZHERperm_amf
#define SYMPERMAMD      ZHERperm_amd
#define SYMPERMMMD      ZHERperm_mmd
#define SYMPERMMETISE   ZHERperm_metis_e
#define SYMPERMMETISN   ZHERperm_metis_n
#define SYMPERMPP       ZHERperm_pp
#define SYMPERMINDSET   ZHERperm_indset

#define SPARTRAN     Zspartran
#define ETREE        Zetree
#define QQSORT       Zqqsort
#define QSORT        Zqsort
#define QQSORT2      Zqqsort2
#define QSORT2       Zqsort2
#define QQSORTS      Zqqsorts
#define QQSORTS2     Zqqsorts2
#define BQSORT       Zbqsort
#define CLEAR        Zclear
#define QSPLIT       Zqsplit
#define QSPLIT2      Zqsplit2
#define QSPLIT3      Zqsplit3
#define BQSPLIT      Zbqsplit
#define QSEPARATE    Zqseparate
#define INDSET       Zindset
#define PQPERM       ZPQpermF
#define PPERM        ZPpermF
#define PPPERM       ZPPpermF
#define WDIAG        ZwDiagRF
#define WPDIAG       ZwPDiagRF
#define SPDWDIAG     ZHPDwDiagRF
#define QSORTR2I     dqsortr2i
#define SWAPM        dswapm
#define FSWAP        zswap
#define BISINIT      Zbisinit
#define STOPBIS      Zstopbis
#define TIDYCG       Ztidycg
#define MGSRO        Zmgsro
#define GIVENS       Zgivens
#define BRKDN        Zbrkdn

#define CSRCSC       Zcsrcsc
#define CSRCSC2      Zcsrcsc2
#define READMTC      Zreadmtc
#define WRITEMTC     Zwritemtc
#define READVECTORS  Zreadvectors
#define WRITEVECTORS Zwritevectors
#define PCG          Zpcg
#define FPCG         Zfpcg
#define BCG          Zbcg
#define SBCGS        ZSYMbcg
#define SBCG         ZHERbcg
#define SQMR         ZHERqmr
#define SQMRS        ZSYMqmr
#define GMRES        Zgmres
#define FGMRES       Zfgmres

#define DISTDOT      zdotc
#define DISTDOTU     zdotu
#define DISTDOT2     zdotc2
#define GEMM         zgemm
#define COPY         zcopy
#define AXPY         zaxpy
#define ROTG         zrotg
#define ROT          zrot
#define NRM          dznrm2
#define I_AMAX       izamax
#define ASUM         dzasum
#define SCAL         zscal

#define MC64I        mc64iz
#define MC64A        mc64az
#define MC64IR       mc64id
#define MC64AR       mc64ad

#define MUMPS_MATCH   zmumps_match
#define MUMPS_MATCHR  dmumps_match

#define FABS(A)      sqrt((double)((A).r*(A).r+(A).i*(A).i))
#define ABS(A)       cdabs(A)
#define FLOATABS(A)  dabs(A)
#define FLOATABSNAME dabs
#define SQRT(A)      dsqrt(A)
#define LOG(A)       dlog(A)
#define LOG10(A)     dlog10(A)
#define IREAL(A)     dble(A)
#define ABSFNAME     cdabs
#define SQRTFNAME    dsqrt
#define LOGFNAME     dlog
#define LOGTENFNAME  dlog10
#define CTOD(A)      dreal(A)
#define CONJ(A)      dconjg(A)
#define FNULL        dcmplx(0.0d0,0.0d0)
#define FONE         dcmplx(1.0d0,0.0d0)
#define RONE         1.0d0
#define RTWO         2.0d0
#define RFOUR        4.0d0
#define RZERO        0.0d0
#define REAL_MAX     1.0d300
#define SIGNUM       signum
#define SIGNUM2      signum2
#define TWICE        2
#endif


#define MALLOC        MAlloc
#define REALLOC       ReAlloc
#define FREE          FRee
#define MC64THRESHOLD 1.0e-2
#define QQSORTI       qqsorti




#endif
