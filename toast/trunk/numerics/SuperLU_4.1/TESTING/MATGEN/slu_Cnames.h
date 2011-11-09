file slu_Cnames.h
 * \brief Macros defining how C routines will be called
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 1, 1997
 *
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this). 
 * </pre>
 */
#ifndef __SUPERLU_CNAMES /* allow multiple inclusions */
#define __SUPERLU_CNAMES



// Force use of TOAST version of BLAS libraries (32-bit integer)
#ifndef WIN32
#define TOAST_BLAS_32INT
#endif


#define ADD_       0
#define ADD__      1
#define NOCHANGE   2
#define UPCASE     3
#define OLD_CRAY   4
#define C_CALL     5

#ifdef UpCase
#define F77_CALL_C  UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C  NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C  ADD_
#endif

#ifdef Add__
#define F77_CALL_C  ADD__
#endif

#ifdef _CRAY
#define F77_CALL_C  OLD_CRAY
#endif

/* Default */
#ifndef F77_CALL_C
#define F77_CALL_C  ADD_
#endif

#ifdef TOAST_BLAS_32INT

/* BLAS */
#define sswap_    toast_sswap_
#define saxpy_    toast_saxpy_
#define sasum_    toast_sasum_
#define isamax_   toast_isamax_
#define scopy_    toast_scopy_
#define sscal_    toast_sscal_
#define sger_     toast_sger_
#define snrm2_    toast_snrm2_
#define ssymv_    toast_ssymv_
#define sdot_     toast_sdot_
#define saxpy_    toast_saxpy_
#define ssyr2_    toast_ssyr2_
#define srot_     toast_srot_
#define sgemv_    toast_sgemv_
#define strsv_    toast_strsv_
#define sgemm_    toast_sgemm_
#define strsm_    toast_strsm_

#define dswap_    toast_dswap_
#define daxpy_    toast_daxpy_
#define dasum_    toast_dasum_
#define idamax_   toast_idamax_
#define dcopy_    toast_dcopy_
#define dscal_    toast_dscal_
#define dger_     toast_dger_
#define dnrm2_    toast_dnrm2_
#define dsymv_    toast_dsymv_
#define ddot_     toast_ddot_
#define dsyr2_    toast_dsyr2_
#define drot_     toast_drot_
#define dgemv_    toast_dgemv_
#define dtrsv_    toast_dtrsv_
#define dgemm_    toast_dgemm_
#define dtrsm_    toast_dtrsm_

#define cswap_    toast_cswap_
#define caxpy_    toast_caxpy_
#define scasum_   toast_scasum_
#define icamax_   toast_icamax_
#define ccopy_    toast_ccopy_
#define cscal_    toast_cscal_
#define scnrm2_   toast_scnrm2_
#define caxpy_    toast_caxpy_
#define cgemv_    toast_cgemv_
#define ctrsv_    toast_ctrsv_
#define cgemm_    toast_cgemm_
#define ctrsm_    toast_ctrsm_
#define cgerc_    toast_cgerc_
#define chemv_    toast_chemv_
#define cher2_    toast_cher2_

#define zswap_    toast_zswap_
#define zaxpy_    toast_zaxpy_
#define dzasum_   toast_dzasum_
#define izamax_   toast_izamax_
#define zcopy_    toast_zcopy_
#define zscal_    toast_zscal_
#define dznrm2_   toast_dznrm2_
#define zaxpy_    toast_zaxpy_
#define zgemv_    toast_zgemv_
#define ztrsv_    toast_ztrsv_
#define zgemm_    toast_zgemm_
#define ztrsm_    toast_ztrsm_
#define zgerc_    toast_zgerc_
#define zhemv_    toast_zhemv_
#define zher2_    toast_zher2_

/* LAPACK */
#define dlamch_   toast_dlamch_
#define slamch_   toast_slamch_
#define xerbla_   toast_xerbla_
#define lsame_    toast_lsame_
#define dlacon_   toast_dlacon_
#define slacon_   toast_slacon_
#define icmax1_   toast_icmax1_
#define scsum1_   toast_scsum1_
#define clacon_   toast_clacon_
#define dzsum1_   toast_dzsum1_
#define izmax1_   toast_izmax1_
#define zlacon_   toast_zlacon_

#else // !TOAST_BLAS_32INT

#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm_(...)
 *
 * This is the default.
 */

#endif

#if (F77_CALL_C == ADD__)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm__(...)
 */
/* BLAS */
#define sswap_    sswap__
#define saxpy_    saxpy__
#define sasum_    sasum__
#define isamax_   isamax__
#define scopy_    scopy__
#define sscal_    sscal__
#define sger_     sger__
#define snrm2_    snrm2__
#define ssymv_    ssymv__
#define sdot_     sdot__
#define saxpy_    saxpy__
#define ssyr2_    ssyr2__
#define srot_     srot__
#define sgemv_    sgemv__
#define strsv_    strsv__
#define sgemm_    sgemm__
#define strsm_    strsm__

#define dswap_    dswap__
#define daxpy_    daxpy__
#define dasum_    dasum__
#define idamax_   idamax__
#define dcopy_    dcopy__
#define dscal_    dscal__
#define dger_     dger__
#define dnrm2_    dnrm2__
#define dsymv_    dsymv__
#define ddot_     ddot__
#define dsyr2_    dsyr2__
#define drot_     drot__
#define dgemv_    dgemv__
#define dtrsv_    dtrsv__
#define dgemm_    dgemm__
#define dtrsm_    dtrsm__

#define cswap_    cswap__
#define caxpy_    caxpy__
#define scasum_   scasum__
#define icamax_   icamax__
#define ccopy_    ccopy__
#define cscal_    cscal__
#define scnrm2_   scnrm2__
#define caxpy_    caxpy__
#define cgemv_    cgemv__
#define ctrsv_    ctrsv__
#define cgemm_    cgemm__
#define ctrsm_    ctrsm__
#define cgerc_    cgerc__
#define chemv_    chemv__
#define cher2_    cher2__

#define zswap_    zswap__
#define zaxpy_    zaxpy__
#define dzasum_   dzasum__
#define izamax_   izamax__
#define zcopy_    zcopy__
#define zscal_    zscal__
#define dznrm2_   dznrm2__
#define zaxpy_    zaxpy__
#define zgemv_    zgemv__
#define ztrsv_    ztrsv__
#define zgemm_    zgemm__
#define ztrsm_    ztrsm__
#define zgerc_    zgerc__
#define zhemv_    zhemv__
#define zher2_    zher2__

/* LAPACK */
#define dlamch_   dlamch__
#define slamch_   slamch__
#define xerbla_   xerbla__
#define lsame_    lsame__
#define dlacon_   dlacon__
#define slacon_   slacon__
#define icmax1_   icmax1__
#define scsum1_   scsum1__
#define clacon_   clacon__
#define dzsum1_   dzsum1__
#define izmax1_   izmax1__
#define zlacon_   zlacon__

/* Fortran interface */
#define c_bridge_dgssv_ c_bridge_dgssv__
#define c_fortran_sgssv_ c_fortran_sgssv__
#define c_fortran_dgssv_ c_fortran_dgssv__
#define c_fortran_cgssv_ c_fortran_cgssv__
#define c_fortran_zgssv_ c_fortran_zgssv__
#endif

#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void DGEMM(...)
 */
/* BLAS */
#define sswap_    SSWAP
#define saxpy_    SAXPY
#define sasum_    SASUM
#define isamax_   ISAMAX
#define scopy_    SCOPY
#define sscal_    SSCAL
#define sger_     SGER
#define snrm2_    SNRM2
#define ssymv_    SSYMV
#define sdot_     SDOT
#define saxpy_    SAXPY
#define ssyr2_    SSYR2
#define srot_     SROT
#define sgemv_    SGEMV
#define strsv_    STRSV
#define sgemm_    SGEMM
#define strsm_    STRSM

#define dswap_    DSWAP
#define daxpy_    DAXPY
#define dasum_    DASUM
#define idamax_   IDAMAX
#define dcopy_    DCOPY
#define dscal_    DSCAL
#define dger_     DGER
#define dnrm2_    DNRM2
#define dsymv_    DSYMV
#define ddot_     DDOT
#define dsyr2_    DSYR2
#define drot_     DROT
#define dgemv_    DGEMV
#define dtrsv_    DTRSV
#define dgemm_    DGEMM
#define dtrsm_    DTRSM

#define cswap_    CSWAP
#define caxpy_    CAXPY
#define scasum_   SCASUM
#define icamax_   ICAMAX
#define ccopy_    CCOPY
#define cscal_    CSCAL
#define scnrm2_   SCNRM2
#define cgemv_    CGEMV
#define ctrsv_    CTRSV
#define cgemm_    CGEMM
#define ctrsm_    CTRSM
#define cgerc_    CGERC
#define chemv_    CHEMV
#define cher2_    CHER2

#define zswap_    ZSWAP
#define zaxpy_    ZAXPY
#define dzasum_   DZASUM
#define izamax_   IZAMAX
#define zcopy_    ZCOPY
#define zscal_    ZSCAL
#define dznrm2_   DZNRM2
#define zgemv_    ZGEMV
#define ztrsv_    ZTRSV
#define zgemm_    ZGEMM
#define ztrsm_    ZTRSM
#define zgerc_    ZGERC
#define zhemv_    ZHEMV
#define zher2_    ZHER2

/* LAPACK */
#define dlamch_   DLAMCH
#define slamch_   SLAMCH
#define xerbla_   XERBLA
#define lsame_    LSAME
#define dlacon_   DLACON
#define slacon_   SLACON
#define icmax1_   ICMAX1
#define scsum1_   SCSUM1
#define clacon_   CLACON
#define dzsum1_   DZSUM1
#define izmax1_   IZMAX1
#define zlacon_   ZLACON

/* Fortran interface */
#define c_bridge_dgssv_ C_BRIDGE_DGSSV
#define c_fortran_sgssv_ C_FORTRAN_SGSSV
#define c_fortran_dgssv_ C_FORTRAN_DGSSV
#define c_fortran_cgssv_ C_FORTRAN_CGSSV
#define c_fortran_zgssv_ C_FORTRAN_ZGSSV
#endif


#if (F77_CALL_C == OLD_CRAY)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void SGEMM(...)
 */
/* BLAS */
#define sswap_    SSWAP
#define saxpy_    SAXPY
#define sasum_    SASUM
#define isamax_   ISAMAX
#define scopy_    SCOPY
#define sscal_    SSCAL
#define sger_     SGER
#define snrm2_    SNRM2
#define ssymv_    SSYMV
#define sdot_     SDOT
#define ssyr2_    SSYR2
#define srot_     SROT
#define sgemv_    SGEMV
#define strsv_    STRSV
#define sgemm_    SGEMM
#define strsm_    STRSM

#define dswap_    SSWAP
#define daxpy_    SAXPY
#define dasum_    SASUM
#define idamax_   ISAMAX
#define dcopy_    SCOPY
#define dscal_    SSCAL
#define dger_     SGER
#define dnrm2_    SNRM2
#define dsymv_    SSYMV
#define ddot_     SDOT
#define dsyr2_    SSYR2
#define drot_     SROT
#define dgemv_    SGEMV
#define dtrsv_    STRSV
#define dgemm_    SGEMM
#define dtrsm_    STRSM

#define cswap_    CSWAP
#define caxpy_    CAXPY
#define scasum_   SCASUM
#define icamax_   ICAMAX
#define ccopy_    CCOPY
#define cscal_    CSCAL
#define scnrm2_   SCNRM2
#define caxpy_    CAXPY
#define cgemv_    CGEMV
#define ctrsv_    CTRSV
#define cgemm_    CGEMM
#define ctrsm_    CTRSM
#define cgerc_    CGERC
#define chemv_    CHEMV
#define cher2_    CHER2

#define zswap_    ZSWAP
#define zaxpy_    ZAXPY
#define dzasum_   DZASUM
#define izamax_   IZAMAX
#define zcopy_    ZCOPY
#define zscal_    ZSCAL
#define dznrm2_   DZNRM2
#define zgemv_    ZGEMV
#define ztrsv_    ZTRSV
#define zgemm_    ZGEMM
#define ztrsm_    ZTRSM
#define zgerc_    ZGERC
#define zhemv_    ZHEMV
#define zher2_    ZHER2

/* LAPACK */
#define dlamch_   DLAMCH
#define slamch_   SLAMCH
#define xerbla_   XERBLA
#define lsame_    LSAME
#define dlacon_   DLACON
#define slacon_   SLACON
#define icmax1_   ICMAX1
#define scsum1_   SCSUM1
#define clacon_   CLACON
#define dzsum1_   DZSUM1
#define izmax1_   IZMAX1
#define zlacon_   ZLACON

/* Fortran interface */
#define c_bridge_dgssv_ C_BRIDGE_DGSSV
#define c_fortran_sgssv_ C_FORTRAN_SGSSV
#define c_fortran_dgssv_ C_FORTRAN_DGSSV
#define c_fortran_cgssv_ C_FORTRAN_CGSSV
#define c_fortran_zgssv_ C_FORTRAN_ZGSSV
#endif


#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm(...)
 */
/* BLAS */
#define sswap_    sswap
#define saxpy_    saxpy
#define sasum_    sasum
#define isamax_   isamax
#define scopy_    scopy
#define sscal_    sscal
#define sger_     sger
#define snrm2_    snrm2
#define ssymv_    ssymv
#define sdot_     sdot
#define saxpy_    saxpy
#define ssyr2_    ssyr2
#define srot_     srot
#define sgemv_    sgemv
#define strsv_    strsv
#define sgemm_    sgemm
#define strsm_    strsm

#define dswap_    dswap
#define daxpy_    daxpy
#define dasum_    dasum
#define idamax_   idamax
#define dcopy_    dcopy
#define dscal_    dscal
#define dger_     dger
#define dnrm2_    dnrm2
#define dsymv_    dsymv
#define ddot_     ddot
#define dsyr2_    dsyr2
#define drot_     drot
#define dgemv_    dgemv
#define dtrsv_    dtrsv
#define dgemm_    dgemm
#define dtrsm_    dtrsm

#define cswap_    cswap
#define caxpy_    caxpy
#define scasum_   scasum
#define icamax_   icamax
#define ccopy_    ccopy
#define cscal_    cscal
#define scnrm2_   scnrm2
#define cgemv_    cgemv
#define ctrsv_    ctrsv
#define cgemm_    cgemm
#define ctrsm_    ctrsm
#define cgerc_    cgerc
#define chemv_    chemv
#define cher2_    cher2

#define zswap_    zswap
#define zaxpy_    zaxpy
#define dzasum_   dzasum
#define izamax_   izamax
#define zcopy_    zcopy
#define zscal_    zscal
#define dznrm2_   dznrm2
#define zgemv_    zgemv
#define ztrsv_    ztrsv
#define zgemm_    zgemm
#define ztrsm_    ztrsm
#define zgerc_    zgerc
#define zhemv_    zhemv
#define zher2_    zher2

/* LAPACK */
#define dlamch_   dlamch
#define slamch_   slamch
#define xerbla_   xerbla
#define lsame_    lsame
#define dlacon_   dlacon
#define slacon_   slacon
#define icmax1_   icmax1
#define scsum1_   scsum1
#define clacon_   clacon
#define dzsum1_   dzsum1
#define izmax1_   izmax1
#define zlacon_   zlacon

/* Fortran interface */
#define c_bridge_dgssv_ c_bridge_dgssv
#define c_fortran_sgssv_ c_fortran_sgssv
#define c_fortran_dgssv_ c_fortran_dgssv
#define c_fortran_cgssv_ c_fortran_cgssv
#define c_fortran_zgssv_ c_fortran_zgssv
#endif

#endif // !TOAST_BLAS_32INT

#endif /* __SUPERLU_CNAMES */
