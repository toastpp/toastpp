#ifndef __BLASNAMES_H
#define __BLASNAMES_H

#if (!defined(_WIN32))&&(!defined(_WIN64))
#define TOASTLOCAL_BLAS_NAMES
#endif

#ifdef TOASTLOCAL_BLAS_NAMES
#define sasum_  toast_sasum_
#define isamax_ toast_isamax_
#define scopy_  toast_scopy_
#define sscal_  toast_sscal_
#define sger_   toast_sger_
#define snrm2_  toast_snrm2_
#define ssymv_  toast_ssymv_
#define sdot_   toast_sdot_
#define saxpy_  toast_saxpy_
#define ssyr2_  toast_ssyr2_
#define srot_   toast_srot_
#define sgemv_  toast_sgemv_
#define strsv_  toast_strsv_
#define sgemm_  toast_sgemm_
#define strsm_  toast_strsm_
#define ssyrk_  toast_ssyrk_

#define dasum_  toast_dasum_
#define idamax_ toast_idamax_
#define dcopy_  toast_dcopy_
#define dscal_  toast_dscal_
#define dger_   toast_dger_
#define dnrm2_  toast_dnrm2_
#define dsymv_  toast_dsymv_
#define ddot_   toast_ddot_
#define daxpy_  toast_daxpy_
#define dsyr2_  toast_dsyr2_
#define drot_   toast_drot_
#define dgemv_  toast_dgemv_
#define dtrsv_  toast_dtrsv_
#define dgemm_  toast_dgemm_
#define dtrsm_  toast_dtrsm_
#define dsyrk_  toast_dsyrk_

#define scasum_ toast_scasum_
#define icamax_ toast_icamax_
#define ccopy_  toast_ccopy_
#define cscal_  toast_cscal_
#define scnrm2_ toast_scnrm2_
#define caxpy_  toast_caxpy_
#define cgemv_  toast_cgemv_
#define ctrsv_  toast_ctrsv_
#define cgemm_  toast_cgemm_
#define ctrsm_  toast_ctrsm_
#define cgerc_  toast_cgerc_
#define chemv_  toast_chemv_
#define cher2_  toast_cher2_
#define csyrk_  toast_csyrk_

#define dzasum_ toast_dzasum_
#define izamax_ toast_izamax_
#define zcopy_  toast_zcopy_
#define zscal_  toast_zscal_
#define dznrm2_ toast_dznrm2_
#define zdotc_  toast_zdotc_
#define zdotu_  toast_zdotu_
#define zaxpy_  toast_zaxpy_
#define zgemv_  toast_zgemv_
#define ztrsv_  toast_ztrsv_
#define zgemm_  toast_zgemm_
#define ztrsm_  toast_ztrsm_
#define zgerc_  toast_zgerc_
#define zhemv_  toast_zhemv_
#define zher2_  toast_zher2_
#define zsyrk_  toast_zsyrk_
#endif // TOASTLOCAL_BLAS_NAMES

#endif // !__BLASNAMES_H
