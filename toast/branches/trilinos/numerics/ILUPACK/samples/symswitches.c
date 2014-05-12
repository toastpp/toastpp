#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
#define CONJG(A)      (A)

#ifdef _DOUBLE_REAL_
#define MYSYMPILUCSOL        DSYMpilucsol
#define MYSYMILUC            DSYMiluc
#define SYMILUPACK           dsymilupack
#define SYMILUPACKFAC        dsymilupackfac
#define SYMILUPACKSOL        dsymilupacksol
#define SYMILUPACKDEL        dsymilupackdel

#define MYGNLSYM             DGNLSYM
#define MYGNLSYMAMGSOLVER    DGNLSYMAMGsolver

#define MYSYMAMGSOLVER       DSYMAMGsolver
#define MYSBCG               DSYMbcg
#define MYSMATVEC            DSYMmatvec
#define MYSYMAMGINIT         DSYMAMGinit
#define MYSYMAMGGETPARAMS    DSYMAMGgetparams
#define MYSYMAMGSETPARAMS    DSYMAMGsetparams
#define MYSYMAMGFACTOR       DSYMAMGfactor

#define MYSYMAMGDELETE       DSYMAMGdelete

#define PERMSMWMAMF          DSYMperm_mwm_amf
#define PERMSMWMAMD          DSYMperm_mwm_amd
#define PERMSMWMMMD          DSYMperm_mwm_mmd
#define PERMSMWMRCM          DSYMperm_mwm_rcm
#define PERMSMWMMETISE       DSYMperm_mwm_metis_e
#define PERMSMWMMETISN       DSYMperm_mwm_metis_n

#define PERMSMC64AMF          DSYMperm_mc64_amf
#define PERMSMC64AMD          DSYMperm_mc64_amd
#define PERMSMC64MMD          DSYMperm_mc64_mmd
#define PERMSMC64RCM          DSYMperm_mc64_rcm
#define PERMSMC64METISE       DSYMperm_mc64_metis_e
#define PERMSMC64METISN       DSYMperm_mc64_metis_n

#define PERMSMWMAMFFCV        DSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        DSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        DSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        DSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     DSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     DSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       DSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       DSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       DSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       DSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    DSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    DSYMperm_mc64_metis_n_fcv

#define PERMSMWMAMFFC        DSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        DSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        DSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        DSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     DSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     DSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       DSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       DSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       DSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       DSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    DSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    DSYMperm_mc64_metis_n_fc

#else

#define MYSYMPILUCSOL        SSYMpilucsol
#define MYSYMILUC            SSYMiluc
#define SYMILUPACK           ssymilupack
#define SYMILUPACKFAC        ssymilupackfac
#define SYMILUPACKSOL        ssymilupacksol
#define SYMILUPACKDEL        ssymilupackdel

#define MYGNLSYM             SGNLSYM
#define MYGNLSYMAMGSOLVER    SGNLSYMAMGsolver

#define MYSYMAMGSOLVER       SSYMAMGsolver
#define MYSBCG               SSYMbcg
#define MYSMATVEC            SSYMmatvec
#define MYSYMAMGINIT         SSYMAMGinit
#define MYSYMAMGGETPARAMS    SSYMAMGgetparams
#define MYSYMAMGSETPARAMS    SSYMAMGsetparams
#define MYSYMAMGFACTOR       SSYMAMGfactor

#define MYSYMAMGDELETE       SSYMAMGdelete

#define PERMSMWMAMF          SSYMperm_mwm_amf
#define PERMSMWMAMD          SSYMperm_mwm_amd
#define PERMSMWMMMD          SSYMperm_mwm_mmd
#define PERMSMWMRCM          SSYMperm_mwm_rcm
#define PERMSMWMMETISE       SSYMperm_mwm_metis_e
#define PERMSMWMMETISN       SSYMperm_mwm_metis_n

#define PERMSMC64AMF          SSYMperm_mc64_amf
#define PERMSMC64AMD          SSYMperm_mc64_amd
#define PERMSMC64MMD          SSYMperm_mc64_mmd
#define PERMSMC64RCM          SSYMperm_mc64_rcm
#define PERMSMC64METISE       SSYMperm_mc64_metis_e
#define PERMSMC64METISN       SSYMperm_mc64_metis_n

#define PERMSMWMAMFFCV        SSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        SSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        SSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        SSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     SSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     SSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       SSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       SSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       SSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       SSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    SSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    SSYMperm_mc64_metis_n_fcv

#define PERMSMWMAMFFC        SSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        SSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        SSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        SSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     SSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     SSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       SSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       SSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       SSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       SSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    SSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    SSYMperm_mc64_metis_n_fc

#endif




#else



#ifdef _COMPLEX_SYMMETRIC_
#define CONJG(A)     (A)

#ifdef _SINGLE_COMPLEX_
#define MYSYMPILUCSOL        CSYMpilucsol
#define MYSYMILUC            CSYMiluc
#define SYMILUPACK           csymilupack
#define SYMILUPACKFAC        csymilupackfac
#define SYMILUPACKSOL        csymilupacksol
#define SYMILUPACKDEL        csymilupackdel

#define MYGNLSYM             CGNLSYM
#define MYGNLSYMAMGSOLVER    CGNLSYMAMGsolver

#define MYSYMAMGSOLVER       CSYMAMGsolver
#define MYSBCG               CSYMbcg
#define MYSMATVEC            CSYMmatvec
#define MYSYMAMGINIT         CSYMAMGinit
#define MYSYMAMGGETPARAMS    CSYMAMGgetparams
#define MYSYMAMGSETPARAMS    CSYMAMGsetparams
#define MYSYMAMGFACTOR       CSYMAMGfactor

#define MYSYMAMGDELETE       CSYMAMGdelete

#define PERMSMWMAMF          CSYMperm_mwm_amf
#define PERMSMWMAMD          CSYMperm_mwm_amd
#define PERMSMWMMMD          CSYMperm_mwm_mmd
#define PERMSMWMRCM          CSYMperm_mwm_rcm
#define PERMSMWMMETISE       CSYMperm_mwm_metis_e
#define PERMSMWMMETISN       CSYMperm_mwm_metis_n

#define PERMSMC64AMF          CSYMperm_mwm_amf
#define PERMSMC64AMD          CSYMperm_mwm_amd
#define PERMSMC64MMD          CSYMperm_mwm_mmd
#define PERMSMC64RCM          CSYMperm_mwm_rcm
#define PERMSMC64METISE       CSYMperm_mwm_metis_e
#define PERMSMC64METISN       CSYMperm_mwm_metis_n


#define PERMSMWMAMFFCV        CSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        CSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        CSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        CSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     CSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     CSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       CSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       CSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       CSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       CSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    CSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    CSYMperm_mc64_metis_n_fcv


#define PERMSMWMAMFFC        CSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        CSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        CSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        CSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     CSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     CSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       CSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       CSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       CSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       CSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    CSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    CSYMperm_mc64_metis_n_fc


#else
#define MYSYMPILUCSOL        ZSYMpilucsol
#define MYSYMILUC            ZSYMiluc
#define SYMILUPACK           zsymilupack
#define SYMILUPACKFAC        zsymilupackfac
#define SYMILUPACKSOL        zsymilupacksol
#define SYMILUPACKDEL        zsymilupackdel

#define MYGNLSYM             ZGNLSYM
#define MYGNLSYMAMGSOLVER    ZGNLSYMAMGsolver

#define MYSYMAMGSOLVER       ZSYMAMGsolver
#define MYSBCG               ZSYMbcg
#define MYSMATVEC            ZSYMmatvec
#define MYSYMAMGINIT         ZSYMAMGinit
#define MYSYMAMGGETPARAMS    ZSYMAMGgetparams
#define MYSYMAMGSETPARAMS    ZSYMAMGsetparams
#define MYSYMAMGFACTOR       ZSYMAMGfactor

#define MYSYMAMGDELETE       ZSYMAMGdelete

#define PERMSMWMAMF          ZSYMperm_mwm_amf
#define PERMSMWMAMD          ZSYMperm_mwm_amd
#define PERMSMWMMMD          ZSYMperm_mwm_mmd
#define PERMSMWMRCM          ZSYMperm_mwm_rcm
#define PERMSMWMMETISE       ZSYMperm_mwm_metis_e
#define PERMSMWMMETISN       ZSYMperm_mwm_metis_n

#define PERMSMC64AMF          ZSYMperm_mc64_amf
#define PERMSMC64AMD          ZSYMperm_mc64_amd
#define PERMSMC64MMD          ZSYMperm_mc64_mmd
#define PERMSMC64RCM          ZSYMperm_mc64_rcm
#define PERMSMC64METISE       ZSYMperm_mc64_metis_e
#define PERMSMC64METISN       ZSYMperm_mc64_metis_n


#define PERMSMWMAMFFCV        ZSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        ZSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        ZSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        ZSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     ZSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     ZSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       ZSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       ZSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       ZSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       ZSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    ZSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    ZSYMperm_mc64_metis_n_fcv


#define PERMSMWMAMFFC        ZSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        ZSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        ZSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        ZSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     ZSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     ZSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       ZSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       ZSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       ZSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       ZSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    ZSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    ZSYMperm_mc64_metis_n_fc

#endif

#else

#ifdef _SINGLE_COMPLEX_
#define MYSYMILUC            CHERiluc
#define MYSYMPILUCSOL        CHERpilucsol
#define SYMILUPACK           cherilupack
#define SYMILUPACKFAC        cherilupackfac
#define SYMILUPACKSOL        cherilupacksol
#define SYMILUPACKDEL        cherilupackdel

#define CONJG(A)     (conjg(A))

#define MYGNLSYM             CGNLHER
#define MYGNLSYMAMGSOLVER    CGNLHERAMGsolver

#define MYSYMAMGSOLVER       CHERAMGsolver
#define MYSBCG               CHERbcg
#define MYSMATVEC            CHERmatvec
#define MYSYMAMGINIT         CHERAMGinit
#define MYSYMAMGGETPARAMS    CHERAMGgetparams
#define MYSYMAMGSETPARAMS    CHERAMGsetparams
#define MYSYMAMGFACTOR       CHERAMGfactor

#define MYSYMAMGDELETE       CHERAMGdelete

#define PERMSMWMAMF          CHERperm_mwm_amf
#define PERMSMWMAMD          CHERperm_mwm_amd
#define PERMSMWMMMD          CHERperm_mwm_mmd
#define PERMSMWMRCM          CHERperm_mwm_rcm
#define PERMSMWMMETISE       CHERperm_mwm_metis_e
#define PERMSMWMMETISN       CHERperm_mwm_metis_n

#define PERMSMC64AMF          CHERperm_mc64_amf
#define PERMSMC64AMD          CHERperm_mc64_amd
#define PERMSMC64MMD          CHERperm_mc64_mmd
#define PERMSMC64RCM          CHERperm_mc64_rcm
#define PERMSMC64METISE       CHERperm_mc64_metis_e
#define PERMSMC64METISN       CHERperm_mc64_metis_n


#define PERMSMWMAMFFCV        CHERperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        CHERperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        CHERperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        CHERperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     CHERperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     CHERperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       CHERperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       CHERperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       CHERperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       CHERperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    CHERperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    CHERperm_mc64_metis_n_fcv


#define PERMSMWMAMFFC        CHERperm_mwm_amf_fc
#define PERMSMWMAMDFC        CHERperm_mwm_amd_fc
#define PERMSMWMMMDFC        CHERperm_mwm_mmd_fc
#define PERMSMWMRCMFC        CHERperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     CHERperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     CHERperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       CHERperm_mc64_amf_fc
#define PERMSMC64AMDFC       CHERperm_mc64_amd_fc
#define PERMSMC64MMDFC       CHERperm_mc64_mmd_fc
#define PERMSMC64RCMFC       CHERperm_mc64_rcm_fc
#define PERMSMC64METISEFC    CHERperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    CHERperm_mc64_metis_n_fc



#else
#define MYSYMILUC            ZHERiluc
#define MYSYMPILUCSOL        ZHERpilucsol
#define SYMILUPACK           zherilupack
#define SYMILUPACKFAC        zherilupackfac
#define SYMILUPACKSOL        zherilupacksol
#define SYMILUPACKDEL        zherilupackdel

#define CONJG(A)     (dconjg(A))

#define MYGNLSYM             ZGNLHER
#define MYGNLSYMAMGSOLVER    ZGNLHERAMGsolver

#define MYSYMAMGSOLVER       ZHERAMGsolver
#define MYSBCG               ZHERbcg
#define MYSMATVEC            ZHERmatvec
#define MYSYMAMGINIT         ZHERAMGinit
#define MYSYMAMGGETPARAMS    ZHERAMGgetparams
#define MYSYMAMGSETPARAMS    ZHERAMGsetparams
#define MYSYMAMGFACTOR       ZHERAMGfactor

#define MYSYMAMGDELETE       ZHERAMGdelete

#define PERMSMWMAMF          ZHERperm_mwm_amf
#define PERMSMWMAMD          ZHERperm_mwm_amd
#define PERMSMWMMMD          ZHERperm_mwm_mmd
#define PERMSMWMRCM          ZHERperm_mwm_rcm
#define PERMSMWMMETISE       ZHERperm_mwm_metis_e
#define PERMSMWMMETISN       ZHERperm_mwm_metis_n

#define PERMSMC64AMF          ZHERperm_mc64_amf
#define PERMSMC64AMD          ZHERperm_mc64_amd
#define PERMSMC64MMD          ZHERperm_mc64_mmd
#define PERMSMC64RCM          ZHERperm_mc64_rcm
#define PERMSMC64METISE       ZHERperm_mc64_metis_e
#define PERMSMC64METISN       ZHERperm_mc64_metis_n


#define PERMSMWMAMFFCV        ZHERperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        ZHERperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        ZHERperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        ZHERperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     ZHERperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     ZHERperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       ZHERperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       ZHERperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       ZHERperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       ZHERperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    ZHERperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    ZHERperm_mc64_metis_n_fcv


#define PERMSMWMAMFFC        ZHERperm_mwm_amf_fc
#define PERMSMWMAMDFC        ZHERperm_mwm_amd_fc
#define PERMSMWMMMDFC        ZHERperm_mwm_mmd_fc
#define PERMSMWMRCMFC        ZHERperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     ZHERperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     ZHERperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       ZHERperm_mc64_amf_fc
#define PERMSMC64AMDFC       ZHERperm_mc64_amd_fc
#define PERMSMC64MMDFC       ZHERperm_mc64_mmd_fc
#define PERMSMC64RCMFC       ZHERperm_mc64_rcm_fc
#define PERMSMC64METISEFC    ZHERperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    ZHERperm_mc64_metis_n_fc

#endif

#endif



#endif
