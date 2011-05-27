include $(TOASTDIR)/Makefile.incl

LIBNAME = libtoastblas

VPATH = $(TOASTDIR)/numerics/BLAS
LIB = $(LIBNAME)$(BLAS_LIB_SUFFIX)

all:: $(LIB)

SSRC = isamax.f sasum.f saxpy.f scasum.f scnrm2.f scopy.f sdot.f sdsdot.f sgbmv.f sgemm.f sgemv.f sger.f snrm2.f srot.f srotg.f srotm.f srotmg.f ssbmv.f sscal.f sspmv.f sspr.f sspr2.f sswap.f ssymm.f ssymv.f ssyr.f ssyr2.f ssyr2k.f ssyrk.f stbmv.f stbsv.f stpmv.f stpsv.f strmm.f strmv.f strsm.f strsv.f

DSRC = idamax.f dasum.f daxpy.f dcabs1.f dcopy.f ddot.f dgbmv.f dgemm.f dgemv.f dger.f dnrm2.f drot.f drotg.f drotm.f drotmg.f dsbmv.f dscal.f dsdot.f dspmv.f dspr.f dspr2.f dswap.f dsymm.f dsymv.f dsyr.f dsyr2.f dsyr2k.f dsyrk.f dtbmv.f dtbsv.f dtpmv.f dtpsv.f dtrmm.f dtrmv.f dtrsm.f dtrsv.f dzasum.f dznrm2.f

CSRC = icamax.f caxpy.f ccopy.f cdotc.f cdotu.f cgbmv.f cgemm.f cgemv.f cgerc.f cgeru.f chbmv.f chemm.f chemv.f cher.f cher2.f cher2k.f cherk.f chpmv.f chpr.f chpr2.f crotg.f cscal.f csrot.f csscal.f cswap.f csymm.f csyr2k.f csyrk.f ctbmv.f ctbsv.f ctpmv.f ctpsv.f ctrmm.f ctrmv.f ctrsm.f ctrsv.f

ZSRC = izamax.f xerbla.f zaxpy.f zcopy.f zdotc.f zdotu.f zdrot.f zdscal.f zgbmv.f zgemm.f zgemv.f zgerc.f zgeru.f zhbmv.f zhemm.f zhemv.f zher.f zher2.f zher2k.f zherk.f zhpmv.f zhpr.f zhpr2.f zrotg.f zscal.f zswap.f zsymm.f zsyr2k.f zsyrk.f ztbmv.f ztbsv.f ztpmv.f ztpsv.f ztrmm.f ztrmv.f ztrsm.f ztrsv.f


TOAST_SSRC = toast_sasum.f toast_isamax.f toast_scopy.f toast_sscal.f toast_sger.f toast_snrm2.f toast_ssymv.f toast_sdot.f toast_saxpy.f toast_ssyr2.f toast_srot.f toast_sgemv.f toast_strsv.f toast_sgemm.f toast_strsm.f toast_ssyrk.f toast_sswap.f

TOAST_DSRC = toast_dasum.f toast_idamax.f toast_dcopy.f toast_dscal.f toast_dger.f toast_dnrm2.f toast_dsymv.f toast_ddot.f toast_daxpy.f toast_dsyr2.f toast_drot.f toast_dgemv.f toast_dtrsv.f toast_dgemm.f toast_dtrsm.f toast_dsyrk.f toast_dswap.f

TOAST_CSRC = toast_scasum.f toast_icamax.f toast_ccopy.f toast_cscal.f toast_scnrm2.f toast_caxpy.f toast_cgemv.f toast_ctrsv.f toast_cgemm.f toast_ctrsm.f toast_cgerc.f toast_chemv.f toast_cher2.f toast_csyrk.f toast_cswap.f

TOAST_ZSRC = toast_dzasum.f toast_izamax.f toast_zcopy.f toast_zscal.f toast_dznrm2.f toast_zdotc.f toast_zdotu.f toast_zaxpy.f toast_zgemv.f toast_ztrsv.f toast_zgemm.f toast_ztrsm.f toast_zgerc.f toast_zhemv.f toast_zher2.f toast_zsyrk.f toast_zswap.f

USRC = lsame.f

SOBJ = ${SSRC:.f=.o}
DOBJ = ${DSRC:.f=.o}
COBJ = ${CSRC:.f=.o}
ZOBJ = ${ZSRC:.f=.o}
UOBJ = ${USRC:.f=.o}

TOAST_SOBJ = ${TOAST_SSRC:.f=.o}
TOAST_DOBJ = ${TOAST_DSRC:.f=.o}
TOAST_COBJ = ${TOAST_CSRC:.f=.o}
TOAST_ZOBJ = ${TOAST_ZSRC:.f=.o}

OBJS = $(SOBJ) $(DOBJ) $(COBJ) $(ZOBJ) $(UOBJ) $(TOAST_SOBJ) $(TOAST_DOBJ) $(TOAST_COBJ) $(TOAST_ZOBJ)

.f.o :
	$(F77) $(FFLAGS) $(SHLIB_CFLAGS) -c $<

# make static version of library
$(LIBNAME)$(STLIB_SUFFIX): $(OBJS)
	@$(RM) $@
	ar cr $@ $^
	$(RANLIB) $@

# make shared version of library
$(LIBNAME)$(SHLIB_SUFFIX): $(OBJS)
	@$(RM) $@
	$(SHLIB_LD) -o $@ $^ $(SHLIB_LD_LIBS) $(FLIBS)

clean::
	$(RM) *.o
	$(RM) $(LIB)

distclean:: clean
	$(RM) Makefile
