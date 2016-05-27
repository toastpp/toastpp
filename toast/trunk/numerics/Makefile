SHELL = /bin/sh
include ../Makefile.incl

NUMERICSDIR = $(TOASTDIR)/numerics
LIBDIR = $(TLIB)
INCDIR = $(TOASTDIR)/include
SUPERLUDIR = SuperLU_5.0

COMPONENT = blas fspblas lapack arpack liblbfgs

ifeq ($(BUILD_SUPERLU),1)
	COMPONENT += superlu
endif

ifeq ($(BUILD_FDOT),1)
	COMPONENT += mesa
endif

ifeq ($(SUPERLU_LOCAL_STATIC),1)
	SUPERLU_NAME = libsuperlu.a
else
	SUPERLU_NAME = libsuperlu$(LIB_SUFFIX)
endif

all:: $(COMPONENT)

blas::
	mkdir -p $(TOASTOBJ)/BLAS; \
	cd $(TOASTOBJ)/BLAS; \
	cp $(NUMERICSDIR)/BLAS/Makefile .; \
	$(MAKE); \
	cp libtoastblas$(BLAS_LIB_SUFFIX) $(TLIB)

fspblas::
	mkdir -p $(TOASTOBJ)/fspblas; \
	cd $(TOASTOBJ)/fspblas; \
	cp $(NUMERICSDIR)/fspblas/Makefile .; \
	$(MAKE); \
	cp *$(LIB_SUFFIX) $(LIBDIR)

lapack::
	mkdir -p $(TOASTOBJ)/LAPACK/SRC; \
	cd $(TOASTOBJ)/LAPACK; \
	cp $(NUMERICSDIR)/LAPACK/Makefile .; \
	cp $(NUMERICSDIR)/LAPACK/make.inc .; \
	cp $(NUMERICSDIR)/LAPACK/SRC/Makefile SRC; \
	$(MAKE) lapacklib; \
	cp SRC/libtoastlapack$(LAPACK_LIB_SUFFIX) $(TLIB)

superlu::lapack
	mkdir -p $(TOASTOBJ)/$(SUPERLUDIR)/SRC; \
	cd $(TOASTOBJ)/$(SUPERLUDIR); \
	cp $(NUMERICSDIR)/$(SUPERLUDIR)/Makefile .; \
	cp $(NUMERICSDIR)/$(SUPERLUDIR)/make.inc .; \
	cp $(NUMERICSDIR)/$(SUPERLUDIR)/SRC/Makefile SRC; \
	$(MAKE) superlulib; \
	cp SRC/$(SUPERLU_NAME) $(LIBDIR)

blzpack::
	mkdir -p $(TOASTOBJ)/blzpack/src/double; \
	cd $(TOASTOBJ)/blzpack; \
	cp $(NUMERICSDIR)/blzpack/creator .; \
	cp $(NUMERICSDIR)/blzpack/cleaner .; \
	cp $(NUMERICSDIR)/blzpack/src/double/Makefile src/double; \
	cp $(NUMERICSDIR)/blzpack/src/*.h src; \
	pushd $(NUMERICSDIR)/blzpack; \
	tar czvf blzpack.tar.gz sys drv etc --exclude "*svn" ; \
	popd; \
	tar xzvf $(NUMERICSDIR)/blzpack/blzpack.tar.gz; \
	rm $(NUMERICSDIR)/blzpack/blzpack.tar.gz; \
	./creator; \
	cp lib/libblzpack$(LIB_SUFFIX) $(LIBDIR)

arpack::
	mkdir -p $(TOASTOBJ)/ARPACK; \
	cd $(TOASTOBJ)/ARPACK; \
	$(MAKE) -f $(NUMERICSDIR)/ARPACK/Makefile.src lib; \
	cp libarpack$(ARPACK_LIB_SUFFIX) $(LIBDIR)

#	mkdir -p $(TOASTOBJ)/ARPACK/SRC; \
#	mkdir -p $(TOASTOBJ)/ARPACK/UTIL; \
#	mkdir -p $(TOASTOBJ)/ARPACK/OBJS; \
#	cd $(TOASTOBJ)/ARPACK; \
#	cp $(NUMERICSDIR)/ARPACK/Makefile .; \
#	cp $(NUMERICSDIR)/ARPACK/ARmake.inc .; \
#	cp $(NUMERICSDIR)/ARPACK/SRC/Makefile SRC; \
#	cp $(NUMERICSDIR)/ARPACK/UTIL/Makefile UTIL; \
#	$(MAKE) lib; \
#	cp OBJS/libarpack* $(LIBDIR)

liblbfgs::
	mkdir -p $(TOASTOBJ)/liblbfgs/lib; \
	mkdir -p $(TOASTOBJ)/liblbfgs/include; \
	mkdir -p $(TOASTOBJ)/liblbfgs/lib; \
	cd $(TOASTOBJ)/liblbfgs; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/configure .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/install-sh .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/ltmain.sh .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/missing .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/depcomp .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/config.sub .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/config.guess .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/config.h.in .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/Makefile.in .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/README .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/INSTALL .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/COPYING .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/AUTHORS .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/ChangeLog .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/NEWS .; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/include/* include; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/lib/* lib; \
	mkdir sample; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/sample/Makefile.in ./sample; \
	cp $(NUMERICSDIR)/liblbfgs-1.10/sample/sample.c ./sample;\
	if test "$(ARCHDIR)" = "darwin32" -o "$(ARCHDIR)" = "darwin64"; then \
		./configure CC=$(CC) CXX=$(CXX) F77=$(F77) FFLAGS='$(FFLAGS)' CFLAGS='$(CFLAGS)' CXXFLAGS='$(CXXFLAGS)' --exec-prefix=/Thisisaveeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeryloooooooooooooooooooooooooooooooongstring; \
	else \
		./configure; \
	fi; \
	make; \
	cp lib/.libs/*$(LIB_SUFFIX) $(LIBDIR); \
	cp include/* $(TOASTDIR)/include; \
	if test "$(ARCHDIR)" = "darwin32" -o "$(ARCHDIR)" = "darwin64"; then \
		install_name_tool -id @rpath/liblbfgs-1.10$(LIB_SUFFIX) $(LIBDIR)/liblbfgs-1.10$(LIB_SUFFIX); \
		install_name_tool -id @rpath/liblbfgs$(LIB_SUFFIX) $(LIBDIR)/liblbfgs$(LIB_SUFFIX); \
	fi

mesa:
	cd $(TOASTOBJ); \
	cp -r $(NUMERICSDIR)/Mesa-7.5.1 . ; \
	cd Mesa-7.5.1; \
	./configure --prefix=$(TOASTDIR) --enable-static --disable-shared --with-driver=osmesa --disable-egl; \
	make; \
	make install; \
	cd ..; \
#	rm -rf Mesa-7.5.1

clean:: blas_clean fspblas_clean lapack_clean superlu_clean blzpack_clean \
	arpack_clean liblbfgs_clean

blas_clean::
	$(RM) -r $(TOASTOBJ)/BLAS

fspblas_clean::
	$(RM) -r $(TOASTOBJ)/fspblas

lapack_clean::
	$(RM) -r $(TOASTOBJ)/LAPACK

superlu_clean::
	$(RM) -r $(TOASTOBJ)/$(SUPERLUDIR)

arpack_clean::
	$(RM) -r $(TOASTOBJ)/ARPACK

blzpack_clean::
	$(RM) -r $(TOASTOBJ)/blzpack

liblbfgs_clean::
	$(RM) -r $(TOASTOBJ)/liblbfgs

distclean:: blas_distclean fspblas_distclean lapack_distclean \
	superlu_distclean blzpack_distclean arpack_distclean liblbfgs_distclean

blas_distclean::blas_clean
	$(RM) $(LIBDIR)/libblas*

fspblas_distclean::fspblas_clean
	$(RM) $(LIBDIR)/libfsbtk*

lapack_distclean::lapack_clean
	$(RM) $(LIBDIR)/liblapack*

superlu_distclean::superlu_clean
	$(RM) $(LIBDIR)/libsuperlu*; \
	$(RM) $(INCDIR)/Cnames.h $(INCDIR)/colamd.h $(INCDIR)/?sp_defs.h \
		$(INCDIR)/?complex.h $(INCDIR)/supermatrix.h \
		$(INCDIR)/util.h

blzpack_distclean::blzpack_clean
	$(RM) -f $(LIBDIR)/libblzpack*

arpack_distclean::arpack_clean
	$(RM) $(LIBDIR)/libarpack*

liblbfgs_distclean::liblbfgs_clean
	$(RM) $(LIBDIR)/liblbfgs*
