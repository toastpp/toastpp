#
#  Top Level Makefile for LAPACK
#  Version 3.0
#  June 30, 1999
#

include make.inc

all: install lib testing blas_testing timing blas_timing

lib: lapacklib tmglib
#lib: blaslib lapacklib tmglib

#clean: cleanlib cleantesting cleantiming
clean: cleanlib

distclean:: clean
	$(RM) make.inc

install:
	( cd INSTALL; $(MAKE); ./testlsame; ./testslamch; \
	  ./testdlamch; ./testsecond; ./testdsecnd; \
	  cp lsame.f ../BLAS/SRC/; cp lsame.f ../SRC; \
	  cp slamch.f ../SRC/; cp dlamch.f ../SRC/; \
	  cp second.f ../SRC/; cp dsecnd.f ../SRC/ )

blaslib:
	( cd BLAS/SRC; $(MAKE) )

lapacklib:
	( cd SRC; $(MAKE) )

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

testing:
	( cd TESTING ; $(MAKE) )

blas_testing:
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 )
	( cd BLAS; ./xblat1s > sblat1.out; \
	           ./xblat1d > dblat1.out; \
	           ./xblat1c > cblat1.out; \
	           ./xblat1z > zblat1.out )

	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 )
	( cd BLAS; ./xblat2s < sblat2.in ; \
	           ./xblat2d < dblat2.in ; \
	           ./xblat2c < cblat2.in ; \
	           ./xblat2z < zblat2.in )

	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 )
	( cd BLAS; ./xblat3s < sblat3.in ; \
	           ./xblat3d < dblat3.in ; \
	           ./xblat3c < cblat3.in ; \
	           ./xblat3z < zblat3.in )

timing:
	( cd TIMING; $(MAKE) )

blas_timing:
	( cd TIMING/LIN; $(MAKE) )
	( cd TIMING; ./xlintims < sblasa.in > sblasa.out ; \
	             ./xlintims < sblasb.in > sblasb.out ; \
	             ./xlintims < sblasc.in > sblasc.out )
	( cd TIMING; ./xlintimd < dblasa.in > dblasa.out ; \
	             ./xlintimd < dblasb.in > dblasb.out ; \
	             ./xlintimd < dblasc.in > dblasc.out )
	( cd TIMING; ./xlintimc < cblasa.in > cblasa.out ; \
	             ./xlintimc < cblasb.in > cblasb.out ; \
	             ./xlintimc < cblasc.in > cblasc.out )
	( cd TIMING; ./xlintimz < zblasa.in > zblasa.out ; \
	             ./xlintimz < zblasb.in > zblasb.out ; \
	             ./xlintimz < zblasc.in > zblasc.out )

cleanlib:
	( cd INSTALL; $(MAKE) clean )
#	( cd BLAS/SRC; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
#	( cd TESTING/MATGEN; $(MAKE) clean )

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd TESTING; rm xlin* xeig* )

cleantiming:
	( cd TIMING/LIN; $(MAKE) clean )
	( cd TIMING/LIN/LINSRC; $(MAKE) clean )
	( cd TIMING/EIG; $(MAKE) clean )
	( cd TIMING/EIG/EIGSRC; $(MAKE) clean )
	( cd TIMING; rm xlin* xeig* )
