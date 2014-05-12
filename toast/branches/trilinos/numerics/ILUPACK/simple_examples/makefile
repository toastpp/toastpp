# makefile using 64 bit long integer



# --------------------------------
# C sample drivers
#
# real SPD driver, single precision
#MAIN=smainspd
# real SPD driver, double precision
#MAIN=dmainspd
# complex HPD driver, single precision
#MAIN=cmainhpd
# complex HPD driver, double precision
#MAIN=zmainhpd
#
# real symmetric driver, single precision
#MAIN=smainsym
# real symmeric driver, double precision
#MAIN=dmainsym
# complex hermitian driver, single precision
#MAIN=cmainher
# complex hermitian driver, double precision
#MAIN=zmainher
# complex symmetric driver, single precision
#MAIN=cmainsym
# complex symmetric driver, double precision
#MAIN=zmainsym
#
# real general driver, single precision
#MAIN=smaingnl
# real general driver, double precision
MAIN=dmaingnl
# complex general driver, single precision
#MAIN=cmaingnl
# complex general driver, double precision
#MAIN=zmaingnl
# --------------------------------




# --------------------------------
# FORTRAN77 sample drivers
#
# real SPD driver, single precision
#MAIN=sfmainspd
# real SPD driver, double precision
#MAIN=dfmainspd
# complex HPD driver, single precision
#MAIN=cfmainhpd
# complex HPD driver, double precision
#MAIN=zfmainhpd
#
# real symmetric driver, single precision
#MAIN=sfmainsym
# real symmeric driver, double precision
#MAIN=dfmainsym
# complex hermitian driver, single precision
#MAIN=cfmainher
# complex hermitian driver, double precision
#MAIN=zfmainher
# complex symmetric driver, single precision
#MAIN=cfmainsym
# complex symmetric driver, double precision
#MAIN=zfmainsym
#
# real general driver, single precision
#MAIN=sfmaingnl
# real general driver, double precision
#MAIN=dfmaingnl
# complex general driver, single precision
#MAIN=cfmaingnl
# complex general driver, double precision
#MAIN=zfmaingnl
# --------------------------------



build: $(MAIN).o
	gfortran -O4 -fPIC -m64 -mcmodel=medium -o$(MAIN).out $(MAIN).o -L ../lib/GNU64_long/ -lilupack -lmumps -lamd -lmetis -lsparspak -llapack -lblaslike -lblas

.f.o:
	gfortran -c -O4 -fPIC -m64 -fdefault-integer-8 -fimplicit-none -mcmodel=medium $(MAIN).f

.c.o:
	gcc -c -O4 -fPIC -m64 -D__UNDERSCORE__ -D_LONG_INTEGER_ -mcmodel=medium -I ../include $(MAIN).c

clean:
	rm -f *.o *.out