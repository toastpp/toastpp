#!/bin/csh
#
# makefile for ILUPACK
# 
STARTDIR=$(PWD)
MYSTARTDIR=$(STARTDIR)

include user.mk
MYPLATFORM=$(PLATFORM)
include makefiles/platforms.mk
include makefiles/default.mk

# ----------------------   program specific settings   -----------------------

DIRSPARSPAK=sparspak/
DIRMETIS=metis/
DIRAMD=amd/
DIRBLAS=blas/
DIRBLASLIKE=blaslike/
DIRLAPACK=lapack/
DIRNOTDISTRIBUTED=notdistributed/

DIRPRECONDITIONERS=preconditioners/
DIRORDERINGS=orderings/
DIRTOOLS=tools/
DIRWRAPS=wraps/
DIRSPARSKIT=sparskit/
DIRSOLVERS=solvers/
DIRMUMPS=mumps/
DIRFORTRAN=fortran/

MYSTARTDIR=$(PWD)

.PHONY: samples/D/D$(MAIN) samples/Z/Z$(MAIN) samples/S/S$(MAIN) samples/C/C$(MAIN)

most: samples/D/D$(MAIN) samples/Z/Z$(MAIN)

all: samples/D/D$(MAIN) samples/Z/Z$(MAIN) samples/S/S$(MAIN) samples/C/C$(MAIN) matlab

samples/D/D$(MAIN):
	@ echo build double real main program
	@ cd samples;\
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=D"
	@ echo sample executable program can be found in samples/D
	@ $(CPLN) -sf samples/D/D$(MAIN) d$(MAIN)

samples/S/S$(MAIN):
	@ echo build single real main program
	@ cd samples;\
        $(MAKE) "ARITHMETIC=-D_SINGLE_REAL_"    "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=S"
	@ echo sample executable program can be found in samples/S
	@ $(CPLN) -sf samples/S/S$(MAIN) s$(MAIN)


samples/Z/Z$(MAIN):
	@ echo build double complex main program
	@ cd samples;  \
        $(MAKE) "ARITHMETIC="                   "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=Z"
	@ echo sample executable program can be found in samples/Z
	@ $(CPLN) -sf samples/Z/Z$(MAIN) z$(MAIN)

samples/C/C$(MAIN):
	@ echo build single complex main program
	@ cd samples;  \
        $(MAKE) "ARITHMETIC=-D_SINGLE_COMPLEX_" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" "FLOAT=C"
	@ echo sample executable program can be found in samples/C
	@ $(CPLN) -sf samples/C/C$(MAIN) c$(MAIN)

matlab:
	@ echo build mex files
	@ cd matlabsrc;  \
        $(MAKE) "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)" 
	@ echo matlab routines can be found in mex








clean:
	rm -rf  samples/?/* mex/*.mex* */*.o




