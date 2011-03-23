#!/bin/csh
#

STARTDIR=$(PWD)/..
MYSTARTDIR=$(STARTDIR)

include ../user.mk
MYPLATFORM=$(PLATFORM)
include ../makefiles/platforms.mk
include ../makefiles/default.mk

MEXDIR=../mex

DOBJECTS=$(MEXDIR)/symilupackmetisn.$(EXT)\
	 $(MEXDIR)/symilupackmetise.$(EXT)\
	 $(MEXDIR)/symmwmilupackmetisn.$(EXT)\
	 $(MEXDIR)/symmwmilupackmetise.$(EXT)\
	 $(MEXDIR)/symmwmilupackamd.$(EXT)\
	 $(MEXDIR)/symmwmilupackmmd.$(EXT)\
	 $(MEXDIR)/symmwmilupackrcm.$(EXT)\
	 $(MEXDIR)/mwmilupackamd.$(EXT)\
	 $(MEXDIR)/mwmilupackmetisn.$(EXT)\
	 $(MEXDIR)/mwmilupackmetise.$(EXT)\
	 $(MEXDIR)/mwmilupackmmd.$(EXT)\
	 $(MEXDIR)/mwmilupackrcm.$(EXT)\
	 $(MEXDIR)/mwmilupacknull.$(EXT)\
         $(MEXDIR)/DSYMilupacksqmr.$(EXT)\
	 $(MEXDIR)/DSPDilupackinit.$(EXT)\
         $(MEXDIR)/DSYMilupackinit.$(EXT)\
         $(MEXDIR)/DGNLilupackinit.$(EXT)\
         $(MEXDIR)/DSPDilupackfactor.$(EXT)\
	 $(MEXDIR)/DSYMilupackfactor.$(EXT)\
	 $(MEXDIR)/DGNLilupackfactor.$(EXT)\
	 $(MEXDIR)/DSPDilupacksolver.$(EXT)\
	 $(MEXDIR)/DSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/DGNLilupacksolver.$(EXT)\
	 $(MEXDIR)/DGNLilupacktsolver.$(EXT)\
	 $(MEXDIR)/DSPDilupackdelete.$(EXT)\
	 $(MEXDIR)/DSYMilupackdelete.$(EXT)\
	 $(MEXDIR)/DGNLilupackdelete.$(EXT)\
	 $(MEXDIR)/DSYMSPDilupacksolver.$(EXT)\
	 $(MEXDIR)/DGNLSPDilupacksolver.$(EXT)\
	 $(MEXDIR)/DGNLSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/DGNLSPDilupacktsolver.$(EXT)\
	 $(MEXDIR)/DGNLSYMilupacktsolver.$(EXT)\
         $(MEXDIR)/Dloadhbo.$(EXT)\
         $(MEXDIR)/DSYMsavehbo.$(EXT)\
         $(MEXDIR)/DGNLsavehbo.$(EXT)\
	 $(MEXDIR)/DSPDilupacksol.$(EXT)\
	 $(MEXDIR)/DSYMilupacksol.$(EXT)\
	 $(MEXDIR)/DGNLilupacksol.$(EXT)\
	 $(MEXDIR)/DGNLilupacktsol.$(EXT)\
	 $(MEXDIR)/DSYMSPDilupackconvert.$(EXT)

ZOBJECTS=$(MEXDIR)/ZSYMilupacksqmr.$(EXT)\
         $(MEXDIR)/ZHERilupacksqmr.$(EXT)\
         $(MEXDIR)/ZHPDilupackinit.$(EXT)\
         $(MEXDIR)/ZSYMilupackinit.$(EXT)\
         $(MEXDIR)/ZHERilupackinit.$(EXT)\
         $(MEXDIR)/ZGNLilupackinit.$(EXT)\
	 $(MEXDIR)/ZHPDilupackfactor.$(EXT)\
	 $(MEXDIR)/ZSYMilupackfactor.$(EXT)\
	 $(MEXDIR)/ZHERilupackfactor.$(EXT)\
	 $(MEXDIR)/ZGNLilupackfactor.$(EXT)\
	 $(MEXDIR)/ZHPDilupacksolver.$(EXT)\
	 $(MEXDIR)/ZSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/ZHERilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZHPDilupackdelete.$(EXT)\
	 $(MEXDIR)/ZSYMilupackdelete.$(EXT)\
	 $(MEXDIR)/ZHERilupackdelete.$(EXT)\
	 $(MEXDIR)/ZGNLilupackdelete.$(EXT)\
	 $(MEXDIR)/ZGNLHPDilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLHPDilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLHPDilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZGNLSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLSYMilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLSYMilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZGNLHERilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLHERilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLHERilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZHERHPDilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLDSPDilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLDSPDilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLDSPDilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZGNLDSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLDSYMilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLDSYMilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZGNLDGNLilupacksolver.$(EXT)\
	 $(MEXDIR)/ZGNLDGNLilupacktsolver.$(EXT)\
	 $(MEXDIR)/ZGNLDGNLilupackhsolver.$(EXT)\
	 $(MEXDIR)/ZHERDSPDilupacksolver.$(EXT)\
	 $(MEXDIR)/ZHERDSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/ZSYMDSPDilupacksolver.$(EXT)\
	 $(MEXDIR)/ZSYMDSYMilupacksolver.$(EXT)\
	 $(MEXDIR)/ZHPDDSPDilupacksolver.$(EXT)\
         $(MEXDIR)/Zloadhbo.$(EXT)\
         $(MEXDIR)/ZGNLsavehbo.$(EXT)\
         $(MEXDIR)/ZHERsavehbo.$(EXT)\
         $(MEXDIR)/ZSYMsavehbo.$(EXT)\
	 $(MEXDIR)/ZHPDilupacksol.$(EXT)\
	 $(MEXDIR)/ZSYMilupacksol.$(EXT)\
	 $(MEXDIR)/ZHERilupacksol.$(EXT)\
	 $(MEXDIR)/ZGNLilupacksol.$(EXT)\
	 $(MEXDIR)/ZGNLilupacktsol.$(EXT)\
	 $(MEXDIR)/ZHERHPDilupackconvert.$(EXT)

OBJECTS=$(DOBJECTS) $(ZOBJECTS)

all: $(OBJECTS)

clean:
	rm -f $(OBJECTS)

.SUFFIXES: .c .$(EXT)

$(DOBJECTS): ../lib/$(PLATFORM)/libilupack.a
$(ZOBJECTS): ../lib/$(PLATFORM)/libilupack.a 
$(MEXDIR)/%.$(EXT): %.c ../lib/$(PLATFORM)/libilupack.a
# old MATLAB releases until 7.2, 
# gcc/gfortran with mexopts_gfortran.sh and 32 bit, gfortran-4.1!!!
#	mex72 -O -DmwSize=int -DmwIndex=int $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $< -L/usr/lib -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis 
# icc/ifort with mexopt_icc.sh and 32 bit
#	mex72 -O -DmwSize=int -DmwIndex=int $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis -lguide -limf -lifcore -lmkl_lapack -lmkl -lpthread



# 
# recent MATLAB releases after MATLAB 72.
#


# 32 BIT ------------------------------
# gcc/gfortran with mexopts77_gfortran.sh
#	mex -O  $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis -lmwlapack -lmwblas
# pgcc/pgf90 with mexopts77_pgi.sh
#	mex -O  $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmwlapack -lmwblas -lmetis -lm -lc
#
# icc/ifort with mexopt77_icc.sh
#	mex -O   $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis -lmwlapack -lmwblas -lm -lc
#	mex -O   $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis -lifcore -lpthread -limf -lguide -lmkl_lapack -lmkl -lm -lc
# 32 BIT ------------------------------


# 64 BIT ------------------------------
# gcc/gfortran with mexopts77_gfortran.sh 
	/opt/matlabR2010a/bin/mex -O -largeArrayDims -D_LONG_INTEGER_ $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmwlapack -lmwblas -lmetis -lm -lc -lgfortran
# pgcc/pgf90 with mexopts77_pgi.sh
#	mex -O -largeArrayDims -D_LONG_INTEGER_ $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmwlapack -lmwblas -lmetis -lm -lc
#
# icc/ifort with mexopt77_icc.sh
#	mex -O -largeArrayDims -D_LONG_INTEGER_ $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis -lmwlapack -lmwblas
##	mex -O -largeArrayDims -D_LONG_INTEGER_ $(MATCHING) $(FORTRANNAMES) -I../include -outdir $(MEXDIR)  $<  -L../lib/$(PLATFORM) -lilupack $(MUMPSDLIB) $(PARDISOLIB) -lamd -lsparspak -lblaslike -lmetis  -lmkl_lapack -lmkl -lguide -lifcore -limf -lpthread
# 64 BIT ------------------------------


$(MEXDIR)/D*.$(EXT): ../lib/$(PLATFORM)/libilupack.a  
$(MEXDIR)/D*.$(EXT): ../lib/$(PLATFORM)/libilupack.a  
$(MEXDIR)/Z*.$(EXT): ../lib/$(PLATFORM)/libilupack.a 
$(MEXDIR)/Z*.$(EXT): ../lib/$(PLATFORM)/libilupack.a 
