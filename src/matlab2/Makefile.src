include $(TOASTDIR)/Makefile.incl
VPATH = $(TSRC)/matlab2

TGTSRC = toastmex.cc
TGTOBJ = ${TGTSRC:.cc=.$(MEXEXT)}

AUXSRC = matlabtoast.cc mexutil.cc util.cc mtMesh.cc mtElement.cc mtQMvec.cc mtBasis.cc mtField.cc mtGradient.cc mtJacobian.cc mtRegul.cc toastlbfgs.cc mtSysmat_basis.cc matlabfdot.cc mtProjector.cc

AUXOBJ = ${AUXSRC:.cc=.o}

all:: $(TGTOBJ)

%.o : %.cc
	$(MEX) $(MEXFLAGS) $(MEXINC) -c $<

%.$(MEXEXT) : %.cc $(AUXOBJ)
	$(MEX) $(MEXINC) $(MEXLIB) $^

clean::
	rm -f *.$(MEXEXT)
	rm -f *.o
