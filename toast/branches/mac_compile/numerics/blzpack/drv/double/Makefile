# 09/28/1999 ###########################################################
#                                                                      #
# This Makefile generates test drivers for BLZPACK. Instructions and   #
# options for the compiler are given in the file MACROS, which will    #
# be included below. Usually MACROS is copied from the subdirectory    #
# blzpack/sys/MACROS, but the user may define his own MACROS, see      #
# script creator.                                                      #
#                                                                      #
# MACROS must define FC, FCOPT and FCLOAD. Depending on the driver   #
# and the installation, BLAS, LAPACK, SOLVER and PARLIB must also be   #
# defined.                                                             #
#                                                                      #
########################################################################

##############
include MACROS
##############

BLZPACK = ../../lib/libblzpack.a
ALLLIBS = $(BLZPACK) $(SOLVER) $(LAPACK) $(BLAS) $(PARLIB)

drvsp: 	drvsp1.x drvsp2.x drvsp3.x drvsp4.x
drvgp:  drvgp1.x drvgp2.x drvgp3.x drvgp4.x
drvmpi: drvmpi.x
drvpvm: drvpvm.x

drvsp1.x: drvsp1.o $(BLZPACK)
	$(FC) $(FCLOAD) drvsp1.o -o $@ $(ALLLIBS)

drvsp2.x: drvsp2.o $(BLZPACK)
	$(FC) $(FCLOAD) drvsp2.o -o $@ $(ALLLIBS)

drvsp3.x: drvsp3.o $(BLZPACK)
	$(FC) $(FCLOAD) drvsp3.o -o $@ $(ALLLIBS)

drvsp4.x: drvsp4.o $(BLZPACK)
	$(FC) $(FCLOAD) drvsp4.o -o $@ $(ALLLIBS)

drvgp1.x: drvgp1.o $(BLZPACK)
	$(FC) $(FCLOAD) drvgp1.o -o $@ $(ALLLIBS)

drvgp2.x: drvgp2.o $(BLZPACK)
	$(FC) $(FCLOAD) drvgp2.o -o $@ $(ALLLIBS)

drvgp3.x: drvgp3.o $(BLZPACK)
	$(FC) $(FCLOAD) drvgp3.o -o $@ $(ALLLIBS)

drvgp4.x: drvgp4.o $(BLZPACK)
	$(FC) $(FCLOAD) drvgp4.o -o $@ $(ALLLIBS)

drvmpi.x: drvmpi.o $(BLZPACK)
	$(FC) $(FCLOAD) drvmpi.o -o $@ $(ALLLIBS)

drvpvm.x: drvpvm.o $(BLZPACK)
	$(FC) $(FCLOAD) drvpvm.o -o $@ $(ALLLIBS)

drvpvm.T3E.x: drvpvm.T3E.o $(BLZPACK)
	$(FC) $(FCLOAD) drvpvm.T3E.o -o $@ $(ALLLIBS)

mytest.x: mytest.o $(BLZPACK)
	$(FC) $(FCLOAD) mytest.o -o $@ $(ALLLIBS)

clean:
	rm -f MACROS core *.o *.x *.out

.f.o:
	$(FC) $(FCOPT) $(FCINCL) -c $?

# end of Makefile ######################################################
