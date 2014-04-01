# 09/28/1999 ###########################################################
#                                                                      #
# This Makefile generates and tests the BLZPACK drivers. It uses the   #
# Makefile in ../single or ../double to generate the drivers.          #
#                                                                      #
# N.B. The MPI and PVM drivers are not automatically tested by this    #
#      makefile.                                                       #
#                                                                      #
########################################################################

###########################
include ../../lib/precision
###########################

BLZPACK = ../../lib/libblzpack.a

# drivers for the standard eigenproblem ################################

drvsp1.x: $(BLZPACK)
	cd ../$(PRECISION); make drvsp1.x; mv drvsp1.x ../data/.
drvsp2.x: $(BLZPACK)
	cd ../$(PRECISION); make drvsp2.x; mv drvsp2.x ../data/.
drvsp3.x: $(BLZPACK)
	cd ../$(PRECISION); make drvsp3.x; mv drvsp3.x ../data/.
drvsp4.x: $(BLZPACK)
	cd ../$(PRECISION); make drvsp4.x; mv drvsp4.x ../data/.

# drivers for the generalized eigenproblem #############################

drvgp1.x: $(BLZPACK)
	cd ../$(PRECISION); make drvgp1.x; mv drvgp1.x ../data/.
drvgp2.x: $(BLZPACK)
	cd ../$(PRECISION); make drvgp2.x; mv drvgp2.x ../data/.
drvgp3.x: $(BLZPACK)
	cd ../$(PRECISION); make drvgp3.x; mv drvgp3.x ../data/.
drvgp4.x: $(BLZPACK)
	cd ../$(PRECISION); make drvgp4.x; mv drvgp4.x ../data/.

# drivers for MPI and PVM ##############################################

drvmpi.x: $(BLZPACK)
	cd ../$(PRECISION); make drvmpi.x; mv drvmpi.x ../data/.
drvpvm.x: $(BLZPACK)
	cd ../$(PRECISION); make drvpvm.x; mv drvpvm.x ../data/.
	cp drvpvm.x $(PVM_ROOT)/bin/$(PVM_ARCH)/.
drvpvm.T3E.x: $(BLZPACK)
	cd ../$(PRECISION); make drvpvm.T3E.x; mv drvpvm.T3E.x ../data/.

# test drivers for the standard eigenproblem ###########################

testdrvsp: drvsp1.x drvsp2.x drvsp3.x drvsp4.x
	drvsp1.x < drvsp1.dat > drvsp1.out
	drvsp2.x < drvsp2.dat > drvsp2.out
	drvsp3.x < drvsp3.dat > drvsp3.out
	drvsp4.x < drvsp4.dat > drvsp4.out

# test drivers for the generalized eigenproblem ########################

testdrvgp: drvgp1.x drvgp2.x drvgp3.x drvgp4.x
	drvgp1.x < drvgp1.dat > drvgp1.out
	drvgp2.x < drvgp2.dat > drvgp2.out
	drvgp3.x < drvgp3.dat > drvgp3.out
	drvgp4.x < drvgp4.dat > drvgp4.out

# cleaning #############################################################

clean:
	rm -f core EIG X *.o *.x *.out ../single/*.o ../double/*.o

# end of Makefile ######################################################
