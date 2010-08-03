TBIN = $(TOASTVER)/bin
RECON = $(TBIN)/supertoast

MESHDIR = ../meshes
DATADIR = ../fwdfem

all :: header test1 test2 test3 test4
ref :: test1_ref test2_ref test3_ref test4_ref

header ::
	@echo ===========================================================
	@echo Testing supertoast: 4 tests
	@echo ===========================================================

test1 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 1 of 4 ...
	@echo Nonlinear conjugate gradient solver
	@echo -----------------------------------------------------------
	@$(RECON) supertoast1.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	cmp iter.dat test1.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out
	@echo "#### Test 1 passed! ####"

test2 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 2 of 4 ...
	@echo Gauss-Newton solver
	@echo -----------------------------------------------------------
	@$(RECON) supertoast2.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	cmp iter.dat test2.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out
	@echo "#### Test 2 passed! ####"

test3 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 3 of 4 ...
	@echo ART solver
	@echo -----------------------------------------------------------
	@$(RECON) supertoast3.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	cmp iter.dat test3.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out
	@echo "#### Test 3 passed! ####"

test4 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 4 of 4 ...
	@echo L-BFGS solver
	@echo -----------------------------------------------------------
	@$(RECON) supertoast4.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	cmp iter.dat test4.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out
	@echo "#### Test 4 passed! ####"



test1_ref ::
	@$(RECON) supertoast1.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test1.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out

test2_ref ::
	@$(RECON) supertoast2.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test2.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out

test3_ref ::
	@$(RECON) supertoast3.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test3.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out

test4_ref ::
	@$(RECON) supertoast4.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test4.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out
