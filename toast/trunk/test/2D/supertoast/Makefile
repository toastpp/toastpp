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
	@sed supertoast.tpl -e s/"#SOLVER#"/"PCG"/ \
		-e s/"#LINSOLVER#"/"BICGSTAB"/ \
		-e s/"#PRIOR#"/"NONE"/ \
		-e s/"#ITMAX#"/"20"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test1.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out tmp.prm
	@echo "#### Test 1 passed! ####"

test2 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 2 of 4 ...
	@echo Gauss-Newton solver
	@sed supertoast.tpl -e s/"#SOLVER#"/"LM"/ \
		-e s/"#LINSOLVER#"/"DIRECT"/ \
		-e s/"#PRIOR#"/"TV"/ \
		-e s/"#ITMAX#"/"10"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test2.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out tmp.prm
	@echo "#### Test 2 passed! ####"

test3 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 3 of 4 ...
	@echo ART solver
	@sed supertoast.tpl -e s/"#SOLVER#"/"ART"/ \
		-e s/"#LINSOLVER#"/"DIRECT"/ \
		-e s/"#PRIOR#"/"TV"/ \
		-e s/"#ITMAX#"/"10"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test3.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out tmp.prm
	@echo "#### Test 3 passed! ####"

test4 ::
	@echo -----------------------------------------------------------
	@echo supertoast test 4 of 4 ...
	@echo L-BFGS solver
	@sed supertoast.tpl -e s/"#SOLVER#"/"LBFGS"/ \
		-e s/"#LINSOLVER#"/"DIRECT"/ \
		-e s/"#PRIOR#"/"TV"/ \
		-e s/"#ITMAX#"/"10"/g > tmp.prm
	@$(RECON) supertoast4.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test4.dat
	@rm gradient*.nim update*.nim recon*.nim iter.dat gridbasis.out tmp.prm
	@echo "#### Test 4 passed! ####"



test1_ref ::
	@sed supertoast.tpl -e s/"#SOLVER#"/"PCG"/ \
		-e s/"#LINSOLVER#"/"BICGSTAB"/ \
		-e s/"#PRIOR#"/"NONE"/ \
		-e s/"#ITMAX#"/"20"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test1.dat
	@rm gradient*.nim update*.nim recon*.nim gridbasis.out tmp.prm

test2_ref ::
	@sed supertoast.tpl -e s/"#SOLVER#"/"LM"/ \
		-e s/"#LINSOLVER#"/"DIRECT"/ \
		-e s/"#PRIOR#"/"TV"/ \
		-e s/"#ITMAX#"/"10"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test2.dat
	@rm gradient*.nim update*.nim recon*.nim gridbasis.out tmp.prm

test3_ref ::
	@sed supertoast.tpl -e s/"#SOLVER#"/"ART"/ \
		-e s/"#LINSOLVER#"/"DIRECT"/ \
		-e s/"#PRIOR#"/"TV"/ \
		-e s/"#ITMAX#"/"10"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test3.dat
	@rm gradient*.nim update*.nim recon*.nim gridbasis.out tmp.prm

test4_ref ::
	@sed supertoast.tpl -e s/"#SOLVER#"/"LBFGS"/ \
		-e s/"#LINSOLVER#"/"DIRECT"/ \
		-e s/"#PRIOR#"/"TV"/ \
		-e s/"#ITMAX#"/"10"/g > tmp.prm
	@$(RECON) tmp.prm >& /dev/null
	@grep Iteration gridbasis.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > test4.dat
	@rm gradient*.nim update*.nim recon*.nim gridbasis.out tmp.prm
