TBIN = $(TOASTVER)/bin

all :: fwdfem1r fwdfem1c fwdfem2r fwdfem2c fwdfem3r fwdfem3c

fwdfem1r ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1: real, direct solver
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/cyl2.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	cmp fmod.fem fmod_test1r.fem
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem1c ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1: complex, direct solver
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/cyl2.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	cmp fmod.fem fmod_test1.fem
	cmp farg.fem farg_test1.fem
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem2r ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 2: real, CG solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 2 1e-12 500 0 ../meshes/cyl2.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	cmp fmod.fem fmod_test2r.fem
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem2c ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 2: complex, BiCGSTAB solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 4 1e-12 500 0 ../meshes/cyl2.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	cmp fmod.fem fmod_test2.fem
	cmp farg.fem farg_test2.fem
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem3r ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 3: real, CG solver, incomplete Choleski
	@echo -----------------------------------------------------------
	@echo "2 2 1e-12 500 3 ../meshes/cyl2.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	cmp fmod.fem fmod_test3r.fem
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem3c ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 3: BiCGSTAB solver, incomplete Choleski
	@echo -----------------------------------------------------------
	@echo "2 4 1e-12 500 3 ../meshes/cyl2.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	cmp fmod.fem fmod_test3.fem
	cmp farg.fem farg_test3.fem
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo
