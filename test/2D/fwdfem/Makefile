TBIN = $(TOASTVER)/bin

all :: fwdfem1r fwdfem1c fwdfem2r fwdfem2c fwdfem3r fwdfem3c

fwdfem1r ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1r: real, direct solver
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 4 ../meshes/tgt_mua_ellips_tri10.nim 4 ../meshes/tgt_mus_ellips_tri10.nim 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_test1r.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem1c ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1c: complex, direct solver
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 4 ../meshes/tgt_mua_ellips_tri10.nim 4 ../meshes/tgt_mus_ellips_tri10.nim 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_test1c.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_test1c.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem2r ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 2r: real, CG solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 2 1e-14 10000 0 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_test2r.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem2c ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 2c: complex, BiCGSTAB solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 4 1e-14 10000 0 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_test2c.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_test2c.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem3r ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 3r: real, CG solver, incomplete Cholesky
	@echo -----------------------------------------------------------
	@echo "2 2 1e-14 10000 3 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_test3r.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem3c ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 3c: complex, BiCGSTAB solver, incomplete Cholesky
	@echo -----------------------------------------------------------
	@echo "2 4 1e-14 10000 3 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_test3c.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_test3c.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

