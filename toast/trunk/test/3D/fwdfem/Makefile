TBIN = $(TOASTVER)/bin

all :: fwdfem1

fwdfem1 ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1 ...
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/cyl3.msh ../meshes/cyl_5ring.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	cmp fmod.fem fmod_cyl3_test1.fem
	cmp farg.fem farg_cyl3_test1.fem
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo
