TBIN = $(TOASTVER)/bin

all :: fwdfem1 fwdfem2

fwdfem1 ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1 ...
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 1 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	cmp fmod.fem fmod_ellips_32x32_100MHz.fem
	cmp farg.fem farg_ellips_32x32_100MHz.fem
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwdfem2 ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 2 ...
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 0 1 2 2.0 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	cmp fmod.fem fmod_ellips_32x32_CW.fem
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo
