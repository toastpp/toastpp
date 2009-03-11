TBIN = $(TOASTDIR)/bin

all :: fwdfem1

fwdfem1 ::
	@echo -----------------------------------------------------------
	@echo fwdfem test 1 ...
	@echo -----------------------------------------------------------
	@echo "0 1" | $(TBIN)/fwdfem fwdfem1.prm >& /dev/null
	cmp fmod.fem fmod_ellips_32x32_100MHz.fem
	cmp farg.fem farg_ellips_32x32_100MHz.fem
	@rm {fmod,farg}.{dat,fem} fwdfem.out
	@echo "#### Test passed! ####"
	@echo
