SHELL = /bin/bash
TBIN = $(TOASTDIR)/bin

all :: real cplx

real :: fwd_real_direct fwd_real_cg fwd_real_cg_ich
cplx :: fwd_cplx_direct fwd_cplx_bicgstab fwd_cplx_bicgstab_ich


fwd_real_direct ::
	@echo -----------------------------------------------------------
	@echo fwd_real_direct: real, direct solver
	@echo -----------------------------------------------------------
	@time (echo "../meshes/cyl2.msh ../meshes/cyl_5ring.qm 2 1 1 2 2 2 2 1 0.01 1 1 1 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null)
	@unzip -q testdata.zip fmod_real_direct.fem
	@$(TBIN)/cmpvec fmod.fem fmod_real_direct.fem -r=1e-8
	@rm fmod.fem fmod_real_direct.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_real_cg ::
	@echo -----------------------------------------------------------
	@echo fwd_real_cg: real, CG solver, no preconditioner
	@echo -----------------------------------------------------------
	@time (echo "../meshes/cyl2.msh ../meshes/cyl_5ring.qm 2 2 1e-15 10000 0 1 2 2 2 2 1 0.01 1 1 1 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null)
	@unzip -q testdata.zip fmod_real_cg.fem
	@$(TBIN)/cmpvec fmod.fem fmod_real_cg.fem -r=1e-8
	@rm fmod.fem fmod_real_cg.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_real_cg_ich ::
	@echo -----------------------------------------------------------
	@echo fwd_real_cg_ich: real, CG solver, incomplete Cholesky
	@echo -----------------------------------------------------------
	@time (echo "../meshes/cyl2.msh ../meshes/cyl_5ring.qm 2 2 1e-15 10000 2 1 2 2 2 2 1 0.01 1 1 1 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null)
	@unzip -q testdata.zip fmod_real_cg_ich.fem
	@$(TBIN)/cmpvec fmod.fem fmod_real_cg_ich.fem -r=1e-8
	@rm fmod.fem fmod_real_cg_ich.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_direct ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_direct: complex, direct solver
	@echo -----------------------------------------------------------
	@time (echo "../meshes/cyl2.msh ../meshes/cyl_5ring.qm 2 1 1 2 2 2 2 100 1 0.01 1 1 1 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null)
	@unzip -q testdata.zip {fmod,farg}_cplx_direct.fem
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_direct.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_direct.fem -r=1e-8
	@rm {fmod,farg}.fem {fmod,farg}_cplx_direct.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_bicgstab ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_bicgstab: complex, BiCGSTAB solver, no preconditioner
	@echo -----------------------------------------------------------
	@time (echo "../meshes/cyl2.msh ../meshes/cyl_5ring.qm 2 4 1e-15 10000 0 1 2 2 2 2 100 1 0.01 1 1 1 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null)
	@unzip -q testdata.zip {fmod,farg}_cplx_bicgstab.fem
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_bicgstab.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_bicgstab.fem -r=1e-8
	@rm {fmod,farg}.fem {fmod,farg}_cplx_bicgstab.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_bicgstab_ich ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_bicgstab_ich: BiCGSTAB solver, incomplete Cholesky
	@echo -----------------------------------------------------------
	@time (echo "../meshes/cyl2.msh ../meshes/cyl_5ring.qm 2 4 1e-15 10000 2 1 2 2 2 2 100 1 0.01 1 1 1 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null)
	@unzip -q testdata.zip {fmod,farg}_cplx_bicgstab_ich.fem
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_bicgstab_ich.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_bicgstab_ich.fem -r=1e-8
	@rm {fmod,farg}.fem {fmod,farg}_cplx_bicgstab_ich.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo
