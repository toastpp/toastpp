TBIN = $(TOASTVER)/bin

# run all tests
all :: fwd_real_direct fwd_cplx_direct fwd_real_cg fwd_cplx_bicgstab fwd_real_cg_ich fwd_cplx_bicgstab_diag fwd_cplx_bicgstab_ich fwd_cplx_bicgstab_dilu fwd_cplx_gmres fwd_cplx_gmres_diag fwd_cplx_gmres_ich fwd_cplx_gmres_dilu

# subgroups problem type (real/complex)
real     :: fwd_real_direct fwd_real_cg fwd_real_cg_diag fwd_real_cg_ich
cplx     :: fwd_cplx_direct fwd_cplx_bicgstab fwd_cplx_bicgstab_diag fwd_cplx_bicgstab_ich fwd_cplx_bicgstab_dilu fwd_cplx_gmres fwd_cplx_gmres_diag fwd_cplx_gmres_ich fwd_cplx_gmres_dilu

# subgroups solver types (direct/cg/bicgstab/gmres)
direct   :: fwd_real_direct fwd_cplx_direct
cg       :: fwd_real_cg fwd_real_cg_diag fwd_real_cg_ich
bicgstab :: fwd_cplx_bicgstab fwd_cplx_bicgstab_diag fwd_cplx_bicgstab_ich fwd_cplx_bicgstab_dilu
gmres    :: fwd_cplx_gmres fwd_cplx_gmres_diag fwd_cplx_gmres_ich fwd_cplx_gmres_dilu

# subgroups preconditioners (nopr/diag/ich/dilu)
nopr     :: fwd_real_cg fwd_cplx_bicgstab fwd_cplx_gmres
diag     :: fwd_real_cg_diag fwd_cplx_bicgstab_diag fwd_cplx_gmres_diag
ich      :: fwd_real_cg_ich fwd_cplx_bicgstab_ich fwd_cplx_gmres_ich
dilu     :: fwd_cplx_bicgstab_dilu fwd_cplx_gmres_dilu


fwd_real_direct ::
	@echo -----------------------------------------------------------
	@echo fwd_real_direct: real, direct solver
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 4 ../meshes/tgt_mua_ellips_tri10.nim 4 ../meshes/tgt_mus_ellips_tri10.nim 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_real_direct.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_real_cg ::
	@echo -----------------------------------------------------------
	@echo fwd_real_cg: real, CG solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 2 1e-14 10000 0 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_real_cg.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_real_cg_diag ::
	@echo -----------------------------------------------------------
	@echo fwd_real_ch_diag: real, CG solver, DIAG preconditioner
	@echo -----------------------------------------------------------
	@echo "2 2 1e-14 10000 1 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_real_cg_diag.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_real_cg_ich ::
	@echo -----------------------------------------------------------
	@echo fwd_real_ch_ich: real, CG solver, ICH preconditioner
	@echo -----------------------------------------------------------
	@echo "2 2 1e-14 10000 2 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem_cw >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_real_cg_ich.fem -r=1e-8
	@rm fmod.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_direct ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_direct: complex, direct solver
	@echo -----------------------------------------------------------
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 4 ../meshes/tgt_mua_ellips_tri10.nim 4 ../meshes/tgt_mus_ellips_tri10.nim 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_direct.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_direct.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_bicgstab ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_bicgstab: complex, BiCGSTAB solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 4 1e-14 10000 0 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_bicgstab.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_bicgstab.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_bicgstab_diag ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_bicgstab_diag: complex, BiCGSTAB solver, DIAG preconditioner
	@echo -----------------------------------------------------------
	@echo "2 4 1e-14 10000 1 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_bicgstab_diag.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_bicgstab_diag.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_bicgstab_ich ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_bicgstab_ich: complex, BiCGSTAB solver, ICH preconditioner
	@echo -----------------------------------------------------------
	@echo "2 4 1e-14 10000 2 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_bicgstab_ich.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_bicgstab_ich.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_bicgstab_dilu ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_bicgstab_dilu: complex, BiCGSTAB solver, DILU preconditioner
	@echo -----------------------------------------------------------
	@echo "2 4 1e-14 10000 3 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_bicgstab_dilu.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_bicgstab_dilu.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_gmres ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_gmres: complex, GMRES solver, no preconditioner
	@echo -----------------------------------------------------------
	@echo "2 5 1e-14 10000 0 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_gmres.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_gmres.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_gmres_diag ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_gmres_diag: complex, GMRES solver, DIAG preconditioner
	@echo -----------------------------------------------------------
	@echo "2 5 1e-14 10000 1 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_gmres_diag.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_gmres_diag.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_gmres_ich ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_gmres_ich: complex, GMRES solver, ICH preconditioner
	@echo -----------------------------------------------------------
	@echo "2 5 1e-14 10000 2 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_gmres_ich.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_gmres_ich.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo

fwd_cplx_gmres_dilu ::
	@echo -----------------------------------------------------------
	@echo fwd_cplx_gmres_dilu: complex, GMRES solver, DILU preconditioner
	@echo -----------------------------------------------------------
	@echo "2 5 1e-14 10000 3 ../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 2 0.01 2 1 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@$(TBIN)/cmpvec fmod.fem fmod_cplx_gmres_dilu.fem -r=1e-8
	@$(TBIN)/cmpvec farg.fem farg_cplx_gmres_dilu.fem -r=1e-8
	@rm {fmod,farg}.fem fwdfem.out
	@echo "#### Test passed! ####"
	@echo
