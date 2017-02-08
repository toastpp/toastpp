SHELL = /bin/bash
TBIN = $(TOASTDIR)/bin

all   :: ncg gn lbfgs fmod_clean

ncg   :: inv_ncg inv_ncg_tk1 inv_ncg_tv inv_ncg_huber
gn    :: inv_gn inv_gn_tk1 inv_gn_tv inv_gn_huber
lbfgs :: inv_lbfgs inv_lbfgs_tk1 inv_lbfgs_tv inv_lbfgs_huber


fmod.fem ::
	@echo "Generating forward data ..."
	@echo "../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 2 1 1 2 2 2 2 100 3 ../meshes/tgt_mua_ellips_tri10.nim 3 ../meshes/tgt_mus_ellips_tri10.nim 1 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@rm fwdfem.out

fmod_clean ::
	@rm fmod.fem farg.fem


inv_ncg :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_ncg: supertoast NCG solver, no regularisation"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_cg.dat
	@cmp iter.dat inv_cg.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_cg.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_ncg_tk1 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_ncg_tk1: supertoast NCG solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_ncg_tk1.dat
	@cmp iter.dat inv_ncg_tk1.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_ncg_tk1.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_ncg_tv :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_ncg_tv: supertoast NCG solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_ncg_tv.dat
	@cmp iter.dat inv_ncg_tv.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_ncg_tv.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_ncg_huber :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_ncg_huber: supertoast NCG solver, Huber reg."
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_ncg_huber.dat
	@cmp iter.dat inv_ncg_huber.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_ncg_huber.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo






inv_gn :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_gn: supertoast GN solver, no regularisation"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_gn.dat
	@cmp iter.dat inv_gn.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_gn.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_gn_tk1 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_gn_tk1: supertoast GN solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_gn_tk1.dat
	@cmp iter.dat inv_gn_tk1.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_gn_tk1.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_gn_tv :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_gn_tv: supertoast GN solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_gn_tv.dat
	@cmp iter.dat inv_gn_tv.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_gn_tv.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_gn_huber :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_gn_huber: supertoast GN solver, Huber reg."
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_gn_huber.dat
	@cmp iter.dat inv_gn_huber.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_gn_huber.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo






inv_lbfgs :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_lbfgs: supertoast LBFGS solver, no regularisation"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_lbfgs.dat
	@cmp iter.dat inv_lbfgs.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_lbfgs.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_lbfgs_tk1 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_lbfgs_tk1: supertoast LBFGS solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_lbfgs_tk1.dat
	@cmp iter.dat inv_lbfgs_tk1.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_lbfgs_tk1.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_lbfgs_tv :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_lbfgs_tv: supertoast LBFGS solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_lbfgs_tv.dat
	@cmp iter.dat inv_lbfgs_tv.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_lbfgs_tv.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

inv_lbfgs_huber :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "inv_lbfgs_huber: supertoast LBFGS solver, Huber reg."
	@echo -----------------------------------------------------------
	@time (echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 1 0.025 1 2 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null)
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@unzip -q testdata.zip inv_lbfgs_huber.dat
	@cmp iter.dat inv_lbfgs_huber.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat inv_lbfgs_huber.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo
