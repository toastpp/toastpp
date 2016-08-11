SHELL = /bin/bash
TBIN = $(TOASTDIR)/bin

all :: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 test12

fmod.fem ::
	@echo "Generating forward data ..."
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 2 1 1 2 2 2 2 100 3 ../meshes/mua_tgt_cyl2.nim 3 ../meshes/mus_tgt_cyl2.nim 1 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@rm fwdfem.out


test1 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 1: NCG solver, no regularisation"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 1 1e-8 1e-8 2 0 0 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test1.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test2 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 2: NCG solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 1 1e-8 1e-8 2 0 0 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test2.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test3 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 3: NCG solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 1 1e-8 1e-8 2 0 0 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test3.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test4 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 4: NCG solver, Huber reg."
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 1 1e-8 1e-8 2 0 0 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test4.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo






test5 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 5: GN solver, no regularisation"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 2 2 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test5.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test6 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 6: GN solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 2 2 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test6.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test7 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 7: GN solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 2 2 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test7.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test8 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 8: GN solver, Huber reg."
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 2 2 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test8.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo






test9 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 9: LBFGS solver, no regularisation"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 5 1e-8 1e-8 5 5 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test9.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test10 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 10: LBFGS solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 5 1e-8 1e-8 5 5 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test10.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test11 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 11: LBFGS solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 5 1e-8 1e-8 5 5 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test11.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo

test12 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 12: LBFGS solver, Huber reg."
	@echo -----------------------------------------------------------
	@echo "../meshes/cyl2.msh ../meshes/cyl_3ring.qm 1 2 2 2 2 40 40 40 n 2 5 1e-8 1e-8 5 5 2 1 1 0.01 1 1 1 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test12.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test passed! ####"
	@echo
