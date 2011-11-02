TBIN = $(TOASTVER)/bin

all :: test1 test2 test3 test4 test5 test6 test7 test8

fmod.fem ::
	@echo "Generating forward data ..."
	@echo "2 1 ../meshes/ellips_tri10.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 100 4 ../meshes/tgt_mua_ellips_tri10.nim 4 ../meshes/tgt_mus_ellips_tri10.nim 2 1.4 0 1" | $(TBIN)/fwdfem >& /dev/null
	@rm fwdfem.out


test1 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 1: NCG solver, no regularisation"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test1.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 1 passed! ####"
	@echo

test2 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 2: NCG solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test2.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 2 passed! ####"
	@echo

test3 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 3: NCG solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test3.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 3 passed! ####"
	@echo

test4 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 4: NCG solver, Huber reg."
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 1 1e-8 1e-8 10 0 0 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test4.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 4 passed! ####"
	@echo






test5 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 5: GN solver, no regularisation"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test5.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 5 passed! ####"
	@echo

test6 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 6: GN solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test6.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 6 passed! ####"
	@echo

test7 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 7: GN solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test7.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 7 passed! ####"
	@echo

test8 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 8: GN solver, Huber reg."
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 2 5 1e-8 5 1e-2 y 0 1 0 1 1 n 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test8.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 8 passed! ####"
	@echo






test9 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 9: LBFGS solver, no regularisation"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 0 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test9.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 9 passed! ####"
	@echo

test10 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 10: LBFGS solver, Tikhonov (TK1)"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 3 1e-4 0" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test10.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 10 passed! ####"
	@echo

test11 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 11: LBFGS solver, total variation (TV)"
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 1 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test11.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 11 passed! ####"
	@echo

test12 :: fmod.fem
	@echo -----------------------------------------------------------
	@echo "supertoast test 12: LBFGS solver, Huber reg."
	@echo -----------------------------------------------------------
	@echo "../meshes/circle25_32.msh ../meshes/circle25_32x32.qm 1 2 2 2 2 80 80 n 2 5 1e-8 1e-8 10 10 2 1 2 0.025 2 2 2 1.4 fmod.fem farg.fem 100 4 1 0 0 0 4 1e-4 0 1e-2" | $(TBIN)/supertoast >& /dev/null
	@grep Iteration supertoast.out | \
		sed -e s/"Iteration [0-9]*  CPU [0-9.]*  OF "//g > iter.dat
	@cmp iter.dat test12.dat
	@rm gradient*.nim update*.nim recon*.raw iter.dat supertoast.out
	@echo "#### Test 12 passed! ####"
	@echo
