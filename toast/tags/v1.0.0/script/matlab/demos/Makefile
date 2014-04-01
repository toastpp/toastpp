# This makes the input data for the demos

all :: demo3 demo6

demo3::
	echo "2 1 ellips_tri10.msh circle25_16x16.qm 1 2 2 2 2 100 4 ../../../test/2D/meshes/tgt_mua_ellips_tri10.nim 4 ../../../test/2D/meshes/tgt_mus_ellips_tri10.nim 2 1.4 0 1" | $(TOASTDIR)/bin/fwdfem
	mv fmod.fem fmod_ellips_16x16_100MHz.fem
	mv farg.fem farg_ellips_16x16_100MHz.fem
	rm fwdfem.out

demo6::
	echo "2 4 1e-10 cyl2.msh cyl_5ring.qm 1 2 2 2 2 100 4 mua_tgt_cyl2.nim 4 mus_tgt_cyl2.nim 2 1.4 0 1" | $(TOASTDIR)/bin/fwdfem
	mv fmod.fem fmod_cyl2_5ring_100MHz.fem
	mv farg.fem farg_cyl2_5ring_100MHz.fem
	rm fwdfem.out
