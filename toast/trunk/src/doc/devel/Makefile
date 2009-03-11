all:: doc

doc::
	mkdir -p libmath
	cp ../../libmath/*.h libmath
	mkdir -p libfe
	cp ../../libfe/*.h libfe
	mkdir -p libstoast
	cp ../../libstoast/*.h libstoast
	mkdir -p libbem
	cp ../../libbem/*.h libbem
	doxygen Doxyfile
	rm -rf libmath
	rm -rf libfe
	rm -rf libstoast
	rm -rf libbem

clean::
	$(RM) -r html
	$(RM) -r latex
