all:: doc

doc::
	mkdir -p libmath
	cp ../../src/libmath/*.h libmath
	mkdir -p libfe
	cp ../../src/libfe/*.h libfe
	mkdir -p libstoast
	cp ../../src/libstoast/*.h libstoast
	doxygen Doxyfile
	rm -rf libmath
	rm -rf libfe
	rm -rf libstoast

clean::
	$(RM) -r html
	$(RM) -r latex
