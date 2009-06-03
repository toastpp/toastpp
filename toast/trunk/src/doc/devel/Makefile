all:: doc

doc::
	cd $(TOASTDIR)/src; \
	doxygen doc/devel/Doxyfile >& doc/devel/doxygen.log
#	mkdir -p libmath
#	cp ../../libmath/*.h libmath
#	mkdir -p libfe
#	cp ../../libfe/*.h libfe
#	mkdir -p libstoast
#	cp ../../libstoast/*.h libstoast
#	mkdir -p libbem
#	cp ../../libbem/*.h libbem
#	doxygen Doxyfile
#	rm -rf libmath
#	rm -rf libfe
#	rm -rf libstoast
#	rm -rf libbem

clean::
	$(RM) -f doxygen.log
	$(RM) -r html
	$(RM) -r latex
