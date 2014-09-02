all:: doc

doc::
	doxygen Doxyfile
	cd $(TOASTDIR)/doc/devel/latex; \
	latex refman.tex; \
	latex refman.tex; \
	dvipdf refman.dvi toastref.pdf

clean::
	$(RM) doxygen.log

distclean::
	$(RM) -r $(TOASTDIR)/doc/devel/html
	$(RM) -r $(TOASTDIR)/doc/devel/latex
