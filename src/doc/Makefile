include $(TOASTDIR)/Makefile.incl

all::devel user

devel::
	cd devel; \
	$(MAKE) all; \
	cd $(TOASTDIR)/doc/devel/latex; \
	$(MAKE) all; \
	dvips -f < refman.dvi > toastref.ps; \
	ps2pdf toastref.ps; \
	rm toastref.ps

user::
	cd user; \
	$(MAKE) all; \
	cd $(TOASTDIR)/doc/user/latex; \
	$(MAKE) all; \
	dvips -f < refman.dvi > toastref.ps; \
	ps2pdf toastref.ps; \
	rm toastref.ps
