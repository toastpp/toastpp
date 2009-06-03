include $(TOASTDIR)/Makefile.incl

all::devel user

devel::
	cd devel; \
	$(MAKE) all; \
	cd $(TOASTDIR)/doc/devel/latex; \
	$(MAKE) all; \
	dvipdf refman

user::
	cd user; \
	$(MAKE) all; \
	cd $(TOASTDIR)/doc/user/latex; \
	$(MAKE) all; \
	dvipdf refman
