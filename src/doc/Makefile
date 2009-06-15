include $(TOASTDIR)/Makefile.incl

all::devel user

clean::devel_clean user_clean

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

devel_clean::
	cd devel; \
	$(MAKE) clean

user_clean::
	cd user; \
	$(MAKE) clean
