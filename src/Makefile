SHELL = /bin/sh
include ../Makefile.incl

#SUBDIRS = libtask libtcl libtoast libctoast libdiffreg libbem displace \
#	toast ctoast fwdfem cfwdfem meshgen viewer tools supertoast
SUBDIRS = libmath libfe libstoast libbemdiff libbem libparam2001 libfdot supertoast bintools

all::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
	if test -d $$i ;\
	then \
		(mkdir -p $(TOASTOBJ)/$$i; \
			cd $(TOASTOBJ)/$$i; \
			cp $(TOASTDIR)/src/$$i/Makefile.src ./Makefile; \
			echo "making" all "in $$i..."; \
			$(MAKE) $(MFLAGS) all); \
			if test -e $(TOASTOBJ)/$$i/$$i$(LIB_SUFFIX); then \
				(cp $(TOASTOBJ)/$$i/$$i$(LIB_SUFFIX) $(TLIB)); \
			fi; \
	fi; \
	done

matlab::
	if test -d matlab ;\
	then \
		(mkdir -p $(TOASTVER)/mex; \
			cd $(TOASTVER)/mex; \
			cp $(TOASTDIR)/src/matlab/Makefile.src ./Makefile; \
			echo "making" all "in matlab..."; \
			$(MAKE) $(MFLAGS) all); \
		(mkdir -p $(TOASTVER)/mex2; \
			cd $(TOASTVER)/mex2; \
			cp $(TOASTDIR)/src/matlab/Makefile2.src ./Makefile; \
			echo "making" all "in matlab..."; \
			$(MAKE) $(MFLAGS) all); \
	fi; \

doc::
	cd doc; \
	$(MAKE) all

clean::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in  $(SUBDIRS) matlab; \
	do \
	if test -d $(TOASTOBJ)/$$i ;\
	then \
		(cd $(TOASTOBJ)/$$i ; echo "making" clean "in $$i..."; \
			$(MAKE) $(MFLAGS)  clean); \
	fi; \
	done
	$(RM) *~

doc_clean::
	cd doc;
	$(MAKE) clean

distclean::
	for i in $(SUBDIRS); \
	do \
	if test -d $$i ;\
	then \
		(cd $(TOASTOBJ)/$$i ; echo "making" distclean "in $$i..."; \
			$(MAKE) $(MFLAGS) distclean); \
	fi; \
	done
