SHELL = /bin/sh
include ../Makefile.incl

SUBDIRS = libmath libfe libstoast libdgfe libbemdiff libbem libparam2001 libfdot supertoast bintools displace liblocalbasisRTE3D supertoast_zoltan

ifdef ILUPACKDIR
SUBDIRS += libvarPN3D
# this library requires ILUPACK support
endif

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

matlab2::
	if test -d matlab2 ;\
	then \
		(mkdir -p $(TOASTVER)/mex2; \
			cd $(TOASTVER)/mex2; \
			cp $(TOASTDIR)/src/matlab2/Makefile.src ./Makefile; \
			echo "making" all "in matlab2..."; \
			$(MAKE) $(MFLAGS) all); \
	fi;

python::
	if test -d python ; then \
		cd python; \
		$(MAKE) python; \
	fi

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
