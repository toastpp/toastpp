include Makefile.incl

SUBDIRS     = $(TSRC)
GUIDIR      = $(TOASTDIR)/gui
NUMERICSDIR = $(TOASTDIR)/numerics

toast::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
		(cd $$i ; echo "making" all "in $$i..."; \
			$(MAKE) $(MFLAGS) all); \
	done

matlab::
	cd src; \
	$(MAKE) matlab;

gui::
	cd $(GUIDIR); \
	echo "making" all "in $(GUIDIR)"; \
	$(MAKE) $(MFLAGS) all;

numerics::
	cd $(NUMERICSDIR); \
	echo "making" all "in $(NUMERICSDIR)"; \
	$(MAKE) $(MFLAGS) all;

doc::
	cd src; \
	$(MAKE) doc

test::
	@cd test; $(MAKE)

clean::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in  $(SUBDIRS);\
	do \
		(cd $$i ; echo "making" clean "in $$i..."; \
			$(MAKE) $(MFLAGS)  clean); \
	done
	$(RM) *~ script/*~ toastrc.d/*~

distclean::
	for i in $(SUBDIRS); \
	do \
		(cd $$i ; echo "making" distclean "in $$i..."; \
			$(MAKE) $(MFLAGS) distclean); \
	done
	$(RM) -f Makefile.incl config.status config.cache config.log
	$(RM) -f *~ script/*~ toastrc.d/*~

uninstall::
	for i in $(SUBDIRS); \
	do \
		(cd $$i ; echo "making" uninstall "in $$i..."; \
			$(MAKE) $(MFLAGS) uninstall); \
	done

gui_clean::
	cd $(GUIDIR); \
	$(MAKE) $(MFLAGS) gui_clean

gui_distclean::
	cd $(GUIDIR); \
	$(MAKE) $(MFLAGS) gui_distclean

numerics_clean::
	cd $(NUMERICSDIR); \
	$(MAKE) $(MFLAGS) clean

numerics_distclean::
	cd $(NUMERICSDIR); \
	$(MAKE) $(MFLAGS) distclean

matlab_clean:
	cd $(TOASTVER)/mex; \
	$(MAKE) clean

doc_clean::
	cd src; \
	$(MAKE) doc_clean

web_distribution_shared::
	rm -f lib
	ln -s lib_shared lib
	rm -f bin
	ln -s bin_shared bin
	tar czvf toast_shared.tar.gz bin/* docs/* scales/* script/* toastrc.d/* lib/* --exclude "lib/*.a"

web_distribution_static::
	rm -f lib
	ln -s lib_static lib
	rm -f bin
	ln -s bin_static bin
	tar czvf toast_static.tar.gz bin/* docs/* scales/* script/* toastrc.d/*

distro_src::
	cd ..; \
	zip -9 -r toast_src.zip toast2008/src -x@toast2008/exclude.lst ; \
	zip -9 -r toast_src.zip toast2008/include -x@toast2008/exclude.lst ; \
	zip -9 -r toast_src.zip toast2008/win32/VS2005 -x@toast2008/exclude.lst ; \
	zip -9 toast_src.zip toast2008/configure toast2008/gpl.txt toast2008/license.html toast2008/Makefile toast2008/mtoast_install.m toast2008/readme.txt toast2008/COMPILE.readme toast2008/numerics.tar.gz; \
	find -H toast2008 -maxdepth 1 -name "*.in" -print | zip toast_src.zip -@ ; \
	find -H toast2008 -maxdepth 1 -name "*.h" -print | zip toast_src.zip -@ ; \
	mv toast_src.zip toast2008

distro_common::
	cd ..; \
	zip -9 -r toast_common.zip toast2008/script -x@toast2008/exclude.lst ;\
	zip -9 -r toast_common.zip toast2008/test -x@toast2008/exclude.lst ; \
	zip -9    toast_common.zip toast2008/mtoast_install.m ; \
	zip -9    toast_common.zip toast2008/readme.txt ; \
	zip -9    toast_common.zip toast2008/license.html toast2008/gpl.txt; \
	mv toast_common.zip toast2008

distro_bin_win32::
	cd ..; \
	zip -9 -r toast_bin_win32.zip toast2008/win32/Release ; \
	mv toast_bin_win32.zip toast2008

distro_bin_linux::
	cd ..; \
	tar czvf toast_bin_linux.tar.gz toast2008/$(ARCHDIR)/bin \
	toast2008/$(ARCHDIR)/lib toast2008/$(ARCHDIR)/mex \
	toast2008/toastenv.csh toast2008/toastenv.sh; \
	mv toast_bin_linux.tar.gz toast2008/$(ARCHDIR)
	@echo "Linux distro tarball created in toast2008/$(ARCHDIR)/toast_bin_linux.tar.gz"
