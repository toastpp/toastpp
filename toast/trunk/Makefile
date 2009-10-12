include Makefile.incl

SUBDIRS     = $(TSRC)
GUIDIR      = $(TOASTDIR)/gui
NUMERICSDIR = $(TOASTDIR)/numerics
TOASTRELDIR = toast2009

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
	zip -9 -r toast_src.zip $(TOASTRELDIR)/src -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/include -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/win32/VS2005 -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 toast_src.zip $(TOASTRELDIR)/configure $(TOASTRELDIR)/gpl.txt $(TOASTRELDIR)/license.html $(TOASTRELDIR)/Makefile $(TOASTRELDIR)/mtoast_install.m $(TOASTRELDIR)/readme.txt $(TOASTRELDIR)/COMPILE.readme $(TOASTRELDIR)/numerics.tar.gz; \
	find -H $(TOASTRELDIR) -maxdepth 1 -name "*.in" -print | zip toast_src.zip -@ ; \
	find -H $(TOASTRELDIR) -maxdepth 1 -name "*.h" -print | zip toast_src.zip -@ ; \
	mv toast_src.zip $(TOASTRELDIR)

distro_common::
	cd ..; \
	ln -T -s trunk $(TOASTRELDIR); \
	zip -9 -r toast_common.zip $(TOASTRELDIR)/script -x@$(TOASTRELDIR)/exclude.lst ;\
	zip -9 -r toast_common.zip $(TOASTRELDIR)/test -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/mtoast_install.m ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/readme.txt ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/license.html $(TOASTRELDIR)/gpl.txt; \
	mv toast_common.zip $(TOASTRELDIR)

distro_bin_win32::
	cd ..; \
	ln -T -s trunk $(TOASTRELDIR); \
	zip -9 -r toast_bin_win32.zip $(TOASTRELDIR)/win32/Release -x@$(TOASTRELDIR)/exclude.lst ; \
	mv toast_bin_win32.zip $(TOASTRELDIR)

distro_bin_linux::
	cd ..; \
	ln -T -s trunk $(TOASTRELDIR); \
	tar czvf toast_bin_linux.tar.gz $(TOASTRELDIR)/$(ARCHDIR)/bin \
	$(TOASTRELDIR)/$(ARCHDIR)/lib $(TOASTRELDIR)/$(ARCHDIR)/mex \
	$(TOASTRELDIR)/toastenv.csh $(TOASTRELDIR)/toastenv.sh; \
	mv toast_bin_linux.tar.gz $(TOASTRELDIR)/toast_bin_$(ARCHDIR).tar.gz
	@echo "Linux distro tarball created in $(TOASTRELDIR)/toast_bin_$(ARCHDIR).tar.gz"
