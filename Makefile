include Makefile.incl

SUBDIRS     = $(TSRC)
GUIDIR      = $(TOASTDIR)/gui
NUMERICSDIR = $(TOASTDIR)/numerics
TOASTRELDIR = toast

# ========================================================
# Making

all::numerics toast

toast::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
		(cd $$i ; echo "making" all "in $$i..."; \
			$(MAKE) $(MFLAGS) all); \
	done

matlab2::
	cd src; \
	$(MAKE) matlab2;

python::
	cd src; \
	$(MAKE) python;

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

# ========================================================
# Cleaning

clean::matlab_clean toast_clean numerics_clean

toast_clean::
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in  $(SUBDIRS);\
	do \
		(cd $$i ; echo "making" clean "in $$i..."; \
			$(MAKE) $(MFLAGS)  clean); \
	done
	$(RM) *~ script/*~ toastrc.d/*~

numerics_clean::
	cd $(NUMERICSDIR); \
	$(MAKE) $(MFLAGS) clean

matlab_clean:
	cd $(TOASTVER)/mex; \
	$(MAKE) clean

matlab2_clean:
	cd $(TOASTVER)/mex2; \
	$(MAKE) clean

doc_clean::
	cd src; \
	$(MAKE) doc_clean

gui_clean::
	cd $(GUIDIR); \
	$(MAKE) $(MFLAGS) gui_clean

# ========================================================
# Dist-cleaning

distclean::matlab_distclean toast_distclean numerics_distclean
	$(RM) -f Makefile.incl config.status config.cache config.log
	$(RM) -f *~ script/*~ toastrc.d/*~

toast_distclean::
	for i in $(SUBDIRS); \
	do \
		(cd $$i ; echo "making" distclean "in $$i..."; \
			$(MAKE) $(MFLAGS) distclean); \
	done

numerics_distclean::
	cd $(NUMERICSDIR); \
	$(MAKE) $(MFLAGS) distclean

matlab_distclean:
	$(RM) -rf $(TOASTVER)/mex
	$(RM) -f $(TOASTVER)mexopts.incl

gui_distclean::
	cd $(GUIDIR); \
	$(MAKE) $(MFLAGS) gui_distclean

# ========================================================
# ???

uninstall::
	for i in $(SUBDIRS); \
	do \
		(cd $$i ; echo "making" uninstall "in $$i..."; \
			$(MAKE) $(MFLAGS) uninstall); \
	done

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
	ln -T -s trunk $(TOASTRELDIR); \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/doc -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/include -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/numerics -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/script -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/src -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 -r toast_src.zip $(TOASTRELDIR)/win32 -x@$(TOASTRELDIR)/exclude.lst ; \
	zip -9 toast_src.zip $(TOASTRELDIR)/configure $(TOASTRELDIR)/gpl.txt $(TOASTRELDIR)/license.html $(TOASTRELDIR)/Makefile $(TOASTRELDIR)/mtoast_install.m $(TOASTRELDIR)/readme.txt; \
	find -H $(TOASTRELDIR) -maxdepth 1 -name "*.in" -print | zip toast_src.zip -@ ; \
	find -H $(TOASTRELDIR) -maxdepth 1 -name "*.h" -print | zip toast_src.zip -@ ; \
	find -H $(TOASTRELDIR) -maxdepth 1 -name "*.readme" -print | zip toast_src.zip -@ ; \
	mv toast_src.zip $(TOASTRELDIR)

# ========================================================
# Distribution packages

distro_common::
	cd ..; \
	ln -T -s trunk $(TOASTRELDIR); \
	zip -9 -r toast_common.zip $(TOASTRELDIR)/script -x@$(TOASTRELDIR)/exclude.lst -x \*/.svn/\* ; \
	zip -9 -r toast_common.zip $(TOASTRELDIR)/test -x@$(TOASTRELDIR)/exclude.lst -x \*/.svn/\* ; \
	zip -9 -r toast_common.zip $(TOASTRELDIR)/examples -x@$(TOASTRELDIR)/exclude.lst -x \*/.svn/\* ; \
	zip -9 -r toast_common.zip $(TOASTRELDIR)/doc -x@$(TOASTRELDIR)/exclude.lst -x \*/.svn/\* ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/mtoast_install.m ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/mtoast2_install.m ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/README.md ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/COPYING ; \
	zip -9    toast_common.zip $(TOASTRELDIR)/INSTALL ; \
	mv toast_common.zip $(TOASTRELDIR)

distro_bin_win64::
	zip -9 -r toast_bin_win64.zip \
	win/x64/Release/bin \
	win/x64/Release/mex2 -x@exclude.lst
	@echo "Windows binary distro created in toast_bin_win64.zip"


distro_bin_linux::
	zip -9 -r toast_bin_$(ARCHDIR).zip \
	$(ARCHDIR)/bin \
	$(ARCHDIR)/lib \
	$(ARCHDIR)/mex2 \
	toastenv.csh toastenv.sh
	@echo "Linux binary distro created in toast_bin_$(ARCHDIR).zip"

distro_bin_darwin::
	ln  -s $(PWD) ../$(TOASTRELDIR); \
	cd ..; \
	zip -9 -r toast_bin_darwin64.zip $(TOASTRELDIR)/$(ARCHDIR)/bin \
		$(TOASTRELDIR)/$(ARCHDIR)/lib \
		$(TOASTRELDIR)/$(ARCHDIR)/mex2 -x "*.DS_Store" ; \
	mv toast_bin_darwin64.zip $(TOASTRELDIR)/toast_bin_$(ARCHDIR).zip
	@echo "Darwin distro tarball created in $(TOASTRELDIR)/toast_bin_$(ARCHDIR).zip"
