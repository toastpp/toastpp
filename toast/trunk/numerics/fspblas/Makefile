include $(TOASTDIR)/Makefile.incl

LIBNAME = libfsbtk

VPATH = $(TOASTDIR)/numerics/fspblas
LIB = $(LIBNAME)$(LIB_SUFFIX)

all:: $(LIB)

SRC = \
	doubblas.f \
        bcomm.f  belmm.f  bsrmm.f  cscmm.f  diamm.f  jadmm.f  vbrmm.f \
        bdimm.f  bscmm.f  coomm.f  csrmm.f  ellmm.f  skymm.f           \
        bdism.f  bscsm.f  cscsm.f  diasm.f  jadsm.f  vbrsm.f           \
        belsm.f  bsrsm.f  csrsm.f  ellsm.f  skysm.f  jadrp.f  xerbla.f \
	dblasz.f csrmmz.f cscmmz.f

OBJS = ${SRC:.f=.o}

.f.o :
	$(F77) $(FFLAGS) $(SHLIB_CFLAGS) -c $<

# make static version of library
$(LIBNAME)$(STLIB_SUFFIX): $(OBJS)
	@$(RM) $@
	ar cr $@ $^
	$(RANLIB) $@

# make shared version of library
$(LIBNAME)$(SHLIB_SUFFIX): $(OBJS)
	@$(RM) $@
	$(SHLIB_LD) -o $@ $^ $(SHLIB_LD_LIBS) $(FLIBS) -L$(TLIB) $(BLASLIB)

clean::
	$(RM) $(OBJS)
	$(RM) $(LIB)

distclean:: clean
	$(RM) Makefile
