# Generated automatically from Makefile.in by configure.
include $(TOASTDIR)/Makefile.incl

LIBNAME = libharmonics

LIB = $(TLIB)/$(LIBNAME).a

all:: $(LIB)

SRCS = 	surfacenet.cc diffusphere.cc optimizzz.cc usefulfan.cc legendre.cc 



OBJS = ${SRCS:.cc=.o}
HDRS = ${SRCS:.cc=.h}

%.cc: %.h
	touch $@

.cc.o :
	c++ -g -DFEM_DEBUG -DNEED_EXPLICIT_INSTANTIATION -fno-implicit-templates -Wno-deprecated  $(INCLUDE) -c $<

$(LIB): $(HDRS) $(OBJS)
	@$(RM) $(LIB)
	ar cr ${LIB} ${OBJS}
	ranlib $(LIB)

clean:: 
	$(RM) *.o *~
	$(RM) $(TLIB)/$(LIBNAME).*

distclean:: clean
	$(RM) Makefile
