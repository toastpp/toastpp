########################################################################
#                                                                      #
# This Makefile generates libblzpack.a. Instructions and options for   #
# the compiler and the command for generating the library are given    #
# in the file MACROS, which will be included below. Usually MACROS     #
# is copied from the subdirectory blzpack/etc, but the user may        #
# define his or her own MACROS, see script creator.                    #
#                                                                      #
# NB. MACROS must define FC, FCOPT, FCNOPT, FCINCL and RANLIB          #
#                                                                      #
########################################################################

.PRECIOUS:libblzpack.a
.SILENT:

##############
include MACROS
##############

SOURCE= blzaux.f blzdrs.f blzeig.f blzend.f blzext.f blzfct.f blzkll.f \
	blzrbx.f blzset.f blzstp.f blzstr.f datchk.f datset.f eigchk.f \
	eigidx.f eigprt.f eigsrt.f identy.f lzchek.f lzcopy.f lzerrs.f \
	lzhead.f lzhist.f lzioop.f lzmmry.f lznrmx.f lzpntr.f lzprt1.f \
	lzprt2.f lzprt3.f lzprt4.f lzprt5.f lzprt6.f lzprt7.f lzprt8.f \
	lzrang.f lzstp1.f lzstp2.f lzstp3.f lzstp4.f lzstp5.f lzstp6.f \
	lzstts.f lztime.f lzwrns.f mgschm.f neigab.f norm2a.f porth.f  \
	purge.f  pythag.f qtbq.f   qtbr.f   randnr.f reorth.f resnrm.f \
	rfactr.f rqalph.f rqbeta.f rvcomp.f rvmngr.f rvrstr.f setdft.f \
	seteps.f setlrm.f setssl.f setto0.f sorth.f  sorth0.f ssback.f \
	ssbext.f ssbnds.f sschek.f ssclsd.f sshist.f ssichk.f ssmovb.f \
	ssneed.f ssorgn.f ssrang.f ssrstr.f sssigl.f sssigr.f ssslct.f \
	ssspec.f sssplt.f ssstat.f sstrsf.f startr.f startx.f tbalph.f \
	tbbeta.f tbbred.f tbconv.f tbeigp.f tbgmin.f tbillc.f tbimql.f \
	upbeta.f upeta.f  uptau.f  $(EIGTB)

OBJECT=	$(SOURCE:.f=.o)

libblzpack.a: $(OBJECT)
	ar -r $@ $(OBJECT) 
	$(RANLIB) $@

seteps.o:
	$(FC) $(FCNOPT) -c seteps.f

.f.o:
	$(FC) $(FCOPT) $(FCINCL) -c $*.f

clean:
	rm -f MACROS blzaux.f *.o

# end of Makefile ######################################################
