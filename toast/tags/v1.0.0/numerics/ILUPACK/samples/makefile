include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk


OBJECTS=$(SAMPLES)

ifeq ($(FLOAT),D)
MUMPSLIB=$(MUMPSDLIB)
else
ifeq ($(FLOAT),Z)
MUMPSLIB=$(MUMPSDLIB)
else
ifeq ($(FLOAT),S)
MUMPSLIB=$(MUMPSSLIB)
else
ifeq ($(FLOAT),C)
MUMPSLIB=$(MUMPSSLIB)
else
endif
endif
endif
endif

$(FLOAT)/$(FLOAT)$(MAIN): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(FLOAT)/$(FLOAT)$(MAIN) $(OBJECTS) -L$(LIBDIR) $(LIBS) 


$(DIRSAMPLES)$(FLOAT)/*.o: ../include/ilupackmacros.h
$(DIRSAMPLES)$(FLOAT)/*.o: ../include/ilupack.h
$(DIRSAMPLES)$(FLOAT)/*.o: ../include/namesilupack.h

$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsyms.o:  $(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym.o
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)fmainsyms.o: $(DIRSAMPLES)$(FLOAT)/$(FLOAT)fmainsym.o
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainilucs.o: $(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainiluc.o
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym3s.o: $(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym3.o

$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym.o:  $(DIRSAMPLES)spdreadmatrix.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym.o:  $(DIRSAMPLES)symmessages.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym.o:  $(DIRSAMPLES)symprintperformance.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym.o:  $(DIRSAMPLES)syminitvectors.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainsym.o:  $(DIRSAMPLES)symfinalres.c

$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainspd.o:  $(DIRSAMPLES)spdreadmatrix.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainspd.o:  $(DIRSAMPLES)spdmessages.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainspd.o:  $(DIRSAMPLES)spdprintperformance.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainspd.o:  $(DIRSAMPLES)spdinitvectors.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)mainspd.o:  $(DIRSAMPLES)spdfinalres.c

$(DIRSAMPLES)$(FLOAT)/$(FLOAT)main.o:     $(DIRSAMPLES)readmatrix.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)main.o:     $(DIRSAMPLES)messages.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)main.o:     $(DIRSAMPLES)printperformance.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)main.o:     $(DIRSAMPLES)initvectors.c
$(DIRSAMPLES)$(FLOAT)/$(FLOAT)main.o:     $(DIRSAMPLES)finalres.c
