include $(TOASTDIR)/Makefile.incl

############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _linux

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = tmglib$(PLAT).a
SUPERLULIB   = superlu$(PLAT).a
BLASDEF      = -DUSE_VENDOR_BLAS

# BLASLIB     = ../blas$(PLAT).a  # defined in Makefile.incl

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
# RANLIB       = ranlib  # defined in Makefile.incl

# CC           = gcc        # defined in Makefile.incl
# CFLAGS       = -Wall -O2  # defined in Makefile.incl
FORTRAN	     = $(F77)
# FFLAGS       = -O         # defined in Makefile.incl
LOADER       = $(CC)
LOADOPTS     = #-xO3

#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/matlab
