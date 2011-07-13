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
#SuperLUroot	= $(HOME)/Codes/SuperLU_4.1
SuperLUroot	= $(TOASTDIR)/numerics/SuperLU_4.1
#SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_4.1.a
#SUPERLULIB	= superlu$(PLAT).a
TMGLIB       	= libtmglib.a

## BLASLIB   	= $(SuperLUroot)/lib/libblas.a

BLASDEF 	= -DUSE_VENDOR_BLAS
#BLASLIB 	= -L$(HOME)/lib/GotoBLAS -lgoto # defined in Makefile.incl

## ATLAS BLAS causes single-precision to fail
#BLASLIB   	= -L/usr/lib/atlas -lblas
## This BLAS causes single-precision to fail the test in SuperLU
#BLASLIB 	= -L/usr/lib -lblas

HAVE_LAPACK     = TRUE
# If TRUE, SuperLU doesn't include the LAPACK auxiliary routines
# dlamch and slamch

#LIBS		= $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
# RANLIB       = ranlib  # defined in Makefile.incl

#CC           = gcc   # defined in Makefile.incl
#CFLAGS       = -DPRNTlevel=0 -O3
NOOPTS       = 
#FORTRAN	     = g77
FORTRAN	     = $(F77)
#FFLAGS       = -O2    # defined in Makefile.incl
#LOADER       = $(FORTRAN)
LOADER       = $(CC)
LOADOPTS     =

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/matlab

