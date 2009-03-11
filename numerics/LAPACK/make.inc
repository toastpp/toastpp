include $(TOASTDIR)/Makefile.incl

####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.0                                             #
#  June 30, 1999                                                  #
####################################################################
#
SHELL = /bin/sh
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _LINUX
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = $(F77)
OPTS     = -funroll-all-loops $(FFLAGS) $(SHLIB_CFLAGS)
DRVOPTS  = $(OPTS)
NOOPT    = $(SHLIB_CFLAGS)
LOADER   = g77
LOADOPTS =

#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr

#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#

# BLASLIB      = ../../blas$(PLAT).a
# included from Makefile.incl

LAPACKLIB    = liblapack.a
TMGLIB       = tmglib.a
EIGSRCLIB    = eigsrc.a
LINSRCLIB    = linsrc.a
