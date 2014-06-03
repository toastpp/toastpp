# $Id: Makefile.am 4 2008-04-25 07:54:10Z naoaki $

EXTRA_DIST = \
	lib.vcproj

liblbfgsincludedir = $(includedir)
liblbfgsinclude_HEADERS = \
	../include/lbfgs.h

lib_LTLIBRARIES = liblbfgs.la

liblbfgs_la_SOURCES = \
	arithmetic_ansi.h \
	arithmetic_sse_double.h \
	arithmetic_sse_float.h \
	../include/lbfgs.h \
	lbfgs.c

liblbfgs_la_LDFLAGS = \
	-no-undefined \
	-release @VERSION@

AM_CFLAGS = @CFLAGS@
INCLUDES = @INCLUDES@
