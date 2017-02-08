#! /bin/bash

dnl Process this file with autoconf to produce a configure script.
AC_INIT(src/toast/toast.cc)

# sanity check
if test "x$TOASTDIR" = x ; then
    echo "TOASTDIR environment variable not set!"
    echo Aborting configure.
    exit
fi

dnl Location of auxiliary configure scripts
AC_CONFIG_AUX_DIR(config)

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
dnl AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_YACC

dnl Switch to C++ tests
AC_LANG_CPLUSPLUS

dnl Checks for header files.
AC_PATH_XTRA
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(fcntl.h limits.h malloc.h sgtty.h sys/file.h sys/ioctl.h sys/time.h termio.h unistd.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_UID_T
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_STRUCT_ST_BLKSIZE
AC_STRUCT_ST_BLOCKS
AC_STRUCT_ST_RDEV
AC_HEADER_TIME
AC_STRUCT_TM

dnl Checks for library functions.
AC_TYPE_GETGROUPS
AC_PROG_GCC_TRADITIONAL
AC_FUNC_MEMCMP
AC_TYPE_SIGNAL
AC_FUNC_UTIME_NULL
AC_FUNC_VFORK
AC_FUNC_VPRINTF
AC_FUNC_WAIT3
AC_CHECK_FUNCS(ftime getcwd gethostname gettimeofday getwd mkdir mkfifo mktime putenv rmdir select socket strcspn strdup strerror strspn strstr strtod strtol strtoul uname)

# ========================================================================
# Check whether X header files work

AC_LANG_CPLUSPLUS

TMP_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS $X_CFLAGS"
NEED_MPI=no
NEED_PTHREAD=no
NEED_PARMETIS=no
NEES_SCOTCH=no

AC_TRY_COMPILE([#include <X11/Xlib.h>], [return 0;], xok=yes, xok=no)

if test "$xok" = "no" -a "$GXX" = "yes" ; then
    CXXFLAGS="$CXXFLAGS -fpermissive"
    AC_TRY_COMPILE([#include <X11/Xlib.h>], [return 0;], xok=yes, xok=no)
    if test "$xok" = "yes" ; then
	X_CFLAGS="$X_CFLAGS -fpermissive"
    fi
fi

# ========================================================================
# Check for xview include files and libraries

XVIEWINC=""
AC_TRY_COMPILE([#include <xview/xview.h>], [return 0;], HAVE_XVIEW=yes, HAVE_XVIEW=no)
if test "$HAVE_XVIEW" = "no" ; then
    XVIEWINC="-I/usr/openwin/include"
    CXXFLAGS="$TMP_CXXFLAGS $XVIEWINC"
    AC_TRY_COMPILE([#include <xview/xview.h>], [return 0;], HAVE_XVIEW=yes, HAVE_XVIEW=no)
fi

CXXFLAGS=$TMP_CXXFLAGS

# ========================================================================
# Check for BLAS

AC_ARG_WITH(blas,
    [  --with-blas=BLAS        Linker command to bind BLAS package, e.g.]
    [                          '-lblas' or '/usr/lib/libblas.a'],
    BLASLIB=$withval, BLASLIB="-ltoastblas")

# ========================================================================
# Check for LAPACK

AC_ARG_WITH(lapack,
    [  --with-lapack=LIB       Linker command to bind LAPACK package, e.g.]
    [                          '-llapack' or '/usr/lib/liblapack.a'],
    LAPACKLIB=$withval, LAPACKLIB="-ltoastlapack")

# ========================================================================
# Check for SuperLU

AC_ARG_WITH(superlu,
    [  --with-superlu=LIB      Linker command to bind SuperLU package, e.g.]
    [                          'libsuperlu.a' or '-lsuperlu (default)'],
    SUPERLULIB=$withval, SUPERLULIB='-lsuperlu')
# Note: Changed this to default to static library to avoid linking problems
# in mex files.
if test "$SUPERLULIB" = '$(TOASTDIR)/lib/libsuperlu.a' -o "$SUPERLULIB" = '-lsuperlu' ; then
    BUILD_SUPERLU=1
    if test "$SUPERLULIB" = '$(TOASTDIR)/lib/libsuperlu.a' ; then
       SUPERLU_LOCAL_STATIC=1
       SUPERLULIB=$TOASTDIR/lib/libsuperlu.a
    fi
else
    BUILD_SUPERLU=0
fi

# ========================================================================
# ========================================================================
# Check for GUI libraries
# ========================================================================
# ========================================================================

GUIFLAGS="$X_CFLAGS"
GUILIBS="$X_LIBS -lX11 $X_EXTRA_LIBS $X_PRE_LIBS"

#make sure any GUI libraries in the toast lib directory are found (last)
LDFLAGS="$LDFLAGS -L$TOASTDIR/lib"

# ========================================================================
# Check for Tcl

AC_ARG_WITH(tcl,
    [  --with-tcl=PATH         Path to Tcl installation directory]
    [                          Default is auto-search for Tcl],
	TCLPATH=$withval, TCLPATH="")

pLDFLAGS=$LDFLAGS
pLIBS=$LIBS
if test "$TCLPATH" != "" ; then
   TCLLIBPATH="-L$TCLPATH/lib"
fi
LDFLAGS="TCLLIBPATH $LDFLAGS"
LIBS=$GUILIBS
AC_CHECK_LIB(tcl, main,
	[GUILIBS="$TCLLIBPATH -ltcl $LIBS" LOCAL_TCL=""],
	[GUILIBS="-ltcl $LIBS" LOCAL_TCL="yes"],
	-lnsl -ldl)
LDFLAGS=$pLDFLAGS
LIBS=$pLIBS
AC_SUBST(LOCAL_TCL)

# ========================================================================
# Check for Tk

AC_ARG_WITH(tk,
    [  --with-tk=PATH          Path to Tk installation directory]
    [                          Default is auto-search for Tk],
	TKPATH=$withval, TKPATH="")

pLDFLAGS=$LDFLAGS
pLIBS=$LIBS
if test "$TKPATH" != "" ; then
   TKLIBPATH="-L$TKPATH/lib"
fi
LDFLAGS="TKLIBPATH $LDFLAGS"
LIBS=$GUILIBS
AC_CHECK_LIB(tk, main,
	[GUILIBS="$TKLIBPATH -ltk $LIBS" LOCAL_TK=""],
	[GUILIBS="-ltk $LIBS" LOCAL_TK="yes"],
	-lnsl -ldl)
LDFLAGS=$pLDFLAGS
LIBS=$pLIBS
AC_SUBST(LOCAL_TK)

# ========================================================================
# Check for Tix
# Note: If we need to build Tix then we also need the Tcl/Tk sources,
# even if Tcl/Tk is already installed

AC_ARG_WITH(tix,
    [  --with-tix=PATH         Path to Tix installation directory]
    [                          Default is auto-search for Tix],
	TIXPATH=$withval, TIXPATH="")

pLDFLAGS=$LDFLAGS
pLIBS=$LIBS
if test "$TIXPATH" != "" ; then
   TIXLIBPATH="-L$TIXPATH/lib"
fi
LDFLAGS="$TIXLIBPATH $LDFLAGS"
LIBS=$GUILIBS
AC_CHECK_LIB(tix, Tix_Init,
	[GUILIBS="$TIXLIBPATH -ltix $LIBS" LOCAL_TIX=""],
	[GUILIBS="-ltix $LIBS" LOCAL_TIX="yes"])
LDFLAGS=$pLDFLAGS
LIBS=$pLIBS
AC_SUBST(LOCAL_TIX)

# ========================================================================
# Check for BLT

AC_ARG_WITH(blt,
    [  --with-blt=PATH         Path to BLT installation directory]
    [                          Default is auto-search for BLT],
       BLTPATH=$withval, BLTPATH="")

pLDFLAGS=$LDFLAGS
pLIBS=$LIBS
if test "$BLTPATH" != "" ; then
   BLTLIBPATH="-L$BLTPATH/lib"
   GUIFLAGS="$GUIFLAGS -I$BLTPATH/include"
fi
LDFLAGS="$BLTLIBPATH $LDFLAGS"
LIBS=$GUILIBS
AC_SEARCH_LIBS(Blt_Init, [BLT24 BLT],
	[GUILIBS="$BLTLIBPATH $LIBS" LOCAL_BLT=""],
	[GUILIBS="-lBLT $LIBS" LOCAL_BLT="yes"])
LDFLAGS=$pLDFLAGS
LIBS=$pLIBS
AC_SUBST(LOCAL_BLT)

# ========================================================================
# Check if ILUPACK support is requested

AC_ARG_WITH(ilu,
    [  --with-ilu=PATH         Path to ILUPACK base directory]
    [                          If not specified, no ILU support is provided],
       ILUPACKDIR=$withval, ILUPACKDIR="")

# ========================================================================
# Check if Mesa custom path is requested

AC_ARG_WITH(mesa,
    [  --with-mesa=PATH        Path to Mesa 3D base directory]
    [                          This should contain the header files in include/GL]
    [                          and the Mesa libraries in lib. You can also use]
    [                          the --with-mesa-incdir and --with-mesa-libdir]
    [                          options.]
    [                          If not specified, default system paths are used],
       MESADIR=$withval, MESADIR="default")
if test "$MESADIR" = "default" ; then
    MESAINCDIR=""
    MESALIBDIR=""
    MESALIB="-lOSMesa"
else
    MESAINCDIR="$MESADIR/include"
    MESALIBDIR="$MESADIR/lib"
    MESALIB="-lOSMesa"
fi

AC_ARG_WITH(mesa-includedir,
    [  --with-mesa-includedir  Path to Mesa 3D include directory]
    [                          This should contain the GL subdirectory],
       MESAINCDIR=$withval,)

AC_ARG_WITH(mesa-libdir,
    [  --with-mesa-libdir      Path to Mesa 3D library directory]
    [                          This should contain libOSMesa.so etc.],
       MESALIBDIR="$withval",)

# ========================================================================
# Check if FMT support is requested

AC_ARG_ENABLE(fdot,
    [  --enable-fdot           Include FMT library support],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    TOAST_FDOT=1
    OGLLIB=$MESALIB
    if test "x$MESALIBDIR" != "x" ; then
	OGLLIB="-L$MESALIBDIR $OGLLIB"
    fi
    if test "x$MESAINCDIR" != "x" ; then
	OGLINC="-I$MESAINCDIR"
    fi
else
    TOAST_FDOT=0
fi

# ========================================================================
# Below are some explicit definitions for system-dependent variables
# mainly concerned with generating shared libraries

# Step 1: set the variable "system" to hold the name and version number
# for the system.  This can usually be done via the "uname" command, but
# there are a few systems, like Next, where this doesn't work.

AC_MSG_CHECKING([system version (for dynamic loading)])
if test -f /usr/lib/NextStep/software_version; then
    system=NEXTSTEP-`awk '/3/,/3/' /usr/lib/NextStep/software_version`
else
    system=`uname -s`-`uname -r`
    if test "$?" -ne 0 ; then
	AC_MSG_RESULT([unknown (can't find uname command)])
	system=unknown
    else
	# Special check for weird MP-RAS system (uname returns weird
	# results, and the version is kept in special file).
    
	if test -r /etc/.relid -a "X`uname -n`" = "X`uname -s`" ; then
	    system=MP-RAS-`awk '{print $3}' /etc/.relid`
	fi
	if test "`uname -s`" = "AIX" ; then
	    system=AIX-`uname -v`.`uname -r`
	fi
        if test "`uname`" = "Darwin" ; then
            system=Darwin
        fi
	AC_MSG_RESULT($system)
    fi
fi

# Step 2: check for existence of -ldl library.  This is needed because
# Linux can use either -ldl or -ldld for dynamic loading.

AC_CHECK_LIB(dl, dlopen, have_dl=yes, have_dl=no)

# Step 3: set configuration options based on system name and version.

# default shared and static library suffixes and flags
SHLIB_SUFFIX=".so"
STLIB_SUFFIX=".a"
DL_LIBS="-ldl"
RPATH=""
VERSIONFLAG=""

case $system in
     Linux*)
	if test "`uname -m`" = "x86_64" ; then
	   ARCHDIR="linux64"
	else
	   ARCHDIR="linux32"
	fi
	;;
     Darwin*)
	if test "`uname -m`" = "x86_64" ; then
	   ARCHDIR="darwin64"
	   VERSIONFLAG="-m64"
	else
	   ARCHDIR="darwin32"
	   VERSIONFLAG="-m32"
	fi
	;;
esac

AC_SUBST(ARCHDIR)

case $system in
    IRIX-5.*|IRIX-6.*|IRIX64-6.5*)
	SHLIB_CFLAGS=""
	SHLIB_LD="ld -n32 -shared -rdata_shared"
	SHLIB_LD_LIBS='${LIBS}'
	DL_LIBS=""
	;;
    IRIX64-6.*)
	SHLIB_CFLAGS=""
	SHLIB_LD="ld -32 -shared -rdata_shared -rpath /usr/local/lib"
	SHLIB_LD_LIBS=""
	DL_LIBS=""
	;;
    Linux*)
	SHLIB_CFLAGS="-fPIC"
	SHLIB_LD_LIBS='${LIBS}'
	if test "$have_dl" = yes; then
	    SHLIB_LD="${CXX} -shared"
	else
	    AC_CHECK_HEADER(dld.h, [
		SHLIB_LD="ld -shared"])
		DL_LIBS="-ldld"
	fi
    	;;
    Darwin*)
	SHLIB_CFLAGS="-fPIC ${VERSIONFLAG}"
	SHLIB_LD_LIBS='${LIBS}'
	SHLIB_SUFFIX=".dylib"
	CFLAGS="${CFLAGS} ${VERSIONFLAG}"
	CXXFLAGS="${CXXFLAGS} ${VERSIONFLAG}"
	LDFLAGS="${LDFLAGS} ${VERSIONFLAG} \$(RPATH)"
	RPATH="-Wl,-rpath,\$(TLIB)"
	if test "$have_dl" = yes; then
	    SHLIB_LD="${CXX} ${VERSIONFLAG} -dynamiclib -install_name @rpath/\$(LIB)"
	else
	    AC_CHECK_HEADER(dld.h, [
		SHLIB_LD="ld ${VERSIONFLAG} -dynamiclib -install_name @rpath/\$(LIB)"])
		DL_LIBS="-ldld"
	fi
    	;;	
    CYGWIN*)
	SHLIB_CFLAGS=""
	SHLIB_LD_LIBS='-Wl,--out-implib=lib${LIB}.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import -Wl,--whole-archive ${LIBS} -Wl,--no-whole-archive ${DEP_LIBS}'
	SHLIB_SUFFIX=".dll"
	SHLIB_LD="${CXX} -shared"
    	;;
    SunOS-4*)
	SHLIB_CFLAGS="-PIC"
	SHLIB_LD="ld"
	SHLIB_LD_LIBS=""
	;;
    SunOS-5*)
	SHLIB_CFLAGS="-KPIC"
	SHLIB_LD="/usr/ccs/bin/ld -G -z text"
	SHLIB_LD_LIBS='${LIBS}'
	;;
esac

if test "x$OGLLIB" != "x" ; then
    case $system in
	Darwin*)
	    OGLLIB="$OGLLIB -framework GLUT -framework OPENGL"
	    ;;
	*)
	    OGLLIB="$OGLLIB -lGLU -lGL"
	    ;;
    esac
fi

# if we are using GNU compiler, set flags accordingly
#if test "$GXX" = "yes" ; then
#    SHLIB_CFLAGS="-fPIC"
#    SHLIB_LD_LIBS='${LIBS}'
#    SHLIB_SUFFIX=".so"
#    if test "$have_dl" = yes; then
#	SHLIB_LD="${CXX} -shared"
#	DL_LIBS="-ldl"
#    else
#	AC_CHECK_HEADER(dld.h, [
#	    SHLIB_LD="ld -shared"])
#	    DL_LIBS="-ldld"
#    fi
#    case $system in
#	SunOS-4*)
#	    SHLIB_LD="$SHLIB_LD -mimpure-text"
#	    ;;
#	SunOS-5*)
#	    SHLIB_LD="$SHLIB_LD -mimpure-text"
#	    ;;
#    esac
#fi

# =========================================================================
# Check whether user wants to build a debugging version

AC_ARG_ENABLE(debug,
    [  --enable-debug          build debug version],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    TOAST_DEBUG=1
    FFLAGS="-g -DFEM_DEBUG"
    CFLAGS="-g -DFEM_DEBUG"
    CXXFLAGS="-g -DFEM_DEBUG"
#   CXXFLAGS="${CXXFLAGS} -DFEM_DEBUG"
else
    TOAST_DEBUG=0
fi

# =========================================================================
# Check whether user wants to build a profiling version

AC_ARG_ENABLE(prof,
    [  --enable-prof           build profiling version],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    TOAST_PROFILE=1
    if test "$GXX" = "yes" ; then
	CXXFLAGS="-pg"
	LDFLAGS="${LDFLAGS} -pg"
    else
        CXXFLAGS="-p"
	LDFLAGS="${LDFLAGS} -p"
    fi
else
    TOAST_PROFILE=0
fi

# =========================================================================
# Check whether user wants to add Zoltan support

AC_ARG_ENABLE(zoltan,
    [  --enable-zoltan         build with Zoltan load-balancing support],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    ENABLE_ZOLTAN=yes
    NEED_PARMETIS=yes
    NEED_SCOTCH=yes
    NEED_PTHREAD=yes
    NEED_MPI=yes

    AC_ARG_WITH(zoltan-dir,
        [  --with-zoltan-dir       Zoltan source installation directory],
        ZOLTANDIR=$withval, ZOLTANDIR=${ZOLTANDIR})

    AC_ARG_WITH(zoltan-incdir,
        [  --with-zoltan-incdir    Zoltan include directory],
        ZOLTANINCDIR=$withval, ZOLTANINCDIR=${ZOLTANDIR}/include)

    AC_ARG_WITH(zoltan-libdir,
        [  --with-zoltan-libdir    Zoltan library directory],
        ZOLTANLIBDIR=$withval, ZOLTANLIBDIR=${ZOLTANDIR}/lib)
fi

AC_SUBST(ENABLE_ZOLTAN)
AC_SUBST(ZOLTANINCDIR)
AC_SUBST(ZOLTANLIBDIR)

# =========================================================================
# Check whether user wants to build a threaded version

AC_ARG_ENABLE(threads,
    [  --enable-threads        build pthreads parallel version],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    CXXFLAGS="${CXXFLAGS} -DTOAST_THREAD"
    MATLAB_TOASTFLAGS="-DTOAST_THREAD"
    NEED_PTHREAD=yes
fi

# =========================================================================
# Check whether user wants to build an MPI version

AC_ARG_ENABLE(mpi,
    [  --enable-mpi            build MPI distributed version],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    NEED_MPI=yes
    CXXFLAGS="${CXXFLAGS} -DTOAST_MPI"
fi

# =========================================================================
# If pthreads support is required, set the pthreads environment

if test "$NEED_PTHREAD" = "yes" ; then
    ENABLE_PTHREAD=yes
    TASKLIB="-lpthread"
fi
AC_SUBST(ENABLE_PTHREAD)

# =========================================================================
# If MPI support is required, set the MPI environment

if test "$NEED_MPI" = "yes" ; then
    ENABLE_MPI=yes
    CXX="mpicxx"
    F77="mpif77"
# Workaround for MPI-2 SEEK bug:
    CXXFLAGS="${CXXFLAGS} -DMPICH_IGNORE_CXX_SEEK"

    AC_ARG_WITH(mpi-dir,
        [  --with-mpi-dir       MPI source installation directory],
        MPIDIR=$withval, MPIDIR=${MPIDIR})

    AC_ARG_WITH(mpi-incdir,
        [  --with-mpi-incdir    MPI include directory],
        MPIINCDIR=$withval, MPIINCDIR=${MPIDIR}/include)

    AC_ARG_WITH(mpi-libdir,
        [  --with-mpi-libdir    MPI library directory],
        MPILIBDIR=$withval, MPILIBDIR=${MPIDIR}/lib)
fi

AC_SUBST(ENABLE_MPI)
AC_SUBST(MPIINCDIR)
AC_SUBST(MPILIBDIR)

# =========================================================================
# If ParMETIS support is required, set the environment

if test "$NEED_PARMETIS" = "yes" ; then
    ENABLE_PARMETIS=yes

    AC_ARG_WITH(parmetis-dir,
        [  --with-parmetis-dir  ParMETIS source installation directory],
        PARMETISDIR=$withval, PARMETISDIR=${PARMETISDIR})

    AC_ARG_WITH(parmetis-incdir,
        [  --with-parmetis-incdir  ParMETIS include directory],
        PARMETISINCDIR=$withval, PARMETISINCDIR=${PARMETISDIR})

    AC_ARG_WITH(parmetis-libdir,
        [  --with-parmetis-libdir  ParMETIS library directory],
        PARMETISLIBDIR=$withval, PARMETISLIBDIR=${PARMETISDIR})
fi

AC_SUBST(ENABLE_PARMETIS)
AC_SUBST(PARMETISINCDIR)
AC_SUBST(PARMETISLIBDIR)

# =========================================================================
# if Scotch support is required, set the environment

if test "$NEED_SCOTCH" = "yes" ; then
    ENABLE_SCOTCH=yes

    AC_ARG_WITH(scotch-dir,
        [  --with-scotch-dir  Scotch source installation directory],
        SCOTCHDIR=$withval, SCOTCHDIR=${SCOTCHDIR})

    AC_ARG_WITH(scotch-incdir,
        [  --with-scotch-incdir  Scotch include directory],
        SCOTCHINCDIR=$withval, SCOTCHINCDIR=${SCOTCHDIR}/include)

    AC_ARG_WITH(scotch-libdir,
        [  --with-scotch-libdir  Scotch library directory],
        SCOTCHLIBDIR=$withval, SCOTCHLIBDIR=${SCOTCHDIR}/lib)
fi

AC_SUBST(ENABLE_SCOTCH)
AC_SUBST(SCOTCHINCDIR)
AC_SUBST(SCOTCHLIBDIR)

# =========================================================================
# Check whether user wants CUDA support

AC_ARG_ENABLE(cuda,
    [  --enable-cuda           enable CUDA GPU support],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    ENABLE_CUDA="yes"
    CXXFLAGS="${CXXFLAGS} -DUSE_CUDA_FLOAT"
fi

# ========================================================================
# Check for CUDA installation directory

AC_ARG_WITH(cuda,
    [  --with-cuda=PATH        Path to CUDA base directory]
    [                          (implies --enable-cuda). If not specified,]
    [                          CUDA_INSTALL_PATH environment variable is used],
    CUDA_INSTALL_PATH=$withval, CUDA_INSTALL_PATH=${CUDA_INSTALL_PATH})

# =========================================================================
# Check whether user wants to disable shared libraries

TOAST_SHARED="--enable-shared"

AC_ARG_ENABLE(shared,
    [  --disable-shared        do not build shared TOAST libraries],
    [isok=yes], [isok=no])
#if test "$isok" = "no" -a "${SHLIB_SUFFIX}" != "" -a "${TOAST_DEBUG}" = 0; then
if test "$isok" = "no" -a "${SHLIB_SUFFIX}" != ""; then
    TOAST_SHARED_BUILD=1
    LIB_SUFFIX=${SHLIB_SUFFIX}
    MAKE_LIB="${SHLIB_LD} -o \${LIB} \${OBJS} ${SHLIB_LD_LIBS} \${FLIBS} -L\${TLIB}"
    RANLIB=":"
# ranlib doesn't work on dynamic libraries
else
    TOAST_SHARED_BUILD=0
    LIB_SUFFIX=${STLIB_SUFFIX}
    MAKE_LIB="ar cr \${LIB} \${OBJS}"
#    SHLIB_CFLAGS=""
    DL_LIBS=""
    TOAST_SHARED="--disable-shared"
fi

# Set static flag if user explicitly wants shared libraries disabled
if test "$isok" = "yes"; then
    LDFLAGS="${LDFLAGS} -static"
fi

# =========================================================================
# Check whether user wants to apply vectorizer

AC_ARG_ENABLE(vectorizer,
    [  --enable-vectorizer     enable vectorizer flags],
    [isok=$enableval], [isok=no])
if test "$isok" = "yes" ; then
    CFLAGS="${CFLAGS} -ftree-vectorize -msse2 -ffast-math -ftree-vectorizer-verbose=2"
    CXXFLAGS="${CXXFLAGS} -ftree-vectorize -msse2 -ffast-math -ftree-vectorizer-verbose=2"
    FFLAGS="${FFLAGS} -ftree-vectorize -msse2 -ffast-math -ftree-vectorizer-verbose=2"
fi

# =========================================================================
# Check if BLAS should be built static

AC_ARG_ENABLE(blas_shared,
    [  --disable-blas-shared   do not build shared BLAS libraries],
    [isok=yes], [isok=no])
if test "$isok" = "no"; then
   BLAS_LIB_SUFFIX=${LIB_SUFFIX}
else
   BLAS_LIB_SUFFIX=${STLIB_SUFFIX}
fi

# =========================================================================
# Check if LAPACK should be built static

AC_ARG_ENABLE(lapack_shared,
    [  --disable-lapack-shared do not build shared LAPACK libraries],
    [isok=yes], [isok=no])
if test "$isok" = "no"; then
   LAPACK_LIB_SUFFIX=${LIB_SUFFIX}
else
   LAPACK_LIB_SUFFIX=${STLIB_SUFFIX}
fi

# =========================================================================
# Check if ARPACK should be built static

AC_ARG_ENABLE(arpack_shared,
    [  --disable-arpack-shared do not build shared ARPACK libraries],
    [isok=yes], [isok=no])
if test "$isok" = "no"; then
   ARPACK_LIB_SUFFIX=${LIB_SUFFIX}
else
   ARPACK_LIB_SUFFIX=${STLIB_SUFFIX}
fi

# =========================================================================
# Check whether we can strip the symbol tables

if test "$TOAST_DEBUG" = "0" -a "$TOAST_PROFILE" = "0" ; then
    LDFLAGS="${LDFLAGS} -s"
fi

# =========================================================================
# Treatment of templates.

TEMPLATEFLAGS="-DNEED_EXPLICIT_INSTANTIATION"

# Turn of implicit instantiation in gcc
#if test "${GXX}" = "yes" ; then
#    TEMPLATEFLAGS="${TEMPLATEFLAGS} -fno-implicit-templates"
#fi

# Add to the compiler flags
CXXFLAGS="${CXXFLAGS} ${TEMPLATEFLAGS}"

#--------------------------------------------------------------------
#	On a few very rare systems, all of the libm.a stuff is
#	already in libc.a.  Set compiler flags accordingly.
#	Also, Linux requires the "ieee" library for math to work
#	right (and it must appear before "-lm").
#--------------------------------------------------------------------

AC_CHECK_FUNC(sin, MATH_LIBS="", MATH_LIBS="-lm")
AC_CHECK_LIB(ieee, main, [MATH_LIBS="-lieee $MATH_LIBS"])

# =========================================================================
# avoid "deprecated header" warnings in gcc
# =========================================================================

if test "${GXX}" = "yes" ; then
    CXXFLAGS="${CXXFLAGS} -Wno-deprecated"
fi


# =========================================================================
# M.S.13.2.07: Check for CC version (don't know if this would work with
# non-gcc compilers)
# =========================================================================

CCVERSION=gcc`$CC -v 2>&1 | sed -n 's/gcc version //p' | sed -n 's/ .*//p'`
AC_SUBST(CCVERSION)

if test "${TOAST_SHARED_BUILD}" = "1" ; then
    OBJDIR="${CCVERSION}_shared"
else
    OBJDIR="${CCVERSION}_static"
fi

if test "${ENABLE_PTHREAD}" = "yes" ; then
    OBJDIR="${OBJDIR}_smp"
fi

if test "${ENABLE_CUDA}" = "yes" ; then
    OBJDIR="${OBJDIR}_cuda"
fi

if test "${TOAST_DEBUG}" = "1" ; then
    OBJDIR="${OBJDIR}_dbg"
fi

if test "$system" = "Darwin" ; then
    OBJDIR="darwin64"
fi

AC_SUBST(OBJDIR)

# =========================================================================
# Check for mex compiler flags
# =========================================================================

MEXFLAGS="-v"
MEXEXT="mexglx"
case $ARCHDIR in
     linux64)
        MEXFLAGS="${MEXFLAGS} -largeArrayDims"
	MEXEXT="mexa64"
	;;
     darwin32)
        MEXEXT="mexmaci"
        ;;
     darwin64)
        MEXFLAGS="${MEXFLAGS} -largeArrayDims"
        MEXEXT="mexmaci64"
        ;;
esac	
if test "${TOAST_DEBUG}" = "1" ; then
   MEXFLAGS="-g ${MEXFLAGS}"
fi

AC_SUBST(MEXFLAGS)
AC_SUBST(MEXEXT)

# =========================================================================
# Check presence of required directories
#if ! test -d "lib" ; then
#    mkdir lib
#    echo "Created directory \$TOASTDIR\lib"
#fi
#if ! test -d "bin" ; then
#    mkdir bin
#    echo "Created directory \$TOASTDIR\bin"
#fi
if ! test -d "include" ; then
    mkdir include
    echo "Created directory \$TOASTDIR\include"
fi

# =========================================================================
# Substitute variables

AC_SUBST(TOAST_SHARED_BUILD)
AC_SUBST(SHLIB_LD)
AC_SUBST(SHLIB_CFLAGS)
AC_SUBST(SHLIB_LD_LIBS)
AC_SUBST(MAKE_LIB)
AC_SUBST(LIB_SUFFIX)
AC_SUBST(SHLIB_SUFFIX)
AC_SUBST(STLIB_SUFFIX)
AC_SUBST(BLAS_LIB_SUFFIX)
AC_SUBST(LAPACK_LIB_SUFFIX)
AC_SUBST(ARPACK_LIB_SUFFIX)
AC_SUBST(MATH_LIBS)
AC_SUBST(DL_LIBS)
AC_SUBST(BLASLIB)
AC_SUBST(LAPACKLIB)
AC_SUBST(SUPERLULIB)
AC_SUBST(BUILD_SUPERLU)
AC_SUBST(SUPERLU_LOCAL_STATIC)
AC_SUBST(TCLPATH)
AC_SUBST(X_CFLAGS)
AC_SUBST(XVIEWINC)
AC_SUBST(GUIFLAGS)
AC_SUBST(GUILIBS)
AC_SUBST(TASKLIB)
#AC_SUBST(TASKDIR)
AC_SUBST(TOAST_SHARED)
AC_SUBST(ILUPACKDIR)
AC_SUBST(RPATH)
AC_SUBST(TOAST_FDOT)
AC_SUBST(OGLLIB)
AC_SUBST(OGLINC)
AC_SUBST(ENABLE_CUDA)
AC_SUBST(CUDA_INSTALL_PATH)
AC_SUBST(MATLAB_TOASTFLAGS)

AC_OUTPUT(Makefile.incl mexopts.incl)
#AC_OUTPUT(Makefile.incl mexopts.incl gui/Makefile gui/Togl/Makefile)

# =========================================================================
# set up directories and symbolic links

echo
if test "${TOAST_SHARED_BUILD}" = "1" ; then
    echo "Using shared library configuration"
else
    echo "Using static library configuration"
fi

# ========================================================================
# Create build directory structure

echo "Setting build directory ${OBJDIR}"
mkdir -p ${OBJDIR}
mkdir -p ${OBJDIR}/bin
mkdir -p ${OBJDIR}/lib
rm -f bin
ln -s ${OBJDIR}/bin bin
rm -f lib
ln -s ${OBJDIR}/lib lib
rm -f ${ARCHDIR}
ln -s ${OBJDIR} ${ARCHDIR}

# ========================================================================
# Prepare TOAST environment file

TOASTENV=${OBJDIR}/toast.env
echo "#!csh" > ${TOASTENV}
echo "# TOAST environment variables generated with configure" >> ${TOASTENV}
echo "# for target ${OBJDIR}" >> ${TOASTENV}
echo "setenv TOASTDIR ${TOASTDIR}" >> ${TOASTENV}
echo "setenv TOASTVER \$TOASTDIR/${OBJDIR}" >> ${TOASTENV}
echo "setenv TOAST_SCRIPT_PATH \$TOASTDIR/script" >> ${TOASTENV}
echo "setenv LD_LIBRARY_PATH \$TOASTVER/lib" >> ${TOASTENV}

##rm -f toast.env
##ln -s ${TOASTENV} toast.env
##ln -f -s ${TOASTDIR}/numerics/cuda/cusp_v0.2.0 ${TOASTDIR}/numerics/cuda/cusp

# write out the distribution configuration shell scripts
AC_SUBST(ARCHDIR)
AC_OUTPUT(toastenv.sh toastenv.csh)

# move toast-mex options file to target directory
##mv -f mexopts.incl ${TOASTDIR}/${OBJDIR}

echo
echo "TOAST environment (csh) stored in ${TOASTENV}"
echo "Please type 'source toast.env' to activate."
