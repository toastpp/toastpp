// -*-C++-*-

#ifndef __MATLAB_TOASTMEX_H
#define __MATLAB_TOASTMEX_H

#include "mathlib.h"
#include "matlabtoast.h"
#include "matlabfdot.h"

// Matlab function indices
#define TOAST_READMESH           0
#define TOAST_MAKEMESH           1
#define TOAST_WRITEMESH          2
#define TOAST_WRITEMESHVTK      82
#define TOAST_CLEARMESH          3
#define TOAST_MESHDATA           4
#define TOAST_SURFDATA           5
#define TOAST_ELEMENTNEIGHBOURS 100 // end
#define TOAST_MARKMESHBOUNDARY  60
#define TOAST_MESHNODECOUNT      6
#define TOAST_MESHELEMENTCOUNT   7
#define TOAST_MESHBB             8
#define TOAST_MESHSIZE          59
#define TOAST_MESHDIMENSION      9
#define TOAST_MESHLIN2QUAD      69
#define TOAST_MESHOPT           70
#define TOAST_MESHREORDER       84
#define TOAST_ELEMENTSIZE       10
#define TOAST_ELEMENTDATA       11
#define TOAST_ELEMENTREGION     85
#define TOAST_FINDELEMENT       12
#define TOAST_SHAPEFUNC         57
#define TOAST_SHAPEGRAD         58
#define TOAST_READQM            13
#define TOAST_SETQM             14
#define TOAST_GETQM             15
#define TOAST_WRITEQM           16
#define TOAST_DATALINKLIST      17
#define TOAST_QVEC              18
#define TOAST_MVEC              19
#define TOAST_QPOS              20
#define TOAST_MPOS              21
#define TOAST_READNIM           22
#define TOAST_WRITENIM          56
#define TOAST_SYSMAT            23
#define TOAST_SYSMATCOMPONENT   74
#define TOAST_SYSMATBASIS       86
#define TOAST_MASSMAT           24
#define TOAST_VOLMAT            65
#define TOAST_BNDMAT            61
#define TOAST_ELMAT             25
#define TOAST_ELDOF             75
#define TOAST_BNDREFLECTIONTERM 62
#define TOAST_SETBASIS          26
#define TOAST_CLEARBASIS        27
#define TOAST_GETBASISSIZE      28
#define TOAST_BASIS_NLEN        71
#define TOAST_BASIS_BLEN        72
#define TOAST_BASIS_SLEN        73
#define TOAST_BASIS_VALUE       81
#define TOAST_BASIS_BUU         87
#define TOAST_BASIS_BVV         88
#define TOAST_BASIS_BUV         89
#define TOAST_BASIS_BVW         93
#define TOAST_BASIS_DUU         98
#define TOAST_BASIS_DVV         99
#define TOAST_BASIS_DUV         97
#define TOAST_BASIS_GUV         94
#define TOAST_BASIS_GVV         95
#define TOAST_BASIS_SUPPORTAREA 90
#define TOAST_BASIS_REFINE      91
#define TOAST_MAPBASIS          29
#define TOAST_MAPMESHTOBASIS    30
#define TOAST_MAPMESHTOGRID     31
#define TOAST_MAPMESHTOSOL      32
#define TOAST_MAPBASISTOMESH    33
#define TOAST_MAPSOLTOMESH      34
#define TOAST_MAPSOLTOBASIS     35
#define TOAST_MAPSOLTOGRID      36
#define TOAST_MAPGRIDTOMESH     37
#define TOAST_MAPGRIDTOBASIS    38
#define TOAST_MAPGRIDTOSOL      39
#define TOAST_BASISTOMESHMATRIX 40
#define TOAST_MESHTOBASISMATRIX 68
#define TOAST_BASISSAMPLE       83
#define TOAST_GRIDELREF         63
#define TOAST_SAMPLEFIELD       66
#define TOAST_IMAGEGRADIENT     64
#define TOAST_SOLUTIONMASK      41
#define TOAST_READVECTOR        42
#define TOAST_WRITEVECTOR       43
#define TOAST_REGUL             44
#define TOAST_CLEARREGUL        45
#define TOAST_REGULVALUE        46
#define TOAST_REGULGRADIENT     47
#define TOAST_REGULHDIAG        48
#define TOAST_REGULHESS         49
#define TOAST_REGULHESS1F       67
#define TOAST_REGULKAPPA        50
#define TOAST_REGULSETLOCALSCALING 96
#define TOAST_FIELDS            51
#define TOAST_GRADIENT          52
#define TOAST_JACOBIAN          53
#define TOAST_JACOBIANCW        54
#define TOAST_KRYLOV            55
#define TOAST_LBFGS             76
#define TOAST_MESHREFINE        77
#define TOAST_SPLITELEMENT      78
#define TOAST_NODENEIGHBOUR     79
#define TOAST_UNWRAPPHASE       80
#define TOAST_SPARSITYSTRUCTURE 92

#define TOAST_SETVERBOSITY    1000
#define TOAST_THREADCOUNT     1001

// Added by SP
#define TOAST_INTFG           4000
#define TOAST_INTGRADFGRADG   4001

// Fluorescence solver interface
#define FDOT_MAKEFWD          2007
#define FDOT_CLEARFWD         2002
#define FDOT_FWDOP            2003
#define FDOT_FWDOPRATIO       2004
#define FDOT_ADJOP            2000
#define FDOT_ADJOPRATIO       2001
#define FDOT_EXCIT            2005
#define FDOT_SYSMAT           2006
#define FDOT_MAKEPROJECTORLIST 2008
#define FDOT_PROJECTTOIMAGE   2009
#define FDOT_PROJECTTOFIELD   2010

#endif // !__MATLAB_TOASTMEX_H
