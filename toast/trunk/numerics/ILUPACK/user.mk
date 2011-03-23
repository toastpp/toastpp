# 1. step. 
# which platform do we use?
# ILUPACK comes along with some standard settings for various compilers and
# platforms. Often enough it suffices to choose one from the list below.
# In case you do not find an acceptable configuration, you can change the
# associated "makefile.include*" configuration in "makefiles/"
# The ILUPACK libraries and some dependencies will be provided in the
# associated sub directory of "libs/"

# GNU-compiler-based options

# 32 BIT gcc/gfortran linux system
# PLATFORM=GNU32

# 64 BIT gcc/gfortran linux system
# PLATFORM=GNU64

# 64 BIT gcc/gfortran linux system with 64 bit integer
 PLATFORM=GNU64_long


# Intel-compiler-based options

# 32 BIT icc/ifort linux system
# PLATFORM=INTEL32

# 64 BIT icc/ifort linux system
# PLATFORM=INTEL64

# 64 BIT icc/ifort linux system with 64 bit integer
# PLATFORM=INTEL64_long


# IBM-AIX
# 64 BIT xlc/xlf aix system
# PLATFORM=aix



# ------------------------------------------


# 2. step 
# which main driver to we want for testing?
# ILUPACK offers you a suite of some sample drivers which show how ILUPACK
# is applied. You find these drivers in "samples/". For testing choose one
# from the following list. Suppose you use e.g. ILUPACK for general real
# matrices. Then select "maingnl". After compiling you can call 
# "dmaingnl <your drop tol> <your condest> <your elbow suggestion> <your matrix>",
# e.g., "dmaingnl 0.01 10 10 lnsp3937.rua".
# If you plan to use ILUPACK, e.g., for complex symmetric matrices, choose
# the complex symmetric driver mainsyms and you can call,e.g.,
# "zmainsyms 0.01 10 10 mplate.csa"


# default general unsymmetric solver
MAIN=maingnl

# default symmetric (Hermitian) positive definite solver
#MAIN=mainspd

# default real symmetric (complex Hermitian) indefinite solver
#MAIN=mainsym

# default complex symmetric indefinite solver
#MAIN=mainsyms

# default real symmetric (complex Hermitian) indefinite solver converted
# to positive definite preconditioner
#MAIN=mainsymspd




# -----------------------------------------------------
# some special solvers

# almost symmetric (Hermitian) matrices preconditioned by its symmetric part 
#MAIN=maingnlsym

# general unsymmetric solver using the transposed matrix
#MAIN=maint

# old-fashioned ILUTP solver from SPARSKIT (revised version for ILUPACK)
#MAIN=mainilutp

# old-fashioned ILUT  solver from SPARSKIT (revised version for ILUPACK)
#MAIN=mainilut

# general inverse-based ILU (one-shot method)
#MAIN=mainiluc

# symmetric positive definite inverse-based ILU (one-shot method)
#MAIN=mainildlc

# symmetric (Hermitian) indefinite inverse-based ILU (one-shot method)
#MAIN=mainsymiluc

# complex symmetric indefinite inverse-based ILU (one-shot method)
#MAIN=mainsymilucs

# solvers that allow use of multiple factorizations (e.g. nonlinear or
# or varying systems)
#MAIN=main3
#MAIN=mainspd3
#MAIN=mainsym3
#MAIN=mainsym3s


# example on Helmholtz problem (real/complex shifts)
#MAIN=helmholtz

# -----------------------------------------------------
# some fortran solvers
# general unsymmetric system
#MAIN=fmain

# real symmetric / complex Hermitian system
#MAIN=fmainsym

# complex symmetric system
#MAIN=fmainsyms

# real symmetric / complex Hermitian positive definite system
#MAIN=fmainspd





