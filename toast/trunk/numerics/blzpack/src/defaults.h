C***********************************************************************
C*                                                                     *
C*    defaults.h is an include file that sets defaults for BLZPACK, it *
C*    is used by the subroutine setdft.f. Changing the defaults below  *
C*    may require an update of the workspace.                          *
C*                                                                     *
C***********************************************************************
C
      INTEGER    FUNIT1,FUNIT2,FUNIT3,FUNIT4
      INTEGER    MAXRUN,MAXSI,MAXVX,MINJT,NVBD
C
C.... default block size ...............................................
C
      PARAMETER (NVBD   = 3)
C
C.... minimum basis size to be computed before changing SIGMA ..........
C
C     used only in the spectrum slicing strategy
C
      PARAMETER (MINJT  = 30)
C
C.... maximum number of subintervals ...................................
C
C     used only in the spectrum slicing strategy
C
      PARAMETER (MAXSI  = 15)
C
C.... maximum number of runs ...........................................
C
      PARAMETER (MAXRUN = MAXSI*3)
C
C.... maximum number of eigenvectors stored in (X) .....................
C
C     MAXVX > 0: The eigenvectors are stored in the array (X).
C
C        In this case, (X) must be declared with dimensions (LNI,LEIG)
C        in the program unit calling BLZDRD/S. The eigenvalues are
C        sorted in ascending order in (EIG) and the corresponding
C        eigenvectors can be printed upon request.
C
C        The generalized eigenvalue problem may need a temporary file
C        opened with direct access and named `blzpack.u.BX', where u
C        is the file unit.
C
C     MAXVX = 0: The eigenvectors are stored in a file.
C
C        In this case, (X) should be declared with dimensions (LNI,1)
C        in the program unit calling BLZDRD/S. The eigenvectors are
C        stored in a file opened with direct access, record length
C        N*8, where N is the dimension of the eigenvalue problem, and
C        named `blzpack.u.X', where u is the file unit. This file is
C        preserved upon succesful exit from BLZDRD/S. The eigenvalues
C        and eigenvectors will neither be sorted nor printed, the
C        eigenvectors are stored in the file accordingly to the
C        order of the eigenvalues stored in EIG.
C
C        The generalized eigenvalue problem may need a temporary file
C        opened with direct access and named `blzpack.u.BX', where u
C        is the file unit.
C
      PARAMETER (MAXVX  = 100000000)
C
C.... file handles .....................................................
C
C     only file units are currently defined
C
      PARAMETER (FUNIT1 = 50)
      PARAMETER (FUNIT2 = 51)
      PARAMETER (FUNIT3 = 52)
      PARAMETER (FUNIT4 = 53)
