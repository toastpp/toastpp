      SUBROUTINE SETDFT (JTMIN,NVB,NRUNMX,NSIMAX,NXMAX,FHNDL,BIGNUM,EPS)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SETDFT sets defaults                                             *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JTMIN  (sio) : minimum basis size before changing SIGMA          *
C*    NVB    (sio) : number of vectors in a block                      *
C*    NRUNMX (sio) : maximum number of runs                            *
C*    NSIMAX (sio) : maximum number of subintervals                    *
C*    NXMAX  (sio) : maximum number of eigenvectors stored in X        *
C*    FHNDL  (aio) : file handle                                       *
C*    BIGNUM (sro) : big number                                        *
C*    EPS    (sro) : roundoff unit                                     *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SETEPS                                                           *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      INCLUDE          'defaults.h'
C 
      DOUBLE PRECISION ONE
      PARAMETER        (ONE=1.0D0)
C   
C==== arguments ========================================================
C
      INTEGER          FHNDL(*),JTMIN,NRUNMX,NSIMAX,NVB,NXMAX
      DOUBLE PRECISION BIGNUM,EPS
C
C**** executable statements ********************************************
C
C.... minimum basis size before changing SIGMA .........................
C
      JTMIN  = MINJT
C
C.... maximum number of subintervals ...................................
C
      NSIMAX = MAXSI
C
C.... maximum number of runs ...........................................
C
      NRUNMX = MAXRUN
C
C.... block size .......................................................
C
      NVB    = NVBD
C
C.... maximum number of eigenvectors stored in (X) .....................
C
      NXMAX  = MAXVX
C 
C.... file handles .....................................................
C
      FHNDL(1) = FUNIT1
      FHNDL(2) = FUNIT2
      FHNDL(3) = FUNIT3
      FHNDL(4) = FUNIT4
C 
C.... machine precision ................................................
C
      CALL SETEPS (EPS)
C
C.... big number .......................................................
C
      BIGNUM = (ONE/EPS)**4
C
      RETURN
C
C**** end of SETDFT ****************************************************
C
      END
