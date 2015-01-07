      SUBROUTINE SETTO0 (N,X,INCX)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SETTO0 fills a vector with zeros                                 *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    N    (sii) : dimension of (X)                                    *
C*    X    (aro) : vector to be filled with zeros                      *
C*    INCX (sii) : increment for the elements of (X)                   *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C 
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          I,INCX,N
      DOUBLE PRECISION X(*)
C
C**** executable statements ********************************************
C
      DO 10 I = 1,N*INCX,INCX
         X(I) = ZERO
   10 CONTINUE
C
      RETURN 
C
C**** end of SETTO0 ****************************************************
C
      END
