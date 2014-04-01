      LOGICAL FUNCTION TBILLC (JT,EPS,THETA)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBILLC checks whether the eigenvalue problem is ill conditioned  *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    EPS   (sri) : roundoff unit                                      *
C*    THETA (ari) : eigenvalues of the block tridiagonal matrix        *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SSRANG                                                           *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MAX                                                              *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      DOUBLE PRECISION EIGHT,ONE,TEN
      PARAMETER        (EIGHT=8.0D0,ONE=1.0D0,TEN=10.0D0)
C
C==== arguments ========================================================
C
      INTEGER          JT
      DOUBLE PRECISION EPS,THETA(JT)
C
C==== local variables ==================================================
C
      DOUBLE PRECISION CONDTB,THTMAX,THTMIN
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SSRANG 
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MAX
C
C**** executable statements ********************************************
C
      THTMIN = MAX(EPS,SSRANG(JT,THETA,'MIN'))
      THTMAX = MAX(EPS,SSRANG(JT,THETA,'MAX'))
C
      CONDTB = THTMAX/THTMIN 
C
      TBILLC = CONDTB*EPS*TEN.GE.ONE .AND. THTMAX.GE.TEN**EIGHT
C
C**** end of TBILLC ****************************************************
C
      END
