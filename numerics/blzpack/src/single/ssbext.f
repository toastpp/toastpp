      SUBROUTINE SSBEXT (NSINT,RITZL,RITZR,SIGMA,SIGMAL,SIGMAR,RSINT)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSBEXT tries to extend the bounds of the current subinterval     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    NSINT  (sii) : number of subintervals                            *
C*    RITZL  (sri) : farthest Ritz value < SIGMA                       *
C*    RITZR  (sri) : farthest Ritz value > SIGMA                       *
C*    SIGMA  (sri) : origin translation                                *
C*    SIGMAL (sri) : origin translation  < SIGMA                       *
C*    SIGMAR (sri) : origin translation  > SIGMA                       *
C*    RSINT  (arb) : lower and upper limits of each subinterval        *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    MAX,MIN                                                          *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          NSINT
      REAL             RITZL,RITZR,SIGMA,SIGMAL,SIGMAR,RSINT(6,*)
C
C==== local variables ==================================================
C
      INTEGER          I
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        MAX,MIN
C
C**** executable statements ********************************************
C
C.... try to extend the bounds .........................................
C
      DO 10 I = 1,NSINT
C
C....... check whether SIGMA is a  lower bound .........................
C
         IF ( RSINT(1,I) .EQ. SIGMA ) THEN
              RSINT(2,I) = MAX(RSINT(2,I),RITZR )
              RSINT(3,I) = MAX(RSINT(3,I),SIGMAR)
         END IF
C
C....... check whether SIGMA is an upper bound .........................
C
         IF ( RSINT(6,I) .EQ. SIGMA ) THEN
              RSINT(5,I) = MIN(RSINT(5,I),RITZL )
              RSINT(4,I) = MIN(RSINT(4,I),SIGMAL)
         END IF
C
   10 CONTINUE
C
      RETURN 
C
C**** end of SSBEXT ****************************************************
C
      END
