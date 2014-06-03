      REAL FUNCTION SIRAND (SEED)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose :                                                        *
C*                                                                     *
C*    SIRAND is an interface for a random number generator             *
C*                                                                     *
C*  - Arguments :                                                      *
C*                                                                     *
C*    SEED (sii) : flag for the random-number generator                *
C*                 - if not zero, it is used as a new seed             *
C*                 - if zero, a random number is returned              *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    REAL                                                             *
C*                                                                     *
C*  - External Functions :                                             *
C*                                                                     *
C*    RANF,RANSET                                                      *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   SEED
C
C==== intrinsic function ===============================================
C
      INTRINSIC REAL
C
C==== external functions ===============================================
C 
      REAL      RANF,RANSET
C
C**** executable statements ********************************************
C
      IF ( SEED .EQ. 0 ) THEN
         SIRAND = REAL(RANF())
      ELSE
         SIRAND = REAL(RANSET(SEED))
      END IF
C
      RETURN
C
C**** end of SIRAND ****************************************************
C
      END
      REAL FUNCTION SITIME (TIME)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose :                                                        *
C*                                                                     *
C*    SITIME is an interface for a CPU time function                   *
C*                                                                     *
C*  - Argument :                                                       *
C*                                                                     *
C*    TIME (sri) : reference time                                      *
C*                                                                     *
C*  - Intrinsic Function :                                             *
C*                                                                     *
C*    REAL                                                             *
C*                                                                     *
C***********************************************************************
C
C==== argument =========================================================
C
      REAL      TIME
C
C==== local variables ==================================================
C 
      INTEGER   COUNT,CRATE,CMAX
C
C==== intrinsic function ===============================================
C
      INTRINSIC REAL
C
C**** executable statements ********************************************
C
      CALL SYSTEM_CLOCK (COUNT,CRATE,CMAX)
C
      SITIME = REAL(COUNT)/REAL(CRATE) - TIME
C
      RETURN
C
C**** end of SITIME ****************************************************
C
      END
