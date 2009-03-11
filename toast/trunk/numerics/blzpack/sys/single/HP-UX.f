      REAL             FUNCTION SIRAND (SEED)
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
C*  - Intrinsic Function :                                             *
C*                                                                     *
C*    INT                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          SEED
C
C==== local variables ==================================================
C 
      REAL             ONEMEP,ONEPEP,S,SCALED,SCALEU,T,XM
      SAVE             ONEMEP,ONEPEP,S,SCALED,SCALEU,T,XM
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        INT
C
C==== initialization of variables ======================================
C
      DATA             SCALEU,XM /2147483648.0,16807.0/
C
C**** executable statements ********************************************
C
      IF ( SEED .GT. 0 ) THEN
         SCALED = 1.0/SCALEU
         ONEPEP = 1.0+SCALED
         ONEMEP = 1.0-SCALED
         T = SEED*SCALED
      END IF
C
      T = T*XM
      S = INT(T*ONEPEP)
      T = T - S*ONEMEP
C
      SIRAND = T*ONEPEP
C
      RETURN
C
C**** end of SIRAND ****************************************************
C
      END
      REAL             FUNCTION SITIME (TIME)
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
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    REAL                                                             *
C*                                                                     *
C*  - External Function :                                              *
C*                                                                     *
C*    ETIME                                                            *
C*                                                                     *
C***********************************************************************
C
C==== argument =========================================================
C
      REAL             TIME
C
C==== local variables ==================================================
C 
      REAL             TARRAY(2)
      REAL             TTOTAL
C
C==== intrinsic function ===============================================
C
      INTRINSIC        REAL
C
C==== external function ================================================
C 
      REAL             ETIME
C
C**** executable statements ********************************************
C
      TTOTAL = REAL(ETIME(TARRAY))
      SITIME = REAL(TARRAY(1))
C
      SITIME = SITIME - TIME
C
      RETURN
C
C**** end of SITIME ****************************************************
C
      END
