      DOUBLE PRECISION FUNCTION SIRAND (SEED)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SIRAND is an interface for a random number generator             *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    SEED (sii) : flag for the random-number generator                *
C*                 - if not zero, it is used as a new seed             *
C*                 - if zero, a random number is returned              *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    DBLE                                                             *
C*                                                                     *
C*  - External Function:                                               *
C*                                                                     *
C*    RAND                                                             *
C*                                                                     *
C***********************************************************************
C
C==== argument =========================================================
C
      INTEGER   SEED
C
C==== intrinsic function ===============================================
C
      INTRINSIC DBLE
C
C==== external function ================================================
C 
      REAL      RAND
C
C**** executable statements ********************************************
C
      SIRAND = DBLE(RAND(SEED))
C
      RETURN
C
C**** end of SIRAND ****************************************************
C
      END
      DOUBLE PRECISION FUNCTION SITIME (TIME)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SITIME is an interface for a CPU time function                   *
C*                                                                     *
C*  - Argument:                                                        *
C*                                                                     *
C*    TIME (sri) : reference time                                      *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    DBLE                                                             *
C*                                                                     *
C*  - External Function:                                               *
C*                                                                     *
C*    ETIME                                                            *
C*                                                                     *
C***********************************************************************
C
C==== argument =========================================================
C
      DOUBLE PRECISION TIME
C
C==== local variables ==================================================
C 
      DOUBLE PRECISION TTOTAL
      REAL             TARRAY(2)
C
C==== intrinsic function ===============================================
C
      INTRINSIC        DBLE
C
C==== external function ================================================
C 
      REAL             ETIME
C
C**** executable statements ********************************************
C
      TTOTAL = DBLE(ETIME(TARRAY))
      SITIME = DBLE(TARRAY(1))
C
      SITIME = SITIME - TIME
C
      RETURN
C
C**** end of SITIME ****************************************************
C
      END
