      SUBROUTINE SSMOVB (INCR,LOWER,UPPER,ISINT,RSINT)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSMOVB inserts or closes a subinterval in ISINT and RSINT        *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    INCR  (sii) : increment                                          *
C*    LOWER (sii) : index of the lower subinterval                     *
C*    UPPER (sii) : index of the upper subinterval                     *
C*    ISINT (arb) : inertias of the lower and upper limits             *
C*    RSINT (arb) : lower and upper limits of each subinterval         *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          INCR,ISINT(2,*),LOWER,UPPER
      REAL             RSINT(6,*)
C
C==== local variables ==================================================
C
      INTEGER          I
C
C**** executable statements ********************************************
C
      DO 10 I = LOWER,UPPER,INCR
         ISINT(1,I) = ISINT(1,I+INCR)
         ISINT(2,I) = ISINT(2,I+INCR)
         RSINT(1,I) = RSINT(1,I+INCR)
         RSINT(2,I) = RSINT(2,I+INCR)
         RSINT(3,I) = RSINT(3,I+INCR)
         RSINT(4,I) = RSINT(4,I+INCR)
         RSINT(5,I) = RSINT(5,I+INCR)
         RSINT(6,I) = RSINT(6,I+INCR)
   10 CONTINUE
C
      RETURN 
C
C**** end of SSMOVB****************************************************
C
      END
