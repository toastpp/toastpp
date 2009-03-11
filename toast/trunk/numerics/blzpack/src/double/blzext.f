      SUBROUTINE BLZEXT (ISTOR,RSTOR,LFLAG)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZEXT records the exit time and flags                           *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aib) : array for integer variables                        *
C*    RSTOR (aib) : array for real variables                           *
C*    LFLAG (sii) : exit flag                                          *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SITIME                                                           *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    ABS                                                              *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*),LFLAG
      DOUBLE PRECISION RSTOR(*)
C
C==== local variable ===================================================
C
      INTEGER          TIME
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
C
C==== intrinsic function ===============================================
C
      INTRINSIC        ABS
C
C**** executable statements ********************************************
C
      ISTOR(EXIT0) = ABS(LFLAG)
C
      IF ( LFLAG .GT. 0 ) THEN
         TIME = ISTOR(IINIT+ITIME-1)
         RSTOR(RINIT+TIME+9) = SITIME(ZERO)
      END IF
C
      RETURN
C
C**** end of BLZEXT ****************************************************
C
      END
