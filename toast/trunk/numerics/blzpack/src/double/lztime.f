      SUBROUTINE LZTIME (LFILE,TIME)           
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZTIME prints the CPU timing for the eigenanalysis               *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE (sii) : file unit for output                               *
C*    TIME  (ari) : time table                                         *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SITIME                                                           *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      DOUBLE PRECISION HNDRD,ZERO
      PARAMETER        (HNDRD=100.0D0,ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          LFILE
      DOUBLE PRECISION TIME(*)  
C
C==== local variables ==================================================
C
      INTEGER          I
      DOUBLE PRECISION TP(7),TTOTAL
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
C
C**** executable statements ********************************************
C
      TTOTAL = SITIME(TIME(13))
C
      IF ( TTOTAL .EQ. ZERO ) THEN 
         TP(1) = ZERO
         TP(2) = ZERO 
         TP(3) = ZERO 
         TP(4) = ZERO 
         TP(5) = ZERO 
         TP(6) = ZERO
         TP(7) = ZERO 
      ELSE
         TP(1) = (TIME(1)/TTOTAL)*HNDRD
         TP(2) = (TIME(2)/TTOTAL)*HNDRD
         TP(3) = (TIME(3)/TTOTAL)*HNDRD
         TP(4) = (TIME(4)/TTOTAL)*HNDRD
         TP(5) = (TIME(5)/TTOTAL)*HNDRD
         TP(6) = (TIME(6)/TTOTAL)*HNDRD
         TP(7) = (TIME(7)/TTOTAL)*HNDRD
      END IF
C
      WRITE (LFILE,1000) (TP(I),TIME(I),I=1,7),TTOTAL
C
      TIME(13) = TTOTAL
C
      RETURN
C
 1000 FORMAT (/,'time table',/,44('-'),
     &        /,'op(A)*vectors .......... (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'op(B)*vectors .......... (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'factorizations ......... (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'vectors generation ..... (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'reorthogonalizations ... (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'reduced problem ........ (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'Ritz vectors ........... (',0P,F5.1,'%) =',1P,E9.2,
     &        /,'total time ...................... =',      1P,E9.2,
     &        /,44('-'))
C
C**** end of LZTIME ****************************************************
C
      END
