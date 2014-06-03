      SUBROUTINE RANDNR (LCOMM,LRERR,NI,NR,NRUN,NVB,R)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RANDNR generates random starting vectors                         *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    NI    (sii) : dimension of the vectors in (R)                    *
C*    NR    (sii) : number of starting vectors available               *
C*    NRUN  (sii) : number of runs                                     *
C*    NVB   (sii) : number of vectors in a block                       *
C*    R     (arb) : starting vectors                                   *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    PIIPID,SETTO0,SIRAND                                             *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LCOMM,LRERR,NI,NR,NRUN,NVB
      DOUBLE PRECISION R(NI,NVB)
C
C==== local variables ==================================================
C
      INTEGER          I,J,PID
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SIRAND
C
C**** executable statements ********************************************
C
      IF ( NR .LT. NVB ) THEN
C
C....... initialize vectors ............................................
C
         CALL SETTO0 (NI*(NVB-NR),R(1,NR+1),1)
C
C....... initialize function, R(1,1) will be overwritten ...............
C
         CALL PIIPID (PID,LCOMM,LRERR)
C
         R(1,1) = SIRAND(NRUN+NR+PID)
C
C....... generate vectors ..............................................
C
         DO 20 I = NR+1,NVB 
            DO 10 J = 1,NI
               R(J,I) = SIRAND(0)
   10       CONTINUE
   20    CONTINUE
C
      END IF
C
      RETURN 
C
C**** end of RANDNR ****************************************************
C
      END
