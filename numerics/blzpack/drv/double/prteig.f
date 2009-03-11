      SUBROUTINE PRTEIG (ISTOR,LEIG,LN,NI,EIG,X,W)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    PRTEIG prints the contents of (EIG) and (X)                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (sii) : BLZPACK array of integer variables                 *
C*    LEIG  (sii) : leading dimension of (EIG)                         *
C*    LN    (sii) : leading dimension of (X)                           *
C*    NI    (sii) : dimension of the vectors in (X)                    *
C*    EIG   (ari) : eigenvalues and associated residuals               *
C*    X     (ari) : eigenvectors                                       *
C*    W     (ari) : workspace                                          *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    ISTORR,PIDGTR                                                    *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*),LEIG,LN,NI
      DOUBLE PRECISION EIG(LEIG,*),X(LN,*),W(*)
C
C==== local variables ==================================================
C
      INTEGER          I,INFO,J,LCOMM,MYPE,N,NTEIG
C
C==== subprogram =======================================================
C
      INTEGER          ISTORR
C
C**** executable statements ********************************************
C
      N     = ISTORR(ISTOR,'N')
      MYPE  = ISTORR(ISTOR,'MYPE')
      LCOMM = ISTORR(ISTOR,'LCOMM')
      NTEIG = ISTORR(ISTOR,'NTEIG')
C
C.... store (EIG) in EIG ...............................................
C
      IF ( MYPE .EQ. 0 ) THEN
         OPEN  (10,FORM='formatted',FILE='EIG')
         WRITE (10,1001) (I,EIG(I,1),EIG(I,2),I=1,NTEIG)
         CLOSE (10)
      END IF
C
C.... store (X) in X, W requires N positions ...........................
C
      IF ( MYPE .EQ. 0 ) OPEN  (11,FORM='formatted',FILE='X')
      DO 10 I = 1,NTEIG
         CALL PIDGTR (NI,X(1,I),W,LCOMM,INFO)
         IF ( INFO .NE. 0 ) STOP '* PRTEIG: global gather error'
         IF ( MYPE .EQ. 0 ) WRITE (11,1002) (W(J),J=1,N)
   10 CONTINUE 
      IF ( MYPE .EQ. 0 ) CLOSE (11)
C
      RETURN
 1001 FORMAT (I6,1P,2E16.8)
 1002 FORMAT (1P,E16.8)
C
C**** end of PRTEIG ****************************************************
C
      END
