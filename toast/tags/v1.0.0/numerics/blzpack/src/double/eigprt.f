      SUBROUTINE EIGPRT (INDXP,INDXS,NPEIG,NPE,LFILE,
     &                   LPRNT,LEIG,LNI,NI,EIG,X)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    EIGPRT prints the eigenpairs                                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    INDXP (sii) : index for the eigenpairs (printing)                *
C*    INDXS (sii) : index for the eigenpairs (algebraic)               *
C*    NPEIG (sii) : number of eigenpairs to be printed                 *
C*    NPE   (sii) : number of processes                                *
C*    LFILE (sii) : file unit for output                               *
C*    LPRNT (sii) : level of printing                                  *
C*    LEIG  (sii) : leading dimension of (EIG)                         *
C*    LNI   (sii) : leading dimension of (X)                           *
C*    NI    (sio) : dimension of the vectors in (X)                    *
C*    EIG   (ari) : eigenvalue approximations and estimated residuals  *
C*    X     (ari) : eigenvector approximations                         *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          INDXP,INDXS,LEIG,LFILE,
     &                 LNI,LPRNT,NI,NPE,NPEIG
      DOUBLE PRECISION EIG(LEIG,2),X(LNI,*)
C
C==== local variables ==================================================
C
      INTEGER          I,J,K,L
C
C==== subprogram =======================================================
C
      LOGICAL          SIBTST
C
C**** executable statements ********************************************
C
      K = INDXS
      L = INDXP
C
C.... print the eigenvalues ............................................
C
      IF ( SIBTST(2,LPRNT) .AND. NPEIG.NE.0 ) THEN
         WRITE (LFILE,1000) (I+K,EIG(I+L,1),EIG(I+L,2),I=0,NPEIG-1)
         WRITE (LFILE,1001)
      END IF
C
C.... print the eigenvectors ...........................................
C
      IF ( SIBTST(5,LPRNT) .AND. NPEIG.NE.0 .AND. NPE.EQ.1 ) THEN
         DO 10 I = 0,NPEIG-1
            WRITE (LFILE,1002) I+K,EIG(I+L,1),(X(J,I+L),J=1,NI) 
   10    CONTINUE
      END IF
C
      IF ( SIBTST(2,LPRNT) .AND. NPEIG.NE.0 ) WRITE (LFILE,1003)
C
      RETURN
C
 1000 FORMAT (/,'eigenvalues',/,35('='),
     &        /,'vector',8X,'value',7X,'residual',/,35('-'),/,
     &        (I6,5X,1P,SP,E12.4,SS,E12.4))
 1001 FORMAT (35('='))
 1002 FORMAT (/,'vector',I5,7X,'value',1P,SP,E12.4,
     &        /,71('='),/,6(E11.4,1X))
 1003 FORMAT (/,71('*'))
C
C**** end of EIGPRT ****************************************************
C
      END
