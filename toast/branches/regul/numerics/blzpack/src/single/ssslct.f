      SUBROUTINE SSSLCT (KEPTLS,KEPTRS,LRMDE,NREIGL,NREIGR,NSLXI,
     &                   NSRXI,NTEIG,INDSI,NSINT,RSINT,ISINT,
     &                   EIG,ORIGIN,SIGMA)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSSLCT selects a subinterval to be examined                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    KEPTLS (sii) : index of the current left  subinterval            *
C*    KEPTRS (sii) : index of the current rigth subinterval            *
C*    LRMDE  (sii) : run mode                                          *
C*    NREIGL (sii) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sii) : number of required eigenvalues greater than EIGL  *
C*    NSLXI  (sii) : number of eigenvalues less    than EIGR           *
C*    NSRXI  (sii) : number of eigenvalues greater than EIGL           *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    INDSI  (sio) : index of the subinterval to be closed             *
C*    NSINT  (sii) : number of subintervals                            *
C*    RSINT  (ari) : lower and upper limits of each subinterval        *
C*    ISINT  (ari) : inertias of the lower and upper limits            *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    ORIGIN (sri) : starting point (first SIGMA)                      *
C*    SIGMA  (sri) : origin translation                                *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    NEIGAB                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          ISINT(2,*),INDSI,KEPTLS,KEPTRS,LRMDE,NREIGL,
     &                 NREIGR,NSINT,NSLXI,NSRXI,NTEIG
      REAL             EIG(*),ORIGIN,SIGMA,RSINT(6,*)
C
C==== local variables ==================================================
C
      INTEGER          I,INDXL,INDXR,J,LOSTL,LOSTR,NERNG,NSRNG
C
C==== subprogram =======================================================
C
      INTEGER          NEIGAB
C
C**** executable statements ********************************************
C
      IF ( LRMDE .EQ. 0 ) RETURN
C
C.... select a subinterval .............................................
C
      IF      ( NSINT .EQ. 1 ) THEN
C
C............ trivial case .............................................
C
              INDSI = 1
C
      ELSE IF ( (SIGMA.LT.ORIGIN) .AND. (KEPTRS.NE.0) ) THEN
C
C............ target subinterval was not closed (right) ................
C
              INDSI = KEPTRS
C
      ELSE IF ( (SIGMA.GT.ORIGIN) .AND. (KEPTLS.NE.0) ) THEN
C
C............ target subinterval was not closed (left) .................
C
              INDSI = KEPTLS
C
      ELSE
C
C............ find the ideal subinterval ...............................
C
              INDXL = 0
              INDXR = 0
C
              DO 10 I = 1,NSINT
                 J = NSINT - I + 1
                 IF ( RSINT(6,I).LE.ORIGIN ) INDXL = I
                 IF ( RSINT(1,J).GE.ORIGIN ) INDXR = J
   10         CONTINUE
C
              IF      ( INDXR .EQ. 0 ) THEN
                      INDSI = INDXL
              ELSE IF ( INDXL .EQ. 0 ) THEN
                      INDSI = INDXR
              ELSE IF ( INDXL.EQ.1 .AND. INDXR.EQ.2 ) THEN
                      IF      ((NREIGL-NSLXI).GE.(NREIGR-NSRXI)) THEN
                              INDSI = 1 
                      ELSE IF ((NREIGR-NSRXI).GE.(NREIGL-NSLXI)) THEN
                              INDSI = 2
                      END IF
              ELSE
                      J = INDXL
                      NSRNG = NEIGAB(NTEIG,RSINT(1,J),RSINT(5,J),EIG)
                      NERNG = ISINT(2,J) - ISINT(1,J)
                      LOSTL = NERNG - NSRNG
                      J = INDXR
                      NSRNG = NEIGAB(NTEIG,RSINT(1,J),RSINT(5,J),EIG)
                      NERNG = ISINT(2,J) - ISINT(1,J)
                      LOSTR = NERNG - NSRNG
                      IF ( LOSTL .GT. LOSTR ) THEN
                         INDSI = INDXL
                      ELSE
                         INDSI = INDXR
                      END IF
              END IF
C
      END IF
C
      RETURN 
C
C**** end of SSSLCT ****************************************************
C
      END
