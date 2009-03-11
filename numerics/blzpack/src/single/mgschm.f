      SUBROUTINE MGSCHM (LBLAS,LCOMM,LRERR,NI,NVB,BR,R,
     &                   BETA,WORK,GNRZD,NPIVT)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    MGSCHM factors (R) as (Q)*(BETA)                                 *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LBLAS (sii) : BLAS level setting                                 *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    NI    (sii) : dimension of the vectors in (BR) and (R)           *
C*    NVB   (sii) : number of vectors in a block                       *
C*    BR    (arb) : (B)*(R) on input, (B)*(Q) on output                *
C*    R     (arb) : (R) on input, (Q) on output                        *
C*    BETA  (aro) : matrix (BETA) in (R)=(Q)*(BETA) at the j-th step   *
C*    WORK  (arw) : workspace                                          *
C*    GNRZD (sli) : problem type flag                                  *
C*    NPIVT (slo) : negative pivot flag                                *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    PISRED,SETTO0                                                    *
C*                                                                     *
C*  - BLAS kernels:                                                    *
C*                                                                     *
C*    SAXPY,SDOT,SGEMV,SGER,SSCAL                                      *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    SQRT                                                             *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      REAL             ONE,ZERO
      PARAMETER        (ONE=1.0E0,ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          LBLAS,LCOMM,LRERR,NI,NVB
      REAL             BETA(NVB,NVB),BR(NI,NVB),R(NI,NVB),WORK(*)
      LOGICAL          GNRZD,NPIVT
C
C==== local variables ==================================================
C
      INTEGER          I,INFO,J
      REAL             DOTP,ZETA
C
C==== BLAS kernel ======================================================
C
      REAL             SDOT
C
C==== intrinsic function ===============================================
C
      INTRINSIC        SQRT
C
C**** executable statements ********************************************
C
      NPIVT = .FALSE.
C
      CALL SETTO0 (NVB*NVB,BETA,1)
C
      DO 30 I = 1,NVB
C
C....... compute the norm of the vector ................................
C
         IF ( GNRZD ) THEN
            DOTP = SDOT(NI,R(1,I),1,BR(1,I),1)
         ELSE
            DOTP = SDOT(NI,R(1,I),1, R(1,I),1)
         END IF
C
         CALL PISRED ('SUM',1,DOTP,ZETA,LCOMM,INFO)
         IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
         IF ( ZETA .LT. ZERO ) THEN
C
C.......... the matrix (B) seems to be indefinite ......................
C
            NPIVT = .TRUE.
            RETURN
C
         ELSE 
C
            ZETA = SQRT(ZETA)
            BETA(I,I) = ZETA
C
            IF ( ZETA .GT. ZERO ) THEN
C
C............. the vector can be normalized ............................
C
               ZETA = ONE/ZETA 
C
               IF ( GNRZD ) THEN
                  CALL SSCAL (NI,ZETA,BR(1,I),1)
                  CALL SSCAL (NI,ZETA, R(1,I),1)
               ELSE
                  CALL SSCAL (NI,ZETA, R(1,I),1)
               END IF
C
               IF      ( LBLAS .EQ. 1 ) THEN
C
C..................... modified Gram-Schmidt using BLAS 1 ..............
C
                       IF ( I .LT. NVB ) THEN
                          DO 10 J = 1,NVB-I
                             IF ( GNRZD ) THEN
                                WORK(J) = SDOT(NI,R(1,I),1,BR(1,I+J),1)
                             ELSE
                                WORK(J) = SDOT(NI,R(1,I),1, R(1,I+J),1)
                             END IF
   10                     CONTINUE
                          CALL PISRED ('SUM',NVB-I,WORK,BETA(I+1,I),
     &                                 LCOMM,INFO)
                          IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
                          DO 20 J = I+1,NVB 
                             IF ( GNRZD ) THEN
                                CALL SAXPY (NI,-BETA(J,I),BR(1,I),
     &                                      1,BR(1,J),1)
                                CALL SAXPY (NI,-BETA(J,I), R(1,I),
     &                                      1, R(1,J),1)
                             ELSE
                                CALL SAXPY (NI,-BETA(J,I), R(1,I),
     &                                      1, R(1,J),1)
                             END IF
   20                     CONTINUE
	               END IF
C
               ELSE IF ( LBLAS .GE. 2 ) THEN
C
C..................... modified Gram-Schmidt using BLAS 2 ..............
C
                       IF ( I .LT. NVB ) THEN
                          IF ( GNRZD ) THEN
	                     CALL SGEMV  ('T',NI,NVB-I,ONE,BR(1,I+1),NI,
     &                                    R(1,I),1,ZERO,WORK,1)
                          ELSE
	                     CALL SGEMV  ('T',NI,NVB-I,ONE, R(1,I+1),NI,
     &                                    R(1,I),1,ZERO,WORK,1)
                          END IF
                          CALL PISRED ('SUM',NVB-I,WORK,BETA(I+1,I),
     &                                 LCOMM,INFO)
                          IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
                          IF ( GNRZD ) THEN
	                     CALL SGER   (NI,NVB-I,-ONE,BR(1,I),1,
     &                                    BETA(I+1,I),1,BR(1,I+1),NI)
	                     CALL SGER   (NI,NVB-I,-ONE, R(1,I),1,
     &                                    BETA(I+1,I),1, R(1,I+1),NI)
                          ELSE
	                     CALL SGER   (NI,NVB-I,-ONE, R(1,I),1,
     &                                    BETA(I+1,I),1, R(1,I+1),NI)
                          END IF
	               END IF
C
	       END IF
C
            END IF
C
	 END IF
C
   30 CONTINUE
C
C.... transpose BETA ...................................................
C
      DO 50 I = 1,NVB
         DO 40 J = I+1,NVB
            BETA(I,J) = BETA(J,I)
            BETA(J,I) = ZERO
   40    CONTINUE
   50 CONTINUE
C
      RETURN
C
C**** end of MGSCHM ****************************************************
C
      END
