      SUBROUTINE LZSTP4 (IR,JL,JTMAX,NVB,NI,FHNDL,LCOMM,LRERR,NQMAX,
     &                   TB,ALPHA,BETAQ,BETAR,R,BASIS,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZSTP4 stores (ALPHA), (BETA) and Lanczos vectors                *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    IR    (aii) : pointers for (R)                                   *
C*    JL    (sii) : number of steps                                    *
C*    JTMAX (sii) : maximum dimension of the block tridiagonal matrix  *
C*    NVB   (sii) : number of vectors in a block                       *
C*    NI    (sii) : dimension of the vectors in (U), (V) and (X)       *
C*    FHNDL (aii) : file handle                                        *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    NQMAX (sii) : maximum number of vectors in (BASIS)               *
C*    TB    (arb) : block tridiagonal matrix                           *
C*    ALPHA (arb) : (Q')*(B)*(R)                                       *
C*    BETAQ (arb) : matrix (BETA) in (R)*(Q)*(BETA) at previous step   *
C*    BETAR (arb) : matrix (BETA) in (R)*(Q)*(BETA) at current  step   *
C*    R     (arb) : work array for Lanczos vectors                     *
C*    BASIS (arb) : Lanczos vectors array                              *
C*    GNRZD (sli) : problem type flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZCOPY,LZIOOP,TBALPH,TBBETA                                      *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    DCOPY                                                            *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          FHNDL(*),IR(4),JL,JTMAX,LCOMM,LRERR,
     &                 NI,NQMAX,NVB
      DOUBLE PRECISION ALPHA(NVB,NVB),BASIS(*),BETAQ(NVB,NVB),
     &                 BETAR(NVB,NVB),R(*),TB(*)
      LOGICAL          GNRZD
C
C==== local variable ===================================================
C
      INTEGER          IQ
C
C**** executable statements ********************************************
C
      IF ( JL .EQ. 0 ) RETURN
C
C.... insert (ALPHA) and (BETAQ) in (TB) ...............................
C
      IF ( JL .GT. 1 ) CALL TBBETA (JL-1,JTMAX,NVB,BETAQ,TB)
                       CALL TBALPH (JL  ,JTMAX,NVB,ALPHA,TB)
C
      CALL DCOPY (NVB*NVB,BETAR,1,BETAQ,1)
C
C.... store the j-th vector ............................................
C
C     R(IR(2)) -> (B)*(Q)
C     R(IR(4)) -> (Q)
C
      IF ( JL .LE. NQMAX ) THEN
C
C....... copy (B)*(Q) and (Q) into (BASIS) .............................
C
         IF (GNRZD) THEN
            IQ = 1 + NI*NVB*(JL-1)*2
            CALL LZCOPY (NI,NI*2,NI,NVB,R(IR(2)),BASIS(IQ+NI))
            CALL LZCOPY (NI,NI*2,NI,NVB,R(IR(4)),BASIS(IQ   ))
         ELSE
            IQ = 1 + NI*NVB*(JL-1)
            CALL LZCOPY (NI,NI  ,NI,NVB,R(IR(4)),BASIS(IQ   ))
         END IF
C
      ELSE 
C
C....... copy (B)*(Q) and (Q) into a file ..............................
C
         IF (GNRZD) THEN
            CALL LZIOOP (FHNDL,LCOMM,LRERR,JL-NQMAX,NI*NVB,
     &                   R(IR(2)),'BQ','PUT')
            IF ( LRERR .NE. 0 ) RETURN
         END IF
         CALL LZIOOP (FHNDL,LCOMM,LRERR,JL-NQMAX,NI*NVB,
     &                R(IR(4)),'Q ','PUT')
         IF ( LRERR .NE. 0 ) RETURN
C
      END IF
C
      RETURN 
C
C**** end of LZSTP4 ****************************************************
C
      END
