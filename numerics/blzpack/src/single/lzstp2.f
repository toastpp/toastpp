      SUBROUTINE LZSTP2 (IR,JL,JLMAX,NVB,LNI,NI,LBLAS,LCOMM,LRERR,
     &                   LRWRN,NULLDQ,NULLDR,V,ANORM,BETAR,BNORM,
     &                   R,TIME,WORK,EPS,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZSTP2 factors (R) as (Q_{j+1})*(BETA_{j+1})                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    IR     (aii) : pointers for (R)                                  *
C*    JL     (sii) : number of steps                                   *
C*    JLMAX  (sii) : maximum number of steps                           *
C*    NVB    (sii) : number of vectors in a block                      *
C*    LNI    (sii) : leading dimension of (U), (V) and (X)             *
C*    NI     (sii) : dimension of the vectors in (U), (V) and (X)      *
C*    LBLAS  (aii) : BLAS level setting                                *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    LRERR  (sio) : code for error messages                           *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NULLDQ (sio) : number of zero diagonal entries in BETAQ          *
C*    NULLDR (sio) : number of zero diagonal entries in BETAR          *
C*    V      (arb) : array for reverse communication, V(LNI,NVB)       *
C*    ANORM  (arb) : extreme singular values of (ALPHA)                *
C*    BETAR  (arb) : matrix (BETA) in (R)*(Q)*(BETA) at current  step  *
C*    BNORM  (arb) : extreme singular values of (BETA)                 *
C*    R      (arb) : work array for Lanczos vectors                    *
C*    TIME   (arb) : time table                                        *
C*    WORK   (arw) : workspace                                         *
C*    EPS    (sri) : roundoff unit                                     *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    IDENTY,LZCOPY,RFACTR,SITIME                                      *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C
      REAL             ZERO
      PARAMETER        (ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          IR(4),JL,JLMAX,LBLAS(*),LCOMM,LNI,LRERR,
     &                 LRWRN,NI,NULLDQ,NULLDR,NVB
      REAL             ANORM(JLMAX,2),BETAR(NVB,NVB),BNORM(JLMAX,2),
     &                 EPS,R(*),TIME(*),V(*),WORK(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          IDELTA,IENORM,IWORK,JLP1
      REAL             ANORMJ,TIME0
C
C==== subprogram =======================================================
C
      REAL             SITIME
C
C**** executable statements ********************************************
C
      IDELTA = 1
      IENORM = IDELTA + NVB*NVB
      IWORK  = IENORM + NVB
      JLP1   = JL + 1
C
C.... the product (BR)=op(B)*(R) is required at this point .............
C
      IF ( JL .GT. 0 ) THEN
         ANORMJ = ANORM(JL,2)
      ELSE
         ANORMJ = ZERO
      END IF
C
C.... deal with an ill conditioned or invariant subspace ...............
C
      IF ( NULLDR .GT. 0 ) THEN
         NULLDQ = NULLDR
      ELSE
         NULLDQ = 0
      END IF
C
C.... factorize (R) and get the next (Q) ...............................
C
      IF ( GNRZD ) THEN
         TIME(2) = TIME(2) + SITIME(TIME(11))
         CALL LZCOPY (LNI,NI,NI,NVB,V,R(IR(3)))
      END IF
C
      TIME0 = SITIME(ZERO)
C
      CALL IDENTY (NVB,BETAR)
C
      CALL RFACTR (LBLAS(4),LCOMM,LRERR,LRWRN,NI,NVB,NULLDR,BETAR,
     &             R(IR(3)),R(IR(1)),WORK(IDELTA),WORK(IENORM),
     &             WORK(IWORK),ANORMJ,BNORM(JLP1,2),
     &             BNORM(JLP1,1),EPS,GNRZD)
C
      TIME(4) = TIME(4) + SITIME(TIME0)
C
      RETURN 
C
C**** end of LZSTP2 ****************************************************
C
      END
