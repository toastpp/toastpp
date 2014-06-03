      SUBROUTINE LZSTP1 (IR,JL,JLMAX,JT,NVB,LNI,NI,LBLAS,LCOMM,LRERR,
     &                   NMOPA,NMOPB,EPS,TIME,ALPHA,ANORM,
     &                   BETAR,R,U,V,WORK,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZSTP1 computes (R)=(V)-(Q_{j})*(ALPHA_{j})-(Q_{j-1})*(BETA_{j}) *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    IR    (aii) : pointers for (R)                                   *
C*    JL    (sib) : number of steps                                    *
C*    JLMAX (sii) : maximum number of steps                            *
C*    JT    (sib) : dimension of the block tridiagonal matrix          *
C*    NVB   (sii) : number of vectors in a block                       *
C*    LNI   (sii) : leading dimension of (U), (V) and (X)              *
C*    NI    (sii) : dimension of the vectors in (U), (V) and (X)       *
C*    LBLAS (aii) : BLAS level setting                                 *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    NMOPA (sib) : number of op(A)*vector performed                   *
C*    NMOPB (sib) : number of op(B)*vector performed                   *
C*    EPS   (sri) : roundoff unit                                      *
C*    TIME  (arb) : time table                                         *
C*    ALPHA (aro) : (Q')*(B)*(R)                                       *
C*    ANORM (arb) : extreme singular values of (ALPHA)                 *
C*    BETAR (ari) : matrix (BETA) in (R)*(Q)*(BETA) at current  step   *
C*    R     (arb) : work array for Lanczos vectors                     *
C*    U     (arb) : array for reverse communication, U(LNI,NVB)        *
C*    V     (arb) : array for reverse communication, V(LNI,NVB)        *
C*    WORK  (arw) : workspace                                          *
C*    GNRZD (sli) : problem type flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZCOPY,NORM2A,QTBR,RQALPH,RQBETA,SITIME                          *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          IR(4),JL,JLMAX,JT,LBLAS(4),LCOMM,
     &                 LNI,LRERR,NI,NMOPA,NMOPB,NVB
      DOUBLE PRECISION ALPHA(NVB,NVB),ANORM(JLMAX,2),BETAR(NVB,NVB),
     &                 EPS,R(*),TIME(*),U(*),V(*),WORK(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      DOUBLE PRECISION DUMMY,TIME0
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
C
C**** executable statements ********************************************
C
C.... the product (R)=op(A)*op(B)*(Q) is required at this point ........
C
      JL = JL + 1
      JT = JL*NVB
      NMOPA = NMOPA + NVB
C
      TIME0 = SITIME(ZERO)
      TIME(1) = TIME(1) + SITIME(TIME(11))
C
      CALL LZCOPY (LNI,NI,NI,NVB,V,R(IR(1))) 
C
C.... orthogonalize (R) against the j-1 th (Q) .........................
C
      IF ( JL .GT. 1 ) CALL RQBETA (R(IR(1)),R(IR(3)),BETAR,
     &                              LBLAS(1),NI,NVB)
C
C.... orthogonalize (R) against the j   th (Q) .........................
C
      CALL QTBR   (NI,R(IR(2)),NI,NVB,R(IR(1)),NI,NVB,ALPHA,
     &             DUMMY,WORK,LBLAS(2),LCOMM,LRERR)
C
      CALL RQALPH (NI,R(IR(1)),NI,NVB,R(IR(4)),NI,NVB,ALPHA,LBLAS(3))
C
      CALL NORM2A (NVB,ALPHA,WORK,EPS,ANORM(JL,2),ANORM(JL,1))
C
C.... set for (BR)=op(B)*(R) if (A)*(x)-eig*(B)*(x)=(0) ................
C
      IF ( GNRZD ) THEN
         CALL LZCOPY (NI,LNI,NI,NVB,R(IR(1)),U)
         NMOPB = NMOPB + NVB
      END IF
C
      TIME(4) = TIME(4) + SITIME(TIME0)
C
      RETURN 
C
C**** end of LZSTP1 ****************************************************
C
      END
