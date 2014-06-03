      SUBROUTINE SSNEED (JT,JTMIN,NVB,LCOMM,NPE,LRMDE,NEWSIG,NSLS,NSRS,
     &                   NSRLS,NSRRS,BIGNUM,EPS,GRNRM,RNORM,THETA,
     &                   THETAL,THETAR,THRES,TIME)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSNEED checks for the need of a new origin translation           *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    JTMIN  (sii) : minimum dimension of the block tridiagonal matrix *
C*    NVB    (sii) : number of vectors in a block                      *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    NPE    (sii) : number of processes                               *
C*    LRMDE  (sii) : run mode                                          *
C*    NEWSIG (sib) : flag for a new starting point                     *
C*    NSLS   (sii) : number of eigenvalues converged < than SIGMA      *
C*    NSRS   (sii) : number of eigenvalues converged > than SIGMA      *
C*    NSRLS  (sii) : number of eigenvalues required  < than SIGMA      *
C*    NSRRS  (sii) : number of eigenvalues required  > than SIGMA      *
C*    BIGNUM (sri) : big number                                        *
C*    EPS    (sri) : roundoff unit                                     *
C*    GRNRM  (srb) : global residual of unconverged eigenpairs         *
C*    RNORM  (ari) : estimated residuals                               *
C*    THETA  (ari) : eigenvalues of the tridiagonal matrix             *
C*    THETAL (sri) : inferior limit to converged eigenvalues           *
C*    THETAR (sri) : superior limit to converged eigenvalues           *
C*    THRES  (sri) : threshold for convergence                         *
C*    TIME   (ari) : time table                                        *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    PIDRED,SITIME,TBILLC                                             *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,LOG10,MAX                                                    *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      DOUBLE PRECISION ONE,THIRD,ZERO
      PARAMETER        (ONE=1.0D0,THIRD=1.0D0/3.0D0,ZERO=0.0D0)
      DOUBLE PRECISION WMAX,WMIN
      PARAMETER        (WMAX=1.40D0,WMIN=0.40D0)
C
C==== arguments ========================================================
C
      INTEGER          JT,JTMIN,LCOMM,LRMDE,NEWSIG,NPE,
     &                 NSLS,NSRS,NSRLS,NSRRS,NVB
      DOUBLE PRECISION BIGNUM,TIME(*),EPS,GRNRM,RNORM(JT),
     &                 THETA(JT),THETAL,THETAR,THRES
C
C==== local variables ==================================================
C
      INTEGER          I,INDEX,INFO,J,NRNRM,NTSOL
      DOUBLE PRECISION C1,C2,W1,W2,W3,RSNRM,T(5),TBASIS,
     &                 TEIGTB,TFACTR,THTMAX,TSUM(5),TPURGE,TTOTAL,
     &                 WEIGHT
C
C==== subprograms ======================================================
C
      DOUBLE PRECISION SITIME
      LOGICAL          TBILLC
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,LOG10,MAX
C
C**** executable statements ********************************************
C
      NRNRM = 0
      RSNRM = ONE
      NTSOL = NSLS + NSRS
C
C.... compute the average time .........................................
C
      TFACTR = TIME( 8)
      TPURGE = TIME( 9)
      TEIGTB = TIME(10) 
      TTOTAL = SITIME(TIME(12))
      TBASIS = TTOTAL - TEIGTB
C
      T(1) = TFACTR
      T(2) = TPURGE
      T(3) = TEIGTB
      T(4) = TTOTAL
      T(5) = TBASIS
C
      CALL PIDRED ('SUM',5,T,TSUM,LCOMM,INFO)
C
      TFACTR = TSUM(1)/NPE
      TPURGE = TSUM(2)/NPE
      TEIGTB = TSUM(3)/NPE
      TTOTAL = TSUM(4)/NPE
      TBASIS = TSUM(5)/NPE
C
C.... global residual of unconverged pairs .............................
C
      DO 10 I = 1,NVB*2
         INDEX = 0
         THTMAX = ZERO 
         DO 20 J = 1,JT
            IF ( (RNORM(J).LT.ZERO ) .AND. 
     &           (ABS(THETA(J)).GT.THTMAX) .AND.
     &           ((THETA(J).LT.THETAL).OR.(THETA(J).GT.THETAR)) ) THEN
                 THTMAX = ABS(THETA(J))
                 INDEX = J
            END IF
   20    CONTINUE
         IF ( (INDEX.NE.0) .AND. (RSNRM*BIGNUM.GE.ONE) ) THEN
            RSNRM = RSNRM*ABS(RNORM(INDEX))
            RNORM(INDEX) = ZERO
            NRNRM = NRNRM + 1
         END IF
   10 CONTINUE
C
C.... extrapolate the global residual ..................................
C
      IF ( NRNRM .NE. 0 ) RSNRM = RSNRM**(ONE/NRNRM)
C
      IF ( GRNRM .GT. RSNRM ) RSNRM = RSNRM*(RSNRM/GRNRM)
C
C.... weight for global residual .......................................
C
      IF      ( RSNRM .GE. ONE ) THEN
C
C             slow convergence: take the maximum weight
C
              W1 = WMAX 
C
      ELSE IF ( RSNRM .LE. THRES ) THEN
C
C             fast convergence: take the minimum weight
C
              W1 = WMIN 
C
      ELSE
C
C             linear model: W1 = f(ABS(LOG10(THRES/RSNRM)) with 
C             W1(0) = WMIN and W1(ABS(LOG10(THRES))) = WMAX
C
              C1 = WMIN
              C2 = (WMAX-WMIN)/ABS(LOG10(THRES))
              W1 = C1 + C2*ABS(LOG10(THRES/RSNRM))
C
      END IF
C
      GRNRM = RSNRM
C
C.... weight for run time ..............................................
C
      IF      ( TEIGTB .GE. TBASIS   ) THEN
C
C             eig(T) is expensive: take the minimum weight
C
              W2 = WMIN
C
      ELSE IF ( TTOTAL .LT. TPURGE*2 ) THEN
C
C             orthogonalization is expensive: take the maximum weight
C
              W2 = WMAX
C
      ELSE IF ( TTOTAL .GE. TPURGE*5 ) THEN
C
C             orthogonalization is cheap: take the minimum weight
C
              W2 = WMIN
C
      ELSE
C
C             linear model: W2 = f(TPURGE/TTOTAL-0.2) with
C             W2(0.2) = WMIN and W2(0.5) = WMAX
C
              C1 = WMIN
              C2 = (WMAX-WMIN)/0.30D0
              W2 = C1 + C2*(TPURGE/TTOTAL-0.20D0)
C
      END IF
C
C.... weight for factorization time ....................................
C
      IF      ( TFACTR .EQ. ZERO   ) THEN
C
C             this is an odd case: take the minimum weight
C
              W3 = WMIN
C
      ELSE IF ( TFACTR .GE. TTOTAL*5 ) THEN
C
C             factorization is expensive: take the minimum weight
C
              W3 = WMIN
C
      ELSE
C
C             linear model: W3 = f(TTOTAL/TFACTR-0.2) with
C             W3(0.2) = WMIN and W3(3.2) = WMAX
C
              C1 = WMIN
              C2 = (WMAX-WMIN)/3.00D0
              W3 = C1 + C2*(TTOTAL/TFACTR-0.20D0)
C
      END IF
C
C.... geometrical average ..............................................
C
      WEIGHT = (W1*W2*W3)**THIRD
C
C.... check for the need of a new origin translation ...................
C
      IF      ( NTSOL.GT.0 .AND. LRMDE.NE.7 .AND. LRMDE.NE.8 .AND.
     &          JTMIN.LT.JT .AND. WEIGHT.GT.ONE ) THEN
                NEWSIG = 2
      ELSE IF ( NTSOL.EQ.0 .AND. LRMDE.NE.7 .AND. LRMDE.NE.8 .AND.
     &          JTMIN*2.LT.JT .AND. WEIGHT.GT.ONE ) THEN
                NEWSIG = 2
      ELSE IF ( NTSOL.GT.0 .AND. TBILLC(JT,EPS,THETA) ) THEN
                NEWSIG = 3
      ELSE IF ( NSLS.GE.NSRLS .AND. NSRS.GE.NSRRS ) THEN
                NEWSIG = 4
      ELSE
                NEWSIG = 1
      END IF
C
      RETURN 
C
C**** end of SSNEED ****************************************************
C
      END
