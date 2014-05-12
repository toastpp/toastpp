      SUBROUTINE TBCONV (JT,NVB,NSLS,NSRS,NSRLS,NSRRS,NULLD,BIGNUM,EPS,
     &                   REPS,BETA,RNORM,S,THETA,THETA0,THETAL,THETAR,
     &                   THRES,THRSH,GNRZD,SLICE,OUTER)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBCONV looks for converged THETA's, the eigenvalues of (TB)      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT     (sii) : dimension of the basis                            *
C*    NVB    (sii) : number of vectors in a block                      *
C*    NSLS   (sio) : number of negative eigenvalues converged          *
C*    NSRS   (sio) : number of positive eigenvalues converged          *
C*    NSRLS  (sii) : number of negative eigenvalues required           *
C*    NSRRS  (sii) : number of positive eigenvalues required           *
C*    NULLD  (sii) : number of eigenvalues equal to zero               *
C*    BIGNUM (sri) : big number                                        *
C*    EPS    (sri) : roundoff unit                                     *
C*    REPS   (sri) : sqrt(EPS)                                         *
C*    BETA   (ari) : matrix (BETA) in (R)=(Q)*(BETA) at the JL-th step *
C*    RNORM  (aro) : estimated residuals of (TB)*(s)-(theta)*(s)       *
C*    S      (ari) : eigenvectors of the tridiagonal matrix            *
C*    THETA  (aro) : eigenvalues  of the tridiagonal matrix            *
C*    THETA0 (sri) : reference point for THETA                         *
C*    THETAL (sri) : inferior limit to converged eigenvalues           *
C*    THETAR (sri) : superior limit to converged eigenvalues           *
C*    THRES  (sro) : threshold for convergence                         *
C*    THRSH  (sri) : threshold for convergence (default)               *
C*    GNRZD  (sli) : problem type flag                                 *
C*    SLICE  (sli) : spectrum slicing flag                             *
C*    OUTER  (slo) : convergence flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    RESNRM,TBGMIN                                                    *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,MAX,MIN                                                      *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*            theta x               .                                  *
C*                  x               . *                                *
C*                  x               .                                  *
C*                  x               .                                  *
C*                  x               .                                  *
C*                  x               .  *                               *
C*                  x               .                                  *
C*                  x               .                                  *
C*                  x               .   *                              *
C*                  x               .                                  *
C*                  x               .    *                             *
C*                  x               .                                  *
C*            zetar x......................*                           *
C*                  |               .      . *                         *
C*                  |               .      .   *                       *
C*                  |               .      .      *                    *
C*                  |	            .	   .          *                *
C*                  |	            .      .		   *           *
C*           theta0 |      lower    .    upper		          *    *
C*    --------------+--------[======+======]-----------------------    *
C*    *             |        .    sigma                        ritz    *
C*           *      |        .      .                                  *
C*                * |        .      .                                  *
C*                  | *      .      .                                  *
C*                  |    *   .      .                                  *
C*                  |      * .      .                                  *
C*            zetal x........*      .                                  *
C*                  x               .                                  *
C*                  x          *    .                                  *
C*                  x               .                                  *
C*                  x           *   .                                  *
C*                  x               .                                  *
C*                  x               .                                  *
C*                  x            *  .                                  *
C*                  x               .                                  *
C*                  x               .                                  *
C*                  x               .                                  *
C*                  x             * .                                  *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      REAL             TWO,ZERO
      PARAMETER        (TWO=2.0E0,ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          JT,NSLS,NSRS,NSRLS,NSRRS,NULLD,NVB
      REAL             BETA(NVB,NVB),BIGNUM,EPS,REPS,RNORM(JT),
     &                 S(JT,JT),THETA(JT),THETA0,THETAL,
     &                 THETAR,THRES,THRSH
      LOGICAL          GNRZD,OUTER,SLICE
C
C==== local variables ==================================================
C
      INTEGER          I,JEND,K,NRLZL,NRRZR,NSLZL,NSRZR
      REAL             BETJ,BETJ2,GAP,RSBND,TBNORM,
     &                 THETA2,TOL,ZETAL,ZETAR
      LOGICAL          IGNORE,OUTERL,OUTERU
C
C==== subprograms ======================================================
C
      REAL             TBGMIN,RESNRM
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,MAX,MIN
C
C**** executable statements ********************************************
C
C.... set convergence parameters and the threshold .....................
C
      NSLZL = 0
      NSRZR = 0
      JEND = JT - NVB + 1
C
      GAP = TBGMIN(JT,REPS,THETA)
      TOL = MIN(TWO**(-18),REPS/TWO)
      TBNORM = MAX(ABS(THETA(1)),ABS(THETA(JT)))
C
      IF ( THRSH .GT. ZERO ) THEN
         THRES = THRSH
      ELSE 
         THRES = MIN(TWO**(-13),TOL*TBNORM)
      END IF
C
C.... compute the residual norms .......................................
C
      DO 10 I = 1,JT
C
	 IF ( NULLD.GT.0 .AND. ABS(THETA(I)).LE.TBNORM*EPS**2 ) THEN
C
C.......... zero eigenvalue corresponding to an invariant subspace .....
C
            THETA(I) = ZERO 
            RNORM(I) = -BIGNUM 
C
         ELSE
C
C.......... residual norms .............................................
C
            BETJ = RESNRM(NVB,BETA,BIGNUM,S(JEND,I))
C
            IF      ( BETJ .LE. THRES ) THEN
                    RSBND  = BETJ
            ELSE IF ( SLICE .AND. BETJ.LE.TOL ) THEN
                    BETJ2  = BETJ**2
                    THETA2 = THETA(I)**2       
                    RSBND  = MIN(BETJ/THETA2,BETJ2/(GAP*THETA2))
            ELSE 
                    RSBND  = BIGNUM 
            END IF 
C
C.......... mark the eigenvalue as converged or not ....................
C
            IF ( RSBND.LE.THRES ) THEN
               RNORM(I) = +MAX(EPS,BETJ)
            ELSE  
               RNORM(I) = -MAX(EPS,BETJ)
            END IF
C
	 END IF
C
   10 CONTINUE
C
C.... set boundaries ...................................................
C
      IF ( THETAR .GE. THETAL ) THEN
         ZETAL = THETAL
         ZETAR = THETAR
         NRLZL = NSRLS
         NRRZR = NSRRS
      ELSE
         ZETAL = THETAR
         ZETAR = THETAL
         NRLZL = NSRRS
         NRRZR = NSRLS
      END IF
C
C.... convergence check for THETA < THETA0 .............................
C
      K = 0
      OUTERL = NRLZL.EQ.0
      IGNORE = OUTERL .AND. GNRZD
C
      DO 20 I = 1,JT,+1
         IF ( THETA(I) .LT. THETA0 ) THEN
            K = K + 1
            IF ( RNORM(I) .GT. ZERO ) THEN
               IF      ( IGNORE ) THEN
                       RNORM(I) = -RNORM(I)
               ELSE IF ( THETA(I) .LT. ZETAL ) THEN
                       NSLZL = NSLZL + 1
               ELSE 
                       IGNORE = .TRUE.
                       IF ( K.EQ.NSLZL+1 .AND. .NOT.SLICE ) THEN
                          NSLZL = NSLZL + 1
                          OUTERL = .TRUE.      
                       ELSE 
                          RNORM(I) = -RNORM(I)
                          OUTERL = SLICE      
                       END IF
               END IF
            END IF
         END IF
   20 CONTINUE
C
C.... convergence check for THETA > THETA0 .............................
C
      K = 0
      OUTERU = NRRZR.EQ.0
      IGNORE = OUTERU .AND. GNRZD
C
      DO 30 I = JT,1,-1
         IF ( THETA(I) .GE. THETA0 ) THEN
            K = K + 1
            IF ( RNORM(I) .GT. ZERO ) THEN
               IF      ( IGNORE ) THEN
                       RNORM(I) = -RNORM(I)
               ELSE IF ( THETA(I) .GT. ZETAR ) THEN
                       NSRZR = NSRZR + 1
               ELSE 
                       IGNORE = .TRUE.
                       IF ( K.EQ.NSRZR+1 .AND. .NOT.SLICE ) THEN
                          NSRZR = NSRZR + 1
                          OUTERU = .TRUE.      
                       ELSE 
                          RNORM(I) = -RNORM(I)
                          OUTERU = SLICE      
                       END IF
               END IF
            END IF
         END IF
   30 CONTINUE
C
C.... set counters .....................................................
C
      IF ( THETAR .GE. THETAL ) THEN
         NSLS = NSLZL
         NSRS = NSRZR
      ELSE
         NSLS = NSRZR
         NSRS = NSLZL
      END IF
C
C.... eigenvalues out of boundaries ....................................
C
      OUTER = OUTERL .AND. OUTERU
C
      RETURN 
C
C**** end of TBCONV ****************************************************
C
      END
