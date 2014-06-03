      SUBROUTINE DATSET (IPSET,RPSET,JLMAX,JTMAX,JTMIN,NVB,LEIG,LTAU,
     &                   LNI,NI,N,LCOMM,MYPE,NPE,LFILE,LPRNT,LOPTS,
     &                   LRERR,LRMDE,NQMAX,NBXMAX,NRUNMX,NSIMAX,
     &                   NSVIN,BIGNUM,EPS,EPS1,REPS,NDEIG,NPEIG,
     &                   NREIG,NREIGL,NREIGR,NTEIG,EIGL,EIGR,
     &                   ENDL,ENDR,SIGMA,THETA0,THRSH,IBUSY,
     &                   RBUSY,AGAIN,STRON)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    DATSET initializes variables                                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    IPSET  (aii) : integer input data                                *
C*    RPSET  (aii) : real input data                                   *
C*    JLMAX  (sio) : maximum number of steps                           *
C*    JTMAX  (sio) : maximum dimension of the block tridiagonal matrix *
C*    JTMIN  (sii) : minimum dimension of the block tridiagonal matrix *
C*    NVB    (sio) : number of vectors in a block                      *
C*    LEIG   (sio) : leading dimension of (EIG)                        *
C*    LTAU   (sio) : leading dimension of (TAU)                        *
C*    LNI    (sio) : leading dimension of (U), (V) and (X)             *
C*    NI     (sio) : dimension of the vectors in (U), (V) and (X)      *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    LCOMM  (sio) : communicator for the parallel version             *
C*    MYPE   (sio) : process rank                                      *
C*    NPE    (sii) : number of processes                               *
C*    LFILE  (sio) : file unit for output                              *
C*    LPRNT  (sio) : level of printing                                 *
C*    LOPTS  (sio) : options for `generalized', `slice' and `purify'   *
C*    LRERR  (sio) : code for error messages                           *
C*    LRMDE  (sio) : run mode                                          *
C*    NQMAX  (sio) : maximum number of vectors in (BASIS)              *
C*    NBXMAX (sio) : maximum number of vectors in (BX)                 *
C*    NRUNMX (sii) : maximum number of runs                            *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    NSVIN  (sio) : number of starting vectors given as input         *
C*    BIGNUM (sri) : big number                                        *
C*    EPS    (sri) : roundoff unit                                     *
C*    EPS1   (sro) : EPS*NVB*sqrt(N)                                   *
C*    REPS   (sro) : sqrt(EPS)                                         *
C*    NDEIG  (sio) : number of eigenpairs required in the run          *
C*    NPEIG  (sio) : number of eigenpairs given as input               *
C*    NREIG  (sio) : number of required eigenpairs                     *
C*    NREIGL (sio) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sio) : number of required eigenvalues greater than EIGL  *
C*    NTEIG  (sio) : number of computed eigenpairs                     *
C*    EIGL   (sro) : inferior bound for eigenvalues                    *
C*    EIGR   (sro) : superior bound for eigenvalues                    *
C*    ENDL   (sro) : inferior bound for Ritz values                    *
C*    ENDR   (sro) : superior bound for Ritz values                    *
C*    SIGMA  (sro) : origin translation                                *
C*    THETA0 (sro) : reference point for THETA                         *
C*    THRSH  (sro) : threshold for convergence                         *
C*    IBUSY  (sib) : number of positions already used in ISTOR         *
C*    RBUSY  (sib) : number of positions already used in RSTOR         *
C*    AGAIN  (slo) : loop control flag                                 *
C*    STRON  (slo) : run starting flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZHEAD,LZMMRY,LZRANG,SETLRM,SIBTST                               *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    INT,MIN,SQRT                                                     *
C*    DBLE                                                             *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      DOUBLE PRECISION ONE,ZERO
      PARAMETER        (ONE=1.0D0,ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          IBUSY,IPSET(*),JLMAX,JTMAX,JTMIN,LCOMM,LEIG,
     &                 LFILE,LNI,LOPTS(*),LPRNT,LRERR,LRMDE,LTAU,
     &                 MYPE,N,NBXMAX,NDEIG,NI,NPEIG,NPE,NQMAX,
     &                 NREIG,NREIGL,NREIGR,NRUNMX,NSIMAX,
     &                 NSVIN,NTEIG,NVB,RBUSY
      DOUBLE PRECISION BIGNUM,EIGL,EIGR,ENDL,ENDR,EPS,EPS1,
     &                 RPSET(*),REPS,SIGMA,THETA0,THRSH
      LOGICAL          AGAIN,STRON
C
C==== local variables ==================================================
C
      INTEGER          LISTOR,LRSTOR,NLSIN,NVBIN
      LOGICAL          EXIT0,GNRZD
C
C==== subprogram =======================================================
C
      LOGICAL          SIBTST
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        INT,MIN,SQRT
      INTRINSIC        DBLE
C
C**** executable statements ********************************************
C
C.... assign input data to the appropriate variables ...................
C
      NI       = IPSET( 1)
      LNI      = IPSET( 2)
      NREIG    = IPSET( 3)
      LEIG     = IPSET( 4)
      NVBIN    = IPSET( 5)
      NLSIN    = IPSET( 6)
      NSVIN    = IPSET( 7)
      NPEIG    = IPSET( 8)
      LOPTS(1) = IPSET( 9)
      LOPTS(2) = IPSET(10)
      LOPTS(3) = IPSET(11)
      LFILE    = IPSET(13)
      LCOMM    = IPSET(14)
      LISTOR   = IPSET(15)
C
      THRSH  = RPSET(3)
      LRSTOR = INT(RPSET(4))
C
      NTEIG  = NPEIG
C
C.... set printing level ...............................................
C
      IF ( LISTOR.GT.0 .AND. LRSTOR.GT.0 ) THEN
         IF ( MYPE.EQ.0 ) THEN
            IF ( IPSET(12) .GT. 0 ) LPRNT = LPRNT + 2**1
            IF ( IPSET(12) .GT. 1 ) LPRNT = LPRNT + 2**2
            IF ( IPSET(12) .GT. 2 ) LPRNT = LPRNT + 2**3
            IF ( IPSET(12) .GT. 3 ) LPRNT = LPRNT + 2**4
            IF ( IPSET(12) .GT. 4 ) LPRNT = LPRNT + 2**5
            IF ( IPSET(12) .GT. 5 ) LPRNT = LPRNT + 2**6
            IF ( IPSET(12) .GT. 6 ) LPRNT = LPRNT + 2**7
         END IF
         IF ( IPSET(12) .GT. 5 ) LPRNT = LPRNT + 2**15
      END IF
C
      IF ( LRERR .GT. 0 ) GO TO 10
C
C.... set the block size and the size of the basis .....................
C
      NVB   = MIN(NVB,N)
      LTAU  = MIN(LEIG,N)
      JTMAX = NLSIN*NVBIN
      GNRZD = LOPTS(1) .GT. 0
C
      IF ( JTMAX .EQ. 0 ) THEN
         JTMAX = MIN(JTMIN*6,N)
      ELSE
         JTMAX = MIN(JTMAX,N)
      END IF
C
C.... check the workspace available ....................................
C
      CALL LZMMRY (JLMAX,NLSIN,JTMAX,NVB,NVBIN,LTAU,NI,N,LCOMM,NPE,
     &             LRERR,NQMAX,NBXMAX,NRUNMX,NSIMAX,LISTOR,
     &             LRSTOR,IBUSY,RBUSY,GNRZD)
C
      IF ( LRERR.EQ.0 .AND. ( LISTOR.LE.0 .OR. LRSTOR.LE.0 ) ) THEN
         IPSET(15) = -1
         RPSET(4) = -1
      END IF
C
      JTMIN = MIN(80,JTMIN+10*NVB,JTMAX)
      EPS1 = EPS*DBLE(NVB)*SQRT(DBLE(N))
      REPS = SQRT(EPS)
C
C.... set the eigenvalue interval ......................................
C
      IF      ( LOPTS(1) .EQ. 1 ) THEN
              THETA0 = ZERO
      ELSE IF ( LOPTS(1) .EQ. 2 ) THEN
              THETA0 = ONE
      ELSE
              LOPTS(2) = 0
              LOPTS(3) = 0
      END IF
C
      CALL LZRANG (LRMDE,NDEIG,NPEIG,NREIG,NREIGL,NREIGR,BIGNUM,
     &             EIGL,EIGR,ENDL,ENDR,RPSET,SIGMA,THETA0,GNRZD)
C
   10 CONTINUE
C
C.... print header .....................................................
C
      IF ( SIBTST(1,LPRNT) ) CALL LZHEAD (JLMAX,NVB,LTAU,N,LFILE,LRERR,
     &                                    NREIG,NTEIG,EIGL,EIGR,SIGMA)
C
C.... set flags accordingly ............................................
C
      EXIT0 = LRERR.EQ.0 .AND. MIN(LISTOR,LRSTOR).LE.0
C
      STRON = LRERR.EQ.0 .AND. .NOT.EXIT0 .AND. .NOT.GNRZD
      AGAIN = STRON
C
      RETURN 
C
C**** end of DATSET ****************************************************
C
      END
