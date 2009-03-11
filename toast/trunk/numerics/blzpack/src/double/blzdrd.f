      SUBROUTINE BLZDRD ( ISTOR ,
     &                    RSTOR ,
     &                    SIGMA ,
     &                    NNEIG ,
     &                    U     ,
     &                    V     ,
     &                    LFLAG ,
     &                    NVOPU ,
     &                    EIG   ,
     &                    X     )
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZDRD is the user interface for a set of subprograms intended   *
C*    for the solution of the standard symmetric eigenvalue problem    *
C*                                                                     *
C*                 (A)*(x) - eig*(I)*(x) = (0)                         *
C*                                                                     *
C*    or the generalized symmetric eigenvalue problem                  *
C*                                                                     *
C*                 (A)*(x) - eig*(B)*(x) = (0)                         *
C*                                                                     *
C*    by means of a block Lanczos type algorithm.                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR  : integer array of dimension 17+LISTOR (defined below).   *
C*             On the first call to BLZDRD, the first 15 entries of    *
C*             ISTOR must contain:                                     *
C*                                                                     *
C*             ISTOR( 1)=NI,     dimension of the vectors in (U), (V)  *
C*                               and (X), global_sum(NI) = N = dim(A)  *
C*                               in the parallel implementation        *
C*             ISTOR( 2)=LNI,    1st dimension of (U), (V) and (X),    *
C*                               LNI >= NI                             *
C*             ISTOR( 3)=LEIG,   maximum number of eigenpairs,         *
C*                               1st dimension of (EIG)                *
C*                               2nd dimension of (X)                  *
C*             ISTOR( 4)=NREIG,  number of required eigenpairs         *
C*             ISTOR( 5)=NVBIN,  number of vectors in a block          *
C*             ISTOR( 6)=NLSIN,  maximum number of steps per run       *
C*             ISTOR( 7)=NSVIN,  number of starting vectors provided   *
C*             ISTOR( 8)=NEPIN,  number of eigenpairs provided         *
C*             ISTOR( 9)=GNRZD,  problem type flag                     *
C*                               : 0, standard    eigenproblem         *
C*                               : 1, generalized eigenproblem         *
C*                               : 2, buckling    eigenproblem         *
C*             ISTOR(10)=SLICE,  spectrum slicing flag                 *
C*                               : 0, slicing off                      *
C*                               : 1, slicing on                       *
C*             ISTOR(11)=PURFY,  eigenvectors purification flag        *
C*                               : 0, purification off                 *
C*                               : 1, purification on                  *
C*             ISTOR(12)=LPRNT,  level of printing                     *
C*                               : 0, nothing is printed               *
C*                               : 1, prints header and exit messages  *
C*                               : 2, the same of 1 plus eigenvalues   *
C*                               : 3, the same of 2 plus statistics    *
C*                               : 4, the same of 3 plus run history   *
C*                               : 5, the same of 4 plus eigenvectors  *
C*                               : 6, the same of 5 plus step history  *
C*             ISTOR(13)=LFILE,  file unit for output                  *
C*             ISTOR(14)=LCOMM,  MPI version: communicator             *
C*                               PVM version: number of PEs            *
C*                               sequential version: not used          *
C*             ISTOR(15)=LISTOR, dimension of (ISTOR), defined below   *
C*                                                                     *
C*    RSTOR  : real array of dimension 5+LRSTOR (define below). On the *
C*             first call to BLZDRD, the first 4 entries of RSTOR must *
C*             contain:                                                *
C*                                                                     *
C*             RSTOR( 1)=EIGL,   inferior limit for eigenvalues        *
C*             RSTOR( 2)=EIGR,   superior limit for eigenvalues        *
C*             RSTOR( 3)=THRSH,  threshold for convergence             *
C*             RSTOR( 4)=LRSTOR, dimension of (RSTOR), defined below   *
C*                                                                     *
C*    SIGMA  : translation of the eigenvalue spectrum                  *
C*                                                                     *
C*    NNEIG  : number of eigenvalues less than SIGMA for the matrix    *
C*             (A)-SIGMA*(B), which can be recovered from (D) in       *
C*             the factorization (A)-SIGMA*(B)=(L)*(D)*(L')            *
C*                                                                     *
C*    U      : real work array, dimension (LNI,NCU)                    *
C*                              NCU=max(NVB,NVBIN)                     *
C*                                                                     *
C*    V      : real work array, dimension (LNI,NCU)                    *
C*                              NCU=max(NVB,NVBIN)                     *
C*                                                                     *
C*    LFLAG  : reverse communication flag                              *
C*                                                                     *
C*             < 0, early or abnormal finalization                     *
C*             = 0, initialization or normal finalization              *
C*             = 1, given (U), compute (V)                             *
C*                  if GNRZD=0, use (V)=(A)*(U)                        *
C*                  if GNRZD>0, use (L)*(D)*(L')*(V)=(B)*(U)           *
C*             = 2, given (U), compute (V)=(B)*(U)                     *
C*             = 3, given SIGMA, compute (A)-SIGMA*(B)=(L)*(D)*(L')    *
C*             = 4, define the starting vectors in (V)                 *
C*             = 5, forces finalization                                *
C*                                                                     *
C*    NVOPU  : number of vectors for (V)=(A)*(U) or (V)=(B)*(U)        *
C*                                                                     *
C*    EIG    : eigenvalue  approximations, dimension EIG(LEIG,2)       *
C*                                                                     *
C*    X      : eigenvector approximations, dimension X(LNI,LEIG)       *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    BLZEIG,BLZEND,BLZEXT,BLZFCT,BLZKLL,BLZRBX,BLZSET,BLZSTP,BLZSTR   *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*    LISTOR .ge. 123 + k1 x 12                                        *
C*                                                                     *
C*    In sequential mode,                                              *
C*                                                                     *
C*    LRSTOR .ge. NI x k2 x 5 + k3 if GNRZD = 0,                       *
C*                NI x ( k2 x 6 + 1 ) + k3 otherwise                   *
C*                                                                     *
C*    and in parallel mode,                                            *
C*                                                                     *
C*    LRSTOR .ge. NI x ( k2 x 4 + k1 ) + k3 if GNRZD = 0,              *
C*                NI x ( k2 x 4 + k1 x 2 + k4 ) + k3 otherwise         *
C*                                                                     *
C*    where                                                            *
C*                                                                     *
C*    k1 = NVBSET x NSTEPS if NVBSET x NSTEPS > 0,                     *
C*         min(N,180) otherwise                                        *
C*                                                                     *
C*    k2 = NVBSET if NVBSET > 0,                                       *
C*         3 otherwise                                                 *
C*                                                                     *
C*    k3 = 484 + k1 x ( 13 + k1 x 2 + k2 + max(18,k2+2) ) +            *
C*         k2 x k2 x 3 + k4 x 2                                        *
C*                                                                     *
C*    k4 = min(LEIG,N)                                                 *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*),LFLAG,NNEIG,NVOPU
      DOUBLE PRECISION EIG(*),RSTOR(*),SIGMA,U(*),V(*),X(*)
C
C==== local variables ==================================================
C
      LOGICAL          ABORT,AGAIN,EIGON,ENDON,KLLON,RBXON,STPON,STRON
C
C**** executable statements ********************************************
C
C.... check flags ......................................................
C
      ABORT = LFLAG.NE.0 .AND. LFLAG.NE.5 .AND. LFLAG.NE.ISTOR(EXIT0)
C
      ENDON = ABORT
      AGAIN = .FALSE.
      EIGON = .FALSE.
      STRON = LFLAG.EQ.4
      KLLON = LFLAG.EQ.5
      STPON = LFLAG.EQ.1 .OR. ( LFLAG.EQ.2 .AND. ISTOR(EXIT1).EQ.0 )
      RBXON = LFLAG.EQ.2 .AND. ISTOR(EXIT1).NE.0
C
C.... loop on the main routines ........................................
C
   10 CONTINUE
C
      IF      ( EIGON ) THEN
C
C............ compute the eigenvectors .................................
C
              CALL BLZEIG (ISTOR(IINIT),RSTOR(RINIT),LFLAG,SIGMA,NNEIG,
     &                     EIG,X,U,AGAIN,ENDON,RBXON,STRON)
C
              EIGON = .FALSE.
              NVOPU = 0
C
      ELSE IF ( ENDON ) THEN            
C
C............ finish ...................................................
C
              CALL BLZEND (ISTOR(IINIT),RSTOR(RINIT),LFLAG,EIG,X,ABORT)
C
              ISTOR(EXIT1) = ISTOR(IINIT+LRWRN-1)
              AGAIN = .FALSE.
C
      ELSE IF ( KLLON ) THEN
C
C............ finalization enforced ....................................
C
              CALL BLZKLL (ISTOR(IINIT),RSTOR(RINIT),SIGMA,EIG,X,U)
C
              KLLON = .FALSE.
              AGAIN = .TRUE.
              ENDON = .TRUE.
C
      ELSE IF ( RBXON ) THEN
C
C............ compute (B)*(X) ..........................................
C
              CALL BLZRBX (ISTOR(IINIT),RSTOR(RINIT),LFLAG,NVOPU,
     &                     X,U,V,AGAIN,ENDON,RBXON,STRON)
C
              ISTOR(EXIT1) = NVOPU
              RBXON = .FALSE.
              STPON = .FALSE.
C
      ELSE IF ( STPON ) THEN
C
C............ perform a block Lanczos step .............................
C
              CALL BLZSTP (ISTOR(IINIT),RSTOR(RINIT),SIGMA,LFLAG,
     &                     NVOPU,EIG,X,U,V,AGAIN,EIGON,ENDON)
C
              ISTOR(EXIT1) = 0
              STPON = .FALSE.
C
      ELSE IF ( STRON ) THEN
C
C............ prepare to start a Lanczos run ...........................
C
              CALL BLZSTR (ISTOR(IINIT),RSTOR(RINIT),SIGMA,LFLAG,
     &                     NVOPU,U,V,EIG,X,AGAIN,ENDON,STPON)
C
              ISTOR(EXIT1) = 0
              STRON = .FALSE.
C
      ELSE IF ( LFLAG .EQ. 0 ) THEN
C
C............ BLZDRD is being called for the first time ................
C
              CALL BLZSET (ISTOR,RSTOR,ISTOR(IINIT),RSTOR(RINIT),
     &                     SIGMA,LFLAG,AGAIN,STRON)
C
      ELSE IF ( LFLAG .EQ. 3 ) THEN
C
C............ factorization performed (update all subintervals) ........
C
              CALL BLZFCT (ISTOR(IINIT),RSTOR(RINIT),SIGMA,NNEIG,
     &                     EIG,AGAIN,ENDON,RBXON,STRON)
C
      END IF
C
      IF ( AGAIN ) GO TO 10
C
C.... record time and exit flag ........................................
C
      CALL BLZEXT (ISTOR,RSTOR,LFLAG)
C
      RETURN 
C
C**** end of BLZDRD ****************************************************
C
      END
