C***********************************************************************
C*                                                                     *
C*    addresses.h is an include file that defines variable addresses   *
C*    in the arrays RSTOR and ISTOR. Variables stored in these arrays  *
C*    can be retrieved by means of calls to the functions RSTORR and   *
C*    ISTORR, respectively; see directory blzpack/etc. Modifications   *
C*    in addresses.h should be followed by appropriate changes in      *
C*    RSTORR and/or ISTORR.                                            *
C*                                                                     *
C***********************************************************************
C
C.... EXIT0 and EXIT1 record exit flags ................................
C
      INTEGER    EXIT0,EXIT1
C
      PARAMETER (EXIT0  = 16)
      PARAMETER (EXIT1  = 17)
C
C.... IINIT and RINIT indicate where the workspaces begin ..............
C
      INTEGER    IINIT,RINIT
C
      PARAMETER (IINIT  = 18)
      PARAMETER (RINIT  =  5)
C
C.... addresses for real scalars (RSTOR) ...............................
C
C     BIGNUM : big number
C     EIGL   : inferior bound for eigenvalues
C     EIGR   : superior bound for eigenvalues
C     ENDL   : inferior bound for Ritz values
C     ENDR   : superior bound for Ritz values
C     EPS    : roundoff unit
C     EPS1   : EPS*NVB*sqrt(N)
C     GRNRM  : global residual of unconverged eigenpairs
C     ORIGIN : starting-point (first SIGMA)
C     RADIUS : radius of convergence around SIGMA
C     REPS   : sqrt(EPS)
C     SFARL  : farthest SIGMA to the left  of ORIGIN
C     SFARR  : farthest SIGMA to the right of ORIGIN
C     THETA0 : reference point for THETA
C     THETAL : inferior limit to eig(TB)               
C     THETAR : superior limit to eig(TB)               
C     THRSH  : threshold for convergence (default)
C     TRUSTL : inferior trust bound
C     TRUSTR : superior trust bound
C
      INTEGER    BIGNUM     , EIGL       , EIGR       , ENDL       , 
     &           ENDR       , EPS        , EPS1       , GRNRM      , 
     &           ORIGIN     , RADIUS     , REPS       , SFARL      ,
     &           SFARR      , THETA0     , THETAL     , THETAR     , 
     &           THRSH      , TRUSTL     , TRUSTR 
C
      PARAMETER (BIGNUM =  1, EIGL   =  2, EIGR   =  3, ENDL   =  4, 
     &           ENDR   =  5, EPS    =  6, EPS1   =  7, GRNRM  =  8, 
     &           ORIGIN =  9, RADIUS = 10, REPS   = 11, SFARL  = 12,
     &           SFARR  = 13, THETA0 = 14, THETAL = 15, THETAR = 16, 
     &           THRSH  = 17, TRUSTL = 18, TRUSTR = 19)
C
C.... RUSED defines the number of positions already used in RSTOR ......
C
      INTEGER    RUSED
      PARAMETER (RUSED = TRUSTR)
C
C.... addresses for integer scalars (ISTOR) ............................
C
C     INDSI  : index of the subinterval to be closed
C     JL     : current number of steps
C     JLMAX  : maximum number of steps
C     JT     : current dimension of the block tridiagonal matrix
C     JTMAX  : maximum dimension of the block tridiagonal matrix
C     JTMIN  : minimum dimension of the block tridiagonal matrix
C     LCOMM  : communicator for the parallel version
C     LEIG   : leading dimension of (EIG)
C     LFILE  : file unit for output
C     LNI    : leading dimension of (U), (V) and (X)
C     LPRNT  : level of printing
C     LRERR  : code for error messages (see subroutine lzerrs.f)
C     LRMDE  : Lanczos run mode
C              = 0 : exit
C              = 1 : origin in EIGL=EIGR
C              = 2 : origin in EIGL, moving to EIGR
C              = 3 : origin in EIGR, moving to EIGL
C              = 4 : splitting a subinterval
C              = 5 : restart without changing SIGMA
C              = 6 : restart with a new SIGMA
C              = 7 : restart with EIGL
C              = 8 : restart with EIGR
C              = 9 : standard problem
C     LRWRN  : code for warning messages (see subroutine lzwrns.f)
C     LTAU   : leading dimension of (TAU)
C     MYPE   : process rank in the parallel version
C     N      : dimension of the eigenvalue problem
C     NBX    : number of vectors in (X) multiplied by (B) 
C     NBXMAX : maximum number of vectors in (BX) 
C     NDEIG  : number of eigenpairs required in the run 
C     NEPIN  : number of eigenpairs given as input
C     NEWSIG : flag for a new SIGMA
C              = 0 : exit
C              = 1 : no need for a new SIGMA
C              = 2 : cost-effectiveness point reached
C              = 3 : potential ill-conditioning detected
C              = 4 : all solutions in subinterval have been found 
C     NFARL  : number of eigenvalues less than SFARL
C     NFARR  : number of eigenvalues less than SFARR
C     NI     : dimension of the vectors in (U), (V) and (X)
C     NMOPA  : number of multiplications op(A)*vector performed
C     NMOPB  : number of multiplications op(B)*vector performed
C     NNSPNT : number of eigenvalues less than ORIGIN
C     NNTRTL : number of eigenvalues less than TRUSTL
C     NNTRTR : number of eigenvalues less than TRUSTR
C     NONEWS : number of runs without any converged solution
C     NPE    : number of processes in the parallel version
C     NPORTH : number of partial reorthogonalizations performed
C     NQMAX  : maximum number of vectors in (BASIS)
C     NREIG  : number of required eigenpairs  
C     NREIGL : number of required eigenvalues less    than EIGR
C     NREIGR : number of required eigenvalues greater than EIGL
C     NRITZ  : number of Ritz values stored in (RITZ)
C     NRUN   : current number of runs
C     NRUNMX : maximum number of runs
C     NSFAIL : number of factorizations discarded
C     NSIGMA : number of origin translations
C     NSIMAX : maximum number of subintervals
C     NSINT  : number of subintervals
C     NSLOG  : number of subintervals recorded in (SSLOG)
C     NSORTH : number of selective orthogonalizations performed
C     NSRLS  : number of eigenvalues required less    than SIGMA
C     NSRRS  : number of eigenvalues required greater than SIGMA
C     NSVIN  : number of starting vectors given as input
C     NTEIG  : number of computed eigenpairs
C     NULLDQ : number of zero diagonal entries in (BETAQ)
C     NULLDR : number of zero diagonal entries in (BETAR)
C     NVB    : number of vectors in a block
C     NWBSY  : memory used
C     NWMAX  : maximum memory needed
C     NWMIN  : minimum memory needed
C     NXMAX  : maximum number of vectors in (BX) and (X)
C
C.... indirect addresses for real arrays ...............................
C
C     ITIME  --> TIME  : time table,
C                        dimension (15)
C                ( 1)  : total time spent with op(A)*vectors
C                ( 2)  : total time spent with op(B)*vectors
C                ( 3)  : total time spent with factorizations
C                ( 4)  : total time spent with vectors generation
C                ( 5)  : total time spent with reorthogonalizations
C                ( 6)  : total time spent with reduced eigenproblem
C                ( 7)  : total time spent computing Ritz vectors
C                ( 8)  : time spent with last factorization
C                ( 9)  : time spent with reorthogonalizations
C                (10)  : time spent with reduced eigenproblem
C                (11)  : records time upon return from blzdr_
C                (12)  : records starting time for current run
C                (13)  : records starting time when LFLAG.EQ.0
C                (14)  : reserved for future use
C                (15)  : reserved for future use
C     IRSINT --> RSINT : lower and upper limits of each subinterval,
C                        dimension (6,NSIMAX)
C                (1,:) : xi_L
C                (2,:) : farthest Ritz value to the right of xi_L
C                (3,:) : new translation     to the right of xi_L
C                (4,:) : new translation     to the left  of xi_R
C                (5,:) : farthest Ritz value to the left  of xi_R
C                (6,:) : xi_R
C     ISSLOG --> SSLOG : spectrum slicing history, 
C                        dimension (8,NRUNMX)
C                (1,:) : SIGMA
C                (2,:) : ENDL 
C                (3,:) : ENDR 
C                (4,:) : time spent in the run
C                (5,:) : number of eigenvalues to the left of SIGMA
C                (6,:) : number of eigenvalues computed in the run
C                (7,:) : number of steps perfomed in the run
C                (8,:) : number of vectors computed in the run
C     IRITZ  --> RITZ  : Ritz values,
C                        dimension (JTMAX,2)
C     ITB    --> TB    : block tridiagonal matrix,
C                        dimension (JTMAX,NVB+1)
C     IALPHA --> ALPHA : matrix (Q')*(B)*(R) at step JL,
C                        dimension (NVB,NVB)
C     IBETAQ --> BETAQ : matrix (BETA) in (R)=(Q)*(BETA) at step JL-1,
C                        dimension (NVB,NVB)
C     IBETAR --> BETAR : matrix (BETA) in (R)=(Q)*(BETA) at step JL,
C                        dimension (NVB,NVB)
C     IANORM --> ANORM : extreme singular values of (ALPHA),
C                        dimension (JLMAX,2)
C     IBNORM --> BNORM : extreme singular values of (BETA),
C                        dimension (JLMAX,2)
C     IETA   --> ETA   : orthog. bounds among (R) and Lanczos vectors,
C                        dimension (JLMAX,2)
C     ITAU   --> TAU   : orthog. bounds among (R) and Ritz    vectors,
C                        dimension (LTAU,2)
C     IR     --> R     : work array for Lanczos vectors,
C                        dimension (NI*NVB,4)
C     ITHETA --> THETA : eigenvalues of (TB) and estimated residuals
C                        dimension (JTMAX,2)
C     IS     --> S     : eigenvectors of (TB) and workspace,
C                        dimension (JTMAX,JTMAX)
C     IBASIS --> BASIS : Lanczos basis,
C                        dimension (NI*NVB,(NQMAX+k)*2) if LOPTS(1) > 0
C                        dimension (NI*NVB,(NQMAX+k) if LOPTS(1) = 0
C                        k = max(0,min(1,JLMAX-NQMAX))
C     IBX    --> BX    : (B)*eigenvectors,
C                        dimension (NI,NBXMAX+k) if LOPTS(1) > 0
C                        dimension (0) if LOPTS(1) = 0
C                        k = max(0,min(1,LTAU-NBXMAX)
C     IRWORK --> RWORK : beginning of the workspace in RSTOR,
C                        dimension (JTMAX*(2+JTMAX+max(NVB+1,18))
C     IIWORK --> IWORK : beginning of the workspace in ISTOR,
C                        dimension (JTMAX*12)
C
C.... addresses for integer arrays .....................................
C
C     INDR  : indirect addresses for (R),
C             dimension (4)
C             (1) : (B)*(Q), index j-1
C             (2) : (B)*(Q), index j
C             (3) : (Q), index j-1
C             (4) : (Q), index j
C     LOPTS : options for the algorithm,
C             dimension (4)
C             (1) : problem type flag 
C             (2) : spectrum slicing flag
C             (3) : eigenvectors purification flag
C             (4) : reserved for future use
C     LBLAS : BLAS levels,
C             dimension (4)
C             (1) : (R):=(R)-(Q)*(BETA')
C             (2) : (ALPHA):=(Q')*(B)*(R)
C             (3) : (R):=(R)-(Q)*(ALPHA)           
C             (4) : (R)=(Q)*(BETA)
C     FHNDL : file handle,
C             dimension (4)
C             (1) : blzpack.__.BQ 
C             (2) : blzpack.__.BX
C             (3) : blzpack.__.Q
C             (4) : blzpack.__.X
C     ISINT : inertia of the lower and upper limits defined in RSINT,
C             dimension (2,NSIMAX)
C             (1,:) : number of eigenvalues to the left of xi_L
C             (2,:) : number of eigenvalues to the left of xi_R
C
      INTEGER    INDSI      , JL         , JLMAX      , JT         ,
     &           JTMAX      , JTMIN      , LCOMM      , LEIG       , 
     &           LFILE      , LNI        , LPRNT      , LRERR      ,
     &           LRMDE      , LRWRN      , LTAU       , MYPE       ,
     &           N          , NBX        , NBXMAX     , NDEIG      , 
     &           NEPIN      , NEWSIG     , NFARL      , NFARR      , 
     &           NI         , NMOPA      , NMOPB      , NNSPNT     , 
     &           NNTRTL     , NNTRTR     , NONEWS     , NPE        ,
     &           NPORTH     , NQMAX      , NREIG      , NREIGL     , 
     &           NREIGR     , NRITZ      , NRUN       , NRUNMX     , 
     &           NSFAIL     , NSIGMA     , NSIMAX     , NSINT      , 
     &           NSLOG      , NSORTH     , NSRLS      , NSRRS      , 
     &           NSVIN      , NTEIG      , NULLDQ     , NULLDR     , 
     &           NVB        , NWBSY      , NWMAX      , NWMIN      , 
     &           NXMAX      , ITIME      , IRSINT     , ISSLOG     , 
     &           IRITZ      , ITB        , IALPHA     , IBETAQ     , 
     &           IBETAR     , IANORM     , IBNORM     , IETA       , 
     &           ITAU       , IR         , ITHETA     , IS         , 
     &           IBASIS     , IBX        , IRWORK     , IIWORK     , 
     &           INDR       , LOPTS      , LBLAS      , FHNDL      , 
     &           ISINT
C
      PARAMETER (INDSI  =  1, JL     =  2, JLMAX  =  3, JT     =  4,
     &           JTMAX  =  5, JTMIN  =  6, LCOMM  =  7, LEIG   =  8, 
     &           LFILE  =  9, LNI    = 10, LPRNT  = 11, LRERR  = 12,
     &           LRMDE  = 13, LRWRN  = 14, LTAU   = 15, MYPE   = 16,
     &           N      = 17, NBX    = 18, NBXMAX = 19, NDEIG  = 20,
     &           NEPIN  = 21, NEWSIG = 22, NFARL  = 23, NFARR  = 24,
     &           NI     = 25, NMOPA  = 26, NMOPB  = 27, NNSPNT = 28,
     &           NNTRTL = 29, NNTRTR = 30, NONEWS = 31, NPE    = 32,
     &           NPORTH = 33, NQMAX  = 34, NREIG  = 35, NREIGL = 36,
     &           NREIGR = 37, NRITZ  = 38, NRUN   = 39, NRUNMX = 40, 
     &           NSFAIL = 41, NSIGMA = 42, NSIMAX = 43, NSINT  = 44, 
     &           NSLOG  = 45, NSORTH = 46, NSRLS  = 47, NSRRS  = 48, 
     &           NSVIN  = 49, NTEIG  = 50, NULLDQ = 51, NULLDR = 52, 
     &           NVB    = 53, NWBSY  = 54, NWMAX  = 55, NWMIN  = 56, 
     &           NXMAX  = 57, ITIME  = 58, IRSINT = 59, ISSLOG = 60, 
     &           IRITZ  = 61, ITB    = 62, IALPHA = 63, IBETAQ = 64, 
     &           IBETAR = 65, IANORM = 66, IBNORM = 67, IETA   = 68, 
     &           ITAU   = 69, IR     = 70, ITHETA = 71, IS     = 72, 
     &           IBASIS = 73, IBX    = 74, IRWORK = 75, IIWORK = 76, 
     &           INDR   = 77, LOPTS  = 81, LBLAS  = 85, FHNDL  = 89, 
     &           ISINT  = 93)
C
C.... IUSED defines the number of positions already used in ISTOR ......
C
      INTEGER    IUSED
      PARAMETER (IUSED = ISINT)
