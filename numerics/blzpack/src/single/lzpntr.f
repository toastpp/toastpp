      SUBROUTINE LZPNTR (OFFSET,IW,JLMAX,JTMAX,NVB,LTAU,NI,
     &                   NQMAX,NBXMAX,NRUNMX,NSIMAX,LOPTS,
     &                   NWBSY,NWMAX,NWMIN)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPNTR defines pointers for the work array                       *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    OFFSET (sii) : initial pointer                                   *
C*    IW     (aio) : pointers for the work array                       *
C*    JLMAX  (sii) : maximum number of steps                           *
C*    JTMAX  (sii) : maximum dimension of the block tridiagonal matrix *
C*    NVB    (sii) : number of vectors in a block                      *
C*    LTAU   (sii) : leading dimension of (TAU)                        *
C*    NI     (sii) : dimension of the vectors in (BASIS) and (X)       *
C*    NQMAX  (sii) : maximum number of vectors in (BASIS)              *
C*    NBXMAX (sii) : maximum number of vectors in (BX)                 *
C*    NRUNMX (sii) : maximum number of runs                            *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    LOPTS  (sii) : options for `generalized', `slice' and `purify'   *
C*    NWBSY  (sio) : memory used                                       *
C*    NWMAX  (sio) : maximum memory needed                             *
C*    NWMIN  (sio) : minimum memory needed                             *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    MAX,MIN                                                          *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   IW(18),JLMAX,JTMAX,LOPTS(*),LTAU,NBXMAX,NI,NQMAX,
     &          NRUNMX,NSIMAX,NWBSY,NWMAX,NWMIN,NVB,OFFSET
C
C==== local variable ===================================================
C
      INTEGER   FACTOR,K,LIW16,LIW17,LIW18
C
C==== intrinsic functions ==============================================
C
      INTRINSIC MAX,MIN
C
C**** executable statements ********************************************
C
      IF ( LOPTS(1) .GT. 0 ) THEN
         FACTOR = 2
      ELSE
         FACTOR = 1
      END IF
C
C.... (TIME)  : timings ................................................
C
      IW( 1) = OFFSET
C                              
C.... (RSINT) : spectrum slicing strategy ..............................
C
      IW( 2) = IW( 1) + 15
C
C.... (SSLOG) : spectrum slicing history ...............................
C
      IW( 3) = IW( 2) + 6*NSIMAX
C
C.... (RITZ)  : eigenvalue approximations ..............................
C
      IW( 4) = IW( 3) + 8*NRUNMX
C
C.... (TB)    : block tridiagonal matrix ...............................
C
      IW( 5) = IW( 4) + JTMAX*2      
C
C.... (ALPHA) : matrix (ALPHA) .........................................
C
      IW( 6) = IW( 5) + JTMAX*(NVB+1)
C
C.... (BETAQ) : matrix (BETA) at Lanczos step i-1 ......................
C
      IW( 7) = IW( 6) + NVB*NVB
C
C.... (BETAR) : matrix (BETA) at Lanczos step i ........................
C
      IW( 8) = IW( 7) + NVB*NVB
C
C.... (ANORM) : extreme singular values of (ALPHA) .....................
C
      IW( 9) = IW( 8) + NVB*NVB
C
C.... (BNORM) : extreme singular values of (BETA) ......................
C
      IW(10) = IW( 9) + JLMAX*2
C
C.... (ETA)   : partial reorthogonalization bounds .....................
C
      IW(11) = IW(10) + JLMAX*2
C
C.... (TAU)   : selective orthogonalization bounds .....................
C
      IW(12) = IW(11) + JLMAX*2
C
C.... (R)     : Lanczos vectors generation .............................
C
      IW(13) = IW(12) + LTAU*2
C
C.... (THETA) : eigenvalues of (TB) ....................................
C
      IW(14) = IW(13) + NI*NVB*4
C
C.... (S)     : eigenvectors of (TB) ...................................
C
      IW(15) = IW(14) + JTMAX*2
C
C.... (BASIS) : Lanczos vectors ........................................
C
      IW(16) = IW(15) + JTMAX*JTMAX
C
C.... (BX)    : (B)*(X) ................................................
C
      K = MAX(0,MIN(1,JLMAX-NQMAX))
C
      IW(17) = IW(16) + NI*NVB*(NQMAX+K)*FACTOR
C
C.... (W)     : scratch area ...........................................
C
      K = MAX(0,MIN(1,LTAU-NBXMAX))
C
      IW(18) = IW(17) + NI*(NBXMAX+K)*(FACTOR-1) 
C
C.... mimimum, maximum and total memory used ...........................
C
      LIW16 = NI*NVB*FACTOR
      LIW17 = NI*(FACTOR-1)
      LIW18 = (2+JTMAX+MAX(NVB+2,18))*JTMAX
C
      NWMIN = ( IW(16)-1 ) + LIW16 + LIW17 + LIW18
      NWMAX = ( IW(16)-1 ) + LIW16*JLMAX + LIW17*LTAU + LIW18
      NWBSY = ( IW(18)-1 ) + LIW18
C
      RETURN
C
C**** end of LZPNTR ****************************************************
C
      END
