      SUBROUTINE TBBETA (JL,JTMAX,NVB,BETA,TB) 
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBBETA inserts (BETA) into (TB)                                  *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL    (sii) : number of steps                                    *
C*    JTMAX (sii) : maximum dimension of the block tridiagonal matrix  *
C*    NVB   (sii) : number of vectors in a block                       *
C*    BETA  (ari) : matrix (BETA) in (R)=(Q)*(BETA) at step JL         *
C*    TB    (arb) : block tridiagonal matrix                           *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*    (TB) is constructed as:                                          *
C*                                                                     *
C*                | ALPHA_1    BETA_2'                         |       *
C*                | BETA_2    ALPHA_2    BETA_3'               |       *
C*                |            BETA_3   ALPHA_3    .           |       *
C*    (TB)      = |                   .         .    .         |       *
C*                |                     .         .    .       |       *
C*                |                       .         .  BETA_j' |       *
C*                |                           BETA_j  ALPHA_j  |       *
C*                                                                     *
C*    where                                                            *
C*                                                                     *
C*                | alpha_1,1  alpha_1,2  ... alpha_1,p |              *
C*                | alpha_2,1  alpha_2,2  ... alpha_2,p |              *
C*    (ALPHA_j) = |     .          .              .     |              *
C*                |     .          .              .     |              *
C*                | alpha_p,1  alpha_p,2  ... alpha_p,p |              *
C*                                                                     *
C*                | beta_1,1   beta_1,2   ... beta_1,p  |              *
C*                |     0      beta_2,2   ... beta_2,p  |              *
C*    (BETA_j)  = |     .          .              .     |              *
C*                |     .          .              .     |              *
C*                |     0          0      ... beta_p,p  |              *
C*                                                                     *
C*    The upper triangle of (TB) is then stored as follows:            *
C*                                                                     *
C*    TB(i,0) := TB(i,i)                                               *
C*    TB(i,1) := TB(i,i+1)                                             *
C*    TB(i,2) := TB(i,i+2)                                             *
C*        .                                                            *
C*        .                                                            *
C*    TB(i,p) := TB(i,i+p)                                             *
C*                                                                     *
C*    such that, at step JL and NVB=3 for example,                     *
C*                                                                     *
C*    |    ...        ...        ...        ...    |                   *
C*    |                                            |                   *
C*    | ALPHA(1,1) ALPHA(1,2) ALPHA(1,3) BETA(1,1) |  row JROW+1       *
C*    | ALPHA(2,2) ALPHA(2,3) BETA(1,2)  BETA(2,2) |  row JROW+2       *
C*    | ALPHA(3,3) BETA(1,3)  BETA(2,3)  BETA(3,3) |  row JROW+3       *
C*    |                                            |                   *
C*    |    ...        ...        ...        ...    |                   *
C*                                                                     *
C*          0          1          2         NVB       columns          *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          JL,JTMAX,NVB
      DOUBLE PRECISION BETA(NVB,NVB),TB(JTMAX,0:NVB)
C
C==== local variables ==================================================
C
      INTEGER          I,J,JCOL,JROW
C
C**** executable statements ********************************************
C
      JCOL = 0
      JROW = (JL-1)*NVB
C
      DO 20 I = 1,NVB
         DO 10 J = I,NVB
            TB(JROW+J,NVB-JCOL) = BETA(J-JCOL,J)
   10    CONTINUE
         JCOL = JCOL + 1
   20 CONTINUE
C
      RETURN
C
C**** end of TBBETA ****************************************************
C
      END
