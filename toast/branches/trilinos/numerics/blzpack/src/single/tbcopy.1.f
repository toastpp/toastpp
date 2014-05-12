      SUBROUTINE TBCOPY (HLFBND,JT,JTMAX,TB,TBREP)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBCOPY copies (TB) into (TBREP) for tridiagonalization           *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    HLFBND (sii) : half bandwidth of (TB)                            *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    JTMAX  (sii) : maximum dimension of the block tridiagonal matrix *
C*    TB     (ari) : block tridiagonal matrix                          *
C*    TBREP  (aro) : copy of (TB)                                      *
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
      INTEGER          HLFBND,JT,JTMAX          
      REAL             TB(JTMAX,HLFBND),TBREP(JT,HLFBND)
C
C==== local variables ==================================================
C
      INTEGER          I,J,K,L
C
C**** executable statements ********************************************
C
C.... this is the format required by TBBRED ............................
C
      K = 0
      L = HLFBND
C
      DO 30 I = 1,HLFBND
         DO 10 J = 1,JT-K
            TBREP(J+K,L) = TB(J,I)
   10    CONTINUE
         DO 20 J = 1,K
            TBREP(J,L) = ZERO
   20    CONTINUE
         K = K + 1
         L = L - 1
   30 CONTINUE
C
      RETURN 
C
C**** end of TBCOPY ****************************************************
C
      END
