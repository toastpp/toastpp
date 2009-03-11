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
C==== arguments ========================================================
C
      INTEGER          HLFBND,JT,JTMAX          
      REAL             TB(JTMAX,HLFBND),TBREP(HLFBND,JT)
C
C==== local variables ==================================================
C
      INTEGER          I,J
C
C**** executable statements ********************************************
C
C.... this is the format required by SSBTRD ............................
C
      DO 20 I = 1,HLFBND
         DO 10 J = 1,JT
            TBREP(I,J) = TB(J,I)
   10    CONTINUE
   20 CONTINUE
C
      RETURN 
C
C**** end of TBCOPY ****************************************************
C
      END
