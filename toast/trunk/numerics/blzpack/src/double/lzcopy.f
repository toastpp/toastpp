      SUBROUTINE LZCOPY (LX,LY,N,NV,X,Y)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZCOPY copies (X) into (Y)                                       *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LX (sii) : leading dimension of (X)                              *
C*    LY (sii) : leading dimension of (Y)                              *
C*    N  (sii) : dimension of the vectors in (X) and (Y)               *
C*    NV (sii) : number of vectors to be copied                        *
C*    X  (ari) : input  vectors                                        *
C*    Y  (aro) : output vectors                                        *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    DCOPY                                                            *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LX,LY,N,NV
      DOUBLE PRECISION X(LX,NV),Y(LY,NV)
C
C==== local variables ==================================================
C
      INTEGER          I
C
C**** executable statements ********************************************
C
      DO 10 I = 1,NV 
         CALL DCOPY (N,X(1,I),1,Y(1,I),1)
   10 CONTINUE
C
      RETURN 
C
C**** end of LZCOPY ****************************************************
C
      END
