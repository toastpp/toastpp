      SUBROUTINE LZPRT5 (JL,JT,NVB,LFILE,ORTH)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT5 prints the basis orthogonality level                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL    (sii) : number of steps                                    *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    NVB   (sii) : number of vectors in a block                       *
C*    LFILE (sii) : file unit for output                               *
C*    ORTH  (arw) : stores (Q')*(B)*(Q)                                *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    NINT                                                             *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          JL,JT,LFILE,NVB
      REAL             ORTH(JT,JT)
C
C==== local variables ==================================================
C
      INTEGER          I,J,K
      CHARACTER        FRMT1*126,FRMT2*17,FRMT3*17,FRMT4*37,FRMT5*43
C
C==== intrinsic function ===============================================
C
      INTRINSIC        NINT
C
C**** executable statements ********************************************
C
C.... set formats ......................................................
C
      WRITE (FRMT1,'(65X,I3,49X,I3)') JT*3,JT*3
C
      FRMT1(  1: 30) = '(/,''basis orthogonality check '
      FRMT1( 31: 65) = '(i=block, j=vector):'',/,6(''-''),''+'','
      FRMT1( 69:104) = '(''-''),/,'' i  j |  |log10(Q''''*B*Q)|'','
      FRMT1(105:117) = '/,6(''-''),''+'','
      FRMT1(121:126) = '(''-''))'
C
      WRITE (FRMT2,'(12X,I2)') JT
      WRITE (FRMT3,'(12X,I2)') JT
C
      FRMT2(  1: 12) = '(I2,''  1 |'','
      FRMT2( 15: 17) = 'I3)'  
      FRMT3(  1: 12) = '(3X,I2,'' |'','
      FRMT3( 15: 17) = 'I3)'  
C
      WRITE (FRMT4,'(12X,I3,17X,I2)') JT*3,JT
      WRITE (FRMT5,'(10X,I2,4X,I2,16X,I3)') JL,(NVB-1)*3,JT*3
C
      FRMT4(  1: 12) = '(6(''-''),''+'','
      FRMT4( 16: 32) = '(''-''),/,4X,''j |'','
      FRMT4( 35: 37) = 'I3)'
      FRMT5(  1: 10) = '(4X,''i |'','
      FRMT5( 13: 16) = '(I3,'
      FRMT5( 19: 34) = 'X),/,6(''-''),''+'','
      FRMT5( 38: 43) = '(''-''))'
C
C.... print the orthogonality level ....................................
C
      WRITE (LFILE,FRMT1)
C
      DO 20 I = 0,JL-1
         WRITE (LFILE,FRMT2) I+1,(NINT(ORTH(I*NVB+1,K)),K=1,JT)
         DO 10 J = 2,NVB
            WRITE (LFILE,FRMT3) J,(NINT(ORTH(I*NVB+J,K)),K=1,JT)
   10    CONTINUE
   20 CONTINUE
C
      WRITE (LFILE,FRMT4) ((I,I=1,NVB),J=1,JL)
      WRITE (LFILE,FRMT5) (J,J=1,JL)
C
      RETURN 
C
C**** end of LZPRT5 ****************************************************
C
      END
