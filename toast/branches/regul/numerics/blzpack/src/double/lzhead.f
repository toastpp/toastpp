      SUBROUTINE LZHEAD (JLMAX,NVB,LEIG,N,LFILE,LRERR,NREIG,
     &                   NTEIG,EIGL,EIGR,SIGMA)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZHEAD prints the BLZPACK header                                 *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JLMAX (sii) : maximum number of steps                            *
C*    NVB   (sii) : number of vectors in a block                       *
C*    LEIG  (sii) : leading dimension of (EIG)                         *
C*    N     (sii) : dimension of the eigenvalue problem                *
C*    LFILE (sii) : file unit for output                               *
C*    LRERR (sii) : code for error messages                            *
C*    NREIG (sii) : number of required eigenpairs                      *
C*    NTEIG (sii) : number of starting eigenpairs                      *
C*    EIGL  (sri) : inferior bound for eigenvalues                     *
C*    EIGR  (sri) : superior bound for eigenvalues                     *
C*    SIGMA (sri) : origin translation                                 *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          JLMAX,LEIG,LFILE,LRERR,N,NREIG,NTEIG,NVB
      DOUBLE PRECISION EIGL,EIGR,SIGMA
C
C**** executable statements ********************************************
C
      WRITE (LFILE,1000)  
C
      IF ( LRERR .EQ. 0 ) WRITE (LFILE,1001) 
     &                    N,NREIG,NTEIG,NVB,LEIG,JLMAX,EIGL,EIGR,SIGMA
C
      RETURN 
C
 1000 FORMAT (  '***********************************************',/
     &          '*                                             *',/
     &          '*                   BLZPACK                   *',/
     &          '*                 ===========                 *',/
     &          '*                                             *',/
     &          '*     Block Lanczos Algorithm Eigensolver     *',/
     &          '*          (real symmetric matrices)          *',/
     &          '*                                             *',/
     &          '*             - release 1999.11 -             *',/
     &          '*                                             *',/
     &          '***********************************************')
 1001 FORMAT (/,'dimension of the problem ........ =',      I08  ,/ 
     &          'number of required eigenpairs ... =',      I08  ,/ 
     &          'number of starting eigenpairs ... =',      I08  ,/ 
     &          'number of vectors in a block .... =',      I08  ,/ 
     &          'maximum number of eigenpairs .... =',      I08  ,/
     &          'maximum number of steps ......... =',      I08  ,/
     &          'lower limit for eigenvalues ..... =',SP,1P,E12.4,/
     &          'upper limit for eigenvalues ..... =',SP,1P,E12.4,/
     &          'starting point .................. =',SP,1P,E12.4) 
C
C**** end of LZHEAD ****************************************************
C
      END
