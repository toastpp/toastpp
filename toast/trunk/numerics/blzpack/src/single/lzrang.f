      SUBROUTINE LZRANG (LRMDE,NDEIG,NPEIG,NREIG,NREIGL,NREIGR,
     &                   BIGNUM,EIGL,EIGR,ENDL,ENDR,RANGE,
     &                   SIGMA,THETA0,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZRANG checks the computational interval set by the user         *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LRMDE  (sio) : run mode                                          *
C*    NDEIG  (sio) : number of required eigenpairs                     *
C*    NPEIG  (sii) : number of eigenpairs given as input               *
C*    NREIG  (sii) : number of required eigenpairs                     *
C*    NREIGL (sio) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sio) : number of required eigenvalues greater than EIGL  *
C*    BIGNUM (sri) : big number                                        *
C*    EIGL   (sro) : inferior bound for eigenvalues                    *
C*    EIGR   (sro) : superior bound for eigenvalues                    *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    RANGE  (ari) : user computational interval                       *
C*    SIGMA  (sro) : origin translation                                *
C*    THETA0 (sri) : reference point for THETA                         *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      REAL             FUDGE,ZERO
      PARAMETER        (FUDGE=0.0010E0,ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          LRMDE,NDEIG,NPEIG,NREIG,NREIGL,NREIGR
      REAL             BIGNUM,EIGL,EIGR,ENDL,ENDR,RANGE(2),SIGMA,THETA0
      LOGICAL          GNRZD
C
C**** executable statements ********************************************
C
C.... eigenvalue range .................................................
C
      ENDL  = RANGE(1)
      ENDR  = RANGE(2)
      EIGL  = BIGNUM*(-1)
      EIGR  = BIGNUM*(+1)
C
C.... check the range for buckling mode ................................
C
      IF ( THETA0 .NE. ZERO ) THEN
         IF ( ENDL .EQ. ZERO ) ENDL = -FUDGE
         IF ( ENDR .EQ. ZERO ) ENDR = -FUDGE
      END IF
C
C.... first SIGMA ......................................................
C
      SIGMA = ENDL
C
      IF      ( .NOT.GNRZD ) THEN
C
C............ standard problem .........................................
C
              SIGMA  = ZERO
              LRMDE = 9
              NREIGL = 0
              NREIGR = 0
C
      ELSE IF ( ENDL .GT. ENDR ) THEN
C
C............ Lanczos run from right to left, interval [ENDR,ENDL] .....
C
              LRMDE = 3
              NREIGL = NREIG + NPEIG
              NREIGR = 0
C
              EIGL = ENDR
              EIGR = SIGMA   
C
      ELSE IF ( ENDL .LT. ENDR ) THEN
C
C............ Lanczos run from left to right, interval [ENDL,ENDR] .....
C
              LRMDE = 2
              NREIGL = 0
              NREIGR = NREIG + NPEIG
C
              EIGL = SIGMA   
              EIGR = ENDR
C
      ELSE 
C
C............ Lanczos run around ENDL (=ENDR) ..........................
C
              LRMDE = 1
              NREIGL = NREIG + NPEIG
              NREIGR = NREIG + NPEIG
C
      END IF
C
      NDEIG = NREIG
      ENDL  = EIGL
      ENDR  = EIGR
C
      RETURN 
C
C**** end of LZRANG ****************************************************
C
      END
