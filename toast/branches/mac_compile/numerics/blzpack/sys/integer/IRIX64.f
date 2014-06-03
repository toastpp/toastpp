      LOGICAL FUNCTION SIBTST (INDEX,INTGR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose :                                                        *
C*                                                                     *
C*    SIBTST is an interface for a a bit test function                 *
C*                                                                     *
C*  - Arguments :                                                      *
C*                                                                     *
C*    INDEX (sii) : the bit index                                      *
C*    INTGR (sii) : the integer number                                 *
C*                                                                     *
C*  - Intrinsic Function :                                             *
C*                                                                     *
C*    BTEST                                                            *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   INDEX,INTGR
C
C==== intrinsic function ===============================================
C 
      INTRINSIC BTEST
C
C**** executable statements ********************************************
C
      SIBTST = BTEST(INTGR,INDEX)
C
      RETURN
C
C**** end of SIBTST ****************************************************
C
      END
      SUBROUTINE SILALG (LBLAS,N,NVB)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose :                                                        *
C*                                                                     *
C*    SILALG sets the level of the basic linear algebra subroutines    *
C*                                                                     *
C*  - Arguments :                                                      *
C*                                                                     *
C*    LBLAS (aio) : BLAS level setting                                 *
C*    N     (sii) : dimension of the problem                           *
C*    NVB   (sii) : number of vectors in a block                       *
C*                                                                     *
C*  - Intrinsic Function :                                             *
C*                                                                     *
C*    MIN                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   LBLAS(*),N,NVB
C
C==== intrinsic function ===============================================
C
      INTRINSIC MIN
C
C**** executable statements ********************************************
C
C.... subroutine RQBETA ................................................
C
      LBLAS(1) = MIN(2,NVB)
C
C.... subroutine QTBR ..................................................
C
      LBLAS(2) = MIN(3,NVB)
C
C.... subroutine RQALPH ................................................
C
      LBLAS(3) = MIN(3,NVB)
C
C.... subroutine MGRAMS ................................................
C
      LBLAS(4) = MIN(2,NVB)
C
      RETURN 
C
C**** end of SILALG ****************************************************
C
      END
