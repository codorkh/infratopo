MODULE PE2D_FD
 IMPLICIT NONE
 CONTAINS
!-----------------------------------------------------------
 INTEGER FUNCTION FCTL(N)
  IMPLICIT NONE
  INTEGER :: I, N
  IF (N == 0) THEN
   FCTL = 1
  ELSE
   FCTL = 1
   DO I = 1,N
    FCTL = FCTL*I
   END DO
  END IF
 END FUNCTION FCTL
!----------------------------------------------------------- 
 FUNCTION COMB(K,N) RESULT(KCN)
  IMPLICIT NONE
  INTEGER :: K, N
  INTEGER :: NUM, DENOM
  REAL :: KCN
  NUM = FCTL(N)
  DENOM = FCTL(K)*FCTL(N-K)
  KCN = NUM/DENOM
 END FUNCTION COMB
!-----------------------------------------------------------
 FUNCTION FD(I) RESULT(COEFFS)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  INTEGER, DIMENSION(0:I) :: COEFFS
  INTEGER :: J, Q, R
  R = MOD(I,2)
  IF (R == 0) THEN
   Q = I/2
  ELSE
   Q = (I-1)/2
  ENDIF
  DO J = 0,I
   COEFFS(J) = ((-1)**(Q-J))*COMB(J,I)
  END DO
 END FUNCTION FD
!-----------------------------------------------------------
END MODULE PE2D_FD
