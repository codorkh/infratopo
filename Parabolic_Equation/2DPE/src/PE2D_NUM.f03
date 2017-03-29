MODULE PE2D_NUM
 USE PE2D_TYPE
 IMPLICIT NONE
 CONTAINS
 !-----------------------------------------------------------
 INTEGER FUNCTION FCTL(N)
 IMPLICIT NONE
 INTEGER(I4B) :: I, N
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
 INTEGER(I4B) :: K, N
 INTEGER(I4B) :: NUM, DENOM
 REAL(DP) :: KCN
 NUM = FCTL(N)
 DENOM = FCTL(K)*FCTL(N-K)
 KCN = NUM/DENOM
 END FUNCTION COMB
 !-----------------------------------------------------------
 FUNCTION FDC(I) RESULT(C)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN) :: I
 INTEGER(I4B) :: J, Q, R
 REAL(DP), DIMENSION(0:I) :: C
 R = MOD(I,2)
 IF (R == 0) THEN
  Q = I/2
 ELSE
  Q = (I-1)/2
 ENDIF
 DO J = 0,I
  C(J) = ((-1)**(Q-J))*COMB(J,I)
 END DO
 END FUNCTION FDC
 !-----------------------------------------------------------
 !FUNCTION FDW(M,N,X0,X) RESULT(T)
 !IMPLICIT NONE
 !REAL(DP), DIMENSION(:), INTENT(IN) :: X
 !REAL(DP), DIMENSION(0:M,0:N,0:N) :: T
 !REAL(DP), INTENT(IN) :: X0
 !INTEGER(I4B) :: M, N, I, J
 !-----------------------------------------------------------
END MODULE PE2D_NUM
