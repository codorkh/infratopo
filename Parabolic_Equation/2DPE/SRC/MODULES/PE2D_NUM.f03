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
 INTEGER FUNCTION COMB(K,N) RESULT(KCN)
  IMPLICIT NONE
  INTEGER(I4B) :: K, N
  INTEGER(I4B) :: NUM, DENOM
  NUM = FCTL(N)
  DENOM = FCTL(K)*FCTL(N-K)
  KCN = NUM/DENOM
 END FUNCTION COMB
 !-----------------------------------------------------------
 ! FUNCTION POLEVAL(P,X,ND) RESULT(Y)
 !  INTEGER(I4B) :: N, I
 !  COMPLEX(DPC) :: P(0:N), X, Y
 !  Y(0) = P(N)
 !  Y(1:ND) = CMPLX(0.0_DP,0.0_DP,KIND=DPC)
 !  DO I = N-1,0,-1
 !   NND = MIN(ND,N+1-I)
 !   DO J = NND,1,-1
 !    Y(J) = Y(J)*X+Y(J-1)
 !   END DO
 !   Y(0) = Y(0)*X+P(I)
 !  END DO
 !  DO I = 0,ND
 !   Y(I)  = FCTL(I)*Y(I)
 !  END DO
 ! END FUNCTION
 !-----------------------------------------------------------
 FUNCTION POLMULT(P,Q) RESULT(R)
  INTEGER(I4B) :: N, M, I, J
  COMPLEX(DPC) :: P(:), Q(:)
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: R
  N = SIZE(P,1)
  M = SIZE(Q,1)
  ALLOCATE(R(0:M+N))
  R = CMPLX(0.0_DP,0.0_DP,KIND=DPC)
  DO I = 0,N
    DO J = 0,M
      R(I+J) = R(I+J)+P(I)*Q(J)
    END DO
  END DO
 END FUNCTION
 !-----------------------------------------------------------
 !FUNCTION POLDIV(P,X,ND) RESULT(DP)
 !
 !END FUNCTION
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
  C(J) = ((-1)**(Q-J))*REAL(COMB(J,I),KIND=DP)
 END DO
 END FUNCTION FDC
 !-----------------------------------------------------------
 FUNCTION FDW(M,N,X0,A) RESULT(D)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN) :: M, N
 REAL(DP), DIMENSION(0:N), INTENT(IN) :: A
 REAL(DP), DIMENSION(0:M,0:N,0:N) :: D        !D(I,J,K) = d^i_{j,k}
 REAL(DP), INTENT(IN) :: X0
 REAL(DP) :: C1, C2, C3
 INTEGER(I4B) :: I, J, K
 D = 0.0_DP
 D(0,0,0) = 1
 C1 = 1
 DO J = 1,N
  C2 = 1
  DO K = 0,J-1
   C3 = A(J)-A(K)
   C2 = C2*C3
   D(0,J,K) = (A(J)-X0)*D(0,J-1,K)/C3
   DO I = 1,MIN(J,M)
    D(I,J,K) = ((A(J)-X0)*D(I,J-1,K)-I*D(I-1,J-1,K))/C3
   END DO
  END DO
  D(0,J,J) = (C1/C2)*(A(J-1)-X0)*D(0,J-1,J-1)
  DO I = 1,MIN(J,M)
   D(I,J,J) = (C1/C2)*(I*D(I-1,J-1,J-1)-(A(J-1)-X0)*D(I,J-1,J-1))
  END DO
  C1 = C2
 END DO
 END FUNCTION
 !-----------------------------------------------------------
END MODULE PE2D_NUM
