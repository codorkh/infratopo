MODULE MATLIB
 USE NRTYPE
 IMPLICIT NONE
 CONTAINS
 !-----------------------------------------------------------
 SUBROUTINE FLIP(X)
 IMPLICIT NONE
 COMPLEX(DP), DIMENSION(:) :: X
 X = X(UBOUND(X,1):LBOUND(X,1):-1)
 END SUBROUTINE
 !-----------------------------------------------------------
 FUNCTION EYE(N) RESULT(M)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 COMPLEX(DP), DIMENSION(N,N) :: M
 INTEGER :: I,J
 DO I = 1,N
  DO J = 1,N
   IF (I.EQ.J) THEN
    M(I,J) = 1._DP
   ELSE
    M(I,J) = 0._DP
   END IF
  END DO
 END DO
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION ONES(N) RESULT(X)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 COMPLEX(DP), DIMENSION(N) :: X
 INTEGER :: I
 DO I = 1,N
  X(I) = 1._DP
 END DO
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION DIAG(X,N) RESULT(M)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 COMPLEX(DP), DIMENSION(:) :: X
 COMPLEX(DP), DIMENSION(N,N) :: M
 INTEGER :: I
 M = 0._DP
 DO I = 1,N
  M(I,I) = X(I)
 END DO
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION BAND(E,N,KL,KU) RESULT(M)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N, KL, KU
 COMPLEX(DP), INTENT(IN), DIMENSION(:) :: E
 COMPLEX(DP), DIMENSION(N,N) :: M
 INTEGER :: I, J
 M = 0._DP
 DO J = 1,N
  DO I = MAX(1,J-KU),MIN(N,J+KL)
   M(I,J) = E(I)
  END DO
 END DO
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION DENSE2BAND(M,N,KL,KU) RESULT(MB)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N, KL, KU
 COMPLEX(DP), INTENT(IN), DIMENSION(:,:) :: M
 COMPLEX(DP), DIMENSION(KL+KU+1,N) :: MB
 INTEGER :: I, J
 DO J = 1,N
  DO I = MAX(1,J-KU),MIN(N,J+KL)
   MB(KU+1+I-J,J) = M(I,J)
  END DO
 END DO
 END FUNCTION
 !-----------------------------------------------------------
 !SUBROUTINE D2SP(M,N)
 !END SUBROUTINE 
 !-----------------------------------------------------------
 FUNCTION INV(A) RESULT(M)
 IMPLICIT NONE
 COMPLEX(DP), DIMENSION(:,:), INTENT(IN) :: A
 COMPLEX(DP), DIMENSION(SIZE(A,1), SIZE(A,2)) :: M
 COMPLEX(DP), DIMENSION(SIZE(A,1)) :: WORK
 COMPLEX(DP), DIMENSION(SIZE(A,1)) :: IPIV
 INTEGER :: N, INFO
 M = A
 N = SIZE(A,1)
 CALL ZGETRF(N,N,M,N,IPIV,INFO)          !LAPACK ROUTINE FOR LU FACTORIZATION
 IF (INFO/= 0) THEN
  STOP "WARNING : Input Matrix is singular !"
 END IF
 CALL ZGETRI(N,M,N,IPIV,WORK,N,INFO)     !INVERSE OF LU
 IF (INFO/= 0) THEN
  STOP "WARNING : Input Matrix is singular !"
 END IF
 END FUNCTION
 !-----------------------------------------------------------
END MODULE
