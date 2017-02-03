MODULE MATLIB
	USE NRTYPE
	IMPLICIT NONE
	CONTAINS
	!-----------------------------------------------------------
	SUBROUTINE FLIP(X)
	IMPLICIT NONE
	COMPLEX(DP), DIMENSION(:) :: X
	X = X(UBOUND(X,1):LBOUND(X,1):-1)
	END SUBROUTINE FLIP
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
	END FUNCTION EYE
	!-----------------------------------------------------------
	FUNCTION ONES(N) RESULT(M)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N
	COMPLEX(DP), DIMENSION(N) :: M
	INTEGER :: I
	DO I = 1,N
	        M(I) = 1._DP
	END DO
	END FUNCTION ONES
	!-----------------------------------------------------------
	FUNCTION DIAG(X) RESULT(M)
	IMPLICIT NONE
	COMPLEX(DP), DIMENSION(:) :: X
	COMPLEX(DP), DIMENSION(SIZE(X,1),SIZE(X,1)) :: M
	INTEGER :: I, N
	N = SIZE(X,1)
	M = 0._DP
	DO I = 1,N
	        M(I,I) = X(I)
	END DO
	END FUNCTION DIAG
	!-----------------------------------------------------------
	SUBROUTINE BAND(E,N,KL,KU,M)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N, KL, KU
	COMPLEX(DP), INTENT(IN), DIMENSION(:) :: E
	COMPLEX(DP), INTENT(OUT), DIMENSION(N,N) :: M
	INTEGER :: I, J, Q, ND
	ND = SIZE(E)
	IF (MOD(Nd,2).EQ.0) THEN
	        STOP 'NDIAG(E,N) : LEN(E) MUST BE ODD'
	END IF
	M = 0.
	Q = (Nd-1)/2          !Nd = 3 gives Q = 2
	DO J = 1,N
	        M(J,J) = E(Q+1)
	END DO
	DO I = 1,Q
	        DO J = 1,N-I
	        	M(J+I,J) = E(Q+1-I)
	        	M(J,J+I) = E(Q+1+I)
	        END DO
	END DO
	END SUBROUTINE
	!-----------------------------------------------------------
	SUBROUTINE D2B(M,N,KL,KU,B)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N, KL, KU
	COMPLEX(DP), INTENT(IN), DIMENSION(:,:) :: M
	COMPLEX(DP), INTENT(OUT), DIMENSION(KL+KU+1,N) :: B
	INTEGER :: I, J
	DO J = 1,N
                DO I = MAX(1,J-KU),MIN(N,J+KL)
			B(KU+1+I-J,J) = M(I,J)
          	END DO
        END DO
	END SUBROUTINE
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
	END FUNCTION INV
	!-----------------------------------------------------------
END MODULE
