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
 	SUBROUTINE APPENDC(X,Y)
  	IMPLICIT NONE
  	INTEGER :: NX, NY
  	COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: BUFR
  	COMPLEX(DP), DIMENSION(:), INTENT(IN) :: Y
  	COMPLEX(DP), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: X
  	NX = SIZE(X)
  	NY = SIZE(Y)
  	ALLOCATE(BUFR(NX+NY))
  	BUFR(1:NX) = X
  	BUFR(NX+1:NX+NY) = Y
  	DEALLOCATE(X)
  	CALL MOVE_ALLOC(BUFR,X)
 	END SUBROUTINE APPENDC
 	!-----------------------------------------------------------
 	SUBROUTINE APPENDI(X,Y)
  	IMPLICIT NONE
  	INTEGER :: NX, NY
  	INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFR
  	INtEGER, DIMENSION(:), INTENT(IN) :: Y
  	INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: X
  	NX = SIZE(X)
  	NY = SIZE(Y)
  	ALLOCATE(BUFR(NX+NY))
  	BUFR(1:NX) = X
  	BUFR(NX+1:NX+NY) = Y
  	DEALLOCATE(X)
  	CALL MOVE_ALLOC(BUFR,X)
 	END SUBROUTINE APPENDI
 	!-----------------------------------------------------------
 	FUNCTION EYE(N) RESULT(M)
  	IMPLICIT NONE
  	INTEGER, INTENT(IN) :: N
  	COMPLEX(DP), DIMENSION(N,N) :: M
  	INTEGER :: I,J
  	DO I = 1,N
   		DO J = 1,N
			IF (I.EQ.J) THEN
    				M(I,J) = 1
			ELSE
				M(I,J) = 0
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
   		M(I) = 1
  	END DO
 	END FUNCTION ONES
 	!-----------------------------------------------------------
 	FUNCTION DIAG(X) RESULT(M)
  	IMPLICIT NONE
  	COMPLEX(DP), DIMENSION(:) :: X
  	COMPLEX(DP), DIMENSION(SIZE(X,1),SIZE(X,1)) :: M
  	INTEGER :: I, N
  	N = SIZE(X,1)
	M = 0.
  	DO I = 1,N
   		M(I,I) = X(I)
  	END DO
 	END FUNCTION DIAG
 	!-----------------------------------------------------------
 	FUNCTION NDIAG(E,N) RESULT(M)
  	IMPLICIT NONE
  	INTEGER, INTENT(IN) :: N
  	COMPLEX(DP), INTENT(IN), DIMENSION(:) :: E
  	COMPLEX(DP), DIMENSION(N,N) :: M
  	INTEGER :: I, J, Q
  	Nd = SIZE(E)
	IF (MOD(Nd,2).EQ.0) THEN
		STOP 'NDIAG(E,N) : LEN(E) MUST BE ODD'
	END IF
	M = 0.
	Q = (Nd-1)/2		!Nd = 3 gives Q = 2
	DO J = 1,N
		M(J,J) = E(Q+1)
	END DO
  	DO I = 1,Q
		DO J = 1,N-I
			M(J+I,J) = E(Q+1-I)
			M(J,J+I) = E(Q+1+I)
  		END DO
	END DO
 	END FUNCTION
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
  	CALL ZGETRF(N,N,M,N,IPIV,INFO)
  	IF (INFO/= 0) THEN
   		STOP "WARNING : Input Matrix is singular !"
   	END IF
  	CALL ZGETRI(N,M,N,IPIV,WORK,N,INFO)
  	IF (INFO/= 0) THEN
   		STOP "WARNING : Input Matrix is singular !"
  	END IF
 	END FUNCTION INV
 	!-----------------------------------------------------------
! 	FUNCTION SPCONV(A) RESULT(SPA)
!  	IMPLICIT NONE
!  	INTEGER :: I, J, NI, NJ, NMAX
!  	COMPLEX(DP), DIMENSION(:,:), INTENT(IN) :: A
!  	TYPE(SPMAT) :: SPA
!  	NI = SIZE(A,1)
!  	NJ = SIZE(A,2)
! 	NMAX = NI*NJ
!  	ALLOCATE(SPA%ROW(0),SPA%COL(0),SPA%VAL(0))
!  	DO I = 1,NI
!   		DO J = 1,NJ
!    			IF (A(I,J) /= 0) THEN
!     				SPA%NELT = SPA%NELT+1
!     				CALL APPEND_INT(SPA%ROW,(/I/))
!     				CALL APPEND_INT(SPA%COL,(/J/))
!    				CALL APPEND_CPLX(SPA%VAL,(/A(I,J)/))
!    			END IF
!   		END DO
!  	END DO
! 	END FUNCTION SPCONV
 	!-----------------------------------------------------------
END MODULE
