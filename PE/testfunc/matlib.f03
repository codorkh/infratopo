MODULE MATLIB
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP1 = KIND(0.d0)
 !-----------------------------------------------------------
 TYPE :: SPMAT
  INTEGER :: NELT = 0
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ROW, COL
  COMPLEX(DP1), ALLOCATABLE, DIMENSION(:) :: VAL
 END TYPE SPMAT
 !-----------------------------------------------------------
 CONTAINS
 !-----------------------------------------------------------
 SUBROUTINE FLIP(X)
  IMPLICIT NONE
  COMPLEX(DP1), DIMENSION(:) :: X
  X = X(UBOUND(X,1):LBOUND(X,1):-1)
 END SUBROUTINE
 !-----------------------------------------------------------
 SUBROUTINE APPEND_CPLX(X,Y)
  IMPLICIT NONE
  INTEGER :: NX, NY
  COMPLEX(DP1), ALLOCATABLE, DIMENSION(:) :: BUFR
  COMPLEX(DP1), DIMENSION(:), INTENT(IN) :: Y
  COMPLEX(DP1), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: X
  NX = SIZE(X)
  NY = SIZE(Y)
  ALLOCATE(BUFR(NX+NY))
  BUFR(1:NX) = X
  BUFR(NX+1:NX+NY) = Y
  DEALLOCATE(X)
  CALL MOVE_ALLOC(BUFR,X) 
 END SUBROUTINE
 !-----------------------------------------------------------
 SUBROUTINE APPEND_INT(X,Y)
  IMPLICIT NONE
  INTEGER :: NX, NY
  INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFR
  INTEGER, DIMENSION(:), INTENT(IN) :: Y
  INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: X
  NX = SIZE(X)
  NY = SIZE(Y)
  ALLOCATE(BUFR(NX+NY))
  BUFR(1:NX) = X
  BUFR(NX+1:NX+NY) = Y
  DEALLOCATE(X)
  CALL MOVE_ALLOC(BUFR,X) 
 END SUBROUTINE
 !-----------------------------------------------------------
 FUNCTION EYE(N) RESULT(eyeN)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  COMPLEX(DP1), DIMENSION(N,N) :: eyeN
  INTEGER :: r, c
  DO r = 1,N
   DO c = 1,N
    eyeN(r,c) = 0
   END DO
   eyeN(r,r) = 1
  END DO
 END FUNCTION EYE
 !-----------------------------------------------------------
 FUNCTION INV(A) RESULT(invA)
  IMPLICIT NONE
  COMPLEX(DP1), DIMENSION(:,:), INTENT(IN) :: A
  COMPLEX(DP1), DIMENSION(SIZE(A,1), SIZE(A,2)) :: invA
  COMPLEX(DP1), DIMENSION(SIZE(A,1)) :: WORK
  COMPLEX(DP1), DIMENSION(SIZE(A,1)) :: IPIV
  INTEGER :: rows, INFO
  invA = A
  rows = SIZE(A,1)
  CALL ZGETRF(rows,rows,invA,rows,IPIV,INFO)
  IF (INFO/= 0) THEN
   STOP "WARNING : Input Matrix is singular !"
   END IF
  CALL ZGETRI(rows,invA,rows,IPIV,WORK,rows,INFO)
  IF (INFO/= 0) THEN
   STOP "WARNING : Input Matrix is singular !"
  END IF
 END FUNCTION INV
 !-----------------------------------------------------------
 FUNCTION SPCONV(A) RESULT(SPA)
  IMPLICIT NONE
  INTEGER :: I, J, NI, NJ
  COMPLEX(DP1), DIMENSION(:,:), INTENT(IN) :: A
  TYPE(SPMAT) :: SPA
  NI = SIZE(A,1)
  NJ = SIZE(A,2)
  ALLOCATE(SPA%ROW(0),SPA%COL(0),SPA%VAL(0))
  DO I = 1,NI
   DO J = 1,NJ
    IF (A(I,J) /= 0) THEN
     SPA%NELT = SPA%NELT+1
     CALL APPEND_INT(SPA%ROW,(/I/))
     CALL APPEND_INT(SPA%COL,(/J/))
     CALL APPEND_CPLX(SPA%VAL,(/A(I,J)/))
    END IF
   END DO
  END DO
 END FUNCTION
END MODULE

