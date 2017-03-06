MODULE PE2D_AUX
 USE PE2D_TYPE
 IMPLICIT NONE
 CONTAINS
 !-----------------------------------------------------------
 SUBROUTINE FLIP(X)
 IMPLICIT NONE
 COMPLEX(DP), DIMENSION(:) :: X
 X = X(UBOUND(X,1):LBOUND(X,1):-1)
 END SUBROUTINE
 !-----------------------------------------------------------
 SUBROUTINE UNIQUE(X,XCNT,XUNI,XDUP,INFO)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(IN) :: X
 REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: XUNI, XCNT, XDUP
 INTEGER :: I,NUM
 LOGICAL, DIMENSION(SIZE(X)) :: MASK, MASK1
 LOGICAL :: INFO
 MASK = .FALSE.
 ALLOCATE(XCNT(SIZE(X)))
 INFO = .FALSE.
 DO I=1,SIZE(X) 
  !COUNT THE NUMBER OF OCCURRENCES OF THIS ELEMENT:
  NUM = COUNT(X(I)==X)
  XCNT(I) = NUM
  IF (NUM==1) THEN
   !THERE IS ONLY ONE, FLAG IT:
   MASK(I) = .TRUE.
   MASK1(I) = .FALSE.
  ELSE
  !FLAG THIS VALUE ONLY IF IT HASN'T ALREADY BEEN FLAGGED:
   INFO = .TRUE.
   IF (.NOT. ANY(X(I)==X.AND.MASK)) MASK(I) = .TRUE.
  END IF
 END DO
 !RETURN ONLY FLAGGED ELEMENTS:
 ALLOCATE(XUNI(COUNT(MASK)),XDUP(COUNT(MASK1)))
 XUNI = PACK(X,MASK)
 XDUP = PACK(X,MASK1)
 !IF YOU ALSO NEED IT SORTED, THEN DO SO.
 ! FOR EXAMPLE, WITH SLATEC ROUTINE:
 !CALL ISORT (VEC_UNIQUE, [0], SIZE(VEC_UNIQUE), 1)
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
 !SUBROUTINE DENSE2SPARSE(M,N) RESULT(MSP)
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
 CALL ZGETRF(N,N,M,N,IPIV,INFO)
 IF (INFO/= 0) THEN
  STOP "WARNING : Input Matrix is singular !"
 END IF
 CALL ZGETRI(N,M,N,IPIV,WORK,N,INFO)
 IF (INFO/= 0) THEN
  STOP "WARNING : Input Matrix is singular !"
 END IF
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION COND(A) RESULT(K)
 IMPLICIT NONE
 COMPLEX(DP), DIMENSION(:,:), INTENT(IN) :: A
 COMPLEX(DP), DIMENSION(2*SIZE(A,1)) :: CWORK
 REAL(DP) :: WORK(SIZE(A,1)), RWORK(2*SIZE(A,1))
 REAL(DP) :: NORM, RCOND, K
 INTEGER(DP) :: N, INFO
 REAL(DP), EXTERNAL :: ZLANGE, ZLANGB
 N = SIZE(A,1)
 NORM = ZLANGE('I',N,N,A,N,WORK)
 CALL ZGECON('I',N,A,N,NORM,RCOND,CWORK,RWORK,INFO)
 IF (INFO/= 0) THEN
  STOP "[COND] WARNING : Input Matrix is singular !"
 END IF
 K = 1./RCOND
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION INTERP1(N,X,NDER,XA,YA) RESULT(YP)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN) :: N, NDER
 REAL(DP), DIMENSION(:), INTENT(IN) :: XA, YA, X
 REAL(DP), DIMENSION(SIZE(X),0:NDER) :: YP
 REAL(DP), DIMENSION(2*SIZE(X)) :: WORK
 REAL(DP), DIMENSION(SIZE(X)) :: C
 REAL(DP) :: YDER(NDER), YY, XX
 INTEGER(I4B) :: I, IERR, M
 M = SIZE(X)
 IF (.NOT.(SIZE(XA).EQ.SIZE(YA))) THEN
  STOP "[INTERP1] ERROR : XA and YA should be same size !"
 END IF
 CALL DPLINT(N,XA,YA,C)
 DO I=1,M
  XX = X(I)
  CALL DPOLVL(NDER,XX,YY,YDER,N,XA,C,WORK,IERR)
  YP(I,0) = YY
  YP(I,1:NDER) = YDER
 END DO
 END FUNCTION
 !-----------------------------------------------------------
 SUBROUTINE STOREARRAY(FILENAME,X,NLINES,NCOL)
 IMPLICIT NONE
 CHARACTER(LEN=*) :: FILENAME
 REAL(DP), DIMENSION(:,:), ALLOCATABLE :: X
 LOGICAL :: FILEEXISTS
 INTEGER(I4B) :: NLINES, NCOL, I, REASON
 INQUIRE(FILE=FILENAME,EXIST=FILEEXISTS)
 IF (.NOT.FILEEXISTS) THEN
  STOP "[STOREARRAY] WARNING : file doesn't exist !"
 END IF
 OPEN(UNIT=100,FILE=FILENAME,ACTION='READ')
 READ(100,*,IOSTAT=REASON) NLINES
 ALLOCATE(X(NLINES,NCOL))
 DO I=1,NLINES
  READ(100,*,IOSTAT=REASON) X(I,:) 
  IF (REASON.GT.0)  THEN
   STOP "[STOREARRAY] WARNING : Read error"
  ELSEIF (REASON.LT.0) THEN
   EXIT
  END IF
 END DO
 END SUBROUTINE
 !-----------------------------------------------------------
END MODULE PE2D_AUX
