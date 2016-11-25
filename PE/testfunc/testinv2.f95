PROGRAM testinv

 IMPLICIT NONE

 INTEGER, PARAMETER :: DP = KIND(0.d0)
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: M, invM
 INTEGER :: i, j, N
 
 PRINT *, "Size"
 READ *, N

 ALLOCATE(M(N,N),invM(N,N))

 PRINT *, "Matrix to invert"
 READ *, ((M(i,j), j=1,N), i=1,N)

 CALL INV(M,invM)
 
 DO i = 1, N
  PRINT *, invM(i,:)
 END DO

CONTAINS

SUBROUTINE INV(A,B)

 IMPLICIT NONE
  
 COMPLEX(DP), DIMENSION(:,:), INTENT(IN) :: A
 COMPLEX(DP), DIMENSION(SIZE(A,1), SIZE(A,2)), INTENT(OUT) :: B
 COMPLEX(DP), DIMENSION(SIZE(A,1)) :: WORK
 COMPLEX(DP), DIMENSION(SIZE(A,1)) :: IPIV
 INTEGER :: rows, INFO
 
 B = A
 rows = SIZE(A,1)

 CALL ZGETRF(rows,rows,B,rows,IPIV,INFO)

 IF (INFO/= 0) THEN
  STOP "WARNING : Input Matrix is singular !"
 END IF

 CALL ZGETRI(rows,B,rows,IPIV,WORK,rows,INFO)

 IF (INFO/= 0) THEN
  STOP "WARNING : Input Matrix is singular !"
 END IF

END SUBROUTINE INV

END PROGRAM testinv
