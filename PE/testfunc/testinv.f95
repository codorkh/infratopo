PROGRAM testinv
 IMPLICIT NONE

 INTEGER, PARAMETER :: DP = KIND(0.d0)
 REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: M, invM
 INTEGER :: i, j, N
 
 PRINT *, "Size"
 READ *, N

 ALLOCATE(M(N,N),invM(N,N))

 PRINT *, "Matrix to invert"
 READ *, ((M(i,j), j=1,N), i=1,N)

 invM = INV(M)
 
 DO i = 1, N
  PRINT *, invM(i,:)
 END DO
	
 CONTAINS

  FUNCTION INV(A) result(invA)

  IMPLICIT NONE

  REAL(DP), DIMENSION(:,:), INTENT(IN) :: A
  REAL(DP), DIMENSION(SIZE(A,1), SIZE(A,2)) :: invA
  REAL(DP), DIMENSION(SIZE(A,1)) :: WORK
  REAL(DP), DIMENSION(SIZE(A,1)) :: IPIV
  INTEGER :: rows, INFO

  invA = A
  rows = SIZE(A,1)

  CALL DGETRF(rows,rows,invA,rows,IPIV,INFO)

   IF (INFO/= 0) THEN
    STOP "WARNING : Input Matrix is singular !"
   END IF

   CALL DGETRI(rows,invA,rows,IPIV,WORK,rows,INFO)

   IF (INFO/= 0) THEN
    STOP "WARNING : Input Matrix is singular !"
   END IF

  END FUNCTION INV

END PROGRAM testinv


