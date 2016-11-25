PROGRAM TestEye
 
 IMPLICIT NONE
 INTEGER :: M
 INTEGER :: k
 REAL, ALLOCATABLE, DIMENSION(:,:) :: I
 
 PRINT *, "Size of the identity matrix :"
 READ *, M
 
 ALLOCATE(I(M,M))
 I = EYE(M)

 DO k = 1,M
  PRINT *, I(k,:)
 END DO
 
CONTAINS
 
 FUNCTION EYE(N)
 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL, DIMENSION(N,N) :: EYE
  INTEGER :: r, c
  
  DO r = 1,N
   DO c = 1,N
    EYE(r,c) = 0
   END DO
   EYE(r,r) = 1
  END DO
 END FUNCTION EYE
 
END PROGRAM TestEye

 


 