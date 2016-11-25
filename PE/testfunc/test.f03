PROGRAM TEST

 IMPLICIT NONE
 COMPLEX, PARAMETER :: i=(0,1)
 COMPLEX, DIMENSION(10) :: vec1, vec2
 INTEGER :: k
 
 vec1 = (/(k*(-1)**k*i+k, k=1,10,1)/)
 
 OPEN(unit=1,file="foo.txt")
 
 DO k=1,10
  WRITE(1,100) vec1(k)
 END DO
 100 FORMAT(F0.0,SP,F0.0,"i")
 
 COMPLEX :: c
 
 c = vec1(3)
 
END PROGRAM TEST