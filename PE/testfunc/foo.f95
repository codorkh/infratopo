PROGRAM FOO

 IMPLICIT NONE
 COMPLEX, DIMENSION(100,100) :: P, Q, R, S
 COMPLEX :: i=(0,1)
 INTEGER :: n, m
 
 OPEN(unit=1,file="foo.txt")
 
 P = 0.
 Q = 0.
 
 P = RESHAPE((/((n+m*i, n=1,100), m=1,100)/),(/100,100/))
 
 S = 0.
 
 DO n = 1,100
  DO m = 1,100
   WRITE(1, 1000, advance='no') P(n,m)
  END DO
  WRITE(1, *) ''
 END DO
 
 1000 FORMAT(' ',F5.0,SP,F0.0,"i",SS,4X)
END PROGRAM FOO