PROGRAM FD_TABLE
 IMPLICIT NONE
 INTEGER :: n_der, n_acc, n_p, n1, n2, m
 REAL :: alphad
 REAL, ALLOCATABLE, DIMENSION(:) :: alpha
 REAL, ALLOCATABLE, DIMENSION(:,:,:) :: delta

 PRINT *, "Enter order of derivative"
 READ *, n_der
 PRINT *, "Enter order of accuracy"
 READ *, n_acc
 ALLOCATE(alpha(0:n_acc),delta(0:n_acc,0:n_acc,0:n_der))
 PRINT *, "Enter grid points x = x(0) ... x(n_acc)"
 READ *, alpha

 alphad = 0
 CALL FDIFF(n_der,n_acc,alphad,alpha,delta)

 OPEN(UNIT=10,FILE="FD_TABLE.dat")
  DO m = 0,n_der
   DO n1 = 0,n_acc
    DO n2 = 0,n_acc
     WRITE(10, 100) delta(n1,n2,m)
    END DO
   END DO
  END DO

  100 FORMAT(F8.3)

CONTAINS

 SUBROUTINE FDIFF(ORDER,N,X0,X,Y)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ORDER, N
  REAL, INTENT(IN) :: X(0:N), X0
  REAL, DIMENSION(0:N,0:N,0:ORDER), INTENT(INOUT) :: Y
  REAL :: c1, c2, c3
  INTEGER :: n1, n2, m
  Y = 0.
  c1 = 1
  DO n1 = 1,N
   c2 = 1
   DO n2 = 0,n1-1
    c3 = X(n1)-X(n2)
    c2 = c2*c3
    !IF (n <= ORDER) THEN
     !Y(n1-1,n2,m) = 0
    !END IF
    DO m = 0,min(n1,ORDER)
     Y(n1,n2,m) = ((X(n1)-X0)*Y(n1-1,n2,m)-m*Y(n1-1,n2,m-1))/c3
    END DO
   END DO
   DO m = 0,min(n1,ORDER)
    Y(n1,n1,m) = (c1/c2)*(m*Y(n1-1,n1-1,m)-(X(n1-1)-X0)*Y(n1-1,n1-1,m))
   END DO
   c1 = c2
  END DO
 END SUBROUTINE FDIFF

END PROGRAM FD_TABLE
