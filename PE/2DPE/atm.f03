MODULE ATM
  USE PARAM
  IMPLICIT NONE
  !-----------------------------------------------------------
  CONTAINS
  !-----------------------------------------------------------
  SUBROUTINE ATMLIN(k,c,c0,f,dz,a)
    IMPLICIT NONE
    REAL, INTENT(IN) :: f, c0, dz, a
    REAL, DIMENSION(:), INTENT(INOUT) :: k, c
    INTEGER :: i, N
    N = SIZE(c,1)
    DO i = 1,N
      c(i) = c0+i*dz*a
      k(i) = 2*PI*f/c(i)
    END DO
  END SUBROUTINE
  !-----------------------------------------------------------
  SUBROUTINE ATMLOG(k,c,c0,f,dz,z0,a)
    IMPLICIT NONE
    REAL, INTENT(IN) :: f, c0, dz, a, z0
    REAL, , DIMENSION(:), INTENT(INOUT) :: k, c
    INTEGER :: i, N
    N = SIZE(c,1)
    DO i = 1,N
      c(i) = c0+a*LOG(1+i*dz/z0)
      k(i) = 2*PI*f/c(i)
    END DO
  END SUBROUTINE
  !-----------------------------------------------------------
END MODULE
