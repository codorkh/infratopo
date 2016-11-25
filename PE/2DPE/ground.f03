MODULE GROUND
 USE PARAM
 IMPLICIT NONE
 !-----------------------------------------------------------
 CONTAINS
 !-----------------------------------------------------------
 FUNCTION HILL(H0,x0,x,s) RESULT(H)
   IMPLICIT NONE
   REAL(DP), INTENT(IN) :: H0, x0, s, x
   REAL(DP) :: H
   H = H0*EXP(-(x-x0)**2/s**2)
 END FUNCTION
 !-----------------------------------------------------------
 FUNCTION PWLIN(alpha,X) RESULT(H)
   IMPLICIT NONE
   REAL(DP), DIMENSION(:), INTENT(IN) :: alpha, X
   REAL(DP) :: H
   INTEGER :: i, N, M
   N = SIZE(X,1)
   M = SIZE(alpha,1)
   IF (N /= M+1) THEN
   ENDIF
 END FUNCTION
 !-----------------------------------------------------------
END MODULE
