MODULE GROUND
 	USE NRTYPE
 	USE NRINT
 	IMPLICIT NONE
 	!-----------------------------------------------------------
 	CONTAINS
 	!-----------------------------------------------------------
 	FUNCTION GHILL(H0,X0,S,X) RESULT(H)
  	IMPLICIT NONE
  	REAL(DP), INTENT(IN) :: H0, x0, s, x
  	REAL(DP) :: H
  	H = H0*EXP(-(x-x0)**2/s**2)
 	END FUNCTION GHILL
 	!-----------------------------------------------------------
 	FUNCTION THILL(H0,X0,D,X) RESULT(H)
 	IMPLICIT NONE
 	REAL(DP), INTENT(IN) :: H0, X0, D, X
 	REAL(DP) :: H
  IF (X.LE.(X0-D)) THEN
    X =
	END FUNCTION THILL
 	!-----------------------------------------------------------
	FUNCTION RHILL
	END FUNCTION RHILL
 	!-----------------------------------------------------------
	FUNCTION REALTERR
	END FUNCTION REALTERR
 	!-----------------------------------------------------------
END MODULE GROUND
