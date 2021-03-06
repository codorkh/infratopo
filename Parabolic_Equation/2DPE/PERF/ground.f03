MODULE GROUND
	!Different ground profiles
 	USE NRTYPE
 	IMPLICIT NONE
 	!-----------------------------------------------------------
 	CONTAINS
 	!-----------------------------------------------------------
	!Gaussian Hill (Default)
 	FUNCTION GHILL(H0,X0,S,X) RESULT(H)
  	IMPLICIT NONE
  	REAL(DP), INTENT(IN) :: H0, x0, s, x
  	REAL(DP) :: H
  	H = H0*EXP(-(x-x0)**2/s**2)
 	END FUNCTION GHILL
 	!-----------------------------------------------------------
	!Symmetric Triangular Hill
 	FUNCTION STHILL(H0,X0,D,X) RESULT(H)
 	IMPLICIT NONE
 	REAL(DP), INTENT(IN) :: H0, X0, D, X
 	REAL(DP) :: H
  	IF ((X.LE.(X0-D)).OR.(X.GE.(X0+D))) THEN
    		H = 0
	ELSEIF (X.LE.X0) THEN
		H = H0*((X-X0)/D+1)
	ELSEIF (X.GT.X0) THEN
		H = H0*((X0-X)/D+1)
	ENDIF
	END FUNCTION STHILL
 	!-----------------------------------------------------------
	!Unsymmetric Triangular Hill
	FUNCTION THILL(H0,X0,D1,D2,X) RESULT(H)
 	IMPLICIT NONE
 	REAL(DP), INTENT(IN) :: H0, X0, D1, D2, X
 	REAL(DP) :: H
  	IF ((X.LE.(X0-D1)).OR.(X.GE.(X0+D2))) THEN
    		H = 0
	ELSEIF (X.LE.X0) THEN
		H = H0*((X-X0)/D1+1)
	ELSEIF (X.GT.X0) THEN
		H = H0*((X0-X)/D2+1)
	ENDIF 
	END FUNCTION THILL
	!-----------------------------------------------------------	
	!Rectangular Hill
	FUNCTION RHILL(H0,X1,X2,X) RESULT(H)
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: H0, X1, X2, X
	REAL(DP) :: H
	IF ((X.LT.X1).OR.(X.GT.X2)) THEN
		H = 0
	ELSE
		H = H0
	ENDIF
	END FUNCTION RHILL
 	!-----------------------------------------------------------
	!Slope
	FUNCTION SLOPE(H1,H2,X1,X2,X) RESULT(H)
 	IMPLICIT NONE
 	REAL(DP), INTENT(IN) :: H1, H2, X1, X2, X
 	REAL(DP) :: H
  	IF ((X.LT.X1).OR.(X.GT.X2)) THEN
		H = 0
	ELSE
		H = ((H2-H1)*X+X2*H1-X1*H2)/(X2-X1)
	ENDIF 
	END FUNCTION SLOPE
	!-----------------------------------------------------------
!	FUNCTION REALTERR
!	END FUNCTION REALTERR
 	!-----------------------------------------------------------
END MODULE GROUND
