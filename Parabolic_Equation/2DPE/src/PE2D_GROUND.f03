MODULE PE2D_GROUND
 USE PE2D_TYPE
 IMPLICIT NONE
 !-----------------------------------------------------------
 CONTAINS
 !-----------------------------------------------------------
 !Gaussian Hill (Default)
 FUNCTION GHILL(H0,X0,S,X) RESULT(GH)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: H0, x0, s, x
 REAL(DP) :: GH(3)
 GH(1) = H0*EXP(-(X-X0)**2/(2*S**2))
 GH(2) = -(H0/S**2)*(X-X0)*EXP(-(X-X0)**2/(2*S**2))
 GH(3) = (H0/S**2)*((X-X0)**2/(S**2)-1)*EXP(-(X-X0)**2/(2*S**2))
 END FUNCTION GHILL
 !-----------------------------------------------------------
 !Symmetric Triangular Hill
 FUNCTION STHILL(H0,X0,D,X) RESULT(STH)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: H0, X0, D, X
 REAL(DP) :: STH
 STH = 0.
 IF ((X.LE.(X0-D)).OR.(X.GE.(X0+D))) THEN
  STH = 0
 ELSEIF (X.LE.X0) THEN
  STH = H0*((X-X0)/D+1)
 ELSEIF (X.GT.X0) THEN
  STH = H0*((X0-X)/D+1)
 ENDIF
 END FUNCTION STHILL
 !-----------------------------------------------------------
 !Unsymmetric Triangular Hill
 FUNCTION THILL(H0,X0,D1,D2,X) RESULT(TH)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: H0, X0, D1, D2, X
 REAL(DP) :: TH
 TH = 0.
 IF ((X.LE.(X0-D1)).OR.(X.GE.(X0+D2))) THEN
  TH = 0
 ELSEIF (X.LE.X0) THEN
  TH = H0*((X-X0)/D1+1)
 ELSEIF (X.GT.X0) THEN
  TH = H0*((X0-X)/D2+1)
 ENDIF 
 END FUNCTION THILL
 !-----------------------------------------------------------     
 !Rectangular Hill
 FUNCTION RHILL(H0,X1,X2,X) RESULT(RH)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: H0, X1, X2, X
 REAL(DP) :: RH
 RH = 0
 IF ((X.LT.X1).OR.(X.GT.X2)) THEN
  RH = 0
 ELSE
  RH = H0
 ENDIF
 END FUNCTION RHILL
 !-----------------------------------------------------------
 !Slope
 FUNCTION SLOPE(H1,H2,X1,X2,X) RESULT(SLP)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: H1, H2, X1, X2, X
 REAL(DP) :: SLP
 SLP = 0
 IF ((X.LT.X1).OR.(X.GT.X2)) THEN
  SLP = 0
 ELSE
  SLP = ((H2-H1)*X+X2*H1-X1*H2)/(X2-X1)
 ENDIF 
 END FUNCTION SLOPE
 !-----------------------------------------------------------
 !FUNCTION REALTERR(ITERR) RESULT(TERR)
 !IMPLICIT NONE
 !INTEGER(I4B), INTENT(IN) :: ITERR
 !INTEGER(I4B) :: M, I
 !REAL(DP), DIMENSION() ::
 !END FUNCTION REALTERR
 !-----------------------------------------------------------
END MODULE PE2D_GROUND
