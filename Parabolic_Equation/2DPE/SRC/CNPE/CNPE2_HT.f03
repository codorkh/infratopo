PROGRAM CNPE2_HT
 USE PE2D_TYPE
 USE PE2D_VAR
 USE PE2D_AUX
 USE PE2D_GROUND
 IMPLICIT NONE
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) "========== Parabolic Equation Package =========="
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "University of Bristol"
 WRITE(*,*) "Mech. Eng. Department"
 WRITE(*,*) "Aerodynamics and Aeroacoustics Research Group"
 WRITE(*,*) ""
 WRITE(*,*) "Codor Khodr, July 2016"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "PE2D_HFR"
 WRITE(*,*) "Dimension  :  2D"
 WRITE(*,*) "Atmosphere : Homogeneous"
 WRITE(*,*) "Boundary  : Flat"
 WRITE(*,*) "Impedance : Rigid"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Enter inputs L, f, zs, As, angle"
 WRITE(*,*) ""
 READ *, L, F, ZS, ABL, ANGLE
 WRITE(*,*) ""
 WRITE(*,*) "Store output ? (Y or N)"
 WRITE(*,*) ""
 READ *, OUTPUT
 !-----------------------------------------------------------
 !----------------------------------------------------------
 C0 = 343.0_DP
 LMBDA0 = C0/F
 K0 = 2*PI*F/C0
 DZ = LMBDA0/10
 DR = LMBDA0/10
 H = 100*LMBDA0
 HA = 30*LMBDA0
 N = FLOOR(H/DZ)
 M = FLOOR(L/DR)
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Size of the system : ", N
 WRITE(*,*) "Number of steps : ", M
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 NS = MAX(1,FLOOR(ZS/DZ))
 NA = FLOOR(HA/DZ)
 A = IM/(2.0_DP*K0*DZ**2)
 KU = 1
 KL = 1
 NDIAG = KU+KL+1
 WIDE = 0
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Memory allocation
 ALLOCATE(K(N), DK2(N), ALT(N), PHI(N,M+1), T(N,N), D(N,N), M1(N,N), &
 M2(N,N), I(N,N), P(N,M), LP(N,M), LP2(N,M), LPG(M), LPG2(M), E(N,N), &
 MB2(NDIAG,N), MB1(NDIAG,N), TEMP(N,1), M1T(N,N), M2T(N,N), IPIV(N), &
 TERR(-1:M+2), TERR1(-1:M+2), TERR2(-1:M+2), LIN(-2:M+2), GHX(3), MCOND(M,2), &
 MB1T(NDIAG,N), MB2T(KL+NDIAG,N), X(-1:M+2))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Atmosphere density and wave velocity
 DO NZ = 1,N
  IF (NZ.GE.(N-NABL)) THEN
   K(NZ) = K0+ABL*IM*(NZ-N+NABL)**2/NABL**2
  ELSE
   K(NZ) = K0
  ENDIF
  DK2(NZ) = IM*(K(NZ)**2-K0**2)/2*K0
  ALT(NZ) = NZ*DZ
 END DO
 TAU1 = 4.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 TAU2 = -1.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 SIGMA1 = 4.0_DP/3.0_DP
 SIGMA2 = -1.0_DP/3.0_DP
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Terrain Boundary condition
 LIN(-2) = 0
 THETA = ANGLE*PI/180.0_DP
 X0 = L/2.0_DP
 S = L/5.0_DP
 H0 = SQRT(EXP(1.0_DP)/2.0_DP)*S*ATAN(THETA)
 DO MX = -1,M+2
  GHX = GHILL(H0,X0,S,MX*DR)
  TERR(MX) = GHX(1)
  !TERR1(MX) = GHX(2)
  !TERR2(MX) = GHX(3)
  X(MX) = MX*DR
  LIN(MX) = LIN(MX-1)+(1+TERR1(MX))*DR
 END DO
 TERR1 = DIFF1(X,TERR)
 TERR2 = DIFF2(X,TERR)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Starting field
 A0 = 1.3717_DP
 A2 = -0.3701_DP
 B = 3.0_DP
 DO NZ = 1,N
  PHI(NZ,1) = SQRT(IM*K0)* &
              ((A0+A2*K0**2*(NZ*DZ-ZS)**2)*EXP(-K0**2*(NZ*DZ-ZS)**2/B) + &
              (A0+A2*K0**2*(NZ*DZ+ZS)**2)*EXP(-K0**2*(NZ*DZ+ZS)**2/B))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Matrix building
 !-----------------------------------------------------------
 !Initialization (Dense Matrices) - T, D and E
 !---------
 T(:,:) = 0.0_DP
 DO NZ = 2,N-1 
  T(NZ,NZ-1) = 1.0_DP
  T(NZ,NZ) = -2.0_DP
  T(NZ,NZ+1) = 1.0_DP
 END DO
 T(1,1) = -2.0_DP+SIGMA1
 T(1,2) = 1.0_DP+SIGMA2
 T(N,N-1) = 1.0_DP+TAU2
 T(N,N) = -2.0_DP+TAU1
 E = DIAG(CMPLX(ALT,0.0_DP),N)
 D = DIAG(DK2,N)
 I = EYE(N)
 !-----------------------------------------------------------
 !System - M1 and M2
 !---------
 ALPHA = CMPLX(1.0_DP,0.0_DP,KIND=DP)
 BETA = CMPLX(0.0_DP,0.0_DP,KIND=DP)
 SELECT CASE (WIDE)
  CASE(0)
   BETA1 = DR/2.0_DP
   BETA2 = -DR/2.0_DP
  CASE(1)
   BETA1 = 1.0_DP/(2.0_DP*IM*K0)+DR/2.0_DP
   BETA2 = 1.0_DP/(2.0_DP*IM*K0)-DR/2.0_DP
  CASE DEFAULT
   STOP "[MAIN] : Enter valid wide-angle case"
 END SELECT
 WRITE(*,*) "Wide-angle order : ", WIDE 
 M1 = I+A*BETA1*T+BETA1*D
 M2 = I+A*BETA2*T+BETA2*D 
 MB1 = DENSE2BAND(M1,N,KL,KU)
 MB2 = DENSE2BAND(M2,N,KL,KU)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Forward-marching procedure
 CALL CPU_TIME(TI)
 !-----------------------------------------------------------
 DO MX = 2,M+1
  WRITE(*,*) "Step :", MX-1, "out of", M
  !-------------------------------
  !Matrix updating
  !-------------------------------
  !TERR1(MX-1) = (TERR(MX)-TERR(MX-1))/DR
  !TERR2(MX-1) = (TERR1(MX)-TERR1(MX-1))/DR       !Varying coefficient to be injected...
  MB1T = MB1
  MB1T(KU+1,:) = MB1(KU+1,:)-BETA1*IM*K0*TERR2(MX-1)*ALT
  MB2T(KL+1:KL+NDIAG,:) = MB2
  MB2T(KL+KU+1,:) = MB2(KU+1,:)-BETA2*IM*K0*TERR2(MX-1)*ALT
  !-------------------------------
  !Right hand side multiplication
  !-------------------------------
  CALL ZGBMV('N',N,N,KL,KU,ALPHA,MB1T,KL+KU+1,PHI(1:N,MX-1),1,BETA,TEMP(:,1),1)
  WRITE(*,*) ".... Right-hand side multiplied"
  !-------------------------------
  !System solving
  !-------------------------------
  !MTEMP2 = M2
  !CALL ZGESV(N,1,M2T,N,IPIV,TEMP,N,INFO)                                                         !General matrix form     
  CALL ZGBSV(N,KL,KU,1,MB2T,2*KL+KU+1,IPIV,TEMP,N,INFO)                                           !Band matrix form
  !CALL ZGTSV(N,1,MT2T(3,1:N-1),MT2T(2,1:N),MT2T(1,2:N),TEMP(:,1),N,INFO)                         !Tridiagonal matrix form
  IF (INFO.NE.0) THEN
   STOP "ERROR : Matrix is singular"
  END IF
  WRITE(*,*) ".... System solved"
  PHI(1:N,MX) = TEMP(:,1)
  P(1:N,MX-1) = EXP(IM*K0*(MX-1)*DR)*PHI(1:N,MX)/SQRT((MX-1)*DR)
  WRITE(*,*)
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 CALL CPU_TIME(TF)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 P0 = ABS(P(NS,1))
 DO MX = 1,M
  DO NZ = 1,N
   LP(NZ,MX) = 20.*LOG10(ABS(P(NZ,MX))/P0)
   IF (LP(NZ,MX).LE.-120) THEN
    LP(NZ,MX) = -120
   END IF
  END DO
 END DO
 LPG = LP(1,1:M)
 LP2 = LP/SQRT(2.0_DP)
 LPG2 = LP2(1,1:M)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 IF (OUTPUT.EQ."Y") THEN
 !Output
  OPEN(UNIT=10,FILE="../results/Topo_LPg.dat")
  OPEN(UNIT=20,FILE="../results/Topo_LP.dat")
  OPEN(UNIT=30,FILE="../results/Topo_P.dat")
  OPEN(UNIT=50,FILE="../results/Topo.dat")
  DO MX = 1,M
   WRITE(10,100) MX*DR, LPG(MX), LPG2(MX)
   WRITE(50,100) MX*DR, TERR(MX), TERR1(MX), TERR2(MX)
   DO NZ = 1,N
    WRITE(20,100) MX*DR, NZ*DZ+TERR(MX), LP(NZ,MX)
    WRITE(30,101) MX*DR, NZ*DZ+TERR(MX), ABS(P(NZ,MX))
   END DO
   WRITE(20,*)
   WRITE(30,*) 
  END DO
  100 FORMAT(3(3X,F12.3))
  101 FORMAT(2(3X,F12.3),3X,F12.3,SP,F12.3,SS,"i")
  !Plotting of results
  CALL SYSTEM('gnuplot -p plot_topo.plt')
  CALL SYSTEM('gnuplot -p plot_topo2.plt')
 END IF
 PRINT *, "Main CPU time (s) :", TF-TI
 PRINT *, "Source pressure P0 (dB) :", 20*LOG10(ABS(P0))
 !-----------------------------------------------------------
 !----------------------------------------------------------
END PROGRAM PE2D_HT_FD
!-----------------------------------------------------------
!-------------      EXTERNAL PROCEDURES      ---------------
!-----------------------------------------------------------

