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
 WRITE(*,*) "GTPE2_HT"
 WRITE(*,*) "Dimension  : 2D"
 WRITE(*,*) "Atmosphere : Homogeneous"
 WRITE(*,*) "Boundary  : Generalized Terrain"
 WRITE(*,*) "Impedance : Rigid"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Enter inputs L, f, zs, As, angle"
 WRITE(*,*) ""
 READ *, L, F, ZS, AS, ANGLE
 WRITE(*,*) ""
 WRITE(*,*) "Store output ? (Y or N)"
 WRITE(*,*) ""
 READ *, OUTPUT
 !-----------------------------------------------------------
 !----------------------------------------------------------
 C0     = 343.0_DP
 LMBDA0 = C0/F
 K0     = 2*PI*F/C0
 DZ     = LMBDA0/10
 DR     = LMBDA0/10
 H      = 100*LMBDA0
 HA     = 30*LMBDA0
 N      = FLOOR(H/DZ)
 M      = FLOOR(L/DR)
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Size of the system : ", N
 WRITE(*,*) "Number of steps : ", M
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 NS     = MAX(1,FLOOR(ZS/DZ))
 NA     = FLOOR(HA/DZ)
 A      = IM/(2.0_DP*K0*DZ**2)
 KU     = 1
 KL     = 1
 NDIAG  = KU+KL+1
 WIDE   = 0
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Memory allocation
 ALLOCATE(K(N), DK2(N), ALT(N), PHI(N,M+1), T(N,N), S(N,N), M1(N,N), &
 M2(N,N), I(N), P(N,M), LP(N,M), LP2(N,M), LPG(M), LPG2(M), E(N,N), &
 MB2(NDIAG,N), MB1(NDIAG,N), TEMP(N,1), M1T(N,N), M2T(N,N), IPIV(N), &
 TERR(-1:M+2), TERR1(-1:M+2), TERR2(-1:M+2), LIN(-2:M+2), GHX(3), MCOND(M,2), &
 MB1T(NDIAG,N), MB2T(KL+NDIAG,N), MT(NDIAG,N), MS(NDIAG,N), RGE(-1:M+2))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Atmosphere density and wave velocity
 DO NZ = 1,N
  K(NZ)   = K0
  DK2(NZ) = K(NZ)**2-K0**2
  ALT(NZ) = NZ*DZ
 END DO
 !Absorbing layer
 DO NZ=N-NABL,N THEN
  K(NZ)  = K0+AS*IM*(NZ-N+NABL)**2/NABL**2
 END
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Terrain Boundary condition
 LIN(-1) = 0
 THETA   = ANGLE*PI/180.0_DP
 X0      = L/2.0_DP
 S       = L/5.0_DP
 H0      = SQRT(EXP(1.0_DP)/2.0_DP)*S*ATAN(THETA)
 DO MX = 0,M
  GHX       = GHILL(H0,X0,S,MX*DR)
  TERR(MX)  = GHX(1)
  TERR1(MX) = GHX(2)
  TERR2(MX) = GHX(3)
  RGE(MX)   = MX*DR
  LIN(MX)   = LIN(MX-1)+(1+TERR1(MX))*DR
  !--------
  SLOPE(MX) = ATAN(TERR1(MX))
  !--------
  ALPHA(MX) = TERR1(MX)**2+1
  BETA(MX)  = 2*IM*K0*TERR1(MX)+TERR2(MX)
  XI(MX)    = TERR2(MX)-2*IM*K0*TERR1(MX)
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Starting field - Initial Condition 
 A0 = 1.3717_DP
 A2 = -0.3701_DP
 B  = 3.0_DP
 DO NZ = 1,N
 PHI(NZ,1) = SQRT(IM*K0)* &
             ((A0+A2*K0**2*(NZ*DZ-ZS)**2)*EXP(-K0**2*(NZ*DZ-ZS)**2/B)+ &
             (A0+A2*K0**2*(NZ*DZ+ZS)**2)*EXP(-K0**2*(NZ*DZ+ZS)**2/B))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Matrix building
 !-----------------------------------------------------------
 !Initialization (Dense Matrices) - T, D and E
 !---------
 FD1 = (/ -0.5_DP, 0.0_DP, 0.5_DP /)
 FD2 = (/ 1.0_DP, -2.0_DP, 1.0_DP /)
 T = BBAND(FD2,N,KL,KU)
 S = BBAND(FD1,N,KL,KU)
 G = DIAG(DK2,N)
 D = D2B(G,N,KL,KU)
 !-----------------------------------------------------------
 !System - M1 and M2
 !---------
 WRITE(*,*) "Wide-angle order : ", WIDE 
 PHI_G(-1:0) = 0.0_DP
 PHI_T(-1:0) = 0.0_DP
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Forward-marching procedure
 CALL CPU_TIME(TI)
 !-----------------------------------------------------------
 DO MX = 2,M+1
  WRITE(*,*) "Step :", MX-1, "out of", M
  !-------------------------------
  !Terrain integration
  !-------------------------------
  A_ALPHA = DR*(ALPHA(MX-2)/3+ALPHA(MX-1)/6)
  B_ALPHA = DR*(ALPHA(MX-2)/6+ALPHA(MX-1)/3)
  A_XI    = DR*(XI(MX-2)/3+XI(MX-1)/6)
  B_XI    = DR*(XI(MX-2)/6+XI(MX-1)/3)
  !-------------------------------
  !GTPE Coefficients updating
  !-------------------------------
  E2 = I+(1/(2*K0))*(1/(2*K0)-IM*DR)*D
  E1 = I+(1/(2*K0))*(1/(2*K0)+IM*DR)*D
  D2 = -(IM*BETA(MX-1)/(2*K0)+2*TERR1(MX-1)-B_XI)/(2*IM*K0*DZ)
  D1 = -(IM*BETA(MX-2)/(2*K0)+2*TERR1(MX-2)+A_XI)/(2*IM*K0*DZ)
  C2 = (IM*ALPHA(MX-1)/(2*K0)+B_ALPHA)/(2*IM*K0*DZ**2)
  C1 = (IM*ALPHA(MX-2)/(2*K0)-A_ALPHA)/(2*IM*K0*DZ**2)
  !-------------------------------
  !Boundary Conditions updating
  !-------------------------------
  EPS2 = 2.0_DP*DZ*COS(SLOPE(MX-1))
  DG2  =          (3.0_DP/EPS2)+(3.0_DP/(2.0_DP*DR)+IM*K0)*SIN(SLOPE(MX-1))
  DT2  = -IM*K(M)+(3.0_DP/EPS2)-(3.0_DP/(2.0_DP*DR)+IM*K0)*SIN(SLOPE(MX-1))
  UG2  = 4.0_DP/(DG2*EPS2)
  UT2  = 4.0_DP/(DT2*EPS2)
  VG2  = -UG2/4.0_DP
  VT2  = -UT2/4.0_DP
  WG2  =  2.0_DP*SIN(SLOPE(MX-1))/(DG2*DR)
  WT2  = -2.0_DP*SIN(SLOPE(MX-1))/(DT2*DR)
  YG2  = -WG2/4.0_DP
  YT2  = -WT2/4.0_DP
  !--------------
  EPS1 = 2.0_DP*DZ*COS(SLOPE(MX-2))
  DG1  =          (3.0_DP/EPS1)+(3.0_DP/(2.0_DP*DR)+IM*K0)*SIN(SLOPE(MX-2))
  DT1  = -IM*K(M)+(3.0_DP/EPS1)-(3.0_DP/(2.0_DP*DR)+IM*K0)*SIN(SLOPE(MX-2))
  UG1  = 4.0_DP/(DG1*EPS1)
  UT1  = 4.0_DP/(DT1*EPS1)
  VG1  = -UG1/4.0_DP
  VT1  = -UT1/4.0_DP
  WG1  =  2.0_DP*SIN(SLOPE(MX-2))/(DG1*DR)
  WT1  = -2.0_DP*SIN(SLOPE(MX-2))/(DT1*DR)
  YG1  = -WG1/4.0_DP
  YT1  = -WT1/4.0_DP
  !--------------
  PHIG(MX-1) = UG2*PHI(1,MX-1)+VG2*PHI(2,MX-1)  +WG2*PHIG(MX-2)+YG2*PHIG(MX-3)
  PHIT(MX-1) = UT2*PHI(N,MX-1)+VT2*PHI(N-1,MX-1)+WT2*PHIT(MX-2)+YT2*PHIT(MX-3)
  !-------------------------------
  !Matrix updating
  !-------------------------------
  !Right-hand side
  T2 = T
  S2 = S
  T2(KL+1,1)      = T2(KL+1,1)+UG2
  T2(KL+1,N)      = T2(KL+1,N)+UT2
  T2(KL,2)        = T2(KL,2)+VG2
  T2(KL+KU+1,N-1) = T2(KL+KU+1,N-1)+VT2
  S2(KL+1,1)      = S2(KL+1,1)-UG2
  S2(KL+1,N)      = S2(KL+1,N)+UT2
  S2(KL,2)        = S2(KL,2)-VG2
  S2(KL+KU+1,N-1) = S2(KL+KU+1,N-1)+VT2
  !Left-hand side
  T1 = T
  S1 = S
  T1(KL+1,1)      = T1(KL+1,1)+UG1
  T1(KL+1,N)      = T1(KL+1,N)+UT1
  T1(KL,2)        = T1(KL,2)+VG1
  T1(KL+KU+1,N-1) = T1(KL+KU+1,N-1)+VT1
  S1(KL+1,1)      = S1(KL+1,1)-UG1
  S1(KL+1,N)      = S1(KL+1,N)+UT1
  S1(KL,2)        = S1(KL,2)-VG1
  S1(KL+KU+1,N-1) = S1(KL+KU+1,N-1)+VT1
  !-------------------------------
  MB2(KL+1:KL+NDIAG,:) = C2*T2+D2*S2+E2
  MB1                  = C1*T1+D1*S1+E1
  !-------------------------------
  !Right hand side multiplication
  !-------------------------------
  CALL ZGBMV('N',N,N,KL,KU,A,MB1,KL+KU+1,PHI(1:N,MX-1),1,B,TEMP,1)
  TEMP(1) = TEMP(1)+C1*(WG1*PHIG(MX-3)+YG1*PHIG(MX-4))        &
                   -D1*(WG1*PHIG(MX-3)+YG1*PHIG(MX-4))/2.0_DP &
                   -C2*(WG2*PHIG(MX-2)+YG2*PHIG(MX-3))        &
                   +D2*(WG2*PHIG(MX-2)+YG2*PHIG(MX-3))/2.0_DP
  TEMP(N) = TEMP(N)+C1*(WG1*PHIT(MX-3)+YG1*PHIT(MX-4))        &
                   -D1*(WG1*PHIT(MX-3)+YG1*PHIT(MX-4))/2.0_DP &
                   -C2*(WG2*PHIT(MX-2)+YG2*PHIT(MX-3))        &
                   +D2*(WG2*PHIT(MX-2)+YG2*PHIT(MX-3))/2.0_DP
  WRITE(*,*) ".... Right-hand side multiplied"
  !-------------------------------
  !System solving
  !-------------------------------
  CALL ZGBSV(N,KL,KU,1,MB2,2*KL+KU+1,IPIV,TEMP,N,INFO)
  IF (INFO.NE.0) THEN
   STOP "ERROR : Matrix is singular"
  END IF
  WRITE(*,*) ".... System solved"
  !-------------------------------
  !Complex pressure
  !-------------------------------
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
 !Output
 IF (OUTPUT.EQ."Y") THEN
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

