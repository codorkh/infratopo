PROGRAM GFPE2_HF
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
 WRITE(*,*) "GFPE2_HT"
 WRITE(*,*) "Dimension 	: 	2D"
 WRITE(*,*) "Atmosphere	:	Homogeneous"
 WRITE(*,*) "Boundary 	:	Flat"
 WRITE(*,*) "Impedance	:	Rigid"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Enter inputs L, H, f, zs, As"
 WRITE(*,*) ""
 READ *, L, F, ZS, AS
 !----------------------------------------------------------
 !-----------------------------------------------------------
 C0     = 343.0_DP
 LMBDA0 = C0/F
 K0     = 2*PI*F/C0
 DZ     = LMBDA0/10
 DR     = LMBDA0/10
 DKZ    = 2*PI/(N*DZ)
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
 BETA   = 0.0_DP
 WIDE   = 0
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Memory allocation
 ALLOCATE(K(N), KZ(N), REF(N), PHASE(N), DK2(N), Z(N), DELTA(N), R(M), &
 BUFR(N), PHI(N,M+1), FPHI(N,M), FPHII(N,M), FPHIR(N,M), FPHIS(M), &
 WA(2,N), CH(2*N), PHI1(N,M), PHI2(N,M), P(N,M), LP(N,M), LPRMS(N,M), &
 P0(N,M), LPG(M), LPRMSG(M))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Atmosphere density and wave velocity
 DO NZ = 1,N
  C(NZ)     = C0
  K(NZ)     = K0
  DK2(NZ)   = K(NZ)**2-K0**2
  ALT(NZ)   = NZ*DZ
  PHASE(NZ) = EXP(IM*DR*(SQRT(K0**2-KZ(NZ)**2)-K0**2))
  DELTA(NZ) = EXP(IM*DR*DK2(NZ)/(2*K0))
 END DO
 !Absorbing layer (PML)
 DO NZ = N-NA,N
  K(NZ)     = K0+AS*IM*(NZ-N+NA)**2/NA**2
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
  GHX      = GHILL(H0,X0,S,MX*DR)
  TERR(MX)  = GHX(1)
  TERR1(MX) = GHX(2)
  TERR2(MX) = GHX(3)
  R(MX)     = MX*DR
  LIN(MX)   = LIN(MX-1)+(1+TERR1(MX))*DR
  SLOPE(MX) = ATAN(TERR1(MX))
 END DO
 !TERR1 = DIFF1(X,TERR)
 !TERR2 = DIFF2(X,TERR)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Starting field
 A0 = 1.3717_DP
 A2 = -0.3701_DP
 B  = 3.0_DP
 DO NZ = 1,N
  PHI(NZ,1) = SQRT(IM*K0)* &
              ((A0+A2*K0**2*(NZ*DZ-ZS)**2)*EXP(-K0**2*(NZ*DZ-ZS)**2/B) + &
              (A0+A2*K0**2*(NZ*DZ+ZS)**2)*EXP(-K0**2*(NZ*DZ+ZS)**2/B))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Forward-marching procedure
 CALL CPU_TIME(TI)
 CALL CFFTI1(N,WA,IFAC)      !IFAC AND WA DEPEND ONLY ON N -> USED IN CFFTF1 AND CFFTB1
 DO MX = 2,M+1
  WRITE(*, *) "STEP :", MX-1, "OUT OF", M
  !-------------------------------
  !-------------------------------
  !Surface term
  DO NZ = 1,N
   FPHIS(MX-1)  = FPHIS(MX-1)+PHI(NZ,MX-1)*DZ*EXP(-J*BETA*NZ*DZ)
  END DO
  !-------------------------------
  !-------------------------------
  !Forward Fourier transform
  TEMP          = PHI(:,MX-1)
  CALL CFFTF1(N,TEMP,CH,WA,IFAC)
  FPHII(:,MX-1) = TEMP                                          !Incident wave
  CALL FLIP(TEMP,FPHIR(:,MX-1))                                 !Reflected wave
  FPHI(:,MX-1)  = (FPHII(:,MX-1)+FPHIR(:,MX-1))*PHASE           !Total wave
  !Inverse Fourier Transform
  TEMP          = FPHI(:,MX-1)
  CALL CFFTB1(N,TEMP,CH,WA,IFAC)
  PHI1(:,MX-1)  = TEMP/N                                        !Normalize Fourier Transform
  PHI2(:,MX-1)  = 2*IM*BETA*FPHIS(MX-1)*EXP(-IM*BETA*Z)*EXP(IM*DR*(SQRT(K0**2-BETA**2)-K0**2))
  PHI(:,MX)     = DELTA*(PHI1(:,MX-1)+PHI2(:,MX-1))*EXP(IM*K0*Z*TERR2(MX-1))
  !-------------------------------
  !-------------------------------
  !Real pressure
  P(:,MX-1)     = EXP(IM*K0*(MX-1)*DR)*PHI(:,MX)*(1/SQRT((MX-1)*DR))
 END DO
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
  OPEN(UNIT=10,FILE="GFPE_LPg.dat")
  OPEN(UNIT=20,FILE="GFPE_LP.dat")
  OPEN(UNIT=30,FILE="GFPE_P.dat")
  OPEN(UNIT=50,FILE="Topo.dat")
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
  CALL SYSTEM('gnuplot -p plot.plt')
  CALL SYSTEM('gnuplot -p plot2.plt')
 END IF
 PRINT *, "Main CPU time (s) :", TF-TI
 PRINT *, "Source pressure P0 (dB) :", 20*LOG10(ABS(P0))
 !-----------------------------------------------------------
 !----------------------------------------------------------
END PROGRAM GFPE2_HF
