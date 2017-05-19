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
 WRITE(*,*) "CNPE2_HT"
 WRITE(*,*) "Dimension  :  	2D"
 WRITE(*,*) "Atmosphere :	Homogeneous"
 WRITE(*,*) "Boundary   : 	Shift Map"
 WRITE(*,*) "Impedance  : 	Rigid"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Enter inputs L, f, zs, As, angle"
 WRITE(*,*) ""
 READ *, L, F, ZS, AS, ANGLE
 WRITE(*,*) ""
 !-----------------------------------------------------------
 !----------------------------------------------------------
 C0     = 343.0_DP
 LMBDA0 = C0/F
 K0     = 2.0_DP*PI*F/C0
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
 !Atmosphere density and wave velocity
 ALLOCATE(C(N),K(N),DK2(N),ALT(N))
 DO NZ = 1,N
  C(NZ)   = C0
  K(NZ)   = K0
  DK2(NZ) = K(NZ)**2-K0**2
  ALT(NZ) = NZ*DZ
 END DO
 !Absorbing layer (PML)
 DO NZ = N-NA,N
  K(NZ)  = K0+AS*IM*(NZ-N+NA)**2/NA**2
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Terrain Boundary condition
 ALLOCATE(LIN(0:M),X(M),TERR(M),TERR1(M),TERR2(M),SLOPX(M))
 LIN(0)  = 0
 THETA   = ANGLE*PI/180.0_DP
 X0      = L/2.0_DP
 S       = L/5.0_DP
 H0      = SQRT(EXP(1.0_DP)/2.0_DP)*S*ATAN(THETA)
 DO MX = 1,M
  GHX      = GHILL(H0,X0,S,MX*DR)
  TERR(MX)  = GHX(1)
  TERR1(MX) = GHX(2)
  TERR2(MX) = GHX(3)
  X(MX)     = MX*DR
  LIN(MX)   = LIN(MX-1)+(1+TERR1(MX))*DR
  SLOPX(MX) = ATAN(TERR1(MX))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Starting field 
 ALLOCATE(PHI(N,M+1),PC(N,M),P(N,M),Q(N,M))
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
 !Matrix updating
 ALLOCATE(T(N,N),D(N,N),E(N,N),I(N,N),DD(N,N))
 TAU1    =  4.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 TAU2    = -1.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 SIGMA1  =  4.0_DP/3.0_DP
 SIGMA2  = -1.0_DP/3.0_DP
 !Matrix Initialization (Dense Matrices) - T, D and E
 T(:,:)  = 0.0_DP
 DO NZ = 2,N-1 
  T(NZ,NZ-1) =  1.0_DP
  T(NZ,NZ)   = -2.0_DP
  T(NZ,NZ+1) =  1.0_DP
  E(NZ,NZ)   = NZ*DZ
  I(NZ,NZ)   =  1.0_DP
  D(NZ,NZ)   = DK2(NZ)
 END DO
 !I(1,1)   =  1.0_DP
 !E(1,1)   =  DZ
 !D(1,1)   =  DK2(1)
 T(1,1)   = -2.0_DP+SIGMA1
 T(1,2)   =  1.0_DP+SIGMA2
 T(N,N-1) =  1.0_DP+TAU2
 T(N,N)   = -2.0_DP+TAU1
 !I(N,N)   =  1.0_DP
 !E(N,N)   =  N*DZ
 !D(N,N)   =  DK2(N)
 E        = DIAG(CMPLX(ALT,0.0_DP,KIND=DPC),N)
 D        = (IM/(2.0_DP*K0))*DIAG(DK2,N)
 I        = EYE(N)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !System - M1 and M2
 ALLOCATE(M1(N,N),M2(N,N),MB2(NDIAG,N),MB1(NDIAG,N))
 C1 = CMPLX(1.0_DP,0.0_DP,KIND=DPC)
 C2  = CMPLX(0.0_DP,0.0_DP,KIND=DPC)
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
 M1  = I+A*BETA1*T+BETA1*D
 M2  = I+A*BETA2*T+BETA2*D 
 MB1 = D2B(M1,N,KL,KU)
 MB2 = D2B(M2,N,KL,KU)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Forward-marching procedure
 ALLOCATE(MB1T(NDIAG,N),MB2T(KL+NDIAG,N),IPIV(N),TEMP(N,1))
 CALL CPU_TIME(TI)
 DO MX = 2,M+1
  WRITE(*,*) "Step :", MX-1, "out of", M
  !Matrix updating
  MB1T                  = MB1
  MB1T(KU+1,:)          = MB1(KU+1,:)-BETA1*IM*K0*TERR2(MX-1)*ALT
  MB2T(KL+1:KL+NDIAG,:) = MB2
  MB2T(KL+KU+1,:)       = MB2(KU+1,:)-BETA2*IM*K0*TERR2(MX-1)*ALT
  !Right hand side multiplication
  !DD = D-IM*K0*TERR2(MX-1)*E
  !M1 = I+(DR/2.0_DP)*(A*T+DD)
  !M2 = I-(DR/2.0_DP)*(A*T+DD)
  !MB2T(KL+1:KL+NDIAG,:) = D2B(M2,N,KL,KU)
  !MB1T                  = D2B(M1,N,KL,KU)
  CALL ZGBMV('N',N,N,KL,KU,C1,MB1T,KL+KU+1,PHI(:,MX-1),1,C2,TEMP,1)   !Band matrix form
  WRITE(*,*) ".... Right-hand side multiplied"
  !System solving
  CALL ZGBSV(N,KL,KU,1,MB2T,2*KL+KU+1,IPIV,TEMP,N,INFO)                           !Band matrix form
  IF (INFO.NE.0) THEN
   STOP "ERROR : Matrix is singular"
  END IF
  WRITE(*,*) ".... System solved"
  !Real pressure
  PHI(1:N,MX)  = TEMP(:,1)
  PC(1:N,MX-1) = EXP(IM*K0*(MX-1)*DR)*PHI(1:N,MX)/SQRT((MX-1)*DR)
  P(1:N,MX-1)  = REALPART(PC(1:N,MX-1))
  WRITE(*,*)
 END DO
 CALL CPU_TIME(TF)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !A posteriori phase correction
 DO MX = 1,M-1
  QC(1:N,MX) = PHI(1:N,MX)+MX*A*(PHI(1:N,MX+1)-2.0_DP*PHI(1:N,MX)+PHI(1:N,MX-1))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 ALLOCATE(LP(N,M),LQ(N,M),LPG(M),LQG(M))
 P0 = ABS(PC(NS,1))                  !Reference pressure (real)
 P0 = ABS(QC(NS,1))                  !Reference corrected pressure (real)
 DO MX = 1,M
  DO NZ = 1,N
   LP(NZ,MX) = 10.0_DP*LOG10(ABS(PC(NZ,MX))**2/(2*P0**2))
   LQ(NZ,MX) = 10.0_DP*LOG10(ABS(QC(NZ,MX))**2/(2*Q0**2)) 
   IF (LP(NZ,MX).LE.-120) THEN
    LP(NZ,MX) = -120
    LQ(NZ,MX) = -120
   END IF
  END DO
 END DO
 LPG = LP(1,1:M)
 LQG = LQ(1,1:M)
 LPR = LP(1,NR)
 LQR = LQ(1,NR)
 PR  = SQRT(2.0_DP)*P0*10**(LPR/20)
 QR  = SQRT(2.0_DP)*Q0*10**(LPR/20)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Output 
 OPEN(UNIT=10,FILE="Gauss_LPg.dat")
 OPEN(UNIT=20,FILE="Gauss_LP.dat")
 OPEN(UNIT=30,FILE="Gauss_PC.dat")
 OPEN(UNIT=40,FILE="Gauss_P.dat")
 OPEN(UNIT=50,FILE="Gauss.dat")
 DO MX = 1,M
  WRITE(10,100) MX*DR, LPG(MX)
  WRITE(50,100) MX*DR, TERR(MX), TERR1(MX), TERR2(MX)
  DO NZ = 1,N
   WRITE(20,100) MX*DR, NZ*DZ+TERR(MX), LP(NZ,MX), LQ(NZ,MX)
   WRITE(30,101) MX*DR, NZ*DZ+TERR(MX), PC(NZ,MX), QC(NZ,MX)
   WRITE(40,100) MX*DR, NZ*DZ+TERR(MX), P(NZ,MX)
  END DO
  WRITE(20,*)
  WRITE(30,*)
  WRITE(40,*)
 END DO
 100 FORMAT(4(3X,F12.6))
 101 FORMAT(2(3X,F12.6),3X,F12.6,SP,F12.6,SS,"i")
 !Plotting of results
 CALL SYSTEM('gnuplot -p plot.plt')
 !CALL SYSTEM('gnuplot -p plot2.plt')
 WRITE(*,*) "Main CPU time (s) :", TF-TI
 !PRINT *, "SPL at the source LP0 (dB) :", 20.0_DP*LOG10(ABS(P0))
 WRITE(*,*) "Real pressure at the source - P0 (Pa) :", P0
 WRITE(*,*) "Real pressure at the target - PR (Pa) :", PR
 WRITE(*,*) "SPL/TL at the target - LPR (dB) :", LPR
 WRITE(*,*) 
 !-----------------------------------------------------------
 !----------------------------------------------------------
END PROGRAM CNPE2_HT

