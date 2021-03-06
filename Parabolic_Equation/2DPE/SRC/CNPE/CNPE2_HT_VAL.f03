PROGRAM CNPE2_HT_VAL
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
 WRITE(*,*) "PE Validation program"
 WRITE(*,*) "Choose source frequency - 1 (10 Hz), 2 (100 Hz), 3 (1000 Hz)"
 WRITE(*,*) ""
 READ *,    FREQ_I
 WRITE(*,*) "Atmospheric refraction - 1 (Isothermal), 2 (+0.1), 3 (-0.1)"
 WRITE(*,*) ""
 READ*,     ATMO_I
 WRITE(*,*) ""
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Ascension Island data driver variables
 WRITE(FREQ,'(I0)') FREQ_I
 WRITE(ATMO,'(I0)') ATMO_I
 FREQ_VAL = (/ 10.0_DP, 100.0_DP, 1000.0_DP /)
 Z_VAL    = (/ CMPLX(38.79_DP,38.41_DP,KIND=DPC), &
               CMPLX(12.81_DP,11.62_DP,KIND=DPC), &
               CMPLX( 5.96_DP, 2.46_DP,KIND=DPC) /)
 ATMO_REF = (/  0.0_DP,   0.1_DP,   -0.1_DP /)
 L        = 10000_DP
 ZS       = 5.000_DP
 ZR       = 1.000_DP
 F	      = FREQ_VAL(FREQ_I)
 Z        = Z_VAL(FREQ_I)
 R        = (Z-1.0_DP)/(Z+1.0_DP)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Discretization of the Atmosphere domain
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
 NR     = MAX(1,FLOOR(ZR/DZ))
 NA     = FLOOR(HA/DZ)
 A      = IM/(2.0_DP*K0*DZ**2)
 KU     = 1
 KL     = 1
 NDIAG  = KU+KL+1
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Absorbing coefficient
 AS = 4.897_DP*(10**(-9))*(F**3) - &
      7.119_DP*(10**(-6))*(F**2) + &
      3.109_DP*(10**(-3))*F + &
      0.113
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Atmosphere density and wave velocity
 ALLOCATE(C(N),K(N),DK2(N),ALT(N))
 DC = ATMO_REF(ATMO_I)
 DO NZ = 1,N
  ALT(NZ) = NZ*DZ
  C(NZ)   = C0+DC*ALT(NZ)
  K(NZ)   = 2.0_DP*PI*F/C(NZ)
  IF (NZ.GT.(N-NA)) THEN
   K(NZ)  = K(NZ)+AS*IM*(NZ-N+NA)**2/NA**2
  END IF
  DK2(NZ) = K(NZ)**2-K0**2
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Terrain Boundary condition
 ALLOCATE(X(M),TERR(M),TERR1(M),TERR2(M),SLOPX(M))
 DO MX = 0,M
  X(MX)     = MX*DR
  TERR(MX)  = 0.000_DP
  TERR1(MX) = 0.000_DP
  TERR2(MX) = 0.000_DP
  SLOPX(MX) = ATAN(TERR1(MX))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Starting field
 ALLOCATE(PHI(N,M+1),PC(N,M),QC(N,M),P(N,M),Q(N,M))
 A0 = 1.3717_DP
 A2 = -0.3701_DP
 B  = 3.0_DP
 DO NZ = 1,N
  PHI(NZ,1) = SQRT(IM*K0)* &
              ((A0+A2*K0**2*(NZ*DZ-ZS)**2)*EXP(-K0**2*(NZ*DZ-ZS)**2/B) + &
               R*(A0+A2*K0**2*(NZ*DZ+ZS)**2)*EXP(-K0**2*(NZ*DZ+ZS)**2/B))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Boundary condition - Finite Difference
 ALLOCATE(T(N,N),D(N,N),E(N,N),I(N,N))
 TAU1    =  4.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 TAU2    = -1.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 SIGMA1  =  4.0_DP*Z/(3.0_DP*Z-2.0_DP*IM*K0*DZ)
 SIGMA2  = -1.0_DP*Z/(3.0_DP*Z-2.0_DP*IM*K0*DZ)
 !Matrix Initialization (Dense Matrices) - T, D and E
 T(:,:)  = 0.0_DP
 DO NZ = 2,N-1
  T(NZ,NZ-1) =  1.0_DP
  T(NZ,NZ)   = -2.0_DP
  T(NZ,NZ+1) =  1.0_DP
 END DO
 T(1,1)   = -2.0_DP+SIGMA1
 T(1,2)   =  1.0_DP+SIGMA2
 T(N,N-1) =  1.0_DP+TAU2
 T(N,N)   = -2.0_DP+TAU1
 E        = DDIAG(CMPLX(ALT,0.0_DP,KIND=DPC),N)
 D        = (IM/(2.0_DP*K0))*DDIAG(DK2,N)
 I        = EYE(N)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !System - M1 and M2
 ALLOCATE(M1(N,N),M2(N,N),MB2(NDIAG,N),MB1(NDIAG,N))
 C1 = CMPLX(1.0_DP,0.0_DP,KIND=DP)
 C2 = CMPLX(0.0_DP,0.0_DP,KIND=DP)
 BETA1 = 1.0_DP/(2.0_DP*IM*K0)+DR/2.0_DP
 BETA2 = 1.0_DP/(2.0_DP*IM*K0)-DR/2.0_DP
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
  CALL ZGBMV('N',N,N,KL,KU,C1,MB1T,KL+KU+1,PHI(1:N,MX-1),1,C2,TEMP(:,1),1)   !Band matrix form
  WRITE(*,*) ".... Right-hand side multiplied"
  !System solving
  CALL ZGBSV(N,KL,KU,1,MB2T,2*KL+KU+1,IPIV,TEMP,N,INFO)                      !Band matrix form
  IF (INFO.NE.0) THEN
   STOP "ERROR : Matrix is singular"
  END IF
  WRITE(*,*) ".... System solved"
  !Complex and real pressure
  PHI(1:N,MX)  = TEMP(:,1)
  PC(1:N,MX-1) = EXP(IM*K0*(MX-1)*DR)*PHI(1:N,MX)/SQRT((MX-1)*DR)
  P(1:N,MX-1)  = REALPART(PC(1:N,MX-1))
  WRITE(*,*)
 END DO
 CALL CPU_TIME(TF)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 ALLOCATE(LP(N,M),LQ(N,M),LPG(M),LQG(M))
 P0 = ABS(PC(NS,1))                   !Reference pressure (real)
 Q0 = ABS(QC(NS,1))                   !Reference corrected pressure (real)
 DO MX = 1,M
  DO NZ = 1,N
   LP(NZ,MX) = 10.0_DP*LOG10(ABS(PC(NZ,MX))**2/(P0**2))
   LQ(NZ,MX) = 10.0_DP*LOG10(ABS(QC(NZ,MX))**2/(Q0**2))
   IF (LP(NZ,MX).LE.-120) THEN
    LP(NZ,MX) = -120
    LQ(NZ,MX) = -120
   END IF
  END DO
 END DO
 LP0 = LP(NS,1)
 LQ0 = LQ(NS,1)
 LPG = LP(1,1:M)
 LQG = LQ(1,1:M)
 LPR = LP(NR,M)
 LQR = LQ(NR,M)
 PR  = SQRT(2.0_DP)*P0*10**(LPR/20)
 QR  = SQRT(2.0_DP)*Q0*10**(LPR/20)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Output
 OUT_EXT  = "W" // WIDE // "_A" // ATMO // "_f" // FREQ
 OUT_LPg  = "LPg_"  // TRIM(OUT_EXT) // ".dat"
 OPEN(UNIT=10,FILE=TRIM(OUT_LPg))
 DO MX = 1,M
  WRITE(10,100) X(MX), LPG(MX)
 END DO
 100 FORMAT(3(3X,F12.3))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 DEALLOCATE(PHI,PC,QC,P,Q,LP,LQ)
 WRITE(*,*) "================================================"
 WRITE(*,*) "Input file  : ", TOPO
 WRITE(*,*) "================================================"
 WRITE(*,*) "Output file : ", TRIM(OUT_LPg)
 WRITE(*,*) "              ", TRIM(OUT_TOPO)
 WRITE(*,*) "================================================"
 WRITE(*,*)
 WRITE(*,*) "System size (N x M) :", N, M
 WRITE(*,*) "Step size :", DR
 WRITE(*,*) "Main CPU time (s) :", TF-TI
 WRITE(*,*) "Real pressure at the source - P0 (Pa) :", P0
 WRITE(*,*) "Real pressure at the target - PR (Pa) :", PR
 WRITE(*,*) "SPL/TL at the target - LPR (dB) :", LPR
 WRITE(*,*)
 !-----------------------------------------------------------
 !----------------------------------------------------------
END PROGRAM CNPE2_HT_VAL
