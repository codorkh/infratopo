PROGRAM CNPE2_HT_SNT
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
 WRITE(*,*) "Ascension Island PE Simulation program"
 WRITE(*,*) "Choose acoustic source - S1 or S2"
 WRITE(*,*) ""
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Read data
 OPEN(UNIT=10,FILE="CNPE2_IN.dat")
 READ(10,*) L          !Domain range
 READ(10,*) F          !Frequency
 READ(10,*) ZS         !Source height
 READ(10,*) ZR         !Target height
 READ(10,*) ZX         !Impedance real part
 READ(10,*) ZY         !Impedance imaginary part
 READ(10,*) ANGLE      !Terrain maximum angle
 READ(10,*) DC         !Atmosphere profile coefficient
 READ(10,*) SRCE_I     !Source number
 READ(10,*) TOPO_I     !Ascension Island topographic profile
 READ(10,*) FREQ_I     !Overtone
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Ascension Island data driver variables
 WRITE(TRACK,'(I0)') TOPO_I
 WRITE(FREQ ,'(I0)') FREQ_I
 WRITE(SRCE ,'(I0)') SRCE_I
 FREQ_NAT  = RESHAPE((/ 1.90_DP, 3.79_DP, 5.69_DP, 7.59_DP, &
                        1.42_DP, 2.84_DP, 4.27_DP, 5.69_DP /), (/2,4/))
 PATH_TOPO = "../../GROUND/SnT_2017/L" // TRACK // "/"
 FILE_TOPO = "asc_L" // TRACK // "_2d_f4_sampled.dat"
 TOPO      = TRIM(PATH_TOPO) // TRIM(FILE_TOPO)
 CALL READ_ARRAY(TRIM(TOPO),GMT_DATA,NLINES,NCOLS)
 ALLOCATE(X(NLINES),TERR(NLINES),TERR1(NLINES),TERR2(NLINES))
 L         = GMT_DATA(NLINES,1)
 X         = GMT_DATA(:,1)
 TERR      = GMT_DATA(:,2)
 TERR1     = GMT_DATA(:,3)
 TERR2     = GMT_DATA(:,4)
 ZS        = 25.000_DP
 ZR        = 0.000_DP
 F	       = FREQ_NAT(SRCE_I,FREQ_I)
 AS        = 0.070_DP
 !-----------------------------------------------------------
 !-----------------------------------------------------------
!Discretization of the Atmosphere domain
 C0     = 343.0_DP
 LMBDA0 = C0/F
 K0     = 2.0_DP*PI*F/C0
 DZ     = LMBDA0/10
 DR     = GMT_DATA(3,1)-GMT_DATA(2,1)
 H      = 100*LMBDA0
 HA     = 30*LMBDA0
 N      = FLOOR(H/DZ)
 M      = NLINES
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
 WIDE_I = 0
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
 DC = 2.5_DP
 Z0 = 0.5_DP
 DO NZ = 1,N
  ALT(NZ) = NZ*DZ
  C(NZ)   = C0+DC*LOG(1+ALT(NZ)/Z0)
  K(NZ)   = 2.0_DP*PI*F/C(NZ)
  IF (NZ.GT.(N-NA)) THEN
   K(NZ)  = K(NZ)+AS*IM*(NZ-N+NA)**2/NA**2
  END IF
  DK2(NZ) = K(NZ)**2-K0**2
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Terrain Boundary condition
 ALLOCATE(LIN(0:M),SLOPX(M))
 DO MX = 0,M
  LIN(MX)   = LIN(MX-1)+(1+TERR1(MX))*DR
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
               (A0+A2*K0**2*(NZ*DZ+ZS)**2)*EXP(-K0**2*(NZ*DZ+ZS)**2/B))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Boundary condition - Finite Difference
 ALLOCATE(T(N,N),D(N,N),E(N,N),I(N,N))
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
 SELECT CASE (WIDE_I)
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
 !A posteriori phase correction
 DO MX = 1,M
  QC(1:N,MX) = (PHI(1:N,MX)+MX*A*(PHI(1:N,MX+1)-2.0_DP*PHI(1:N,MX)+PHI(1:N,MX-1)))/SQRT(MX*DR)
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 ALLOCATE(LP(N,M),LQ(N,M),LPG(M),LQG(M))
 P0 = ABS(PC(NS,1))                   !Reference pressure (real)
 Q0 = ABS(QC(NS,1))                   !Reference corrected pressure (real)
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
 OUT_EXT  = "L" // TRACK // "_f" // FREQ
 OUT_LP   = "LP_"   // TRIM(OUT_EXT) // ".dat"
 OUT_LPg  = "LPg_"  // TRIM(OUT_EXT) // ".dat"
 OUT_P    = "P_"    // TRIM(OUT_EXT) // ".dat"
 OUT_PC   = "PC_"   // TRIM(OUT_EXT) // ".dat"
 OUT_TOPO = "TOPO_" // TRIM(OUT_EXT) // ".dat"
 OPEN(UNIT=10,FILE=TRIM(OUT_LPg))
 OPEN(UNIT=20,FILE=TRIM(OUT_LP))
 OPEN(UNIT=30,FILE=TRIM(OUT_PC))
 OPEN(UNIT=40,FILE=TRIM(OUT_P))
 OPEN(UNIT=50,FILE=TRIM(OUT_TOPO))
 DO MX = 1,M
  WRITE(10,100) X(MX), LPG(MX), LQG(MX)
  WRITE(50,100) X(MX), TERR(MX), TERR1(MX), TERR2(MX), SLOPX(MX)
  DO NZ = 1,N
   WRITE(20,100) X(MX), NZ*DZ+TERR(MX), LP(NZ,MX), LQ(NZ,MX)
   WRITE(30,101) MX*DR, NZ*DZ+TERR(MX), PC(NZ,MX), QC(NZ,MX)
   WRITE(40,100) MX*DR, NZ*DZ+TERR(MX), P(NZ,MX)
  END DO
  WRITE(20,*)
  WRITE(30,*)
  WRITE(40,*)
 END DO
 100 FORMAT(5(3X,F12.3))
 101 FORMAT(2(3X,F12.3),3X,F12.3,SP,F12.3,SS,"i")
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 DEALLOCATE(PHI,PC,QC,P,Q,LP,LQ)
 WRITE(*,*) "================================================"
 WRITE(*,*) "Input file  : ", TOPO
 WRITE(*,*) "================================================"
 WRITE(*,*) "Output file : ", TRIM(OUT_LPg)
 WRITE(*,*) "              ", TRIM(OUT_LP)
 WRITE(*,*) "              ", TRIM(OUT_PC)
 WRITE(*,*) "              ", TRIM(OUT_P)
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
END PROGRAM CNPE2_HT_SNT
