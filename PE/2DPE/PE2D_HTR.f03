PROGRAM PE2D_HTR
 USE NRTYPE
 USE MATLIB
 USE PARAM
 USE GROUND
 IMPLICIT NONE
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) "========== Parabolic Equation Package =========="
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "University of Bristol"
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
 !-----------------------------------------------------------
 !----------------------------------------------------------
 C0 = 343.0_DP
 LMBDA0 = C0/F
 K0 = 2*PI*F/C0
 DZ = LMBDA0/10
 DR = LMBDA0/10
 H = 100*LMBDA0
 HABL = 30*LMBDA0
 N = FLOOR(H/DZ)
 M = FLOOR(L/DR)
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Size of the system : ", N
 WRITE(*,*) "Number of steps : ", M
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 NS = MAX(1,FLOOR(ZS/DZ))
 NABL = FLOOR(HABL/DZ)
 ALPHA = IM/(2.0_DP*K0*DZ**2)
 KU = 1
 KL = 1
 NDIAG = KU+KL+1
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Memory allocation
 ALLOCATE(K(N), DK2(N), ALT(N), PHI(N,M+1), T(N,N), D(N,N), M1(N,N), &
 M2(N,N), I(N,N), P(N,M), LP(N,M), LP2(N,M), LPG(M), LPG2(M), E(N,N), &
 MB2(NDIAG,N), MB1(NDIAG,N), TEMP(N,1), M1T(N,N), M2T(N,N), IPIV(N), &
 TERR(1:M+1), TERR1(1:M+1), TERR2(1:M+1), LIN(0:M+1), GHX(3), COND(M,2), &
 MB1T(NDIAG,N), MB2T(KL+NDIAG,N))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Atmosphere density and wave velocity
 DO NZ = 1,N
  IF (NZ >= N-NABL) THEN
   K(NZ) = K0+ABL*IM*(NZ-N+NABL)**2/NABL**2
  ELSE
   K(NZ) = K0
  ENDIF
  DK2(NZ) = IM*(K(NZ)**2-K0**2)/(2*K0)
  ALT(NZ) = NZ*DZ
 END DO
 TAU1 = 4.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 TAU2 = -1.0_DP/(3.0_DP-2.0_DP*IM*K0*DZ)
 SIGMA1 = 4.0_DP/3.0_DP
 SIGMA2 = -1.0_DP/3.0_DP
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Terrain Boundary condition
 LIN(0) = 0
 THETA = ANGLE*PI/180.0_DP
 X0 = L/2.0_DP
 S = L/5.0_DP
 H0 = SQRT(EXP(1.0_DP)/2.0_DP)*S*ATAN(THETA)
 DO MX = 1,M
  GHX = GHILL(H0,X0,S,MX*DR)
  TERR(MX) = GHX(1)
  TERR1(MX) = GHX(2)
  TERR2(MX) = GHX(3)
  LIN(MX) = LIN(MX-1)+(1+TERR1(MX))*DR
 END DO
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
 E = DIAG(ALT,N)
 D = DIAG(DK2,N)
 I = EYE(N)
 !-----------------------------------------------------------
 !System - M1 and M2
 !---------
 BETA1 = 1.0_DP/(2.0_DP*IM*K0)+DR/2.0_DP
 BETA2 = 1.0_DP/(2.0_DP*IM*K0)-DR/2.0_DP
 M1 = I+ALPHA*BETA1*T+BETA1*D
 M2 = I+ALPHA*BETA2*T+BETA2*D 
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
  M1T = M1-BETA1*IM*K0*TERR2(MX-1)*E
  !M2T = M2-BETA2*IM*K0*TERR2(MX-1)*E
  MB1T = MB1
  MB1T(KU+1,:) = MB1(KU+1,:)-BETA1*IM*K0*TERR2(MX-1)*ALT
  MB2T(KL+1:KL+NDIAG,:) = MB2
  MB2T(KL+KU+1,:) = MB2(KU+1,:)-BETA2*IM*K0*TERR2(MX-1)*ALT
  !-------------------------------
  !Right hand side multiplication
  !-------------------------------
  !TERR2(MX-1) = (TERR1(MX)-TERR1(MX-1))/DR  			!Varying coefficient to be injected...  
  !TEMP(:,1) = MATMUL(M1T,PHI(1:N,MX-1))
  !zgbmv (TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  CALL ZGBMV('N',N,N,KL,KU,CMPLX(1.0_DP,0.0_DP,KIND=DP),MB1T,KL+KU+1,PHI(1:N,MX-1),1,CMPLX(0.0_DP,0.0_DP,KIND=DP),TEMP(:,1),1)
    !---------
    !---------
    !ALLOCATE(WORK(N),RWORK(2*N),CWORK(2*N))
    !NORM = ZLANGE('I',N,N,M1T,N,WORK)
    WRITE(*,*) "||M1|| = ", NORM 
    CALL ZGECON('I',N,M1T,N,NORM,RCOND,CWORK,RWORK,INFO)
    WRITE(*,*) "K(M1) = ", 1.0_DP/RCOND
    COND(MX-1,1) = 1.0_DP/RCOND
    DEALLOCATE(WORK,RWORK,CWORK)
    !---------
    !---------
  WRITE(*,*) ".... Right-hand side multiplied"
  WRITE(*,*)
  !-------------------------------
  !System solving
  !-------------------------------
  !MTEMP2 = M2
  !CALL ZGESV(N,1,M2T,N,IPIV,TEMP,N,INFO)                                                         !General matrix form     
  CALL ZGBSV(N,KL,KU,1,MB2T,2*KL+KU+1,IPIV,TEMP,N,INFO)                                           !Band matrix form
  !CALL ZGTSV(N,1,MT2T(3,1:N-1),MT2T(2,1:N),MT2T(1,2:N),TEMP(:,1),N,INFO)                         !Tridiagonal matrix form
  IF (INFO.NE.0) THEN
   WRITE(*, *) INFO
   STOP "ERROR : Matrix is singular"
  END IF
    !---------
    !---------
    !ALLOCATE(WORK(N),RWORK(N),CWORK(2*N))
    !NORM = ZLANGB('I',N,KL,KU,MB2T(KL+1:2*KL+KU+1,:),KL+KU+1,WORK)
    !WRITE(*,*) "||M2|| = ", NORM 
    !CALL ZGBCON('I',N,KL,KU,MB2T,2*KL+KU+1,IPIV,NORM,RCOND,CWORK,RWORK,INFO)
    !WRITE(*,*) "K(M2) = ", 1.0_DP/RCOND 
    !COND(MX-1,2) = 1.0_DP/RCOND
    !DEALLOCATE(WORK,RWORK,CWORK)
    !---------
    !---------
  WRITE(*,*) ".... System solved"
  WRITE(*,*)
  PHI(1:N,MX) = TEMP(:,1)
  P(1:N,MX-1) = EXP(IM*K0*(MX-1)*DR)*PHI(1:N,MX)*(1/SQRT((MX-1)*DR))
 END DO
 !-----------------------------------------------------------
 !ELSEIF (SP == 1) THEN
 !Sparse SLATEC
  !DO MX = 2,M+1
  ! WRITE(*, *) "Step :", MX-1, "out of", M
  ! TERR2(MX-1) = (TERR1(MX)-TERR1(MX-1))/Dr
  ! D = D-im*K0*TERR2(MX-1)*DIAG(ALT)
  ! WRITE(*, *) "1"
  ! M1 = I+(Dr/2)*(aa*T+D)+(1/(2*im*K0))*(aa*T+D)
  ! M2 = I-(Dr/2)*(aa*T+D)+(1/(2*im*K0))*(aa*T+D)
  ! !MAT = MATMUL(INV(M2),M1)
  ! WRITE(*, *) "2"
  ! !MAT_SP = SPCONV(MAT) ! SLOW PROCEDURE !
  ! CALL SPRSIN_DP(M2,0,SPM2)
  ! CALL SPRSIN_DP(M1,0,SPM1)
  ! NELT = MAT_SP%NELT
  ! ALLOCATE(IMAT(NELT),JMAT(NELT))
  ! IMAT = MAT_SP%ROW
  ! JMAT = MAT_SP%COL
  ! VALMAT = MAT_SP%VAL
  ! CALL DSMV(N,PHI(1:N,MX-1),PHI(1:N,MX),NELT,IMAT,JMAT,VALMAT,0)
  ! P(1:N,mx-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
  !END DO
 !-----------------------------------------------------------
 CALL CPU_TIME(TF)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 P0 = P(NS,1)
 DO MX = 1,M
  DO NZ = 1,N
   IF (ABS(P(NZ,MX)).LE.10**(-6)) THEN
    P(NZ,MX) = 10**(-6)
   END IF
   LP(NZ,MX) = 20.*LOG10(ABS(P(NZ,MX)/P0))
  END DO
 END DO
 LPG = LP(1,1:M)
 LP2 = LP/SQRT(2.0_DP)
 LPG2 = LP2(1,1:M)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 WRITE(*,*)
 WRITE(*,*) "================================================"
 WRITE(*,*)
 WRITE(*,*) " T = "
 DO NZ = 1,4
  WRITE(*,101) T(NZ,1:4)
 END DO
 WRITE(*,*)
 WRITE(*,*) "================================================"
 WRITE(*,*)
 WRITE(*,*) " E = "
 DO NZ = 1,4
  WRITE(*,101) E(NZ,1:4)
 END DO
 WRITE(*,*)
 WRITE(*,*) "================================================"
 WRITE(*,*)
 WRITE(*,*) " D = "
 DO NZ = N-3,N
  WRITE(*,101) D(NZ,N-3:N)
 END DO
 WRITE(*,*)
 WRITE(*,*) "================================================"
 WRITE(*,*)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Output
 OPEN(UNIT=10,FILE="results/Test_LPg.dat")
 OPEN(UNIT=20,FILE="results/Test_LP.dat")
 OPEN(UNIT=30,FILE="results/Test_P.dat")
 !OPEN(UNIT=40,FILE="results/Test_Cond.dat")
 OPEN(UNIT=50,FILE="results/Test_Terr.dat")
 DO MX = 1,M
  WRITE(10,100) MX*DR, LPG(MX), LPG2(MX)
  !WRITE(40,100) COND(MX,1), COND(MX,2)
  WRITE(50,100) MX*DR, TERR(MX), TERR1(MX), TERR2(MX)
  DO NZ = 1,M
   WRITE(20,100) MX*DR, NZ*DZ, LP(NZ,MX)
   WRITE(30,101) MX*DR, NZ*DZ, P(NZ,MX)
  END DO 
 END DO
 100 FORMAT(3(3X,F12.3))
 101 FORMAT(2(3X,F12.3),3X,F12.3,SP,F12.3,SS,"i")
 PRINT *, "Main CPU time (s) :", TF-TI
 PRINT *, "Source pressure P0 (dB) :", 20*LOG10(ABS(P0))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Plotting of results
 CALL SYSTEM('gnuplot -p plot_topo.plt')
 !----------------------------------------------------------
END PROGRAM PE2D_HTR
!-----------------------------------------------------------
!-------------      EXTERNAL PROCEDURES      ---------------
!-----------------------------------------------------------

