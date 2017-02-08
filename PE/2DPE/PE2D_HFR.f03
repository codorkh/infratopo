PROGRAM PE2D_HFR
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
 WRITE(*,*) "Dimension : 2D"
 WRITE(*,*) "Atmosphere : Homogeneous"
 WRITE(*,*) "Boundary : Flat"
 WRITE(*,*) "Impedance : Rigid"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Enter inputs L, f, zs, As"
 WRITE(*,*) ""
 READ *, L, F, ZS, ABL
 WRITE(*,*) ""
 !-----------------------------------------------------------
 !-----------------------------------------------------------
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
 MB2(NDIAG,N), MB1(NDIAG,N), WORK(N), RWORK(2*N), CWORK(2*N), TEMP(N,1), &
 MTEMP(2*KL+KU+1,N), IPIV(N), TEMP(N))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 DO NZ = 1,N
  IF (NZ.GE.N-NABL) THEN
   K(NZ) = K0+ABL*IM*(NZ-N+NABL)**2/NABL**2
  ELSE
   K(NZ) = K0
  ENDIF
  DK2(NZ) = IM*(K(NZ)**2-K0**2)/(2*K0)
  ALT(NZ) = NZ*DZ
 END DO
 TAU1 = 4./(3.-2.*IM*K0*DZ)
 TAU2 = -1./(3.-2.*IM*K0*DZ)
 SIGMA1 = 4./3.
 SIGMA2 = -1./3.
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 A0 = 1.3717
 A2 = -0.3701
 B = 3.
 DO NZ = 1,N
  PHI(NZ,1) = SQRT(IM*K0)* &
              ((A0+A2*K0**2*(NZ*DZ-ZS)**2)*EXP(-K0**2*(NZ*DZ-ZS)**2/B) + &
              (A0+A2*K0**2*(NZ*DZ+ZS)**2)*EXP(-K0**2*(NZ*DZ+ZS)**2/B))
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Initialization (Dense Matrices)
 !FDW = [1 -2 1]  !Finite Difference Weights
 !---------
 T(:,:) = 0._DP
 DO NZ = 2,N-1 
  T(NZ,NZ-1) = T(NZ,NZ-1)+1.
  T(NZ,NZ) = T(NZ,NZ)-2.
  T(NZ,NZ+1) = T(NZ,NZ+1)+1.
 END DO
 T(1,1) = T(1,1)-2.+SIGMA1
 T(1,2) = T(1,2)+1.+SIGMA2
 T(N,N-1) = T(N,N-1)+1.+TAU2
 T(N,N) = T(N,N)-2.+TAU1
 E = DIAG(ALT)
 D = DIAG(DK2)
 I = EYE(N)
 !-----------------------------------------------------------
 !System - M1 and M2
 !---------
 BETA1 = 1./(2.*IM*K0)+DR/2.
 BETA2 = 1./(2.*IM*K0)-DR/2.
 M1 = I+ALPHA*BETA1*T+BETA1*D
 M2 = I+ALPHA*BETA2*T+BETA2*D 
 !CALL DENSE2BAND(M1,N,KL,KU,MB1)
 MB2 = DENSE2BAND(M2,N,KL,KU)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Forward-marching procedure
 CALL CPU_TIME(TI)
 !-----------------------------------------------------------
 DO MX = 2,M+1
  WRITE(*, *) "Step :", MX-1, "out of", M 
  TEMP(:,1) = MATMUL(M1,PHI(1:N,MX-1))
  !CALL ZGBMV('N',N,N,KL,KU,1.,MB1,KL+KU+1,PHI(1:N,MX-1),1,0.,TEMP(:,1),1)  
  WRITE(*, *) ".... Right-hand side multiplied"
  !MTEMP = M2
  MTEMP(KL+1:2*KL+KU+1,:) = MB2  
  !CALL ZGESV(N,1,MTEMP,N,IPIV,TEMP,N,INFO)
  CALL ZGBSV(N,KL,KU,1,MTEMP,2*KL+KU+1,IPIV,TEMP(:,1),N,INFO)
  !CALL ZGTSV(N,1,MTEMP(2*KL+KU+1,1:N-1),MTEMP(KL+KU+1,1:N),MTEMP(KL+1,2:N),TEMP(:,1),N,INFO)
  IF (INFO.NE.0) THEN
   WRITE(*,*) INFO
   STOP "ERROR : Matrix is singular"
  END IF
  WRITE(*, *) ".... System solved"
  PHI(1:N,MX) = TEMP(:,1)
  P(1:N,MX-1) = EXP(IM*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
 END DO
! !-----------------------------------------------------------
! IF (SP == 0) THEN
! !-----------------------------------------------------------
!  MAT = MATMUL(INV(M2),M1)
!  DO MX = 2,M+1
!   WRITE(*, *) "Step :", MX-1, "out of", M
!   PHI(1:N,MX) = MATMUL(MAT,PHI(1:N,mx-1))
!   P(1:N,MX-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
!  END DO
! !-----------------------------------------------------------
! ELSEIF (SP == 1) THEN
! !-----------------------------------------------------------
!   MAT_SP = SPCONV(MAT) ! SLOW PROCEDURE !
!   NELT = MAT_SP%NELT
!   ALLOCATE(IMAT(NELT),JMAT(NELT))
!   IMAT = MAT_SP%ROW
!   JMAT = MAT_SP%COL
!   VALMAT = MAT_SP%VAL
!   DO MX = 2,M+1
!    WRITE(*, *) "Step :", mx-1, "out of", M
!    CALL DSMV(N,PHI(1:N,MX-1),PHI(1:N,MX),NELT,IMAT,JMAT,VALMAT,0)
!    P(1:N,mx-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
!   END DO
! !-----------------------------------------------------------
! ELSE
! !-----------------------------------------------------------
!  ALLOCATE(SA1(N**3),SA2(N**3),IJA(N**3))
!  THRESH = 0
!  TOL = 1.0e-4_dp
!  CALL SPRSIN_DP(M1,THRESH,SP1)
!  CALL SPRSIN_DP(M2,THRESH,SP2)
!  DO MX = 2,M+1
!   WRITE(*, *) "Step :", mx-1, "out of", M
!   CALL SPRSAX_DP(SP1,PHI(1:N,MX-1),TEMP)
!   CALL LINBCG(SP2,TEMP,PHI(1:N,MX),1,TOL,10,ITER,ERR)
!   P(1:N,mx-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
!  END DO
! !-----------------------------------------------------------
! END IF
! !-----------------------------------------------------------
 CALL CPU_TIME(TF)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 P0 = P(NS,1)
 DO MX = 1,M
  DO NZ = 1,N
   LP(NZ,MX) = 20.*LOG10(ABS(P(NZ,MX)/P0))
  END DO
 END DO
 LPG = LP(1,1:M)
 LP2 = LP/SQRT(2.)
 LPG2 = LP2(1,1:M)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Output
 OPEN(UNIT=10,FILE="results/PE2D_HFR_LPg.dat")
 OPEN(UNIT=20,FILE="results/PE2D_HFR_LP.dat")
 OPEN(UNIT=30,FILE="results/PE2D_HFR_P.dat")
 DO NZ = 1,N
  WRITE(20, 102) LP(NZ,1:M)
  WRITE(30, 101) P(NZ,1:M)
 END DO
 DO MX = 1,M
  WRITE(10, 100) MX*Dr, LPG(MX), LPG2(MX)
 END DO
 100 FORMAT(3(3X,F12.3))
 101 FORMAT(3X,F12.3,SP,F12.3,SS,"i")
 102 FORMAT(100000(3X,F12.3)) 
 PRINT *, "Main CPU time (s) :", TF-TI
 PRINT *, "Source pressure P0 (dB) :", 20*LOG10(ABS(P0))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Checking Matrix conditioning
 NORM = ZLANGE('I',N,N,M1,N,WORK)
 WRITE(*,*) "||M1|| = ", NORM 
 CALL ZGECON('I',N,M1,N,NORM,RCOND,CWORK,RWORK,INFO)
 WRITE(*,*) "K(M1) = ", RCOND, "COND(M1) = ", 1./RCOND
 WRITE(*,*)
 NORM = ZLANGE('I',N,N,M2,N,WORK)
 WRITE(*,*) "||M2|| = ", NORM 
 CALL ZGECON('I',N,M2,N,NORM,RCOND,CWORK,RWORK,INFO)
 WRITE(*,*) "K(M2) = ", RCOND, "COND(M2) = ", 1./RCOND
 WRITE(*,*)
 NORM = ZLANGB('I',N,KL,KU,MB2,KL+KU+1,WORK)
 WRITE(*,*) "||M2|| = ", NORM 
 CALL ZGBCON('I',N,KL,KU,MTEMP,2*KL+KU+1,IPIV,NORM,RCOND,CWORK,RWORK,INFO)
 WRITE(*,*) "K(M2) = ", RCOND, "COND(M2) = ", 1./RCOND
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Plotting of results
 CALL SYSTEM('gnuplot -p plot_flat.plt')
 !----------------------------------------------------------
END PROGRAM PE2D_HFR
!-----------------------------------------------------------
!-------------      EXTERNAL PROCEDURES      ---------------
!-----------------------------------------------------------
