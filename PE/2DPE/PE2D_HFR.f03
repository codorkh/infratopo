PROGRAM PE2D_HFR
	USE NRTYPE
	USE NRUTIL
	USE NRSP77
	USE PARAM
	USE MATLIB
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
	WRITE(*,*) "Dimension 	: 	2D"
	WRITE(*,*) "Atmosphere	:	Homogeneous"
	WRITE(*,*) "Boundary 	:	Flat"
	WRITE(*,*) "Impedance	:	Rigid"
	WRITE(*,*) ""
	WRITE(*,*) "================================================"
	WRITE(*,*) ""
	WRITE(*,*) "Enter inputs L, f, zs, As"
	WRITE(*,*) ""
	READ *, L, F, ZS, ABL
	WRITE(*,*) ""
	WRITE(*,*) "Matrix Storage Method"
	READ *, SP
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	lmbda0 = C0/F
	K0 = 2*PI*F/C0
	Dz = lmbda0/10
	Dr = lmbda0/10
	HABL = 30*lmbda0
	H = 100*lmbda0
	N = FLOOR(H/dz)
	M = FLOOR(L/dr)
	NS = MAX(1,FLOOR(ZS/Dz))
	NABL = FLOOR(HABL/Dz)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Memory allocation
	ALLOCATE(K(N), DK2(N), ALT(N), PHI(N,M+1), T(N,N), D(N,N), M1(N,N), &
	M2(N,N), MAT(N,N), I(N,N), P(N,M), LP(N,M), LP2(N,M), LPG(M), LPG2(M))
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	DO NZ = 1,N
		IF (NZ >= N-NABL) THEN
			K(NZ) = K0+ABL*im*(NZ-N+NABL)**2/NABL**2
		ELSE
			K(NZ) = K0
		ENDIF
		DK2(NZ) = im*(K(NZ)**2-(K0**2))/(2*K0)
	END DO
	TAU1 = 4/(3-2*im*K0*Dz)
	TAU2 = -1/(3-2*im*K0*Dz)
	SIGMA1 = 4/3
	SIGMA2 = -1/3
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	A0 = 1.3717
	A2 = -0.3701
	B = 3
	DO NZ = 1,N
		PHI(NZ,1) = SQRT(im*K0)* &
		((A0+A2*K0**2*(NZ*Dz-ZS)**2)*EXP(-K0**2*(NZ*Dz-ZS)**2/B) &
		+ (A0+A2*K0**2*(NZ*Dz+ZS)**2)*EXP(-K0**2*(NZ*Dz+ZS)**2/B))
	END DO
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	I = EYE(N)
	aa = im/(2*K0*Dz**2)
	T = 0.
	D = 0.
	T(1,1) = T(1,1) - 2 + SIGMA1
	T(1,2) = T(1,2) + 1 + SIGMA2
	D(1,1) = dk2(1)
	DO nz = 2,N-1
		T(NZ,NZ-1) = T(NZ,NZ-1)+1
		T(NZ,NZ) = T(NZ,NZ)-2
		T(NZ,NZ+1) = T(NZ,NZ+1)+1
		D(NZ,NZ) = DK2(NZ)
	END DO
	T(N,N-1) = T(N,N-1) + 1 + TAU2
	T(N,N) = T(N,N) - 2 + TAU1
	D(N,N) = DK2(N)
	M1 = I+(dr/2)*(aa*T+D)+(1/(2*im*K0))*(aa*T+D)
	M2 = I-(dr/2)*(aa*T+D)+(1/(2*im*K0))*(aa*T+D)
	!MAT = MATMUL(INV(M2),M1)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Forward-marching procedure
	CALL CPU_TIME(TI)
	!-----------------------------------------------------------
	IF (SP == 0) THEN
	!-----------------------------------------------------------
		MAT = MATMUL(INV(M2),M1)
		DO MX = 2,M+1
			WRITE(*, *) "Step :", MX-1, "out of", M
			PHI(1:N,MX) = MATMUL(MAT,PHI(1:N,mx-1))
			P(1:N,MX-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
		END DO
	!-----------------------------------------------------------
	ELSEIF (SP == 1) THEN
	!-----------------------------------------------------------
	  	MAT_SP = SPCONV(MAT) ! SLOW PROCEDURE !
	  	NELT = MAT_SP%NELT
	  	ALLOCATE(IMAT(NELT),JMAT(NELT))
	  	IMAT = MAT_SP%ROW
	  	JMAT = MAT_SP%COL
	  	VALMAT = MAT_SP%VAL
	  	DO MX = 2,M+1
	   		WRITE(*, *) "Step :", mx-1, "out of", M
	   		CALL DSMV(N,PHI(1:N,MX-1),PHI(1:N,MX),NELT,IMAT,JMAT,VALMAT,0)
	   		P(1:N,mx-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
	  	END DO
	!-----------------------------------------------------------
	ELSE
	!-----------------------------------------------------------
		ALLOCATE(SA1(N**3),SA2(N**3),IJA(N**3))
		THRESH = 0
		TOL = 1.0e-4_dp
		CALL SPRSIN_DP(M1,THRESH,SP1)
		CALL SPRSIN_DP(M2,THRESH,SP2)
		DO MX = 2,M+1
			WRITE(*, *) "Step :", mx-1, "out of", M
			CALL SPRSAX_DP(SP1,PHI(1:N,MX-1),TEMP)
			CALL LINBCG(SP2,TEMP,PHI(1:N,MX),1,TOL,10,ITER,ERR)
			P(1:N,mx-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
		END DO
	!-----------------------------------------------------------
	END IF
	!-----------------------------------------------------------
	CALL CPU_TIME(TF)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Transmission loss conversion
	P0 = P(NS,1)
	DO MX = 1,M
		DO NZ = 1,N
			LP(NZ,MX) = 20*LOG10(ABS(P(NZ,MX-1)/P0))
		END DO
	END DO
	LPG = LP(1,1:M)
	LP2 = LP/SQRT(2.)
	LPG2 = LP2(1,1:M)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Output
	OPEN(UNIT=10,FILE="PE2D_HFR_LPg.dat")
	OPEN(UNIT=20,FILE="PE2D_HFR_LP.dat")
	OPEN(UNIT=30,FILE="PE2D_HFR_P.dat")
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
	!----------------------------------------------------------
END PROGRAM PE2D_HFR
!-----------------------------------------------------------
!-------------      EXTERNAL PROCEDURES      ---------------
!-----------------------------------------------------------
