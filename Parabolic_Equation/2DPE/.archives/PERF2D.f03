PROGRAM PERF2D
	USE NRTYPE
	USE MATLIB
	USE PARAM
	USE GROUND
	IMPLICIT NONE
	!-----------------------------------------------------------
	!Compile instructions in infratopo/PE/2DPE:
	!
	!gfortran -g -fbacktrace -Wtabs -o PERF2D PERF2D.f03 nrtype.f03 param.f03 ground.f03 matlib.f03 -llapack -lblas
	!*** nrtype, param, ground, matlib : module files ***
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
	WRITE(*,*) "PERF2D - Performance Test Function for PE"
	WRITE(*,*) "Dimension 	: 	2D"
	WRITE(*,*) "Atmosphere	:	Homogeneous"
	WRITE(*,*) "Boundary 	:	Flat"
	WRITE(*,*) "Impedance	:	Rigid"
	WRITE(*,*) ""
	WRITE(*,*) "================================================"
	WRITE(*,*) ""
	WRITE(*,*) "propagation range = 10000 m"
	WRITE(*,*) "frequency = 10Hz , wavelength = 34.3 m"
	WRITE(*,*) "source height = 5 m"
	WRITE(*,*) "max slope of terrain = 40 degrees"
	WRITE(*,*) ""
	!-----------------------------------------------------------
	!----------------------------------------------------------
	L = 5000
	F = 10
	ZS = 5
	ABL = 0.07
	ANGLE = 40
	C0 = 343
	lmbda0 = C0/F
	K0 = 2*PI*F/C0
	Dz = lmbda0/10
	Dr = lmbda0/10
	H = 100*lmbda0
	HABL = 30*lmbda0
	N = FLOOR(H/dz)
	M = FLOOR(L/dr)
	WRITE(*,*) "================================================"
	WRITE(*,*) ""
	WRITE(*,*) "Size of the system : ", N
	WRITE(*,*) "Number of loops : ", M
	WRITE(*,*) ""
	WRITE(*,*) "================================================"
	NS = MAX(1,FLOOR(ZS/DZ))
	NABL = FLOOR(HABL/dz)
	ALPHA = IM/(2*K0*DZ**2)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Memory allocation
	ALLOCATE(K(N), DK2(N), ALT(N), PHI(N,M+1), T(N,N), D(N,N), M1(N,N), &
	M2(N,N), I(N,N), P(N,M), LP(N,M), LP2(N,M), LPG(M), LPG2(M), &
	TERR2(1:M), LIN(0:M+1), TERR1(1:M+1), E(N,N))
	!-----------------------------------------------------------
	!Atmosphere
	!-----------------------------------------------------------
	DO NZ = 1,N
		IF (NZ >= N-NABL) THEN
			K(NZ) = K0+ABL*im*(NZ-N+NABL)**2/NABL**2
		ELSE
			K(NZ) = K0
		ENDIF
		DK2(NZ) = im*(K(NZ)**2-(K0**2))/(2*K0)
		ALT(NZ) = NZ*Dz
	END DO
	TAU1 = 4./(3.-2.*im*K0*Dz)
	TAU2 = -1./(3.-2.*im*K0*Dz)
	SIGMA1 = 4./3.
	SIGMA2 = -1./3.
	!-----------------------------------------------------------
	!Terrain
	!-----------------------------------------------------------
	LIN(0) = 0
	THETA = ANGLE*PI/180
	H0 = SQRT(EXP(1.)/2)*L*ATAN(THETA)
	S = 3*L/10
	DO MX = 1,M+1
		XL = GHILL(H0,X0,(MX-1)*DR,S)
		XR = GHILL(H0,X0,MX*DR,S)
		TERR1(MX) = (XR-XL)/Dr
		LIN(MX) = LIN(MX-1)+(1+TERR1(MX))*Dr
	END DO
	!-----------------------------------------------------------
	!Starting field
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
	!Matrix building
	!-----------------------------------------------------------
	!Initialization (Dense Matrices)
	!---------
	T(:,:) = 0.
	D(:,:) = 0.
	E(:,:) = 0.
	I = EYE(N)
	!---------
	DO NZ = 2,N !Out of bound elements ignored
		T(NZ,NZ-1) = T(NZ,NZ-1)+1
		T(NZ,NZ) = T(NZ,NZ)-2
		T(NZ,NZ+1) = T(NZ,NZ+1)+1
		D(NZ,NZ) = DK2(NZ)
		E(NZ,NZ) = NZ*DZ
	END DO
	T(1,1) = T(1,1)-2+SIGMA1
	T(1,2) = T(1,2)+1+SIGMA2
	T(N,N-1) = T(N,N-1)+1+TAU2
	T(N,N) = T(N,N)-2+TAU1
	!-----------------------------------------------------------
	!Forward-marching procedure
	CALL CPU_TIME(TI)
	!-----------------------------------------------------------
	DO MX = 2,M+1
		WRITE(*, *) "Step :", MX-1, "out of", M
		TERR2(MX-1) = (TERR1(MX)-TERR1(MX-1))/Dr		!Varying coefficient to be injected...
		D = D-IM*K0*TERR2(MX-1)*E				!... inside the D matrix
		WRITE(*, *) ".... D Updated"
		M1 = I+(Dr/2)*(ALPHA*T+D)+(1/(2*im*K0))*(ALPHA*T+D)	!Now the matrices M1 and M2 vary
		M2 = I-(Dr/2)*(ALPHA*T+D)+(1/(2*im*K0))*(ALPHA*T+D)	!from step to step because of D
		ALLOCATE(TEMP(N),MTEMP(N,N),IPIV(N))		
		TEMP = MATMUL(M1,PHI(1:N,MX-1))
		WRITE(*, *) ".... Right-hand side multiplied"
		MTEMP = M2			
		CALL ZGESV(N,1,MTEMP,N,IPIV,TEMP,N,INFO)
		WRITE(*, *) ".... System solved"
		PHI(1:N,MX) = TEMP
		P(1:N,MX-1) = EXP(IM*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
		DEALLOCATE(IPIV,TEMP,MTEMP)
	END DO
	!-----------------------------------------------------------
	CALL CPU_TIME(TF)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Transmission loss conversion
	P0 = P(NS,1)
	DO MX = 1,M
		DO NZ = 1,N
			LP(NZ,MX) = 20.*LOG10(ABS(P(NZ,MX-1)/P0))
		END DO
	END DO
	LPG = LP(1,1:M)
	LP2 = LP/SQRT(2.)
	LPG2 = LP2(1,1:M)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Output
	OPEN(UNIT=10,FILE="GROUNDLVL_TL.dat")
	OPEN(UNIT=20,FILE="FULLFIELD_TL.dat")
	OPEN(UNIT=30,FILE="FULLFIELD_P.dat")
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
	END PROGRAM PERF2D
	!-----------------------------------------------------------
	!-------------      EXTERNAL PROCEDURES      ---------------

