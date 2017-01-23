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
	WRITE(*,*) "Dimension 	: 	2D"
	WRITE(*,*) "Atmosphere	:	Homogeneous"
	WRITE(*,*) "Boundary 	:	Flat"
	WRITE(*,*) "Impedance	:	Rigid"
	WRITE(*,*) ""
	WRITE(*,*) "================================================"
	WRITE(*,*) ""
	WRITE(*,*) "Enter inputs L, f, zs, As, angle"
	WRITE(*,*) ""
	READ *, L, F, ZS, ABL, ANGLE
	WRITE(*,*) ""
	!-----------------------------------------------------------
	!----------------------------------------------------------
	C0 = CMPLX(343.0_DP,0.0_DP,KIND=DP)
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
	WRITE(*,*) "Number of loops : ", M
	WRITE(*,*) ""
	WRITE(*,*) "================================================"
	NS = MAX(1,FLOOR(ZS/DZ))
	NABL = FLOOR(HABL/dz)
	ALPHA = IM/(2.*K0*DZ**2)
	KU = 1
	KL = 1
	NDIAG = KU+KL+1
	NDIAGX = KU+2*KL+1
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Memory allocation
	ALLOCATE(K(N), DK2(N), ALT(N), PHI(N,M+1), T(N,N), D(N,N), I(N,N), M1(N,N), &
	M2(N,N), MB1(NDIAG,N), MB2(NDIAG,N), P(N,M), LP(N,M), LP2(N,M), LPG(M), LPG2(M), &
	TERR2(1:M), LIN(0:M+1), TERR1(1:M+1), MB1X(NDIAG,N), MB2X(NDIAG,N), M2X(NDIAGX,N))
	!-----------------------------------------------------------
	!Atmosphere
	!-----------------------------------------------------------
	!
	DO NZ = 1,N
		IF (NZ >= N-NABL) THEN
			K(NZ) = K0+ABL*IM*(NZ-N+NABL)**2/NABL**2
		ELSE
			K(NZ) = K0
		ENDIF
		DK2(NZ) = IM*(K(NZ)**2-(K0**2))/(2*K0)
		ALT(NZ) = NZ*DZ
	END DO
	TAU1 = 4./(3.-2.*IM*K0*DZ)
	TAU2 = -1./(3.-2.*IM*K0*DZ)
	SIGMA1 = 4./3.
	SIGMA2 = -1./3.
	!-----------------------------------------------------------
	!Terrain
	!-----------------------------------------------------------
	LIN(0) = 0
	THETA = ANGLE*PI/180.
	S = 3*L/10
	H0 = SQRT(EXP(1.)/2.)*S*ATAN(THETA)
	DO MX = 1,M+1
		XL = GHILL(H0,X0,S,(MX-1)*Dr)
		XR = GHILL(H0,X0,S,MX*DR)
		TERR1(MX) = (XR-XL)/Dr
		LIN(MX) = LIN(MX-1)+(1+TERR1(MX))*Dr
	END DO
	!-----------------------------------------------------------
	!Starting field
	!-----------------------------------------------------------
	A0 = 1.3717
	A2 = -0.3701
	B = 3.
	DO NZ = 1,N
		PHI(NZ,1) = SQRT(IM*K0)* &
		((A0+A2*K0**2*(NZ*Dz-ZS)**2)*EXP(-K0**2*(NZ*Dz-ZS)**2/B) &
		+ (A0+A2*K0**2*(NZ*Dz+ZS)**2)*EXP(-K0**2*(NZ*Dz+ZS)**2/B))
	END DO
	!-----------------------------------------------------------
	!Matrix building
	!-----------------------------------------------------------
	!Initialization
	!---------
!	T(:,:) = 0.
!	D(:,:) = 0.
!	I = EYE(N)
!	!---------
!	T(1,1) = T(1,1)-2+SIGMA1
!	T(1,2) = T(1,2)+1+SIGMA2
!	D(1,1) = DK2(1)
!	DO NZ = 2,N-1
!		T(NZ,NZ-1) = T(NZ,NZ-1)+1
!		T(NZ,NZ) = T(NZ,NZ)-2
!		T(NZ,NZ+1) = T(NZ,NZ+1)+1
!		D(NZ,NZ) = DK2(NZ)
!	END DO
!	T(N,N-1) = T(N,N-1)+1+TAU2
!	T(N,N) = T(N,N)-2+TAU1
!	D(N,N) = DK2(N)
!	!---------
!	M1 = I+(DR/2.)*(ALPHA*T+D)+(1./(2.*IM*K0))*(ALPHA*T+D)
!	M2 = I-(DR/2.)*(ALPHA*T+D)+(1./(2.*IM*K0))*(ALPHA*T+D)
!	!---------
!	DO NZ = 1,N
!                NN = KU+1-NZ
!                DO MX = MAX(1,NZ-KU),MIN(N,NZ+KL)
!                      	MB1(NN+MX,NZ) = M1(MX,NZ)
!			MB2(NN+MX,NZ) = M2(MX,NZ)
!          	END DO
!        END DO
!	MB1X = MB1
!	MB2X = MB2
	!Tridiagonal FD matrix
	T(1,2:N) = ONES(N-1)		!Supdiag
	T(2,1:N) = -2*ONES(N)		!Diag
	T(3,1:N-1) = ONES(N-1)		!Subdiag
	T(1,2) = T(1,2)+SIGMA2
	T(2,1) = T(2,1)+SIGMA1
	T(2,N) = T(2,N)+TAU1
	T(3,N-1) = T(3,N-1)+TAU2
	DO I = 1,KU
		M1(I,:) = (1./(2.*IM*K0)+DR/2.)*ALPHA*T(I,:)
		M2(I,:) = (1./(2.*IM*K0)-DR/2.)*ALPHA*T(I,:)
	END DO
	DO I = NDIAG,KU+1,1
		M1(I,:) = (1./(2.*IM*K0)+DR/2.)*ALPHA*T(I,:)
		M2(I,:) = (1./(2.*IM*K0)-DR/2.)*ALPHA*T(I,:)
	END DO
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Forward-marching procedure
	CALL CPU_TIME(TI)
	!-----------------------------------------------------------
	DO MX = 2,M+1
		WRITE(*, *) "Step :", MX-1, "out of", M
		TERR2(MX-1) = (TERR1(MX)-TERR1(MX-1))/Dr
		WRITE(*, *) "1"
		M1(2,:) = ONES(N)+(1./(2.*IM*K0)+DR/2.)*ALPHA*T(2,:)+(1./(2.*IM*K0)+DR/2.)*(DK2(N)-IM*K0*TERR2(MX-1)*ALT(:))
		M2(2,:) = ONES(N)+(1./(2.*IM*K0)-DR/2.)*ALPHA*T(2,:)+(1./(2.*IM*K0)-DR/2.)*(DK2(N)-IM*K0*TERR2(MX-1)*ALT(:))
		!MB1X(KU+1,:) = MB1(KU+1,:)-(1./(2.*IM*K0)+DR/2.)*(IM*K0*TERR2(MX-1)*ALT(:))
		!MB2X(KU+1,:) = MB2(KU+1,:)-(1./(2.*IM*K0)-DR/2.)*(IM*K0*TERR2(MX-1)*ALT(:))
		ALLOCATE(AA(N-1),BB(N),CC(N-1),TEMP(N,1))
		!CGBMV (TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
		CALL ZGBMV('N',N,N,KL,KU,1.,MB1X,NDIAG,PHI(1:N,MX-1),1,0.,TEMP(1:N,1),1)
		AA = MB2X(3,1:N-1)  
		BB = MB2X(2,:)
		CC = MB2X(1,2:N)
		
		!AA, BB, CC destroyed after CGTSL
		CALL ZGBSV(N,KL,KU,AA,BB,CC,TEMP,N,INFO)
		IF (INFO.NE.0) THEN
			STOP "WARNING : Input Matrix is singular !"
		END IF
		PHI(1:N,MX) = TEMP(1:N,1)
		DEALLOCATE(AA,BB,CC,TEMP)
		WRITE(*, *) "2"
		P(1:N,MX-1) = EXP(IM*K0*(MX-1)*DR)*PHI(1:N,MX)*(1/SQRT((MX-1)*DR))
	END DO
	!-----------------------------------------------------------
	!ELSEIF (SP == 1) THEN
	!Sparse SLATEC
		!DO MX = 2,M+1
		!	WRITE(*, *) "Step :", MX-1, "out of", M
		!	TERR2(MX-1) = (TERR1(MX)-TERR1(MX-1))/Dr
		!	D = D-im*K0*TERR2(MX-1)*DIAG(ALT)
		!	WRITE(*, *) "1"
		!	M1 = I+(Dr/2)*(aa*T+D)+(1/(2*im*K0))*(aa*T+D)
		!	M2 = I-(Dr/2)*(aa*T+D)+(1/(2*im*K0))*(aa*T+D)
		!	!MAT = MATMUL(INV(M2),M1)
		!	WRITE(*, *) "2"
		!	!MAT_SP = SPCONV(MAT) ! SLOW PROCEDURE !
		!	CALL SPRSIN_DP(M2,0,SPM2)
		!	CALL SPRSIN_DP(M1,0,SPM1)
		!	NELT = MAT_SP%NELT
		!	ALLOCATE(IMAT(NELT),JMAT(NELT))
		!	IMAT = MAT_SP%ROW
		!	JMAT = MAT_SP%COL
		!	VALMAT = MAT_SP%VAL
		!	CALL DSMV(N,PHI(1:N,MX-1),PHI(1:N,MX),NELT,IMAT,JMAT,VALMAT,0)
		!	P(1:N,mx-1) = EXP(im*K0*(MX-1)*Dr)*PHI(1:N,MX)*(1/SQRT((MX-1)*Dr))
		!END DO
	!-----------------------------------------------------------
	CALL CPU_TIME(TF)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Transmission loss conversion
	P0 = ABS(P(NS,2))
	DO MX = 1,M
		DO NZ = 1,N
			LP(NZ,MX) = 20.*LOG10(ABS(P(NZ,MX-1)))
		END DO
	END DO
	LPG = LP(1,1:M)
	LP2 = LP/SQRT(2.)
	LPG2 = LP2(1,1:M)
	!-----------------------------------------------------------
	!-----------------------------------------------------------
	!Output
	OPEN(UNIT=10,FILE="PE2D_HTR_LPg.dat")
	OPEN(UNIT=20,FILE="PE2D_HTR_LP.dat")
	OPEN(UNIT=30,FILE="PE2D_HTR_P.dat")
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
END PROGRAM PE2D_HTR
!-----------------------------------------------------------
!-------------      EXTERNAL PROCEDURES      ---------------
!-----------------------------------------------------------

