PROGRAM PE2D_HT_SSF
 USE PE2D_TYPE
 USE PE2D_VAR
 USE PE2D_AUX
 USE PE2D_GROUND
 IMPLICIT NONE
 !-----------------------------------------------------------
 !Variable statements
 !-----------------------------------------------------------
 !Constant values / parameters
 REAL(DP), PARAMETER :: c0 = 343, PI = 3.14159
 COMPLEX(DP), PARAMETER :: j = (0,1)
 !Inputs
 REAL(DP) :: L, H, f, zs, As
 !Propagation domain discretization
 INTEGER :: N, Ns, Na, M, nz, mx, nn
 REAL(DP) :: k0, lmbda0, Ha, dz, dr, dkz
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: r, z
 COMPLEX(DP) :: beta
 !Atmosphere variables
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: k, dk2, kz, Phase, Ref, DELTA
 !Terrain coordinate transform
 REAL(DP) :: H0, x0, s, xl, xr
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: alti, linear, alpha, dalpha
 !z-domain Pressure field
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: P, PHI, PHI1, PHI2, Preal, Pimag
 COMPLEX(DP) :: A0, A2, B
 REAL(DP) :: Ps
 !kz-domain Pressure field
 REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: WA, LP, LPrms, P0
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: CH, LPg, LPrmsg
 !FFT Work array
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: FPHI, FPHIi, FPHIr
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: BUFR, FPHIs
 INTEGER, DIMENSION(100) :: IFAC
 !CPU Time
 REAL :: start, finish
 !-----------------------------------------------------------
 !User interface
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
 WRITE(*,*) "GF2D_HFR"
 WRITE(*,*) "Dimension 	: 	2D"
 WRITE(*,*) "Atmosphere	:	Homogeneous"
 WRITE(*,*) "Boundary 	:	Flat"
 WRITE(*,*) "Impedance	:	Rigid"
 WRITE(*,*) ""
 WRITE(*,*) "================================================"
 WRITE(*,*) ""
 WRITE(*,*) "Enter inputs L, H, f, zs, As"
 WRITE(*,*) ""
 READ *, L, H, f, zs, As
 !----------------------------------------------------------
 !Problem parameters
 !-----------------------------------------------------------
 lmbda0 = c0/f
 k0 = 2*PI*f/c0
 dz = lmbda0/10
 dr = lmbda0/10
 Ha = 30*lmbda0
 N = floor(H/dz)
 M = floor(L/dr)
 Ns = MAX(1,floor(zs/dz))
 Na = floor(Ha/dz)
 !-----------------------------------------------------------
 !Memory allocation
 !-----------------------------------------------------------
 ALLOCATE(k(N), kz(N), Ref(N), Phase(N), dk2(N), z(N), DELTA(N), r(M), &
 BUFR(N), PHI(N,M+1), FPHI(N,M), FPHIi(N,M), FPHIr(N,M), FPHIs(M), &
 WA(2,N), CH(2*N), PHI1(N,M), PHI2(N,M), P(N,M), LP(N,M), LPrms(N,M), &
 P0(N,M), LPg(M), LPrmsg(M))
 !-----------------------------------------------------------
 !Discretization
 !-----------------------------------------------------------
 !z-direction wavenumber and kz-domain discretization
 z = (/ (nz*dz, nz = 1,N) /)
 dkz = 2*PI/(N*dz)
 DO nz = 1,N
  IF (nz >= N-Na) THEN
   k(nz) = k0+As*j*(nz-N+Na)**2/Na**2
  ELSE
   k(nz) = k0
  ENDIF
  kz(nz) = nz*dkz
  dk2(nz) = j*(k(nz)**2-(k0**2))/(2*k0)
  Phase(nz) = EXP(j*dr*(SQRT(k0**2-kz(nz)**2)-k0**2))
 END DO
 DELTA = EXP(j*dr*dk2/(2*k0))
 !-----------------------------------------------------------
 ALLOCATE(dalpha(1:M),r(0:M+1),linear(0:M+1),alpha(1:M+1))
 r = (/ (mx*dr, mx = 0,M+1) /)
 linear(0) = 0
 DO mx = 1,M+1
   xl = HILL(H0,x0,r(mx-1),s)
   xr = HILL(H0,x0,r(mx),s)
   alpha(mx) = (xr - xl)/dr
   linear(mx) = linear(mx-1) + (1 + alpha(mx))*dr
 END DO
 !-----------------------------------------------------------
 !Starting field
 !-----------------------------------------------------------
 A0 = 1.3717
 A2 = -0.3701
 B = 3
 DO nz = 1,N
  PHI(nz,1) = SQRT(j*k0)* &
  ((A0+A2*k0**2*(nz*dz-zs)**2)*EXP(-k0**2*(nz*dz-zs)**2/B) &
  + (A0+A2*k0**2*(nz*dz+zs)**2)*EXP(-k0**2*(nz*dz+zs)**2/B))
 END DO
 !-----------------------------------------------------------
 !Main Loop
 !-----------------------------------------------------------
 !Start of Forward-marching procedure
 CALL CPU_TIME(START)
 !Initialization of the FFT subroutines :
 !IFAC and WA depend only on N -> used in CFFTF1 and CFFTB1
 CALL CFFTI1(N,WA,IFAC)
 DO mx = 2,M+1
  WRITE(*, *) "Step :", mx-1, "out of", M
  !Surface term
  dalpha(mx-1) = (alpha(mx)-alpha(mx-1))/dr
  DO nz = 1,N
   FPHIs(mx-1) = FPHIs(mx-1)+PHI(nz,mx-1)*dz*EXP(-j*beta*nz*dz)
  END DO
  BUFR = PHI(:,mx-1)
  CALL CFFTF1(N,BUFR,CH,WA,IFAC)
  FPHIi(:,mx-1) = BUFR
  CALL FLIP(BUFR)
  FPHIr(:,mx-1) = BUFR
  FPHI(:,mx-1) = (FPHIi(:,mx-1)+FPHIr(:,mx-1))*Phase
  BUFR = FPHI(:,mx-1)
  CALL CFFTB1(N,BUFR,CH,WA,IFAC)
  PHI1(:,mx-1) = BUFR
  PHI2(:,mx-1) = 2*j*beta*FPHIs(mx-1)*EXP(-j*beta*z)*EXP(j*dr*(SQRT(k0**2-beta**2)-k0**2))
  PHI(:,mx) = DELTA*(PHI1(:,mx-1)+PHI2(:,mx-1))*EXP(j*k0*z*dalpha(mx-1))
  P(:,mx-1) = EXP(j*k0*(mx-1)*dr)*PHI(:,mx)*(1/SQRT((mx-1)*dr))
 END DO
 !End of Forward-marching procedure
 CALL CPU_TIME(FINISH)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Transmission loss conversion
 Ps = ABS(P(Ns,1))
 Preal = REAL(REAL(P))
 Pimag = REAL(AIMAG(P))
 DO mx = 1,M
  DO nz = 1,N
   LP(nz,mx) = 20*LOG10(ABS(P(nz,mx-1))/Ps)
  END DO
 END DO
 LPg = LP(1,1:M)
 LPrms = LP/SQRT(2.)
 LPrmsg = LPrms(1,1:M)
 !-----------------------------------------------------------
 !Output
 !-----------------------------------------------------------
 OPEN(UNIT=10,FILE="GF2D_HFR_LPg.dat")
 OPEN(UNIT=20,FILE="GF2D_HFR_LP.dat")
 OPEN(UNIT=30,FILE="GF2D_HFR_P.dat")
 DO nz = 1,N
   WRITE(20, 100, advance="no") LP(nz,1:M)
   WRITE(30, 101, advance="no") P(nz,1:M)
 END DO
 DO mx = 1,M
  WRITE(10, 100) mx*dr, LPg(mx), LPrmsg(mx)
 END DO
 100 FORMAT(3X,F8.3,$)
 101 FORMAT(3X,F8.3,SP,F8.3,SS,"i",$)
 PRINT *, "Main CPU time (s) :", FINISH-START
 PRINT *, "Source pressure P0 (dB) :", 20*LOG10(Ps)
 !-----------------------------------------------------------
END PROGRAM PE2D_HT_SSF
