PROGRAM PE2D_HFR
 USE PARAM
 USE MATLIB
 USE FDLIB
 USE GROUND
 IMPLICIT NONE
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Inputs
 REAL(DP) :: L, H, f, zs, As
 !Propagation domain discretization
 INTEGER :: N, Ns, Na, M, nz, mx, nn
 REAL(DP) :: k0, lmbda0, Ha, dz, dr
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: r, z
 COMPLEX(DP) :: aa
 !Absorbing layer
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: k, dk2
 !Boundary conditions
 COMPLEX(DP) :: tau1, tau2, sigma1, sigma2
 !Terrain coordinate transform
 REAL(DP) :: H0, x0, s, xl, xr, angle
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: r, z, alti, linear, alpha, dalpha
 !Initial pressure field
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: P, PHI
 COMPLEX(DP) :: A0, A2, B
 REAL(DP) :: Ps
 !Matrices
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: T, D, M1, M2, MAT, I
 REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: LP, LPrms, Preal, Pimag, P0
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: LPg, LPrmsg
 ! CPU Time
 REAL :: start, finish
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
 WRITE(*,*) "Enter inputs L, H, f, zs, As"
 WRITE(*,*) ""
 READ *, L, H, f, zs, As
 WRITE(*,*) ""
 WRITE(*,*) "Enter terrain - Guassian hill maximum angle (degrees)"
 WRITE(*,*) ""
 READ *, angle
 !-----------------------------------------------------------
 !----------------------------------------------------------
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
 !-----------------------------------------------------------
 !Memory allocation
 ALLOCATE(k(N), dk2(N), z(N), r(M), PHI(N,M+1), T(N,N), D(N,N) ,M1(N,N), &
 M2(N,N), MAT(N,N), I(N,N), P(N,M), LP(N,M), LPrms(N,M), Preal(N,M), &
 Pimag(N,M), P0(N,M), LPg(M), LPrmsg(M))
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 DO nz = 1,N
  IF (nz >= N-Na) THEN
   k(nz) = k0+As*j*(nz-N+Na)**2/Na**2
  ELSE
   k(nz) = k0
  ENDIF
  dk2(nz) = j*(k(nz)**2-(k0**2))/(2*k0)
 END DO
 tau1 = 4/(3-2*j*k0*dz)
 tau2 = -1/(3-2*j*k0*dz)
 sigma1 = 4/3
 sigma2 = -1/3
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 ALLOCATE(dalpha(1:M),r(0:M+1),linear(0:M+1),alpha(1:M+1))
 r = (\(mx*dr, mx = 0,M+1)\)
 linear(0) = 0
 theta = angle*PI/180
 H0 = SQRT(EXP(1)/2)*L*theta
 DO mx = 1,M+1
   xl = HILL(H0,x0,r(mx-1),s)
   xr = HILL(H0,x0,r(mx),s)
   alpha(mx) = (xr - xl)/dr
   linear(mx) = linear(mx-1) + (1 + alpha(mx))*dr
 END DO
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 A0 = 1.3717
 A2 = -0.3701
 B = 3
 DO nz = 1,N
  PHI(nz,1) = SQRT(j*k0)* &
  ((A0+A2*k0**2*(nz*dz-zs)**2)*EXP(-k0**2*(nz*dz-zs)**2/B) &
  + (A0+A2*k0**2*(nz*dz+zs)**2)*EXP(-k0**2*(nz*dz+zs)**2/B))
 END DO
 I = EYE(N)
 aa = j/(2*k0*dz**2)
 T = 0.
 D = 0.
 T(1,1) = T(1,1) - 2 + sigma1
 T(1,2) = T(1,2) + 1 + sigma2
 D(1,1) = dk2(1)
 DO nz = 2,N-1
  T(nz,nz-1) = T(nz,nz-1)+1
  T(nz,nz) = T(nz,nz)-2
  T(nz,nz+1) = T(nz,nz+1)+1
  D(nz,nz) = dk2(nz)
 END DO
 T(N,N-1) = T(N,N-1) + 1 + tau2
 T(N,N) = T(N,N) - 2 + tau1
 D(N,N) = dk2(N)
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 !Forward-marching procedure
 CALL CPU_TIME(START)
 DO mx = 2,M+1
  WRITE(*, *) "Step :", mx-1, "out of", M
  dalpha(mx-1) = (alpha(mx)-alpha(mx-1))/dr
  D = D - j*k0*dalpha(mx-1)*MATMUL(I,z) !SLOW!
  M1 = I+(dr/2)*(aa*T+D)
  M2 = I-(dr/2)*(aa*T+D)
  MAT = MATMUL(INV(M2),M1)
  PHI(1:N,mx) = MATMUL(MAT,PHI(1:N,mx-1))
  P(1:N,mx-1) = EXP(j*k0*(mx-1)*dr)*PHI(1:N,mx)*(1/SQRT((mx-1)*dr))
 END DO
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
 !-----------------------------------------------------------
 !Output
 OPEN(UNIT=10,FILE="PE2D_HFR_LPg.dat")
 OPEN(UNIT=20,FILE="PE2D_HFR_LP.dat")
 OPEN(UNIT=30,FILE="PE2D_HFR_P.dat")
 DO nz = 1,N
   WRITE(20, 100) LP(nz,1:M)
   WRITE(30, 101) P(nz,1:M)
 END DO
 DO mx = 1,M
  WRITE(10, 100) mx*dr, LPg(mx), LPrmsg(mx)
 END DO
 100 FORMAT(*(3X,F8.3))
 101 FORMAT(*(3X,F8.3,SP,F8.3,SS,"i"))
 PRINT *, "Main CPU time (s) :", FINISH-START
 PRINT *, "Source pressure P0 (Pa) :", Ps
 !-----------------------------------------------------------
END PROGRAM PE2D_HTR
!-----------------------------------------------------------
!-------------      EXTERNAL PROCEDURES      ---------------
!-----------------------------------------------------------
