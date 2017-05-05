MODULE PE2D_VAR
 USE PE2D_TYPE
 IMPLICIT NONE
 !-----------------------------------------------------------
 COMPLEX(DPC), PARAMETER :: IM = CMPLX(0.0_DP,1.0_DP,KIND=DP)
 CHARACTER :: OUTPUT
 INTEGER :: N, NT, NSYS, M, NS, NA, NZ, MX, NELT, KU, KL, INFO, NDIAG, WIDE
 INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV
 REAL(DP) :: C0, L, LX, LY, H, HA, F, ZS, AS, K0, LMBDA0, DKZ, DZ, DR, A0, A2, B, TI, TF, P0, NORM, RCOND
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: ALT, X, LPG, LPG2, WORK, RWORK
 REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: LP, LP2, MCOND
 REAL(DP) :: H0, X0, S, XL, XR, ANGLE, THETA
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: GHX, LIN, TERR, TERR1, TERR2
 COMPLEX(DPC) :: A, ALPHA, BETA, BETA1, BETA2, TAU1, TAU2, SIGMA1, SIGMA2, IMP, R
 COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: C, K, DK2, CWORK
 COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: PHI, P, T, D, I, E, M1, M2, MB1, MB2, M1T, M2T, MB1T, MB2T, TEMP
 TYPE(US), SAVE :: SPMAT, SP1, SP2
 !-----------------------------------------------------------
END MODULE PE2D_VAR

! !Propagation domain discretization
!### INTEGER :: N, Ns, Na, M, nz, mx
!### REAL(DP) :: k0, lmbda0, Ha, dz, dr, dkz
! REAL(DP), ALLOCATABLE, DIMENSION(:) :: r, z
! COMPLEX(DP) :: beta
! !Atmosphere variables
! COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: k, dk2, kz, Phase, Ref, DELTA
! !Terrain coordinate transform
! REAL(DP) :: H0, x0, s, xl, xr
! REAL(DP), ALLOCATABLE, DIMENSION(:) :: alti, linear, alpha, dalpha
! !z-domain Pressure field
! COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: P, PHI, PHI1, PHI2, Preal, Pimag
! COMPLEX(DP) :: A0, A2, B
! REAL(DP) :: Ps
! !kz-domain Pressure field
! REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: WA, LP, LPrms, P0
! REAL(DP), ALLOCATABLE, DIMENSION(:) :: CH, LPg, LPrmsg
! !FFT Work array
! COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: FPHI, FPHIi, FPHIr
! COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: BUFR, FPHIs
! INTEGER, DIMENSION(100) :: IFAC
! !CPU Time
! REAL :: start, finish
