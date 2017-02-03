MODULE PARAM
	USE NRTYPE
	IMPLICIT NONE
	!-----------------------------------------------------------
	COMPLEX(DP), PARAMETER :: IM = CMPLX(0.0_DP,1.0_DP,KIND=DP)
	INTEGER :: N, M, NS, NABL, NZ, MX, NELT, ITER, KU, KL, INFO, NDIAG, NN
	INTEGER, ALLOCATABLE, DIMENSION(:) :: IMAT, JMAT, IJA, IPIV
	REAL(DP) :: C0, L, H, HABL, F, ZS, ABL, K0, LMBDA0, DZ, DR, A0, A2, B, TI, TF, P0, NORM, RCOND
	REAL(DP), EXTERNAL :: ZLANGE, ZLANGB
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: LPG, LPG2, WORK, RWORK, C1, C2
	REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: LP, LP2
	REAL(DP) :: H0, X0, S, XL, XR, ANGLE, THETA
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: GHX, LIN, TERR, TERR1, TERR2
	COMPLEX(DP) :: ALPHA, BETA1, BETA2, TAU1, TAU2, SIGMA1, SIGMA2, IMP, R
	COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: ALT, C, K, DK2, VALMAT, AA, BB, CC, CWORK
	COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI, P, T, D, I, E, M1, M2, MTEMP, MB1, MB2, TEMP, MTEMP1
	TYPE(SPRS2_DP), SAVE :: SPMAT, SP1, SP2
	!-----------------------------------------------------------
END MODULE PARAM
