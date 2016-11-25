PROGRAM PE2D_HFR

 IMPLICIT NONE
 
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 
 !Precision : double precision
 INTEGER, PARAMETER :: DP = KIND(0.d0)
 
 !Constant values / parameters
 REAL(DP), PARAMETER :: c0 = 343, PI = 3.14159
 COMPLEX(DP), PARAMETER :: j = (0,1)
 
 !Inputs 
 REAL(DP) :: L, H, f, zs, As
 
 !Propagation domain discretization
 INTEGER :: N, Ns, Na, M, nz, mx
 REAL(DP) :: k0, lmbda0, Ha, dz, dr
 COMPLEX(DP) :: aa
 
 !Absorbing layer
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: c, k, dk2
 
 !Boundary conditions
 COMPLEX(DP) :: tau1, tau2, sigma1, sigma2
 
 !Initial pressure field
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: PHI
 COMPLEX(DP) :: A0, A2, B, P0
 
 !Matrices
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: T, D, M1, M2, MAT, I
 COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: P, LP
 
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 
 WRITE(*,*) "======================================================="
 WRITE(*,*) "==============        PE2D_HFR        ================="
 WRITE(*,*) "======================================================="
 WRITE(*,*) ""
 WRITE(*,*) "2D Parabolic Equation, Homogeneous, Flat rigid boundary"
 WRITE(*,*) "University of Bristol A2RG - AWE Blacknest"
 WRITE(*,*) ""
 WRITE(*,*) ""
 WRITE(*,*) "Codor Khodr, July 2016"
 WRITE(*,*) ""
 WRITE(*,*) "======================================================="
 
 WRITE(*,*) "Enter inputs L, H, f, zs, As"
 WRITE(*,*) ""
 READ *, L, H, f, zs, As
 
 lmbda0 = c0/f
 k0 = 2*PI*f/c0
 dz = 1/k0
 dr = 1/k0
 Ha = 30*lmbda0
 
 N = floor(H/dz)
 M = floor(L/dr)
 Ns = floor(zs/dz)
 Na = floor(Ha/dz)
 
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 
 ALLOCATE(c(N),k(N))
 ALLOCATE(PHI(N,M+1))
 ALLOCATE(T(N,N), D(N,N) ,M1(N,N), M2(N,N), MAT(N,N), I(N,N))
 ALLOCATE(P(N,M),LP(M,M))
 
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 
 c = c0*1.
 k = k0*1.
 DO nz = N-Na,N
  k(N) = k(N)+As*j*(nz-N+Na)**2/Na**2
 END DO
 dk2 = j*(k**2-(k0**2)*1.)/(2*k0)
 
 tau1 = 4/(3-2*j*k0*dz)
 tau2 = -1/(3-2*j*k0*dz)
 sigma1 = 4/3
 sigma2 = -1/3
 
 A0 = 1.3717
 A2 = -0.3701
 B = 3
 PHI = 0.
 DO nz = 1,N
  PHI(n,1) = SQRT(j*k0)*(A0+A2*k0**2*(n*dz-zs)**2)*EXP(-k0**2*(n*dz-zs)**2/B) &
   +SQRT(j*k0)*(A0+A2*k0**2*(n*dz+zs)**2)*EXP(-k0**2*(n*dz+zs)**2/B)
 END DO

 I = EYE(N)
 aa = j/(2*k0*dz**2)
 T = 0.
 D = 0.
 
 T(1,1) = sigma1
 T(1,2) = sigma2
 D(1,1) = dk2(1)
 DO nz = 2,N-1
  T(nz,nz-1) = 1
  T(nz,nz) = -2
  T(nz,nz+1) = 1
  D(nz,nz) = dk2(nz)
 END DO
 T(N,N-1) = tau2
 T(N,N) = tau1
 D(N,N) = dk2(N)
 
 M1 = I+(dr/2)*(aa*T+D)+(1/(2*j*k0))*(aa*T+D)
 M2 = I-(dr/2)*(aa*T+D)+(1/(2*j*k0))*(aa*T+D)
 MAT = MATMUL(INV(M2),M1)

 !-----------------------------------------------------------
 !-----------------------------------------------------------

 DO mx = 2,M+1
  PHI(:,mx) = MATMUL(MAT,PHI(:,mx-1))
  P(:,mx-1) = EXP(j*k0*(mx-1)*dr)*PHI(:,mx)*(1/SQRT((mx-1)*dr))
  IF (mx==2) P0 = P(Ns,mx-1)
  LP(:,mx-1) = ABS(P(:,mx-1)/P0)
 END DO
 
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 
 !Output format
 
 OPEN(unit=1,file="Hom_2D_flat_P.txt")
 OPEN(unit=2,file="Hom_2D_flat_LP.txt")
 DO nz = 1,N
  WRITE(1, 100) (P(nz,mx), mx = 1,M,1)
  WRITE(2, 100) (LP(nz,mx), mx = 1,M,1)
 END DO
 100 FORMAT(F0.0,SP,F0.0,"j")
 
 !-----------------------------------------------------------
 !-----------------------------------------------------------
 
CONTAINS
 
 REAL FUNCTION INFINITY()
 
  IMPLICIT NONE
  CHARACTER(len=3) :: inf = "INF"
  READ(inf,*) INFINITY
  
 END FUNCTION INFINITY

 FUNCTION EYE(Ni) RESULT(eyeN)
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Ni
  REAL, DIMENSION(Ni,Ni) :: eyeN
  INTEGER :: r, c
  
  DO r = 1,Ni
   DO c = 1,Ni
    eyeN(r,c) = 0
   END DO
   eyeN(r,r) = 1
  END DO
 END FUNCTION EYE
 
 FUNCTION INV(A) RESULT(invA)

  IMPLICIT NONE

  COMPLEX(DP), DIMENSION(:,:), INTENT(IN) :: A
  COMPLEX(DP), DIMENSION(SIZE(A,1), SIZE(A,2)) :: invA
  COMPLEX(DP), DIMENSION(SIZE(A,1)) :: WORK
  COMPLEX(DP), DIMENSION(SIZE(A,1)) :: IPIV
  INTEGER :: rows, INFO

  invA = A
  rows = SIZE(A,1)

  CALL DGETRF(rows,rows,invA,rows,IPIV,INFO)

  IF (INFO/= 0) THEN
   STOP "WARNING : Input Matrix is singular !"
   END IF

  CALL DGETRI(rows,invA,rows,IPIV,WORK,rows,INFO)

  IF (INFO/= 0) THEN
   STOP "WARNING : Input Matrix is singular !"
  END IF

 END FUNCTION INV
 
END PROGRAM PE2D_HFR

 
 
 
 
