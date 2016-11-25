PROGRAM PE2DFH

 IMPLICIT NONE
 
 PRINT *, ""
 PRINT *, "2D Parabolic Equation, Homogeneous, Flat rigid boundary"
 PRINT *, "University of Bristol A2RG - AWE Blacknest"
 PRINT *, "_______________________________________________________"
 PRINT *, ""
 PRINT *, "Codor Khodr, May 2016"
 PRINT *, ""
 PRINT *, ""
 
 !Global numerical values
 
 REAL, PARAMETER :: c0 = 343, PI = 3.14159
 COMPLEX, PARAMETER :: i = (0,1)
 !INTEGER, PARAMETER :: DP = KIND(0.d0)
 
 !Inputs 
 
 REAL :: L, H, f, zs, As
 COMPLEX :: Z
 
 PRINT *, "Enter inputs L, H, f, zs, As, Z"
 PRINT *, ""
 READ *, L, H, f, zs, As
 
 !Propagation domain discretization
 
 REAL :: k0, lmbda0, h, dz, dr
 INTEGER :: N, Ns, Na, M
 INTEGER :: nz, mx
 
 lmbda0 = c0/f
 k0 = 2*PI*f/c0
 dz = 1/k0
 dr = 1/k0
 
 N = floor(H/dz)
 M = floor(L/dx)
 Ns = floor(zs/dz)
 
 !Absorbing layer
 
 COMPLEX, ALLOCATABLE, DIMENSION(:) :: c, k, dk2
 ALLOCATE(c(N),k(N))
 
 c = c0*1.
 k = k0*1.
 DO nz = N-Na,N
  k(N) = k(N)+As*i*(nz-N+Na)**2/Na**2
 END DO
 dk2 = i*(k**2-(k0**2)*1.)/(2*k0)

 !Boundary conditions
 
 COMPLEX :: tau1, tau2, sigma1, sigma2
 
 tau1 = 4/(3-2*i*k0*dz)
 tau2 = -1/(3-2*i*k0*dz)
 sigma1 = 4/3
 sigma2 = -1/3
 
 !Initial pressure field
 
 COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: PHI
 COMPLEX :: A0, A2, B
 ALLOCATE(PHI(N,M+1))
 
 A0 = 1.3717
 A2 = -0.3701
 B = 3
 PHI = 0.
 DO nz = 1,N
  PHI(n,1) = SQRT(i*k0)*(A0+A2*k0**2*(n*dz-zs)**2)*EXP(-k0**2*(n*dz-zs)**2/B) &
   +SQRT(i*k0)*(A0+A2*k0**2*(n*dz+zs)**2)*EXP(-k0**2*(n*dz+zs)**2/B)
 END DO

 !Matrices allocation
 
 COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: T, D, M1, M2, I
 COMPLEX :: aa
 ALLOCATE(T(N,N), D(N,N) ,M1(N,N), M2(N,N), I(N,N), b(N))
 
 I = EYE(N)
 aa = i/(2*k0*dz**2)
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
 
 M1 = I+(dr/2)*(aa*T+D)+(A/(2*i*k0))*(aa*T+D)
 M2 = I-(dr/2)*(aa*T+D)+(A/(2*i*k0))*(aa*T+D)
 M = MATMUL(INV(M2),M1)

 !Forward-Marching procedure
 
 COMPLEX :: P0
 COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: P, LP
 ALLOCATE(P(N,M),LP(M,M))
 
 DO mx = 2,M+1
  PHI(:,mx) = MATMUL(M,PHI(:,mx-1))
  P(:,mx-1) = EXP(i*k0*(mx-1)*dr)*PHI(:,mx)*(1/SQRT((mx-1)*dr))
  IF (mx==2) P0 = P(Ns,mx-1)
  LP(:,mx-1) = ABS(P(:,mx-1)/P0)
 END DO
 
 !Output format
 
 OPEN(unit=1,file="Hom_2D_flat_P.txt")
 OPEN(unit=2,file="Hom_2D_flat_LP.txt")
 DO nz = 1,N
  WRITE(1, 100) (P(nz,mx), mx = 1,M,1)
  WRITE(2, 100) (LP(nz,mx), mx = 1,M,1)
 END DO
 100 FORMAT(F0.0,SP,F0.0,"i")
 
CONTAINS

 !Subfunctions
 
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

  REAL(DP), DIMENSION(:,:), INTENT(IN) :: A
  REAL(DP), DIMENSION(SIZE(A,1), SIZE(A,2)) :: invA
  REAL(DP), DIMENSION(SIZE(A,1)) :: WORK
  REAL(DP), DIMENSION(SIZE(A,1)) :: IPIV
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
 
END PROGRAM PE2FH

 
 
 
 
