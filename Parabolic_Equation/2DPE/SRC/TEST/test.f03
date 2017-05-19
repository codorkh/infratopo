PROGRAM TEST
 USE PE2D_TYPE
 USE PE2D_AUX
 USE PE2D_GROUND
 USE PE2D_ATMOS
 IMPLICIT NONE 
 INTEGER(I4B) :: I, NLINES, NCOLS
 CHARACTER :: TRACK, FREQ
 CHARACTER(LEN=100) :: PATH_TOPO, FILE_TOPO, TOPO
 REAL(DP), DIMENSION(:,:), ALLOCATABLE :: S
 !REAL(DP), DIMENSION(:,:), ALLOCATABLE :: YP
 !REAL(DP), DIMENSION(:), ALLOCATABLE :: X, XA, YA
 REAL(DP) :: F, C0, L, L1, DR, DR1, LAMBDA
 !INTEGER(I4B) :: M, NDER
 !--------
 !Data extraction
 TRACK     = "2"
 FREQ      = "2"
 PATH_TOPO = "../../GROUND/SnT_2017/L" // TRACK // "/"
 FILE_TOPO = "asc_L" // TRACK // "_2d_f" // FREQ // "_sampled.dat"
 TOPO      = TRIM(PATH_TOPO) // TRIM(FILE_TOPO)
 WRITE(*,*) TRIM(TOPO)
 CALL READ_ARRAY(TRIM(TOPO),S,NLINES,NCOLS)
 WRITE(*,*) NLINES, NCOLS
 OPEN(UNIT=999,FILE='read_data.dat')
 DO I=1,NLINES
  WRITE(999,100) S(I,1:NCOLS)
 END DO
 F = 3.79_DP
 C0 = 343.0_DP
 LAMBDA = C0/F
 DR = LAMBDA/10
 DR1 = S(2,1)-S(1,1)
 L = NLINES*DR
 L1 = S(NLINES,1)
 WRITE(*,*) DR, DR1
 WRITE(*,*) L, L1
 100 FORMAT(4(F16.5,X))
 !--------
 !Interpolation
 !ALLOCATE(XA(NLINES),YA(NLINES))
 !XA = S(1,:)
 !YA = S(2,:)
 !NDER = 2
 !M = 3000
 !L = XA(N)
 !DR = L/M
 !ALLOCATE(X(M),YP(M,0:NDER))
 !X = (/ (DR*I, I=1,M) /)
 !YP = POLINT1(N,X,NDER,XA,YA)
 !OPEN(UNIT=998,FILE='test_polint.dat')
 !DO I=1,M
 ! WRITE(998,101) X(I), YP(I,:)  
 !END DO
 !101 FORMAT(4(F16.5,X))
 !--------
 !CALL GET_ENVIRONMENT_VARIABLE("nlines", FILENAME)
 !WRITE(*,*) TRIM(FILENAME)
END PROGRAM
