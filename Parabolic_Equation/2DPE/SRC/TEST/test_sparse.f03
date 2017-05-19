PROGRAM TEST_SPARSE
 USE PE2D_TYPE
 USE PE2D_AUX
 USE PE2D_GROUND
 USE PE2D_ATMOS
 USE PE2D_NUM
 IMPLICIT NONE
 !--------
 INTEGER(I4B) :: I, N, NELT, ISYM, KL, KU
 INTEGER(I4B) :: ITOL, ITMAX, ITER, IERR, IUNIT, LENW, LENIW
 INTEGER(I4B), ALLOCATABLE, DIMENSION(:) :: IA, JA, U, V, W, IWORK
 REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: A
 REAL(DP), ALLOCATABLE, DIMENSION(:) :: X, X0, X1, Y, B, VAL, R, Z, P, RR, ZZ, PP, DZ, AA, RWORK, H
 COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: FD, K
 COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: M, MB
 REAL(DP) :: THRESH
 REAL(DP) :: TOL, ERR
 TYPE(US) :: SPA
 !EXTERNAL :: MATVEC, MTTVEC, MSOLVE, MTSOLV
 !--------
 N = 6
 ALLOCATE(AA(N**2), A(N,N), B(N), Y(N), X(N), X0(N), X1(N), R(N), Z(N), P(N), RR(N), ZZ(N), PP(N), DZ(N))
 THRESH = 0.0_DP
 AA = (/ 1.0_DP,  0.0_DP, -1.0_DP,  0.0_DP,  1.0_DP,  5.0_DP, &
        -1.0_DP,  1.0_DP,  0.0_DP,  2.0_DP, -1.0_DP,  0.0_DP, &
         0.0_DP,  1.0_DP, -2.0_DP,  0.0_DP, -3.0_DP,  0.0_DP, &
         0.0_DP,  1.0_DP, -2.0_DP,  0.0_DP,  1.0_DP, -2.0_DP, &
         0.0_DP,  0.0_DP,  0.0_DP,  1.0_DP,  5.0_DP,  2.0_DP, &
        -2.0_DP,  0.0_DP,  1.0_DP, -5.0_DP,  0.0_DP,  4.0_DP /)
 A = RESHAPE(AA, (/N,N/))
 X = (/ 1.0_DP, 0.0_DP, -1.0_DP, 2.0_DP, -1.0_DP, 3.0_DP /)
 DO I = 1,N
  WRITE(*,1000) A(I,:)
 END DO
 WRITE(*,*)
 CALL D2US(A,N,SPA,THRESH)
 !--------
 NELT = SPA%NELT
 ALLOCATE(IA(NELT),JA(NELT),VAL(NELT))
 IA = SPA%IROW
 JA = SPA%JCOL
 VAL = SPA%VAL
 DO I = 1,NELT
  WRITE(*,'(I5,X)',ADVANCE='no') IA(I)
 END DO
 WRITE(*,*)
 DO I = 1,NELT
  WRITE(*,'(I5,X)',ADVANCE='no') JA(I)
 END DO
 WRITE(*,*)
 DO I = 1,NELT
  WRITE(*,1001,ADVANCE='no') VAL(I)
 END DO
 WRITE(*,*)
 WRITE(*,*)
 !--------
 ISYM = 0
 CALL DS2Y(N,NELT,IA,JA,VAL,ISYM)
 CALL DSMV(N,X,Y,NELT,IA,JA,VAL,ISYM)
 !-------- 
 B = MATMUL(A,X)
 !--------
 CALL SPARSE_MV(SPA,Z,X)
 !--------
 DO I = 1,N
  WRITE(*,1000) Y(I), B(I), Z(I)
 END DO
 WRITE(*,*)
 !--------
 ITOL=1
 TOL=0.0005_DP
 ITMAX=1000
 IUNIT=20
 LENW=NELT+8*N
 LENIW=NELT+4*N+12
 ALLOCATE(IWORK(LENIW),RWORK(LENW))
 X0 = 0.0_DP
 X1 = 0.0_DP
 CALL DSLUBC(N,B,X0,NELT,IA,JA,VAL,ISYM,ITOL,TOL,ITMAX,ITER,ERR,IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
 CALL SPARSE_SOLVE(SPA,X1,B)
 !--------
 DO I=1,N
  WRITE(*,1000) X1(I), X0(I), X(I)
 END DO
 WRITE(*,*)
 !--------
 KL = 1
 KU = 1
 ALLOCATE(FD(KL+KU+1), M(N,N))
 FD = (/1.0_DP, -2.0_DP, 1.0_DP/)
 M = DBAND(FD,N,KL,KU)
 MB = BBAND(FD,N,KL,KU)
 DO I=1,N
  WRITE(*,999) M(I,:) 
 END DO
 WRITE(*,*)
 DO I=1,KL+KU+1
  WRITE(*,999) MB(I,:) 
 END DO
 WRITE(*,*)
 !--------
 ALLOCATE(H(N),K(N))
 H = (/ (I*2.0_DP, I=1,N) /)
 K = CMPLX(H,0.0_DP,KIND=DPC)
 WRITE(*,*)
 DO I=1,N
  WRITE(*,*) K(I)
 END DO
 !--------
 ALLOCATE(U(10),V(10))
 U = (/ 1, 2, 3, 8, 9, 0, 0, 0, 5, 6 /)
 V = (/ 1, 1, 2, 4, 5, 6, 0, 2, 7, 7 /)
 CALL INTERSECT(U,V,W)
 WRITE(*,*)
 WRITE(*,*) "U =", U
 WRITE(*,*) "V =", V
 WRITE(*,*) "W =", W
 WRITE(*,*)
 DO I = 1,SIZE(Y)
  WRITE(*,1001) Y(I)
 END DO
 !--------
 
 !--------
 999 FORMAT(6(F5.2,SP,F5.2,SS,"i",X))
 1000 FORMAT(6(F5.1,X))
 1001 FORMAT(F6.1,X)
END PROGRAM
