PROGRAM TEST3
 USE PE2D_TYPE
 USE PE2D_AUX
 USE PE2D_GROUND
 USE PE2D_ATMOS
 USE PE2D_NUM
 IMPLICIT NONE
 !--------
 INTEGER(I4B) :: I, J, N
 REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: C
 REAL(DP), DIMENSION(:,:), ALLOCATABLE :: YA, YVAL1, YVAL2, YCOMP
 REAL(DP), DIMENSION(:), ALLOCATABLE :: XVAL, XA, A
 REAL(DP) :: L, DR, X0, H0, S, THETA, ANGLE, X(3,3), Y(3,3)
 INTEGER(I4B) :: M, NDER
 COMPLEX(DPC) :: Z(3,3), IM
 !--------
 WRITE(*,*) "Enter number of data points N and number of interpolation points M"
 READ(*,*) N, M
 !--------
 OPEN(UNIT=10,FILE='test/test_f_analytic.dat')
 OPEN(UNIT=20,FILE='test/test_f_numeric.dat')
 OPEN(UNIT=30,FILE='test/test_f_numeric2.dat')
 OPEN(UNIT=40,FILE='test/test_f_comp.dat')
 !--------
 NDER = 2
 L = 5000.0_DP
 ANGLE = 40.0_DP
 THETA = ANGLE*PI/180.0_DP
 X0 = L/2.0_DP
 S = L/5.0_DP
 H0 = SQRT(EXP(1.0_DP)/2.0_DP)*S*ATAN(THETA)
 DR = L/N
 ALLOCATE(XA(N+1), YA(N+1,0:NDER))
 XA = (/ (I*DR, I=0,N) /)
 DO I=1,N+1
  YA(I,:) = GHILL(H0,X0,S,XA(I))
  WRITE(10,100) XA(I), YA(I,:)
 END DO
 !--------
 DR = L/M
 ALLOCATE(XVAL(M+1),YVAL1(M+1,0:NDER),YVAL2(M+1,0:NDER),YCOMP(M+1,1:NDER))
 XVAL = (/ (I*DR, I=0,M) /)
 YVAL1 = POLINT1(XVAL,NDER,XA,YA)
 YVAL2(:,0) = INTERP1(XA,YA,XVAL,.TRUE.,0.0_DP)
 YVAL2(:,1) = DIFF1(XA,YA(:,0))
 YVAL2(:,2) = DIFF2(XA,YA(:,0))
 DO I=3,N-1
  DO J=1,NDER
   YCOMP(I,J) = YVAL2(I,J)/YA(I,J)
  END DO
  WRITE(20,100) XVAL(I), YVAL1(I,:)
  WRITE(30,100) XVAL(I), YVAL2(I,:)
  WRITE(40,100) XVAL(I), YCOMP(I,:)
 END DO
 !--------
 100 FORMAT(4(F16.5,X))
 !--------
 ALLOCATE(C(0:4,0:8,0:8),A(0:8))
 A = (/0,1,-1,2,-2,3,-3,4,-4/)
 !A = (/-4,-3,-2,-1,0,1,2,3,4/)
 C = FDW(4,8,0.0_DP,A)
 DO J = 0,4
  WRITE(*,*) "DERIVATIVE ORDER = ", J
  WRITE(*,*) "-----------------------"
  DO I = 1,4
   WRITE(*,101) 2*I, C(J,2*I,:)
  END DO
  WRITE(*,*) "-----------------------"
  WRITE(*,*)
 END DO
 !--------
 101 FORMAT(I1,X,9(F10.3,X))
 !--------
 IM = CMPLX(0.0_DP,1.0_DP,KIND=DP)
 DO I = 1,3
  DO J = 1,3
   Z(I,J) = I+J*IM
  END DO
 END DO
 X = REAL(Z)
 Y = AIMAG(Z)
 WRITE(*,*) "Z = ", Z
 WRITE(*,*) "RE(Z) = ", X
 WRITE(*,*) "IM(Z) = ", Y
 !--------
END PROGRAM
