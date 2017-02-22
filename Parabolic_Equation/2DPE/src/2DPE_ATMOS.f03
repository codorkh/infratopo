MODULE 2DPE_ATMOS
 USE NRTYPE
 USE NRINT
 IMPLICIT NONE
 !-----------------------------------------------------------
 CONTAINS
 !-----------------------------------------------------------
 FUNCTION ATMLIN(C0,A,Z) RESULT(C)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: C0, A, Z
 REAL(DP) :: C
 C(I) = C0+I*DZ*A
 END END FUNCTION ATMLIN
 !-----------------------------------------------------------
 FUNCTION ATMLOG(C0,A,Z0,Z) RESULT(C)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: C0, A, Z0, Z
 REAL(DP) :: C
 INTEGER :: I
 C(I) = C0+A*LOG(1+I*DZ/Z0)
 END FUNCTION ATMLOG
 !-----------------------------------------------------------   
 ! SUBROUTINE ATMREAL(N,K,C,NATM)
 ! IMPLICIT NONE
 ! INTEGER :: I, N, NATM
 ! REAL, DIMENSION(:), INTENT(INOUT) :: K(N), C(N)
 ! REAL, ALLOCATABLE, DIMENSION(:) :: CRAW, RANGE
 ! I = 0
 ! IF (NATM==1) THEN
 !     OPEN(UNIT=200,FILE="atmos_2d_adiabatic_sound_sp.dat",ACTION="READ")
 !     DO
 !         I = I+1
 !         READ(UNIT=200,FMT=*,IOSTAT=IOS) RANGE(I), C(I)
 !         IF (IOS/=0) THEN
 !     END DO
 ! ELSE

 ! END SUBROUTINE ATMREAL
 !-----------------------------------------------------------
END MODULE 2DPE_ATMOS
