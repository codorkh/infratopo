MODULE PE2D_ATMOS
 USE PE2D_TYPE
 IMPLICIT NONE
 !-----------------------------------------------------------
 CONTAINS
 !-----------------------------------------------------------
 FUNCTION ATMLIN(C0,A,Z) RESULT(C)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: C0, A, Z
 REAL(DP) :: C
 C = C0+Z*A
 END FUNCTION ATMLIN
 !-----------------------------------------------------------
 FUNCTION ATMLOG(C0,A,Z0,Z) RESULT(C)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: C0, A, Z0, Z
 REAL(DP) :: C
 C = C0+A*LOG(1+Z/Z0)
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
END MODULE PE2D_ATMOS
