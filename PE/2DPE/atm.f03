MODULE ATM
    USE NRTYPE
    USE NRINT
    IMPLICIT NONE
    !-----------------------------------------------------------
    CONTAINS
    !-----------------------------------------------------------
    SUBROUTINE ATMLIN(N,K,C,C0,F,DZ,A)
    IMPLICIT NONE
    INTEGER :: N
    REAL, INTENT(IN) :: F, C0, DZ, A
    REAL, DIMENSION(:), INTENT(INOUT) :: K(N), C(N)
    INTEGER :: I
    DO I = 1,N
        C(I) = C0+I*DZ*A
        K(I) = 2*PI*F/C(I)
    END DO
    END SUBROUTINE ATMLIN
    !-----------------------------------------------------------
    SUBROUTINE ATMLOG(N,K,C,C0,F,DZ,Z0,A)
    IMPLICIT NONE
    INTEGER :: N
    REAL, INTENT(IN) :: F, C0, DZ, A, Z0
    REAL, , DIMENSION(:), INTENT(INOUT) :: K(N), C(N)
    INTEGER :: I
    DO I = 1,N
        C(I) = C0+A*LOG(1+I*DZ/Z0)
        K(I) = 2*PI*F/C(I)
    END DO
    END SUBROUTINE ATMLOG
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
END MODULE ATM
