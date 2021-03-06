MODULE NRTYPE
	IMPLICIT NONE
 	!-----------------------------------------------------------
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_SP
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_SP
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_SP
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_SP
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_SP
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_DP
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_DP
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_DP
	TYPE SPRS2_SP
		INTEGER(I4B) :: N,NELT
		REAL(SP), DIMENSION(:), POINTER :: VAL
		INTEGER(I4B), DIMENSION(:), POINTER :: IROW
		INTEGER(I4B), DIMENSION(:), POINTER :: JCOL
	END TYPE SPRS2_SP
	TYPE SPRS2_DP
		INTEGER(I4B) :: N,NELT
		REAL(DP), DIMENSION(:), POINTER :: VAL
		INTEGER(I4B), DIMENSION(:), POINTER :: IROW
		INTEGER(I4B), DIMENSION(:), POINTER :: JCOL
	END TYPE SPRS2_DP
	 !-----------------------------------------------------------
END MODULE NRTYPE

