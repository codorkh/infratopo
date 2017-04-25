#========================================
#
#	 Include Makefile Flags
#
#========================================

FC=		gfortran
FCFLAGS=	-g -ftraceback -fcheck=all -Wall
#FCFLAGS=	-O2 -funroll-loops
LDFLAGS=	-Wl,-rpath="/usr/local/lib/gcc" -L$(HOME)/lib \
		-lslatec -lsparseblas -llapack -lblas

