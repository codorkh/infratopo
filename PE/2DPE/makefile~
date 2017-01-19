FC = gfortran
FFLAGS = -g -O3 -fbacktrace -Wall
LDFLAGS = -L$HOME/slatec/src/static -lslatec -llapack -lblas

all : PE2D_HFR

PE2D_HFR : PE2D_HFR.o matlib.o fdlib.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o PE2D_HFR

PE2D_HFR.o : PE2D_HFR.f03 matlib.mod fdlib.mod param.mod
	$(FC) $(FFLAGS) $(LDFLAGS) -c PE2D_HFR.f03

param.mod param.o : param.f03
	$(FC) $(FFLAGS) -c param.f03

matlib.mod matlib.o : matlib.f03 param.mod
	$(FC) $(FFLAGS) -llapack -lblas -c matlib.f03

fdlib.mod fdlib.o : fdlib.f03 param.mod
        $(FC) $(FFLAGS) -c fdlib.f03

clean :
	rm matlib.o fdlib.o param.o matlib.mod fdlib.mod param.mod PE2D_HFR

