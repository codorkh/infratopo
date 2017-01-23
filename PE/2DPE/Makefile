# Modified by Anton Shterenlikht

FC=		gfortran6
FFLAGS=		-g -O3 -fbacktrace -Wall -fcheck=all
LDFLAGS=	-Wl,-rpath="/usr/local/lib/gcc6" -llapack -lblas \
		-L$(HOME)/slatec/src/static -lslatec -llapack -lblas

MOD_SRC=	nrtype.f03 matlib.f03 ground.f03 param.f03
MOD_MOD=	$(MOD_SRC:.f03=.mod)
MOD_OBJ=	$(MOD_SRC:.f03=.o)
CLEAN+=		$(MOD_MOD) $(MOD_OBJ)

PROG_SRC=	PERF2D.f03
PROG_OBJ=	$(PROG_SRC:.f03=.o)
PROG_EXE=	$(PROG_SRC:.f03=.x)
CLEAN+=		$(PROG_OBJ) $(PROG_EXE)

.SUFFIXES:
.SUFFIXES:	.f03 .mod .o

all:		$(MOD_MOD) $(MOD_OBJ) $(PROG_OBJ) $(PROG_EXE)

.f03.mod:
	$(FC) -c $(FFLAGS) $<
	touch $@

.f03.o:
	$(FC) -c $(FFLAGS) $<

$(PROG_OBJ):	$(PROG_SRC) $(MOD_MOD)
	$(FC) -c $(FFLAGS) $(PROG_SRC)

$(PROG_EXE):	$(PROG_OBJ) $(MOD_OBJ)
	$(FC) -o $@ $(FFLAGS) $(PROG_OBJ) $(MOD_OBJ) $(LDFLAGS)

clean :
	rm $(CLEAN)

