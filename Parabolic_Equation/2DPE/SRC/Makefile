# Modified by Anton Shterenlikht

FC=		gfortran
FFLAGS=		-g -O3 -fbacktrace -Wall -fcheck=all
LDFLAGS=	-Wl,-rpath="/usr/local/lib/gcc" -llapack -lblas \
		-L$(HOME)/slatec/src/static -lslatec -llapack -lblas

MOD_SRC=	PE2D_TYPE.f03 PE2D_VAR.f03 PE2D_GROUND.f03 PE2D_ATMOS.f03 PE2D_AUX.f03
MOD_MOD=	$(MOD_SRC:.f03=.mod)
MOD_OBJ=	$(MOD_SRC:.f03=.o)
CLEAN+=		$(MOD_MOD) $(MOD_OBJ)

PROG1_SRC=	PE2D_HT_FD.f03
PROG1_OBJ=	$(PROG1_SRC:.f03=.o)
PROG1_EXE=	$(PROG1_SRC:.f03=.x)
CLEAN+=		$(PROG1_OBJ) $(PROG1_EXE)

PROG2_SRC=	PE2D_HF_FD.f03
PROG2_OBJ=	$(PROG2_SRC:.f03=.o)
PROG2_EXE=	$(PROG2_SRC:.f03=.x)
CLEAN+=		$(PROG2_OBJ) $(PROG2_EXE)

.SUFFIXES:	.f03 .mod .o

all:		$(PROG1_EXE) $(PROG2_EXE)

$(MOD_MOD) $(MOD_OBJ): $(MOD_SRC)
$(PROG1_OBJ) $(PROG2_OBJ): $(MOD_MOD)

.f03.mod:
	$(FC) -c $(FFLAGS) $<
	touch $@

.f03.o:
	$(FC) -c $(FFLAGS) $<

$(PROG1_EXE):	$(PROG1_OBJ) $(MOD_OBJ)
	$(FC) -o $@ $(FFLAGS) $(PROG1_OBJ) $(MOD_OBJ) $(LDFLAGS)

$(PROG2_EXE):	$(PROG2_OBJ) $(MOD_OBJ)
	$(FC) -o $@ $(FFLAGS) $(PROG2_OBJ) $(MOD_OBJ) $(LDFLAGS)

clean :
	rm $(CLEAN)

