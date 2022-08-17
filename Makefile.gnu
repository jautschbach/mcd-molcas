#
# This program may need lapack and blas routines
#
# For MPI support, we obviously also need the MPI libs and include files
#
# Note: we define two compiler variables: $FC for fixed format F77 or F90 and $FC90 for
# free format F90


### Fortran compiler
FC     = gfortran
FC90   = gfortran

LINKER = gfortran

### Disable/Enable OpenMP support
OMP = 

### Compilation/Linking options
FOPT = -g -Wall -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow # gfortran 

### Extra libraries
FLIB=-llapack

FFLAGS = $(OMP)

.SUFFIXES: .f .F .F90 .f90

BIN = mcd-c-molcas mchd-c-molcas mcd-a-molcas \
      mcd-b-molcas transition-dip-rot plot-mcdspectrum 

OBJ = 

OBJ90 = definitions.o \
        namelist-module.o \
        constants-parameters.o \
        shared-variables.o \
        read-options.o \
        read-data-files.o \
        print-rec-matrix.o \
        diagonalize-matrix.o \
        diagonalize-magdip-gs.o \
        print-constants.o \
        vector-product.o \
        process-energies.o


.F.o:
	$(FC) $(FOPT) -c $< -o $@
.F90.o: 
	$(FC90) $(FOPT) -c $< -o $@
.f90.o: 
	$(FC90) $(FOPT) -c $< -o $@

all: $(BIN)

mcd-c-molcas: $(OBJ) $(OBJ90) mcd-c-molcas.o
	$(LINKER) -o mcd-c-molcas mcd-c-molcas.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

mcd-a-molcas: $(OBJ) $(OBJ90) mcd-a-molcas.o
	$(LINKER) -o mcd-a-molcas mcd-a-molcas.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

mcd-b-molcas: $(OBJ) $(OBJ90) mcd-b-molcas.o
	$(LINKER) -o mcd-b-molcas mcd-b-molcas.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

mchd-c-molcas: $(OBJ) $(OBJ90) mchd-c-molcas.o
	$(LINKER) -o mchd-c-molcas mchd-c-molcas.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

transition-dip-rot: $(OBJ) $(OBJ90) transition-dip-rot.o
	$(LINKER) -o transition-dip-rot transition-dip-rot.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

plot-mcdspectrum: plot-mcdspectrum.o
	$(LINKER) -o plot-mcdspectrum plot-mcdspectrum.o $(FFLAGS) $(FLIB)

clean:
	rm -f $(BIN) *.o *.mod



