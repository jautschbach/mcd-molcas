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

BIN = mcd-c-molcas mcd-a-molcas mcd-b-molcas plot-spectrum

OBJ = 

OBJ90 = definitions.o \
        print-rec-matrix.o \
        diagonalize-matrix.o \
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

plot-spectrum: plot-spectrum.o
	$(LINKER) -o plot-spectrum plot-spectrum.o $(FFLAGS) $(FLIB)

clean:
	rm -f $(BIN) *.o *.mod



