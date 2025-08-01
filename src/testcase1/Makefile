FC=gfortran
FFLAGS=-O2

all: shallow_water_test1

shallow_water_test1: shallow_water_test1.f90 ../cost_functions/cost_module.f90
	$(FC) $(FFLAGS) -I../cost_functions -o $@ ../cost_functions/cost_module.f90 shallow_water_test1.f90

clean:
	rm -f shallow_water_test1 error.dat snapshot_*.bin snapshot_*.png

.PHONY: all clean
