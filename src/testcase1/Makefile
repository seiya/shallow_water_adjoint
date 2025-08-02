FC=gfortran
FFLAGS=-O2

COMMON=../common/constants_module.f90 ../common/variables_module.f90 \
       ../common/equations_module.f90 ../common/rk4_module.f90 \
       ../common/io_module.f90 ../common/cost_functions.f90

all: shallow_water_test1

shallow_water_test1: shallow_water_test1.f90 $(COMMON)
	$(FC) $(FFLAGS) -I../common -o $@ \
	$(COMMON) shallow_water_test1.f90

clean:
	rm -f shallow_water_test1 error.dat cost.log grid_params.txt snapshot_*.bin snapshot_*.png

.PHONY: all clean
