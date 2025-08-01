FC=gfortran
FFLAGS=-O2

all: shallow_water_test1

shallow_water_test1: shallow_water_test1.f90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f shallow_water_test1 error.dat

.PHONY: all clean
