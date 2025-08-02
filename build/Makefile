.PHONY: all ad clean ensure_fautodiff
.RECIPEPREFIX := >

FC = gfortran
FFLAGS ?= -O2 -ffree-line-length-none

all: shallow_water_test1 shallow_water_forward shallow_water_reverse

shallow_water_test1:
>$(FC) $(FFLAGS) ../src/common/constants_module.f90 ../src/common/variables_module.f90 ../src/common/equations_module.f90 ../src/common/rk4_module.f90 ../src/common/io_module.f90 ../src/common/cost_functions.f90 ../src/testcase1/shallow_water_test1.f90 -o $@

shallow_water_forward: ad ../src/testcase1/shallow_water_forward.f90
>$(FC) $(FFLAGS) -c fautodiff_stack.f90
>$(FC) $(FFLAGS) -c ../src/common/constants_module.f90
>$(FC) $(FFLAGS) -c ../src/common/cost_functions.f90
>$(FC) $(FFLAGS) -c constants_module_ad.f90
>$(FC) $(FFLAGS) -c cost_module_ad.f90
>$(FC) $(FFLAGS) -I. fautodiff_stack.o constants_module.o cost_functions.o constants_module_ad.o cost_module_ad.o ../src/testcase1/shallow_water_forward.f90 -o $@

shallow_water_reverse: ad ../src/testcase1/shallow_water_reverse.f90
>$(FC) $(FFLAGS) -c fautodiff_stack.f90
>$(FC) $(FFLAGS) -c ../src/common/constants_module.f90
>$(FC) $(FFLAGS) -c ../src/common/cost_functions.f90
>$(FC) $(FFLAGS) -c constants_module_ad.f90
>$(FC) $(FFLAGS) -c cost_module_ad.f90
>$(FC) $(FFLAGS) -I. fautodiff_stack.o constants_module.o cost_functions.o constants_module_ad.o cost_module_ad.o ../src/testcase1/shallow_water_reverse.f90 -o $@

ad: ensure_fautodiff
>../scripts/copy_fautodiff_modules.sh
>../scripts/fautodiff/bin/fautodiff ../src/common/constants_module.f90 -M . -o constants_module_ad.f90
>../scripts/fautodiff/bin/fautodiff ../src/common/cost_functions.f90 -M . -I . -o cost_module_ad.f90

ensure_fautodiff:
>../scripts/setup_fautodiff.sh

clean:
>rm -f *.o *.mod shallow_water_test1 shallow_water_forward shallow_water_reverse *_ad.f90 *.fadmod fautodiff_stack.f90
