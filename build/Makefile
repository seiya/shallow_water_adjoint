.PHONY: all ad clean update_fautodiff
.RECIPEPREFIX := >

FC = gfortran
FFLAGS ?= -O2 -ffree-line-length-none
FFLAGS += -ffpe-trap=invalid,zero,overflow -fbounds-check -finit-real=nan -g -fbacktrace
FFLAGS += -fopenmp

COMMON_OBJS = constants_module.o \
              variables_module.o \
			  equations_module.o \
			  rk4_module.o \
			  io_module.o \
			  cost_module.o

COMMON_AD_OBJS = $(patsubst %.o,%_ad.o,$(COMMON_OBJS))

FAUTODIFF_OBJS = fautodiff_stack.o
FAUTODIFF_SRC  = fautodiff_stack.f90

VPATH = ../src/common

all: update_fautodiff \
     shallow_water_test1.out shallow_water_test1_forward.out shallow_water_test1_reverse.out \
     shallow_water_test2.out shallow_water_test2_forward.out shallow_water_test2_reverse.out \
     shallow_water_test5.out shallow_water_test5_forward.out shallow_water_test5_reverse.out

shallow_water_test1.out: ../src/testcase1/shallow_water_test1.f90 $(COMMON_OBJS)
>$(FC) $(FFLAGS) -o $@ $^

shallow_water_test1_forward.out: ../src/testcase1/shallow_water_test1_forward.f90 $(COMMON_OBJS) $(COMMON_AD_OBJS) $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -I. -o $@ $^

shallow_water_test1_reverse.out: ../src/testcase1/shallow_water_test1_reverse.f90 $(COMMON_OBJS) $(COMMON_AD_OBJS) $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -I. -o $@ $^

shallow_water_test2.out: ../src/testcase2/shallow_water_test2.f90 $(COMMON_OBJS)
>$(FC) $(FFLAGS) -o $@ $^

shallow_water_test2_forward.out: ../src/testcase2/shallow_water_test2_forward.f90 $(COMMON_OBJS) $(COMMON_AD_OBJS) $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -I. -o $@ $^

shallow_water_test2_reverse.out: ../src/testcase2/shallow_water_test2_reverse.f90 $(COMMON_OBJS) $(COMMON_AD_OBJS) $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -I. -o $@ $^

shallow_water_test5.out: ../src/testcase5/shallow_water_test5.f90 $(COMMON_OBJS)
>$(FC) $(FFLAGS) -o $@ $^

shallow_water_test5_forward.out: ../src/testcase5/shallow_water_test5_forward.f90 $(COMMON_OBJS) $(COMMON_AD_OBJS) $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -I. -o $@ $^

shallow_water_test5_reverse.out: ../src/testcase5/shallow_water_test5_reverse.f90 $(COMMON_OBJS) $(COMMON_AD_OBJS) $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -I. -o $@ $^

%_ad.o: %_ad.f90 $(FAUTODIFF_OBJS)
>$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
>$(FC) $(FFLAGS) -c $< -o $@

%_ad.f90: %.f90 $(FAUTODIFF_SRC)
>../scripts/fautodiff/bin/fautodiff $< -I . -M . -o $@

ad: update_fautodiff

update_fautodiff:
>../scripts/setup_fautodiff.sh
>../scripts/copy_fautodiff_modules.sh

$(FAUTODIFF_SRC):
>../scripts/copy_fautodiff_modules.sh

clean:
>rm -f *.o *.mod *.out *_ad.f90 *.fadmod fautodiff_stack.f90

.PRECIOUS: %_ad.f90
