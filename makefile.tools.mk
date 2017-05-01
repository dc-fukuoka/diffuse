CC       = gcc
FC       = gfortran
CFLAGS   = -g -O3 -mavx
CPPFLAGS = 
FPPFLAGS = 
LDFLAGS  = 
LIBS     = -lrt
OPENMP   = -fopenmp

ALL: cg cg_cr expl check

dclock.o: dclock.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

cg: cg.F90 dclock.o
	$(FC) $(FPPFLAGS) $(CFLAGS) $(OPENMP) $(LDFLAGS) $(LIBS) $< dclock.o -o $@
cg_cr: cg.F90 dclock.o
	$(FC) $(FPPFLAGS) -D_CR $(CFLAGS) $(OPENMP) $(LDFLAGS) $(LIBS) $< dclock.o -o $@
expl: expl.F90 dclock.o
	$(FC) $(FPPFLAGS) $(CFLAGS) $(OPENMP) $(LDFLAGS) $(LIBS) $< dclock.o -o $@
check: check.f90
	$(FC) $(FPPFLAGS) $(CFLAGS) $(OPENMP) $< -o $@
clean:
	rm -f cg cg_cr expl check *.o *.mod *~
