CC       = icc
FC       = ifort
CFLAGS   = -g -O3 -mavx
CPPFLAGS = 
FPPFLAGS = 
LDFLAGS  = 
LIBS     = -lrt

ALL: cg expl check

dclock.o: dclock.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

cg: cg.F90 dclock.o
	$(FC) $(FPPFLAGS) $(CFLAGS) $(OPENMP) $(LDFLAGS) $(LIBS) $< dclock.o -o $@
expl: expl.F90 dclock.o
	$(FC) $(FPPFLAGS) $(CFLAGS) $(OPENMP) $(LDFLAGS) $(LIBS) $< dclock.o -o $@
check: check.f90
	$(FC) $(FPPFLAGS) $(CFLAGS) $(OPENMP) $< -o $@
clean:
	rm -f cg expl check *.o *.mod *~
