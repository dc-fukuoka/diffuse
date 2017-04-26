fc       = ifort
mpifc    = mpiifort
fppflags = -fpp -D_OVERLAP_MPI
fflags   = -g -O3 -mavx
openmp   = -fopenmp
ldflags  =
libs     =

all: create_input create_output diffuse
create_input: create_input.F90
	$(fc) $(openmp) $(fflags) $^ -o $@
create_output: create_output.F90
	$(fc) $(openmp) $(fflags) $^ -o $@
diffuse: diffuse.F90
	$(mpifc) $(fppflags) $(fflags) $^ -o $@
clean:
	rm -f f create_input create_output diffuse *.mod *~ core.* data_in data_out
	find . -type f -name "fort.*"|grep -v fort.11|xargs rm	
