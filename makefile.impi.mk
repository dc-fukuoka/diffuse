fc       = ifort
mpifc    = mpiifort
fppflags = -fpp -D_OVERLAP_MPI
fflags   = -g -O3 -mavx
openmp   = -fopenmp
ldflags  =
libs     =

all: create_input create_output diffuse diffuse_hyb
create_input: create_input.F90
	$(fc) $(openmp) $(fflags) $^ -o $@
create_output: create_output.F90
	$(fc) $(openmp) $(fflags) $^ -o $@
diffuse: diffuse.F90
	$(mpifc) $(fppflags) $(fflags) $^ -o $@
diffuse_hyb: diffuse.F90
	$(mpifc) $(fppflags) $(openmp) $(fflags) $^ -o $@
clean:
	rm -f create_input create_output diffuse diffuse_hyb diffuse.gif *.mod *~ core.* data_in data_out fort.100 fort.599 fort.600 fort.777 fort.2222 fort.9999
