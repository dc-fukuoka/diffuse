fc       = ifort
mpifc    = mpiifort
fppflags = -fpp
fflags   = -g -O3 -mavx
ldflags  =
libs     =

all: create_input create_output diffuse
create_input: create_input.F90
	$(fc) $(fflags) $^ -o $@
create_output: create_output.F90
	$(fc) $(fflags) $^ -o $@
diffuse: diffuse.F90
	$(mpifc) $(fflags) $^ -o $@
clean:
	rm -f f create_input create_output diffuse *.mod *~ core.* data_in data_out
	find . -type f -name "fort.*"|grep -v fort.11|xargs rm	
