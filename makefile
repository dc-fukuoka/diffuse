all: create_input create_output diffuse
create_input: create_input.F90
	ifort -g -O3 -mavx create_input.F90 -o create_input
create_output: create_output.F90
	ifort -g -O3 -mavx create_output.F90 -o create_output
diffuse: diffuse.F90
	mpif90 -g -O3 -mavx diffuse.F90 -o diffuse
clean:
	rm -f f create_input create_output diffuse *.mod *~ core.* data_in data_out
	find . -type f -name "fort.*"|grep -v fort.11|xargs rm	
