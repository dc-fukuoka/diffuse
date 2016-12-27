all: create_input create_output d
create_input: create_input.F90
	ifort -g -traceback create_input.F90 -o create_input
create_output: create_output.F90
	ifort -g -traceback create_output.F90 -o create_output
d: d.F90
	mpif90 -fpe0 -D_DEBUG -g -traceback d.F90 -o d
clean:
	rm -f f create_input create_output d *.mod *~ core.* data_in data_out
	find . -type f -name "fort.*"|grep -v fort.11|xargs rm	
