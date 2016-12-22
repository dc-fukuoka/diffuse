all: create_input d
create_input: create_input.F90
	ifort -g -traceback create_input.F90 -o create_input
d: d.F90
	mpiifort -D_DEBUG -g -traceback d.F90 -o d
clean:
	rm -f fort.20* create_input d *.mod *~ data_in
