all: create_input d.impi
create_input: create_input.F90
	ifort -g -traceback create_input.F90 -o create_input
d.impi: d.F90
	mpiifort -D_DEBUG -g -traceback d.F90 -o d.impi
clean:
	rm -f fort.20* fort.30* fort.1000 create_input d.impi *.mod *~ core.* data_in
