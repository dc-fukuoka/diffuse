*diffuse* - 3D diffusion equation solver by using CG method(no precondition) with MPI.
===============	  
the result was compared to serial version(cg.F90)  
  
input file: fort.11

    $ cat fort.11 
    &param0
    dt        = 1.0d-6  ! stride of time
    dx        = 1.0d-2  ! stride of x, y, z direction
    tol       = 1.0d-40 ! convergence tolerance of CG method
    diff_coef = 1.0d2   ! diffusion coefficient
    /
    &param1
    imax = 64 ! global size for x direction
    jmax = 64 ! global size for y direction
    kmax = 64 ! global size for z direction
    /
    &param2
    idiv = 1 ! # of division for x direction
    jdiv = 1 ! # of division for y direction
    kdiv = 4 ! # of division for z direction
    ! idiv*jdiv*kdiv must equal to total # of processes(np)
    ! imax/idiv, jmax/jdiv and kmax/kdiv must be an integer
    /
    
    &param3
    iter_max   = 262144 ! maximum CG method iteration
    tstep_max  = 100    ! maximum timestep
    freq_write = 100    ! frequency of writing the result, write each tstep_max/freq_write time
    /


how to run:  
    
    $ make # if you have intel compiler and intelmpi, try "make -f makefile.impi.mk"  
    $ vi fort.11 # adjust the parameters  
    $ ./create_input  
    $ mpirun -np $NP ./diffuse # or mpirun -np $NP ./diffuse_hyb where $NP = idiv*jdiv*kdiv  
    $ ./create_output  
    $ ./create_anime.sh  
  
performance comparison:  
cpu: Intel(R) Xeon(R) CPU E5-2450 0 @ 2.10GHz  
problem size: 128x128x128  
tstep_max   : 100  
freq_write  :  1  
with 1   core           : 381.2 s  
with 32 cores(flat MPI) : 20.2  s  
with 32 cores(hbyrid)   : 15.8  s  

compare the results between implicit method and explicit method
1.  
    dt        = 1.0d-6  ! stride of time
    dx        = 1.0d-2  ! stride of x, y, z direction
    tol       = 1.0d-40 ! convergence tolerance of CG method
    diff_coef = 1.0d2   ! diffusion coefficient


implicit method  
![Alt text](gifs/diffuse.1.impl.gif?raw=true "implicit method 1")
explicit method  
![Alt text](gifs/diffuse.1.expl.gif?raw=true "explicit method 1")


2.  
    dt        = 1.0d-6  ! stride of time
    dx        = 1.0d-3  ! stride of x, y, z direction
    tol       = 1.0d-40 ! convergence tolerance of CG method
    diff_coef = 1.0d2   ! diffusion coefficient


implicit method  
![Alt text](gifs/diffuse.2.impl.gif?raw=true "implicit method 2")
explicit method  
![Alt text](gifs/diffuse.2.expl.gif?raw=true "explicit method 2")
